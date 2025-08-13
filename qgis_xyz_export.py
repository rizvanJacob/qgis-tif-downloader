"""
qgis_xyz_export.py — Export XYZ tiles and build DEFLATE-compressed GeoTIFF mosaics (one-step).

Load in QGIS Python Console:
    exec(open(r"F:\qgis_scripts\qgis_xyz_export.py").read())

Single-step export → mosaics:
    export_bbox_to_geotiff_tiles(
        layer_name="Google Satellite Images",
        top_left=(-117, 44), bottom_right=(-112, 40),
        zoom=12,
        out_dir=r"F:\TIF\mosaic",
        # Mosaic controls:
        tiles_per_side=128,          # try 128 first; raise to 144/160 if you want bigger files
        mosaic_internal_tiled=False,  # strongly recommended for multi-GB files
        mosaic_blocksize=512,        # 256 or 512
        # Cache controls:`
        tile_cache_dir=None,         # default: <out_dir>\_tilecache
        keep_tile_cache=False,       # set True to keep per-tile GeoTIFFs
        # Perf/codec:
        zlevel=9,                    # DEFLATE level for both per-tile and mosaic
        chunk_size=20, pause_secs=0.7,
        verbose=False
    )

Notes:
- Mosaics are written directly into out_dir as files like:
      z19_x<start>-<end>_y<start>-<end>.tif
- CRS is EPSG:3857 (Web Mercator). We pass WKT directly to avoid EPSG lookup issues.
- Per-tile GeoTIFFs in cache are **strip-based** (no internal tiling).
- Mosaic GeoTIFFs are **internally tiled** by default (faster IO for big files).
"""

import os, math, time
from qgis.core import (
    QgsProject, QgsCoordinateReferenceSystem, QgsRectangle,
    QgsMapSettings, QgsMapRendererCustomPainterJob, QgsApplication, Qgis
)
from qgis.PyQt.QtGui import QImage, QPainter, QColor
from qgis.PyQt.QtCore import QSize
from osgeo import gdal

# -------------------------------------------
# Globals / defaults
# -------------------------------------------
DRY_RUN_MODE = False
DEFAULT_ZOOM_FALLBACK = 19
DEFAULT_TILE_SIZE = 256  # set 512 if your XYZ source serves 512px tiles

WEBMERC = QgsCoordinateReferenceSystem("EPSG:3857")
R = 6378137.0
ORIGIN_SHIFT = math.pi * R

if not hasattr(QImage, 'Format_ARGB32'):
    try:
        # PyQt6 style
        QImage.Format_ARGB32 = QImage.Format.Format_ARGB32
    except AttributeError:
        raise RuntimeError("Could not find QImage.Format_ARGB32 in this PyQt version.")

# -------------------------------------------
# PROJ/GDAL environment hardening (Windows conflicts)
# -------------------------------------------
def _qgis_osgeo4w_roots():
    prefix = QgsApplication.prefixPath()  # e.g. D:\Non_Vital_Programs\OSGeo4W\apps\qgis
    root = os.path.abspath(os.path.join(prefix, r"..", r".."))
    share = os.path.join(root, "share")
    proj  = os.path.join(share, "proj")
    gdaldata = os.path.join(share, "gdal")
    return root, share, proj, gdaldata

def _ensure_gdal_proj_env():
    _, _, proj, gdaldata = _qgis_osgeo4w_roots()
    os.environ["PROJ_LIB"]  = proj
    os.environ["PROJ_DATA"] = proj
    os.environ["GDAL_DATA"] = gdaldata
    gdal.SetConfigOption("PROJ_LIB",  proj)
    gdal.SetConfigOption("PROJ_DATA", proj)
    gdal.SetConfigOption("GDAL_DATA", gdaldata)

_ensure_gdal_proj_env()

def _webmerc_wkt():
    try:
        return WEBMERC.toWkt(Qgis.WktVariant.WKT2_2018)
    except Exception:
        return WEBMERC.toWkt()

# -------------------------------------------
# Helpers — XYZ / WebMercator math & layer utils
# -------------------------------------------
def clamp_lat(lat):  return max(min(lat, 85.05112878), -85.05112878)

def lonlat_to_tilexy(lon, lat, z):
    lat = clamp_lat(lat)
    n = 2 ** z
    xtile = (lon + 180.0) / 360.0 * n
    ytile = (1.0 - math.log(math.tan(math.radians(lat)) + 1/math.cos(math.radians(lat))) / math.pi) / 2.0 * n
    return xtile, ytile

def merc_bounds_for_tile(x, y, z):
    n = 2 ** z
    minx = -ORIGIN_SHIFT + (x / n) * (2 * ORIGIN_SHIFT)
    maxx = -ORIGIN_SHIFT + ((x + 1) / n) * (2 * ORIGIN_SHIFT)
    miny = ORIGIN_SHIFT - ((y + 1) / n) * (2 * ORIGIN_SHIFT)
    maxy = ORIGIN_SHIFT - (y / n) * (2 * ORIGIN_SHIFT)
    return (minx, miny, maxx, maxy)

def find_layer_by_name(layer_name):
    for lyr in QgsProject.instance().mapLayers().values():
        if lyr.name() == layer_name:
            return lyr
    return None

def guess_layer_max_zoom(rlayer, fallback=DEFAULT_ZOOM_FALLBACK):
    for k in ("zmax", "maxZoom", "maximumZoomLevel"):
        v = rlayer.customProperty(k, None)
        if v is not None:
            try: return int(v)
            except Exception: pass
    return fallback

def bbox_to_tile_range(top_left, bottom_right, zoom):
    min_lon, min_lat = bottom_right
    max_lon, max_lat = top_left
    if min_lon > max_lon: min_lon, max_lon = max_lon, min_lon
    if min_lat > max_lat: min_lat, max_lat = max_lat, min_lat
    tlx, tly = lonlat_to_tilexy(min_lon, max_lat, zoom)
    brx, bry = lonlat_to_tilexy(max_lon, min_lat, zoom)
    x_min = math.floor(min(tlx, brx)); x_max = math.floor(max(tlx, brx))
    y_min = math.floor(min(tly, bry)); y_max = math.floor(max(tly, bry))
    return x_min, x_max, y_min, y_max

# -------------------------------------------
# Single-tile renderer (strip-based GeoTIFF into cache)
# -------------------------------------------
def _render_one_tile(rlayer, z, x, y, out_dir, zlevel=9, verbose=True):
    minx, miny, maxx, maxy = merc_bounds_for_tile(x, y, z)
    tile_dir = os.path.join(out_dir, str(z), str(x))
    os.makedirs(tile_dir, exist_ok=True)
    tmp_png = os.path.join(tile_dir, f"{y}_tmp.png")
    out_tif = os.path.join(tile_dir, f"{y}.tif")

    if os.path.exists(out_tif):
        if verbose: print(f"[tile skip] exists -> {out_tif}")
        return out_tif

    width = height = DEFAULT_TILE_SIZE
    ms = QgsMapSettings()
    ms.setLayers([rlayer])
    ms.setDestinationCrs(WEBMERC)
    ms.setExtent(QgsRectangle(minx, miny, maxx, maxy))
    ms.setOutputSize(QSize(width, height))

    img = QImage(width, height, QImage.Format_ARGB32)
    img.fill(QColor("white"))
    painter = QPainter(img)
    job = QgsMapRendererCustomPainterJob(ms, painter)
    job.start(); job.waitForFinished(); painter.end()

    if verbose: print(f"[render] -> {tmp_png}")
    if not img.save(tmp_png, "PNG"):
        raise RuntimeError(f"Failed to save temp PNG: {tmp_png}")

    wkt = _webmerc_wkt()
    translate_opts = gdal.TranslateOptions(
        options = [
            "-a_srs", wkt,
            "-a_ullr", str(minx), str(maxy), str(maxx), str(miny),
            "-co", "COMPRESS=DEFLATE",
            "-co", "PREDICTOR=2",
            "-co", "BIGTIFF=IF_SAFER",
            "-co", f"ZLEVEL={zlevel}",
        ],
        format = "GTiff"
    )
    if verbose: print(f"[compress] gdal.Translate -> {out_tif}")
    ds_out = gdal.Translate(out_tif, tmp_png, options=translate_opts)
    if ds_out is None:
        raise RuntimeError("gdal.Translate returned None.")
    ds_out = None

    try:
        os.remove(tmp_png)
        aux = tmp_png + ".aux.xml"
        if os.path.exists(aux): os.remove(aux)
    except Exception:
        pass

    return out_tif

# -------------------------------------------
# Mosaic helpers (VRT -> GTiff)
# -------------------------------------------
def _collect_tile_paths(cache_dir, z, x_min, x_max, y_min, y_max):
    paths = []
    for x in range(x_min, x_max + 1):
        x_dir = os.path.join(cache_dir, str(z), str(x))
        if not os.path.isdir(x_dir): continue
        for y in range(y_min, y_max + 1):
            tif = os.path.join(x_dir, f"{y}.tif")
            if os.path.exists(tif):
                paths.append(tif)
    return paths

def _build_mosaic(tile_paths, mosaic_path, internal_tiled=True, blocksize=512,
                  compress="DEFLATE", zlevel=9, predictor=2, bigtiff="YES"):
    if len(tile_paths) == 0:
        return None
    vrt_opts = gdal.BuildVRTOptions(resolution="highest")
    vrt_ds = gdal.BuildVRT("", tile_paths, options=vrt_opts)
    if vrt_ds is None:
        raise RuntimeError("gdal.BuildVRT failed.")

    co = [
        f"COMPRESS={compress}",
        f"PREDICTOR={predictor}",
        f"BIGTIFF={bigtiff}",
        f"ZLEVEL={zlevel}",
    ]
    if internal_tiled:
        co += [ "TILED=YES", f"BLOCKXSIZE={blocksize}", f"BLOCKYSIZE={blocksize}" ]

    trans_opts = gdal.TranslateOptions(format="GTiff", creationOptions=co)
    ds_out = gdal.Translate(mosaic_path, vrt_ds, options=trans_opts)
    if ds_out is None:
        raise RuntimeError("gdal.Translate (mosaic) failed.")
    ds_out = None
    vrt_ds = None
    return mosaic_path

# -------------------------------------------
# One-step export: tiles (to cache) → mosaics (to out_dir)
# -------------------------------------------
def export_bbox_to_geotiff_tiles(
    layer_name,
    top_left, bottom_right,
    zoom=None,
    out_dir=".",
    # Mosaic configuration:
    tiles_per_side=128,             # choose 128 first; 144/160 for larger files
    mosaic_internal_tiled=True,
    mosaic_blocksize=512,
    # Cache configuration:
    tile_cache_dir=None,            # default: <out_dir>/_tilecache
    keep_tile_cache=False,
    # Codec & performance:
    zlevel=9,
    chunk_size=20,
    pause_secs=1.5,
    dry_run=None,
    verbose=True
):
    """
    Single-step pipeline:
      1) Render per-XYZ tiles for bbox into a temp cache (strip-based, DEFLATE).
      2) Build mosaics of size (tiles_per_side x tiles_per_side) into out_dir.
      3) Optionally delete the tile cache when done.

    Output mosaics:
      <out_dir>/z<zoom>_x<x0>-<x1>_y<y0>-<y1>.tif
    """
    if dry_run is None:
        dry_run = DRY_RUN_MODE

    if not os.path.isdir(out_dir):
        raise RuntimeError(f"Output directory does not exist: {out_dir}")

    rlayer = find_layer_by_name(layer_name)
    if rlayer is None:
        raise RuntimeError(f"Layer '{layer_name}' not found.")

    if zoom is None:
        zoom = guess_layer_max_zoom(rlayer)

    # Determine tile index range
    x_min, x_max, y_min, y_max = bbox_to_tile_range(top_left, bottom_right, zoom)
    total_tiles = (x_max - x_min + 1) * (y_max - y_min + 1)

    # Cache dir (where per-tile GeoTIFFs go)
    if tile_cache_dir is None:
        tile_cache_dir = os.path.join(out_dir, "_tilecache")
    os.makedirs(tile_cache_dir, exist_ok=True)

    print(f"[run] layer='{layer_name}' zoom={zoom} tiles={total_tiles} dry_run={dry_run}")
    print(f"[cache] {tile_cache_dir}")

    # 1) Render tiles into cache (resume-friendly)
    if not dry_run:
        processed = skipped = batch = 0
        for i, x in enumerate(range(x_min, x_max + 1), 1):
            for y in range(y_min, y_max + 1):
                out_tif = os.path.join(tile_cache_dir, str(zoom), str(x), f"{y}.tif")
                if os.path.exists(out_tif):
                    skipped += 1
                    if verbose: print(f"[tile skip] {zoom}/{x}/{y}")
                else:
                    _render_one_tile(rlayer, zoom, x, y, tile_cache_dir, zlevel=zlevel, verbose=verbose)
                    processed += 1
                    batch += 1
                    if batch >= chunk_size:
                        done = processed + skipped
                        print(f"[pause] tiles done={done}/{total_tiles}; sleeping {pause_secs}s...")
                        time.sleep(pause_secs)
                        batch = 0
        if verbose: print(f"[tiles] new={processed} skipped={skipped} total={total_tiles}")
    else:
        print("[dry_run] Skipping tile rendering.")

    # 2) Build mosaics into out_dir
    print("[mosaic] building mosaics...")
    made = []
    for x0 in range(x_min, x_max + 1, tiles_per_side):
        x1 = min(x0 + tiles_per_side - 1, x_max)
        for y0 in range(y_min, y_max + 1, tiles_per_side):
            y1 = min(y0 + tiles_per_side - 1, y_max)
            tile_paths = _collect_tile_paths(tile_cache_dir, zoom, x0, x1, y0, y1)
            if not tile_paths:
                continue
            mosaic_name = f"z{zoom}_x{x0}-{x1}_y{y0}-{y1}.tif"
            mosaic_path = os.path.join(out_dir, mosaic_name)
            if os.path.exists(mosaic_path):
                print(f"[mosaic skip] exists -> {mosaic_path}")
                made.append(mosaic_path)
                continue
            print(f"[mosaic] {mosaic_name} from {len(tile_paths)} tiles")
            if not dry_run:
                _build_mosaic(tile_paths, mosaic_path,
                              internal_tiled=mosaic_internal_tiled,
                              blocksize=mosaic_blocksize,
                              compress="DEFLATE", zlevel=zlevel, predictor=2, bigtiff="YES")
            made.append(mosaic_path)

    print(f"[mosaic done] wrote {len(made)} mosaics into {out_dir}")

    # 3) Optionally delete cache
    if not dry_run and not keep_tile_cache:
        try:
            import shutil
            print(f"[cleanup] removing tile cache: {tile_cache_dir}")
            shutil.rmtree(tile_cache_dir, ignore_errors=True)
        except Exception as e:
            print(f"[warn] cache cleanup failed: {e}")

    return made
