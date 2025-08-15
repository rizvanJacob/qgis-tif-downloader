"""
qgis_xyz_export.py — One-step: export XYZ tiles to DEFLATE GeoTIFF mosaics (Web Mercator WKT1).

Add XYZ Tile Source as Layer

Load in QGIS Python Console (Windows example):
    exec(open(r"<dir>\qgis_xyz_export.py").read())

Step to build cache:
    build_tile_cache(
        layer_name="<Name of XYZ Tile Layer>",
        top_left=(1, 1), bottom_right=(-1, -1),
        zoom=15,                         # detect max zoom from layer (fallback 19)
        out_dir=r"F:\TIF\mosaic",
        # Cache configuration:
        tile_cache_dir=None,             # default: <out_dir>\_tilecache
        # Codec & performance:
        zlevel=9,                        # DEFLATE level for tiles & mosaics
        chunk_size=50, pause_secs=0.2,   # throttling for tile rendering
        dry_run=False,
        verbose=False
    )    

Single-step run example:
    export_bbox_to_geotiff_tiles(
        layer_name="<Name of XYZ Tile Layer>",
        top_left=(1, 1), bottom_right=(-1, -1),
        zoom=15,                         # detect max zoom from layer (fallback 19)
        out_dir=r"F:\TIF\mosaic",
        # Mosaic controls:
        tiles_per_side=128,              # try 128; bump to 144/160 for larger files
        mosaic_internal_tiled=False,     # recommended for multi-GB files
        mosaic_blocksize=512,            # 256 or 512
        # Cache controls:
        tile_cache_dir=None,             # default: <out_dir>\_tilecache
        keep_tile_cache=False,           # True = keep per-tile GeoTIFFs
        # Perf/codec:
        zlevel=9,                        # DEFLATE level for tiles & mosaics
        chunk_size=50, pause_secs=0.2,   # throttling for tile rendering
        dry_run=False,
        verbose=False
    )

Outputs:
- Per-tile GeoTIFFs (strip-based) staged in: <out_dir>\_tilecache\<z>\<x>\<y>.tif
- Mosaics (internally tiled by default) in out_dir:
      z<z>_x<x0>-<x1>_y<y0>-<y1>.tif

Both tiles & mosaics embed Web Mercator (WKT1/GDAL), avoiding EPSG DB lookups.
"""

import os, math, time
from qgis.core import (
    QgsProject, QgsCoordinateReferenceSystem, QgsRectangle,
    QgsMapSettings, QgsMapRendererCustomPainterJob, QgsApplication, Qgis
)
from qgis.PyQt.QtGui import QImage, QPainter, QColor
from qgis.PyQt.QtCore import QSize
from osgeo import gdal, osr

# -------------------------------------------
# Globals / defaults
# -------------------------------------------
DRY_RUN_MODE = False
DEFAULT_ZOOM_FALLBACK = 19
DEFAULT_TILE_SIZE = 256   # set 512 if your XYZ serves 512px tiles

WEBMERC = QgsCoordinateReferenceSystem("EPSG:3857")
WEBMERC_WKT1_EPSG_3857 = (
    'PROJCS["WGS 84 / Pseudo-Mercator",'
    'GEOGCS["WGS 84",'
        'DATUM["WGS_1984",'
            'SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],'
            'AUTHORITY["EPSG","6326"]],'
        'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
        'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
        'AUTHORITY["EPSG","4326"]],'
    'PROJECTION["Mercator_1SP"],'
    'PARAMETER["central_meridian",0],'
    'PARAMETER["scale_factor",1],'
    'PARAMETER["false_easting",0],'
    'PARAMETER["false_northing",0],'
    'UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
    'AXIS["X",EAST],AXIS["Y",NORTH],'
    'AUTHORITY["EPSG","3857"]]'
)

R = 6378137.0
ORIGIN_SHIFT = math.pi * R

# Qt5/Qt6 compatibility for QImage format enums
try:
    _FMT_ARGB32 = QImage.Format_ARGB32            # Qt5
except AttributeError:
    try:
        _FMT_ARGB32 = QImage.Format.Format_ARGB32 # Qt6
    except AttributeError:
        try:
            _FMT_ARGB32 = QImage.Format_ARGB32_Premultiplied
        except AttributeError:
            _FMT_ARGB32 = getattr(getattr(QImage, "Format", QImage), "Format_ARGB32_Premultiplied", getattr(QImage, "Format_RGB32", None))
        if _FMT_ARGB32 is None:
            raise RuntimeError("Cannot resolve a usable QImage ARGB32 format on this Qt build.")

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

# -------------------------------------------
# SRS helpers (WKT1/GDAL for max compatibility; avoids EPSG lookups)
# -------------------------------------------
def _webmerc_wkt1():
    # Prefer classic WKT1 (GDAL flavor)
    try:
        return WEBMERC.toWkt(Qgis.WktVariant.WKT1_GDAL)
    except Exception:
        return WEBMERC.toWkt()

# -------------------------------------------
# XYZ / WebMercator math & layer utils
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
            try:
                return int(v)
            except Exception:
                pass
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
# Single-tile renderer (strip-based DEFLATE GeoTIFF into cache)
# -------------------------------------------
def _render_one_tile(rlayer, z, x, y, cache_dir, zlevel=9, verbose=True):
    # Compute exact Web Mercator bounds for this XYZ tile
    minx, miny, maxx, maxy = merc_bounds_for_tile(x, y, z)

    # Output paths in the cache
    tile_dir = os.path.join(cache_dir, str(z), str(x))
    os.makedirs(tile_dir, exist_ok=True)
    tmp_png = os.path.join(tile_dir, f"{y}_tmp.png")
    out_tif = os.path.join(tile_dir, f"{y}.tif")

    # Skip if already done
    if os.path.exists(out_tif):
        if verbose: print(f"[tile skip] {z}/{x}/{y}")
        return out_tif

    # --- Render the tile (force opaque RGB) ---
    width = height = DEFAULT_TILE_SIZE

    ms = QgsMapSettings()
    ms.setLayers([rlayer])
    ms.setDestinationCrs(WEBMERC)
    ms.setExtent(QgsRectangle(minx, miny, maxx, maxy))
    ms.setOutputSize(QSize(width, height))

    # Use RGB32 (alpha byte is present but always 0xFF -> opaque)
    try:
        fmt_rgb = QImage.Format_RGB32
    except AttributeError:
        # Qt6 path
        fmt_rgb = getattr(getattr(QImage, "Format", QImage), "Format_RGB32", _FMT_ARGB32)

    img = QImage(width, height, fmt_rgb)
    img.fill(QColor("white"))  # ensure any empty pixels are opaque white

    painter = QPainter(img)
    job = QgsMapRendererCustomPainterJob(ms, painter)
    job.start()
    job.waitForFinished()
    painter.end()

    if verbose: print(f"[render] -> {tmp_png}")
    if not img.save(tmp_png, "PNG"):
        raise RuntimeError(f"Failed to save temp PNG: {tmp_png}")

    # --- PNG -> GeoTIFF (force 3-band RGB, no alpha), embed Web Mercator (WKT1) ---
    translate_opts = gdal.TranslateOptions(
        options = [
            "-a_srs", WEBMERC_WKT1_EPSG_3857,
            "-a_ullr", str(minx), str(maxy), str(maxx), str(miny),
            "-b", "1", "-b", "2", "-b", "3",          # drop any alpha from PNG just in case
            "-co", "PHOTOMETRIC=RGB",                 # make intent explicit
            "-co", "COMPRESS=DEFLATE",
            "-co", "PREDICTOR=2",
            "-co", "BIGTIFF=IF_SAFER",
            "-co", f"ZLEVEL={zlevel}",
        ],
        format="GTiff",
    )

    try:
        if verbose: print(f"[compress] gdal.Translate -> {out_tif}")
        ds_out = gdal.Translate(out_tif, tmp_png, options=translate_opts)
    except Exception as e:
        # Reassert env and retry once (guards against external env churn)
        _ensure_gdal_proj_env()
        if verbose: print(f"[warn] Translate failed ({e}); retrying...")
        ds_out = gdal.Translate(out_tif, tmp_png, options=translate_opts)

    if ds_out is None:
        raise RuntimeError("gdal.Translate returned None.")
    ds_out = None

    # Cleanup temp files
    try:
        os.remove(tmp_png)
        aux = tmp_png + ".aux.xml"
        if os.path.exists(aux):
            os.remove(aux)
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

    # Explicitly stamp SRS on mosaic (WKT1/GDAL)
    trans_opts = gdal.TranslateOptions(
        format="GTiff",
        creationOptions=co,
        outputSRS=WEBMERC_WKT1_EPSG_3857   
    )
    ds_out = gdal.Translate(mosaic_path, vrt_ds, options=trans_opts)
    if ds_out is None:
        raise RuntimeError("gdal.Translate (mosaic) failed.")
    ds_out = None
    vrt_ds = None
    return mosaic_path

def build_tile_cache(
    layer_name,
    top_left, bottom_right,
    zoom=None,
    out_dir=".",
    # Cache configuration:
    tile_cache_dir=None,
    # Codec & performance:
    zlevel=9,
    chunk_size=20,
    pause_secs=1.0,
    dry_run=None,
    verbose=True
):
    """
    Render per-XYZ tiles for bbox into a resume-friendly cache (strip-based, DEFLATE).

    Returns:
      {
        "tile_cache_dir": <resolved cache dir>,
        "zoom": <zoom used>,
        "x_min": int, "x_max": int,
        "y_min": int, "y_max": int
      }
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

    x_min, x_max, y_min, y_max = bbox_to_tile_range(top_left, bottom_right, zoom)
    total_tiles = (x_max - x_min + 1) * (y_max - y_min + 1)

    if tile_cache_dir is None:
        tile_cache_dir = os.path.join(out_dir, "_tilecache")
    os.makedirs(tile_cache_dir, exist_ok=True)

    print(f"[run] layer='{layer_name}' zoom={zoom} tiles={total_tiles} dry_run={dry_run}")
    print(f"[cache] {tile_cache_dir}")

    # Render tiles to cache
    if dry_run:
        print("[dry_run] Skipping tile rendering.")
    else:
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
                    # Refresh env periodically to resist external churn
                    if processed % 1000 == 0:
                        _ensure_gdal_proj_env()
                    if batch >= chunk_size:
                        if verbose:
                            done = processed + skipped
                            print(f"[pause] tiles done={done}/{total_tiles}; sleeping {pause_secs}s...")
                        time.sleep(pause_secs)
                        batch = 0
            print(f"[tiles] row {x}: new={processed} skipped={skipped} total={total_tiles}")

    return {
        "tile_cache_dir": tile_cache_dir,
        "zoom": zoom,
        "x_min": x_min, "x_max": x_max,
        "y_min": y_min, "y_max": y_max
    }


# -------------------------------------------
# One-step export: tiles (cache) → mosaics (out_dir)
# -------------------------------------------
def export_bbox_to_geotiff_tiles(
    layer_name,
    top_left, bottom_right,
    zoom=None,
    out_dir=".",
    # Mosaic configuration:
    tiles_per_side=128,
    mosaic_internal_tiled=True,
    mosaic_blocksize=512,
    # Cache configuration:
    tile_cache_dir=None,
    keep_tile_cache=False,
    # Codec & performance:
    zlevel=9,
    chunk_size=20,
    pause_secs=1.0,
    dry_run=None,
    verbose=True
):
    # 1) Build Tile Cache (or skip if dry_run)
    info = build_tile_cache(
        layer_name=layer_name,
        top_left=top_left, bottom_right=bottom_right,
        zoom=zoom,
        out_dir=out_dir,
        tile_cache_dir=tile_cache_dir,
        zlevel=zlevel,
        chunk_size=chunk_size,
        pause_secs=pause_secs,
        dry_run=dry_run,
        verbose=verbose
    )

    # Unpack what we need for mosaicking
    tile_cache_dir = info["tile_cache_dir"]
    zoom = info["zoom"]
    x_min, x_max = info["x_min"], info["x_max"]
    y_min, y_max = info["y_min"], info["y_max"]

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
                _build_mosaic(
                    tile_paths, mosaic_path,
                    internal_tiled=mosaic_internal_tiled,
                    blocksize=mosaic_blocksize,
                    compress="DEFLATE", zlevel=zlevel, predictor=2, bigtiff="YES"
                )
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

