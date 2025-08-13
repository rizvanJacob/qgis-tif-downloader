"""
qgis_xyz_export.py — Export XYZ tiles to DEFLATE-compressed GeoTIFFs from QGIS.

Load in QGIS Python Console:
    exec(open(r"F:/qgis_scripts/qgis_xyz_export.py").read())   # Windows example
    # or
    exec(open("/path/to/qgis_xyz_export.py").read())           # macOS/Linux

Quick start:
    preview_tiles("Google Satellite",
                  top_left=(103.82, 1.31),
                  bottom_right=(103.86, 1.27),
                  zoom=None)

    resume_report("Google Satellite",
                  top_left=(103.82, 1.31),
                  bottom_right=(103.86, 1.27),
                  out_dir=r"F:\geo_output")

    export_bbox_to_geotiff_tiles(
        layer_name="Google Satellite",
        top_left=(103.82, 1.31),
        bottom_right=(103.86, 1.27),
        zoom=None,  # detect from layer; fallback DEFAULT_ZOOM_FALLBACK
        out_dir=r"F:\geo_output",
        dry_run=False,
        first_n=30, chunk_size=10, pause_secs=1.5, zlevel=9
    )

Notes:
- Outputs are one GeoTIFF per XYZ tile at chosen zoom:
      <out_dir>/<z>/<x>/<y>.tif
- GeoTIFFs are lossless DEFLATE-compressed (PREDICTOR=2, ZLEVEL configurable),
  with NO internal TIFF tiling (strip-based).
- CRS is EPSG:3857 (Web Mercator), matching XYZ services.
- Resume-friendly: existing tiles are skipped.
"""

import os, math, time
from qgis.core import (
    QgsProject, QgsCoordinateReferenceSystem, QgsRectangle,
    QgsMapSettings, QgsMapRendererCustomPainterJob,
    QgsApplication, Qgis
)
from qgis.PyQt.QtGui import QImage, QPainter, QColor
from qgis.PyQt.QtCore import QSize
from osgeo import gdal

# -------------------------------------------
# Globals / defaults
# -------------------------------------------
DRY_RUN_MODE = False            # Default dry-run toggle
DEFAULT_ZOOM_FALLBACK = 19      # Used if layer doesn't expose max zoom
DEFAULT_TILE_SIZE = 256         # 256 for most XYZ; set 512 if your source is 512px

# Web Mercator constants
WEBMERC = QgsCoordinateReferenceSystem("EPSG:3857")
R = 6378137.0
ORIGIN_SHIFT = math.pi * R

def _qgis_osgeo4w_roots():
    """
    Returns (osgeo4w_root, share_dir, proj_dir, gdal_dir) for the running QGIS.
    Works for OSGeo4W layout like: <root>\apps\qgis\  and <root>\share\{proj,gdal}
    """
    prefix = QgsApplication.prefixPath()  # e.g. D:\Non_Vital_Programs\OSGeo4W\apps\qgis
    # go up two levels to OSGeo4W root
    root = os.path.abspath(os.path.join(prefix, r"..", r".."))
    share = os.path.join(root, "share")
    proj  = os.path.join(share, "proj")
    gdaldata = os.path.join(share, "gdal")
    return root, share, proj, gdaldata

def _ensure_gdal_proj_env():
    _, _, proj, gdaldata = _qgis_osgeo4w_roots()
    # Hint both legacy and modern env var names so whichever GDAL/PROJ binding is loaded will find them.
    os.environ["PROJ_LIB"]  = proj
    os.environ["PROJ_DATA"] = proj
    os.environ["GDAL_DATA"] = gdaldata
    # Also tell GDAL via its own config API (works even if env was ignored)
    gdal.SetConfigOption("PROJ_LIB",  proj)
    gdal.SetConfigOption("PROJ_DATA", proj)
    gdal.SetConfigOption("GDAL_DATA", gdaldata)

# Call this once on import:
_ensure_gdal_proj_env()

def _webmerc_wkt():
    # Provide WKT directly to avoid EPSG DB lookup entirely.
    try:
        return WEBMERC.toWkt(Qgis.WktVariant.WKT2_2018)  # QGIS 3.22+
    except Exception:
        return WEBMERC.toWkt()  # fallback

# -------------------------------------------
# Helpers — XYZ / WebMercator math & layer utils
# -------------------------------------------
def clamp_lat(lat):
    "Clamp latitude to WebMercator supported range."
    return max(min(lat, 85.05112878), -85.05112878)

def lonlat_to_tilexy(lon, lat, z):
    "Convert lon/lat at zoom z into fractional XYZ tile indices."
    lat = clamp_lat(lat)
    n = 2 ** z
    xtile = (lon + 180.0) / 360.0 * n
    ytile = (1.0 - math.log(math.tan(math.radians(lat)) + 1/math.cos(math.radians(lat))) / math.pi) / 2.0 * n
    return xtile, ytile

def merc_bounds_for_tile(x, y, z):
    "Get EPSG:3857 (minx, miny, maxx, maxy) bounds for XYZ tile x,y at zoom z."
    n = 2 ** z
    minx = -ORIGIN_SHIFT + (x / n) * (2 * ORIGIN_SHIFT)
    maxx = -ORIGIN_SHIFT + ((x + 1) / n) * (2 * ORIGIN_SHIFT)
    miny = ORIGIN_SHIFT - ((y + 1) / n) * (2 * ORIGIN_SHIFT)
    maxy = ORIGIN_SHIFT - (y / n) * (2 * ORIGIN_SHIFT)
    return (minx, miny, maxx, maxy)

def find_layer_by_name(layer_name):
    "Return the first layer with matching name from the current project."
    for lyr in QgsProject.instance().mapLayers().values():
        if lyr.name() == layer_name:
            return lyr
    return None

def guess_layer_max_zoom(rlayer, fallback=DEFAULT_ZOOM_FALLBACK):
    "Try to read stored max zoom from XYZ layer properties; fallback if not present."
    for k in ("zmax", "maxZoom", "maximumZoomLevel"):
        v = rlayer.customProperty(k, None)
        if v is not None:
            try:
                return int(v)
            except Exception:
                pass
    return fallback

# -------------------------------------------
# Tile enumeration for a bbox
# -------------------------------------------
def tiles_for_bbox(bbox_wgs84, zoom):
    """
    bbox_wgs84 = (min_lon, min_lat, max_lon, max_lat)  [bottom-left, top-right]
    returns list of (z, x, y)
    """
    (min_lon, min_lat, max_lon, max_lat) = bbox_wgs84
    tlx, tly = lonlat_to_tilexy(min_lon, max_lat, zoom)  # top-left
    brx, bry = lonlat_to_tilexy(max_lon, min_lat, zoom)  # bottom-right

    xmin = math.floor(min(tlx, brx))
    xmax = math.floor(max(tlx, brx))
    ymin = math.floor(min(tly, bry))
    ymax = math.floor(max(tly, bry))

    tiles = []
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            tiles.append((zoom, x, y))
    return tiles

def preview_tiles(layer_name, top_left, bottom_right, zoom=None, max_preview=10):
    """
    Print summary of tiles for a bbox at given zoom.
    top_left=(lon,lat), bottom_right=(lon,lat)
    """
    lyr = find_layer_by_name(layer_name)
    if not lyr:
        print(f"[preview] Layer '{layer_name}' not found.")
        return

    if zoom is None:
        zoom = guess_layer_max_zoom(lyr)

    # normalize bbox
    min_lon, min_lat = bottom_right[0], bottom_right[1]
    max_lon, max_lat = top_left[0], top_left[1]
    if min_lon > max_lon: min_lon, max_lon = max_lon, min_lon
    if min_lat > max_lat: min_lat, max_lat = max_lat, min_lat

    tiles = tiles_for_bbox((min_lon, min_lat, max_lon, max_lat), zoom)
    print(f"[preview] Zoom={zoom}, total tiles={len(tiles)}")
    for t in tiles[:max_preview]:
        print("  sample tile:", t)
    if len(tiles) > max_preview:
        print(f"  ... +{len(tiles)-max_preview} more")

# -------------------------------------------
# Render one tile -> compressed GeoTIFF
# -------------------------------------------
def render_one_tile_to_geotiff(rlayer, zxy, out_dir, zlevel=9, keep_tmp=False, verbose=True):
    if not os.path.isdir(out_dir):
        raise RuntimeError(f"Output directory does not exist: {out_dir}")

    z, x, y = zxy
    minx, miny, maxx, maxy = merc_bounds_for_tile(x, y, z)

    # Resume-friendly directories
    tile_dir = os.path.join(out_dir, str(z), str(x))
    os.makedirs(tile_dir, exist_ok=True)
    tmp_png = os.path.join(tile_dir, f"{y}_tmp.png")
    out_tif = os.path.join(tile_dir, f"{y}.tif")

    if os.path.exists(out_tif):
        if verbose: print(f"[skip] exists -> {out_tif}")
        return out_tif

    # Render map
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
    job.start()
    job.waitForFinished()
    painter.end()

    if verbose: print(f"[render] -> {tmp_png}")
    if not img.save(tmp_png, "PNG"):
        raise RuntimeError(f"Failed to save temp PNG: {tmp_png}")

    webmerc_wkt = _webmerc_wkt()
    translate_opts = gdal.TranslateOptions(
        options = [
            "-a_srs", webmerc_wkt,           # <-- use WKT, not 'EPSG:3857'
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
        raise RuntimeError("gdal.Translate failed.")
    ds_out = None

    if not keep_tmp:
        try:
            os.remove(tmp_png)
            aux = tmp_png + ".aux.xml"
            if os.path.exists(aux):
                os.remove(aux)
        except Exception as e:
            if verbose: print(f"[warn] temp cleanup: {e}")

    return out_tif

# -------------------------------------------
# Resume helpers and bbox exporter
# -------------------------------------------
def _tile_out_paths(out_dir, z, x, y):
    tile_dir = os.path.join(out_dir, str(z), str(x))
    tmp_png  = os.path.join(tile_dir, f"{y}_tmp.png")
    out_tif  = os.path.join(tile_dir, f"{y}.tif")
    return tile_dir, tmp_png, out_tif

def resume_report(layer_name, top_left, bottom_right, zoom=None, out_dir="."):
    lyr = find_layer_by_name(layer_name)
    if lyr is None:
        raise RuntimeError(f"Layer '{layer_name}' not found.")
    if zoom is None:
        zoom = guess_layer_max_zoom(lyr)

    min_lon, min_lat = bottom_right[0], bottom_right[1]
    max_lon, max_lat = top_left[0], top_left[1]
    if min_lon > max_lon: min_lon, max_lon = max_lon, min_lon
    if min_lat > max_lat: min_lat, max_lat = max_lat, min_lat

    tiles = tiles_for_bbox((min_lon, min_lat, max_lon, max_lat), zoom)
    total = len(tiles)
    exist = 0
    for (z, x, y) in tiles:
        _, _, out_tif = _tile_out_paths(out_dir, z, x, y)
        if os.path.exists(out_tif): exist += 1
    miss = total - exist
    print(f"[resume] zoom={zoom} total={total} exist={exist} missing={miss}")
    return {"zoom": zoom, "total": total, "exist": exist, "missing": miss}

def export_bbox_to_geotiff_tiles(
    layer_name,
    top_left, bottom_right,
    zoom=None,
    out_dir=".",
    zlevel=9,
    dry_run=None,
    first_n=None,
    chunk_size=20,
    pause_secs=2.0,
    verbose=True
):
    if dry_run is None:
        dry_run = DRY_RUN_MODE

    if not os.path.isdir(out_dir):
        raise RuntimeError(f"Output directory does not exist: {out_dir}")

    rlayer = find_layer_by_name(layer_name)
    if rlayer is None:
        raise RuntimeError(f"Layer '{layer_name}' not found.")

    if zoom is None:
        zoom = guess_layer_max_zoom(rlayer)

    min_lon, min_lat = bottom_right[0], bottom_right[1]
    max_lon, max_lat = top_left[0], top_left[1]
    if min_lon > max_lon: min_lon, max_lon = max_lon, min_lon
    if min_lat > max_lat: min_lat, max_lat = max_lat, min_lat

    tiles = tiles_for_bbox((min_lon, min_lat, max_lon, max_lat), zoom)
    if first_n is not None:
        tiles = tiles[:first_n]

    total = len(tiles)
    print(f"[run] layer='{layer_name}' zoom={zoom} tiles={total} dry_run={dry_run}")

    processed = 0
    skipped   = 0
    batch     = 0

    for i, (z, x, y) in enumerate(tiles, 1):
        tile_dir, tmp_png, out_tif = _tile_out_paths(out_dir, z, x, y)

        if os.path.exists(out_tif):
            skipped += 1
            if verbose:
                print(f"[skip] {z}/{x}/{y} -> {out_tif}")
        else:
            if dry_run:
                if verbose:
                    print(f"[would] {z}/{x}/{y} -> {out_tif}")
            else:
                os.makedirs(tile_dir, exist_ok=True)
                render_one_tile_to_geotiff(rlayer, (z, x, y), out_dir, zlevel=zlevel, verbose=verbose)

            processed += 1
            batch += 1

            if batch >= chunk_size and i < total:
                if verbose:
                    done = processed + skipped
                    print(f"[pause] done={done}/{total}; sleeping {pause_secs}s...")
                time.sleep(pause_secs)
                batch = 0

    print(f"[done] total={total} new={processed} skipped={skipped} dry_run={dry_run}")
