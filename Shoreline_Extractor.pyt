# -*- coding: utf-8 -*-
"""
Shoreline Extraction Toolbox (SE+ v2.0.0) — ArcGIS Pro .pyt

Release 2.0.0 (2026-10-11): first public release. It consolidates the development builds v2.6–v2.8.5
listed below (history kept for reference; see CHANGELOG.md).

v2.8.5
- NEW: the tool warns on the form, and stops the DEM step with a clear message, when the DEM does not overlap the image
  (e.g., Chenier Plain lidar used with a Barataria scene).

v2.8.4
- Every classified output (segmentation, DEM zones, DEM overrides, marsh mask) now stores its class NAMES (CLASS field)
  and COLOURS (colormap) inside the file, so the meaning of each value shows in ArcGIS (legend, Identify, attribute
  table) without consulting the README. Map layers get descriptive names; the marsh mask is added (hidden).

v2.8.3
- FIX: shoreline segments that trace the edge of the image data (NoData margins inside the raster extent) are removed
  when "Remove study-area/extent border segments" is on. These appeared as straight lines in the Gulf.
- FIX: segmentation colours are embedded in the file (colormap); the layer is no longer restyled by the tool.

v2.8.2
- FIX: marsh-shoreline test now uses a majority vote of the LAND pixels within a search radius (default 2.5 px)
  instead of only the nearest pixel. At 10 m the edge pixel of a marsh pond is mixed water/grass with low NDVI, so
  v2.8.1 often missed marsh pond edges and kept them. Beaches a few pixels wide still win the vote and are kept.
- FIX: Segmentation colours now also match when ArcGIS keys the classes by the CLASS field.

v2.8.1 (speed + diagnostics)
- FIX (speed): marsh-shoreline exclusion now tests each segment midpoint directly against the marsh pixels.
  v2.8 converted the whole marsh mask to polygons and buffered it, which could run for hours on large scenes.
- FIX: a flat lidar surface only counts as water when it is not vegetated (NDVI <= marsh threshold); lidar DEMs are
  often interpolated (near-flat) under dense marsh, which was being mistaken for hydro-flattened water.
- NEW: <prefix>_dem_overrides.tif (1 MNDWI water→land by elevation, 2 land→water by low elevation,
  3 land→water by flat surface, 4 marsh fill) and a message splitting land→water by cause.
- Segments carry IS_MARSH (0/1) before being split into kept and excluded outputs.

v2.8 (exclude marsh shorelines)
- NEW (default on): "Exclude marsh shorelines". Marsh is identified as low-lying vegetation: pixels in the DEM-uncertain
  zone (near MHW) with NDVI above the marsh threshold. Shoreline segments that border marsh are removed from the main
  output, so the final shoreline covers only beach / sandy / hard (non-marsh) coast. The removed segments are saved
  separately as <prefix>_marsh_shoreline_segments.shp for reference.
- Marsh is still handled internally (marsh = land option) so that MNDWI false water inside the marsh does not create
  spurious lines elsewhere; the marsh-edge lines are then dropped from the deliverable.
- NEW outputs: <prefix>_marsh_mask.tif, <prefix>_marsh_shoreline_segments(_projected).shp;
  metadata reports kept vs. excluded shoreline length.

v2.7 (marsh = land)
- Goal: the shoreline is the boundary between LAND (including vegetated marsh) and open WATER. MNDWI + Otsu often
  labels wet, vegetated marsh as water, producing false shoreline segments inside the marsh. v2.7 corrects this.
- NEW (default on): "Treat vegetated marsh as land". In the DEM's UNCERTAIN zone (marsh platform near MHW), pixels
  whose NDVI exceeds a threshold are classified as land, overriding MNDWI. Open water (low NDVI) is unaffected, and the
  DEM's water-sure zone still keeps flattened/deep water as water even where floating vegetation is present.
- NEW (optional, default off): "Treat lidar ground at/above MHW as land". In the uncertain zone, pixels where the
  (vegetation-corrected) lidar surface is at or above MHW and is not a flattened water surface become land. Useful for
  sparsely vegetated marsh and mudflat platforms that NDVI misses.
- NEW zone codes in <prefix>_dem_zones.tif: 3 = marsh→land by vegetation, 4 = marsh→land by lidar ground ≥ MHW.
- NEW statistics: pixels changed water→land by each marsh rule (metadata JSON + messages).
- NDVI is now computed whenever Red + NIR are given (no longer only when vegetation bias > 0).

v2.6 (DEM / NOAA lidar integration)
- NEW (optional): lidar DEM constraint. Every pixel is zoned from elevation relative to Mean High Water (MHW):
    * LAND-SURE   : elevation > MHW + land margin      -> forced to land (ridges, levees, cheniers, spoil banks, upland)
    * WATER-SURE  : elevation < MHW - water margin, or a hydro-flattened water surface -> forced to water
    * UNCERTAIN   : everything in between (marsh platform, edges) -> MNDWI/Otsu decides, as before
  Zones are applied after Otsu AND re-applied after morphology, so cleanup cannot undo them.
- NEW (optional): vegetation bias correction. Lidar over dense marsh reads high; where NDVI > threshold the DEM is
  lowered by a user-set bias before zoning (needs Red + NIR bands).
- NEW (optional): DEM-assisted ocean connectivity. Water bodies joined to the Gulf through sub-pixel tidal channels
  (seen in the lidar as low cells not crossed by high ground) are kept, instead of only the single largest polygon.
- NEW: DEM is clipped, projected, aggregated (mean/min/max, preserving narrow channels and levee crests) and snapped
  to the composite grid; unit conversion for metres / US survey feet / international feet.
- FIX: bands are now read on one explicit reference grid (honours Study Area extent/mask + snap raster). In v2.5
  bands were read at full-composite extent but written with the clipped extent's lower-left corner, which could
  shift outputs when a Study Area was used.
- FIX: band NoData is now carried as NaN instead of being treated as reflectance 0.
- All v2.5 behaviour is unchanged when no DEM is supplied (parameters 0–7 keep the same order).

Requires: ArcGIS Pro 3.x, Spatial Analyst; scikit-image (scipy is bundled with ArcGIS Pro)
Developed by Curtis Amo Dwira, MSc in Civil Engineering, Louisiana State University.
With support from Dr. Ahmmed Abdullah, PhD, Assistant Professor, Louisiana State University.
Prof. Bernard Kumi-Boateng, PhD, Professor, University of Mines and Technology, Ghana.
"""

import arcpy
import numpy as np
import os
import json
import datetime
import warnings

from arcpy.sa import *
from skimage.filters import threshold_otsu, gaussian
from skimage.morphology import disk, opening, closing, remove_small_objects, remove_small_holes

try:
    from scipy import ndimage as _ndi
except Exception:  # scipy ships with ArcGIS Pro; skimage fallback just in case
    _ndi = None

warnings.filterwarnings("ignore")

NODATA_F = -9999.0
UNIT_FACTORS = {
    "Meters": 1.0,
    "US Survey Feet": 1200.0 / 3937.0,
    "International Feet": 0.3048,
}


# -----------------------
# Helpers (v2.5)
# -----------------------

def _resolve_band(raster_path, user_key):
    """Resolve a band by exact child name, fuzzy child match, or numeric index."""
    # 1) Exact child name/path
    try:
        cand = arcpy.Raster(f"{raster_path}/{user_key}")
        _ = cand.meanCellWidth
        return cand
    except Exception:
        pass
    # 2) Fuzzy child match
    try:
        desc = arcpy.Describe(raster_path)
        children = getattr(desc, "children", []) or []
        for ch in children:
            if user_key.lower() in ch.name.lower():
                return arcpy.Raster(ch.catalogPath)
    except Exception:
        pass
    # 3) Numeric index
    try:
        idx = int(user_key)
        cand = arcpy.Raster(raster_path + f"/Band_{idx}")
        _ = cand.meanCellWidth
        return cand
    except Exception:
        raise RuntimeError(
            f"Could not locate band '{user_key}' in composite: {raster_path}. "
            f"Provide an exact sublayer name (e.g., SR_B2) or a numeric index."
        )


def _autoscale(arr):
    """If array looks like 0..10000 scaled reflectance, normalize ~0..1."""
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return arr
    p99 = np.nanpercentile(finite, 99)
    return arr / 10000.0 if p99 > 10 else arr


def _define_projection_safely(path, spatial_ref, messages):
    try:
        arcpy.DefineProjection_management(path, spatial_ref)
    except Exception as e:
        messages.addWarningMessage(f"DefineProjection failed for {os.path.basename(path)}: {e}")


def _clip_percentiles(a, lo=0.5, hi=99.5):
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return a, None, None
    low, high = np.nanpercentile(finite, [lo, hi])
    b = np.clip(a, low, high)
    return b, float(low), float(high)


# -----------------------
# Helpers (v2.6) — grid I/O
# -----------------------

def _make_grid(ref_raster):
    """Reference grid (lower-left, size, cell) taken from a raster evaluated under the current env."""
    ext = ref_raster.extent
    return {
        "ll": arcpy.Point(ext.XMin, ext.YMin),
        "ncols": int(ref_raster.width),
        "nrows": int(ref_raster.height),
        "cs": float(ref_raster.meanCellWidth),
        "extent": ext,
    }


def _read_on_grid(ras, grid):
    """Read a FLOAT raster (path or Raster) onto the reference grid. NoData -> NaN."""
    arr = arcpy.RasterToNumPyArray(ras, grid["ll"], grid["ncols"], grid["nrows"], nodata_to_value=np.nan)
    return arr.astype(np.float64)


def _write_float(arr, grid, path, sr, messages):
    out = np.where(np.isfinite(arr), arr, NODATA_F).astype(np.float32)
    r = arcpy.NumPyArrayToRaster(out, grid["ll"], grid["cs"], grid["cs"], value_to_nodata=NODATA_F)
    r.save(path)
    _define_projection_safely(path, sr, messages)


SEG_CLASSES = {
    1: ("Open water (shoreline-forming)", [0, 92, 230]),
    2: ("Inland / unconnected water", [115, 223, 255]),
    3: ("Land", [205, 170, 102]),
    4: ("Marsh (excluded)", [76, 168, 0]),
}

ZONE_CLASSES = {
    0: ("Near MHW, not classed as marsh - MNDWI decides", [190, 190, 190]),
    1: ("High ground - land (lidar)", [230, 200, 60]),
    2: ("Low or flat water surface - water (lidar)", [0, 112, 255]),
    3: ("Marsh - vegetated, near MHW", [76, 168, 0]),
    4: ("Bare ground at/above MHW - land (lidar)", [168, 112, 0]),
}

OVERRIDE_CLASSES = {
    0: ("No change from MNDWI", [235, 235, 235]),
    1: ("MNDWI water -> land (high ground)", [255, 170, 0]),
    2: ("MNDWI land -> water (low elevation)", [0, 112, 255]),
    3: ("MNDWI land -> water (flat surface)", [132, 0, 168]),
    4: ("MNDWI water -> land (marsh fill)", [76, 168, 0]),
}

MARSH_CLASSES = {
    0: ("Not marsh", [235, 235, 235]),
    1: ("Marsh", [76, 168, 0]),
}


def _label_classified(path, classes, messages=None):
    """
    Store colours (colormap) and class names (CLASS field in the attribute table) inside a classified raster,
    so the class meanings show in ArcGIS (legend, Identify, attribute table) without consulting the README.
    """
    try:
        clr = os.path.splitext(path)[0] + ".clr"
        with open(clr, "w") as f:
            for k, (_, rgb) in classes.items():
                f.write(f"{k} {rgb[0]} {rgb[1]} {rgb[2]}\n")
        arcpy.management.AddColormap(path, "#", clr)
    except Exception as e:
        if messages:
            messages.addWarningMessage(f"Could not embed colours in {os.path.basename(path)}: {e}")
    try:
        arcpy.management.BuildRasterAttributeTable(path, "Overwrite")
        if "CLASS" not in [f.name.upper() for f in arcpy.ListFields(path)]:
            arcpy.management.AddField(path, "CLASS", "TEXT", field_length=60)
        with arcpy.da.UpdateCursor(path, ["Value", "CLASS"]) as u:
            for r in u:
                r[1] = classes.get(int(r[0]), ("", None))[0]
                u.updateRow(r)
    except Exception as e:
        if messages:
            messages.addWarningMessage(f"Could not add class names to {os.path.basename(path)}: {e}")


def _segmentation(valid, water_all, ocean, marsh=None):
    """Display classes: 1 open water, 2 inland water, 3 land, 4 marsh, 255 NoData."""
    seg = np.full(valid.shape, 255, dtype=np.uint8)
    seg[valid] = 3
    if marsh is not None:
        seg[valid & marsh & ~water_all] = 4
    seg[valid & water_all & ~ocean] = 2
    seg[valid & ocean] = 1
    return seg


def _style_unique_raster(lyr, classes):
    """Unique-value colours + labels for a classified raster layer (best effort)."""
    sym = lyr.symbology
    if not hasattr(sym, "colorizer"):
        return
    sym.updateColorizer("RasterUniqueValueColorizer")
    try:
        sym.colorizer.field = "Value"
    except Exception:
        pass
    by_name = {name.lower(): (name, rgb) for name, rgb in classes.values()}
    for grp in sym.colorizer.groups:
        for itm in grp.items:
            key = str(itm.values[0]).strip() if itm.values else ""
            match = None
            try:
                match = classes.get(int(float(key)))          # item keyed by numeric Value
            except Exception:
                match = by_name.get(key.lower()) or by_name.get(str(itm.label).strip().lower())  # keyed by CLASS
            if match:
                itm.label = match[0]
                itm.color = {"RGB": list(match[1]) + [100]}
    lyr.symbology = sym


def _style_line(lyr, rgb, width):
    sym = lyr.symbology
    if not hasattr(sym, "renderer"):
        return
    sym.updateRenderer("SimpleRenderer")
    sym.renderer.symbol.color = {"RGB": rgb + [100]}
    sym.renderer.symbol.width = width
    lyr.symbology = sym


def _near_marsh(xs, ys, marsh, grid, tol_px):
    """
    Vectorised test: is any marsh pixel within tol_px * cell size of each point (x, y)?
    Distance is measured from the point to the nearest edge of each pixel square.
    """
    n = len(xs)
    out = np.zeros(n, dtype=bool)
    if n == 0 or marsh is None:
        return out
    cs = grid["cs"]
    llx, lly = grid["ll"].X, grid["ll"].Y
    nrows, ncols = marsh.shape
    ymax = lly + nrows * cs
    r_tol = tol_px * cs
    col0 = np.floor((xs - llx) / cs).astype(np.int64)
    row0 = np.floor((ymax - ys) / cs).astype(np.int64)
    k = int(np.ceil(tol_px)) + 1
    for di in range(-k, k + 1):
        for dj in range(-k, k + 1):
            r = row0 + di
            c = col0 + dj
            inside = (r >= 0) & (r < nrows) & (c >= 0) & (c < ncols)
            if not inside.any():
                continue
            x0 = llx + c * cs
            y1 = ymax - r * cs
            dx = np.maximum(np.maximum(x0 - xs, 0.0), xs - (x0 + cs))
            dy = np.maximum(np.maximum((y1 - cs) - ys, 0.0), ys - y1)
            close = inside & (np.hypot(dx, dy) <= r_tol + 1e-9)
            idx = np.where(close)[0]
            if idx.size:
                out[idx] |= marsh[r[idx], c[idx]]
    return out


def _near_data_edge(xs, ys, valid, grid, tol_px=1.0):
    """True where a point lies within tol_px of a NoData pixel or of the raster's outer edge."""
    near_nodata = _near_marsh(xs, ys, ~valid, grid, tol_px)
    cs = grid["cs"]
    nrows, ncols = valid.shape
    x0, y0 = grid["ll"].X, grid["ll"].Y
    x1, y1 = x0 + ncols * cs, y0 + nrows * cs
    d = tol_px * cs + 1e-9
    near_frame = (xs - x0 <= d) | (x1 - xs <= d) | (ys - y0 <= d) | (y1 - ys <= d)
    return near_nodata | near_frame


def _marsh_majority(xs, ys, marsh, land, grid, radius_px, min_frac=0.5):
    """
    v2.8.2: a segment borders marsh when, among the LAND pixels within radius_px of its midpoint, at least min_frac
    are marsh. Looking a few pixels inland avoids the mixed water/vegetation edge pixel (low NDVI at 10 m) that made
    the nearest-pixel test miss marsh pond edges, while beaches several pixels wide still win the vote.
    """
    n = len(xs)
    out = np.zeros(n, dtype=bool)
    if n == 0 or marsh is None:
        return out
    cs = grid["cs"]
    llx, lly = grid["ll"].X, grid["ll"].Y
    nrows, ncols = marsh.shape
    ymax = lly + nrows * cs
    r_tol = radius_px * cs
    col0 = np.floor((xs - llx) / cs).astype(np.int64)
    row0 = np.floor((ymax - ys) / cs).astype(np.int64)
    n_land = np.zeros(n, dtype=np.int32)
    n_marsh = np.zeros(n, dtype=np.int32)
    k = int(np.ceil(radius_px)) + 1
    for di in range(-k, k + 1):
        for dj in range(-k, k + 1):
            r = row0 + di
            c = col0 + dj
            inside = (r >= 0) & (r < nrows) & (c >= 0) & (c < ncols)
            if not inside.any():
                continue
            x0 = llx + c * cs
            y1 = ymax - r * cs
            dx = np.maximum(np.maximum(x0 - xs, 0.0), xs - (x0 + cs))
            dy = np.maximum(np.maximum((y1 - cs) - ys, 0.0), ys - y1)
            idx = np.where(inside & (np.hypot(dx, dy) <= r_tol + 1e-9))[0]
            if idx.size:
                is_land = land[r[idx], c[idx]]
                n_land[idx] += is_land
                n_marsh[idx] += is_land & marsh[r[idx], c[idx]]
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (n_land > 0) & (n_marsh >= min_frac * n_land)
    return out


def _sum_length(fc):
    total = 0.0
    with arcpy.da.SearchCursor(fc, ["SHAPE@LENGTH"]) as cur:
        for (l,) in cur:
            total += l or 0.0
    return total


def _write_uint8(arr, grid, path, sr, messages):
    r = arcpy.NumPyArrayToRaster(arr.astype(np.uint8), grid["ll"], grid["cs"], grid["cs"], value_to_nodata=255)
    r.save(path)
    _define_projection_safely(path, sr, messages)


# -----------------------
# Helpers (v2.6) — DEM preparation (arcpy)
# -----------------------

def _prepare_dem(dem_in, units, grid, ref_sr, snap_path, work_dir, prefix, messages):
    """
    Clip -> project (native res) -> convert to metres -> flat-surface mask (hydro-flattening) ->
    aggregate to composite cell size (MEAN / MINIMUM / MAXIMUM) -> snap/resample onto the reference grid.
    Returns dict of NumPy arrays on the reference grid (metres, DEM vertical datum).
    """
    os.makedirs(work_dir, exist_ok=True)
    f = UNIT_FACTORS.get(units, 1.0)

    d = arcpy.Describe(dem_in)
    dem_path = d.catalogPath
    dem_sr = d.spatialReference
    dem_cs = float(arcpy.Raster(dem_path).meanCellWidth)

    both_projected = (getattr(dem_sr, "type", "") == "Projected" and getattr(ref_sr, "type", "") == "Projected")
    if both_projected:
        dem_cs_t = dem_cs * dem_sr.metersPerUnit / ref_sr.metersPerUnit   # DEM cell in composite map units
    else:
        dem_cs_t = dem_cs
    messages.addMessage(f"DEM: {dem_path} | SR: {dem_sr.name} | cell {dem_cs:g} ({dem_cs_t:g} in composite units) | z-units: {units}")

    neutral = dict(snapRaster=None, cellSize="MAXOF", mask=None, extent="MAXOF")

    # 1) Clip to the processing extent (+3 cells) in the DEM's own SR
    ext = grid["extent"]
    buf = 3.0 * grid["cs"]
    ring = arcpy.Array([
        arcpy.Point(ext.XMin - buf, ext.YMin - buf), arcpy.Point(ext.XMin - buf, ext.YMax + buf),
        arcpy.Point(ext.XMax + buf, ext.YMax + buf), arcpy.Point(ext.XMax + buf, ext.YMin - buf),
        arcpy.Point(ext.XMin - buf, ext.YMin - buf)])
    e = arcpy.Polygon(ring, ref_sr).projectAs(dem_sr).extent
    de = d.extent
    if e.XMax < de.XMin or e.XMin > de.XMax or e.YMax < de.YMin or e.YMin > de.YMax:
        raise RuntimeError(
            "The DEM does not overlap the image. Check that the lidar covers this study site "
            f"(DEM: {os.path.basename(dem_path)}).")
    clipped = os.path.join(work_dir, f"{prefix}_dem_clip.tif")
    with arcpy.EnvManager(**neutral):
        arcpy.management.Clip(dem_path, f"{e.XMin} {e.YMin} {e.XMax} {e.YMax}", clipped,
                              "#", "#", "NONE", "NO_MAINTAIN_EXTENT")

        # 2) Project at native resolution if needed
        if dem_sr.name != ref_sr.name:
            native = os.path.join(work_dir, f"{prefix}_dem_native.tif")
            arcpy.management.ProjectRaster(clipped, native, ref_sr, "BILINEAR", str(dem_cs_t))
            messages.addMessage("DEM projected to composite coordinate system (native resolution).")
        else:
            native = clipped

        # 3) Metres + flat-surface mask at native resolution (hydro-flattened water = constant elevation)
        dem_m = Float(Raster(native)) * f
        rng = FocalStatistics(dem_m, NbrRectangle(3, 3, "CELL"), "RANGE", "DATA")
        flat = Con(rng <= 0.002, 1.0, 0.0)

        # 4) Aggregate to the composite cell size
        factor = int(round(grid["cs"] / dem_cs_t)) if dem_cs_t > 0 else 1
        if both_projected and factor >= 2:
            r_mean = Aggregate(dem_m, factor, "MEAN", "EXPAND", "DATA")
            r_min = Aggregate(dem_m, factor, "MINIMUM", "EXPAND", "DATA")
            r_max = Aggregate(dem_m, factor, "MAXIMUM", "EXPAND", "DATA")
            r_flat = Aggregate(flat, factor, "MEAN", "EXPAND", "DATA")
            messages.addMessage(f"DEM aggregated x{factor} (mean / min / max / flat-fraction).")
        else:
            r_mean = r_min = r_max = dem_m
            r_flat = flat
            factor = 1

        tmp = {}
        for key, ras in (("mean", r_mean), ("min", r_min), ("max", r_max), ("flat", r_flat)):
            p = os.path.join(work_dir, f"{prefix}_dem_{key}_agg.tif")
            Float(ras).save(p)
            tmp[key] = p

    # 5) Snap/resample onto the composite grid and read as arrays
    out = {}
    with arcpy.EnvManager(snapRaster=snap_path, cellSize=grid["cs"], mask=None, extent=grid["extent"]):
        for key, p in tmp.items():
            rp = os.path.join(work_dir, f"{prefix}_dem_{key}_grid.tif")
            method = "NEAREST" if key in ("min", "max") else "BILINEAR"
            arcpy.management.Resample(p, rp, f"{grid['cs']} {grid['cs']}", method)
            out[key] = _read_on_grid(rp, grid)

    out["factor"] = factor
    out["native_cell"] = dem_cs_t
    return out


# -----------------------
# Helpers (v2.6) — pure NumPy logic (testable outside ArcGIS)
# -----------------------

def _dem_zones(dem_mean, flat_frac, mhw, land_margin, water_margin, use_flat=True,
               ndvi=None, ndvi_thr=0.3, veg_bias=0.0, flat_ndvi=None, flat_ndvi_max=None):
    """
    Zone pixels from elevation (metres, same vertical datum as mhw).
    Returns land_sure, water_sure (bool), zone (uint8: 0 uncertain, 1 land, 2 water, 255 no DEM), veg mask, z.
    v2.8.1: if flat_ndvi/flat_ndvi_max are given, a flat (hydro-flattened-looking) surface only counts as water when it
    is NOT vegetated (NDVI <= flat_ndvi_max). Lidar DEMs are often interpolated (smooth, near-flat) under dense marsh,
    which otherwise looks like flattened water.
    """
    z = np.array(dem_mean, dtype=np.float64, copy=True)
    veg = np.zeros(z.shape, dtype=bool)
    if ndvi is not None and veg_bias:
        veg = np.isfinite(ndvi) & (ndvi > ndvi_thr) & np.isfinite(z)
        z[veg] -= float(veg_bias)

    valid = np.isfinite(z)
    with np.errstate(invalid="ignore"):
        land = valid & (z > mhw + land_margin)
        water = valid & (z < mhw - water_margin)
        if use_flat and flat_frac is not None:
            flat_water = valid & np.isfinite(flat_frac) & (flat_frac >= 0.5) & (dem_mean <= mhw + land_margin)
            if flat_ndvi is not None and flat_ndvi_max is not None:
                flat_water &= ~(np.isfinite(flat_ndvi) & (flat_ndvi > flat_ndvi_max))
            water |= flat_water
    water &= ~land

    zone = np.full(z.shape, 255, dtype=np.uint8)
    zone[valid] = 0
    zone[land] = 1
    zone[water] = 2
    return land, water, zone, veg, z


def _marsh_to_land(zone, ndvi=None, marsh_ndvi_thr=0.3, use_veg_rule=True,
                   dem_z=None, flat_frac=None, mhw=None, use_ground_rule=False):
    """
    v2.7 marsh = land rules, applied ONLY in the DEM-uncertain zone (zone == 0).
      veg_land    : NDVI > marsh_ndvi_thr                       (vegetated marsh)
      ground_land : lidar surface (veg-corrected) >= MHW and not a flattened water surface (bare/sparse marsh platform)
    Returns (veg_land, ground_land) boolean masks.
    """
    uncertain = (zone == 0)
    veg_land = np.zeros(zone.shape, dtype=bool)
    ground_land = np.zeros(zone.shape, dtype=bool)
    with np.errstate(invalid="ignore"):
        if use_veg_rule and ndvi is not None:
            veg_land = uncertain & np.isfinite(ndvi) & (ndvi > marsh_ndvi_thr)
        if use_ground_rule and dem_z is not None and mhw is not None:
            not_flat = np.ones(zone.shape, dtype=bool)
            if flat_frac is not None:
                not_flat = ~(np.isfinite(flat_frac) & (flat_frac >= 0.5))
            ground_land = uncertain & np.isfinite(dem_z) & (dem_z >= mhw) & not_flat & ~veg_land
    return veg_land, ground_land


def _apply_zones(binary, land, water):
    out = np.array(binary, dtype=bool, copy=True)
    out[land] = False
    out[water] = True
    return out


def _label(mask):
    if _ndi is not None:
        lab, n = _ndi.label(mask, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
        return lab, int(n)
    from skimage.measure import label as sk_label
    lab = sk_label(mask, connectivity=1)
    return lab, int(lab.max())


def _ocean_connected(binary, valid, dem_min=None, dem_max=None, mhw=None, land_margin=None, water_margin=None):
    """
    Keep water connected to the 'ocean' seed (largest water body touching the processing-area edge; largest overall
    if none touch). Low DEM cells (min < MHW - water margin) not crossed by high ground (max < MHW + land margin)
    act as sub-pixel channel connectors — they link water bodies but are NOT added to the output water.
    4-connectivity is used so water cannot leak diagonally through a levee.
    Returns (ocean_mask, n_water_bodies_kept).
    """
    water = np.asarray(binary, dtype=bool) & valid
    lab_w, n_w = _label(water)
    if n_w == 0:
        return water, 0

    sizes = np.bincount(lab_w.ravel())
    sizes[0] = 0
    # edge of the valid region (array border or next to NoData/outside study area)
    if _ndi is not None:
        interior = _ndi.binary_erosion(valid, structure=np.ones((3, 3), bool), border_value=0)
    else:
        from skimage.morphology import binary_erosion
        padded = np.pad(valid, 1, constant_values=False)
        interior = binary_erosion(padded, np.ones((3, 3), bool))[1:-1, 1:-1]
    edge = valid & ~interior
    touching = np.unique(lab_w[edge & water])
    touching = touching[touching > 0]
    if touching.size:
        seed = touching[np.argmax(sizes[touching])]
    else:
        seed = int(np.argmax(sizes))

    connector = water.copy()
    if dem_min is not None and dem_max is not None and mhw is not None:
        with np.errstate(invalid="ignore"):
            channel = (np.isfinite(dem_min) & np.isfinite(dem_max)
                       & (dem_min < mhw - water_margin) & (dem_max < mhw + land_margin))
        connector |= channel & valid

    lab_c, _ = _label(connector)
    ocean_label = lab_c[lab_w == seed][0]
    ocean = water & (lab_c == ocean_label)
    kept = int(np.unique(lab_w[ocean]).size)
    return ocean, kept


# -----------------------
# Toolbox
# -----------------------

class Toolbox(object):
    def __init__(self):
        self.label = "Shoreline Extraction Toolbox (Segments)"
        self.alias = "shoreline_toolbox_segments"
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        self.label = "Shoreline Extractor (SE+ v2.0)"
        self.description = (
            "Semi-automated shoreline extraction using MNDWI + Otsu with adaptive smoothing, morphology, "
            "optional lidar-DEM land/water constraint (NOAA lidar) with vegetated marsh classified as land, "
            "dissolved-water boundary → shoreline segments, "
            "optional border-segment removal, and projection."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        def P(**kw):
            return arcpy.Parameter(**kw)

        params = [
            # ---- v2.5 parameters (order unchanged: 0–7) ----
            P(displayName="Composite Raster", name="composite_raster", datatype="GPRasterLayer",
              parameterType="Required", direction="Input"),
            P(displayName="Green Band (name or index, e.g., SR_B2 or 2)", name="green_band_name",
              datatype="GPString", parameterType="Required", direction="Input"),
            P(displayName="SWIR Band (name or index, e.g., SR_B5 or 5)", name="swir_band_name",
              datatype="GPString", parameterType="Required", direction="Input"),
            P(displayName="Study Area (optional)", name="study_area", datatype="GPFeatureLayer",
              parameterType="Optional", direction="Input"),
            P(displayName="Output Folder", name="output_folder", datatype="DEFolder",
              parameterType="Required", direction="Input"),
            P(displayName="Output File Prefix", name="output_prefix", datatype="GPString",
              parameterType="Required", direction="Input"),
            P(displayName="Output Coordinate System", name="output_crs", datatype="GPCoordinateSystem",
              parameterType="Required", direction="Input"),
            P(displayName="Remove study-area/extent border segments", name="drop_border",
              datatype="GPBoolean", parameterType="Optional", direction="Input"),

            # ---- v2.6 DEM parameters (8–18) ----
            P(displayName="Lidar DEM (bare-earth, e.g., NOAA Digital Coast)", name="dem_raster",
              datatype="GPRasterLayer", parameterType="Optional", direction="Input"),
            P(displayName="DEM vertical units", name="dem_units", datatype="GPString",
              parameterType="Optional", direction="Input"),
            P(displayName="MHW elevation in the DEM's vertical datum (metres, e.g., NAVD88)", name="mhw_elev",
              datatype="GPDouble", parameterType="Optional", direction="Input"),
            P(displayName="Land margin above MHW (m) — higher = land for certain", name="land_margin",
              datatype="GPDouble", parameterType="Optional", direction="Input"),
            P(displayName="Water margin below MHW (m) — lower = water for certain", name="water_margin",
              datatype="GPDouble", parameterType="Optional", direction="Input"),
            P(displayName="Treat hydro-flattened (flat) surfaces as water", name="use_hydroflat",
              datatype="GPBoolean", parameterType="Optional", direction="Input"),
            P(displayName="Keep water connected to the Gulf through DEM channels", name="dem_connectivity",
              datatype="GPBoolean", parameterType="Optional", direction="Input"),
            P(displayName="Red Band (for marsh detection & vegetation correction, e.g., B4)", name="red_band_name",
              datatype="GPString", parameterType="Optional", direction="Input"),
            P(displayName="NIR Band (for marsh detection & vegetation correction, e.g., B8)", name="nir_band_name",
              datatype="GPString", parameterType="Optional", direction="Input"),
            P(displayName="Lidar vegetation bias to remove (m)", name="veg_bias",
              datatype="GPDouble", parameterType="Optional", direction="Input"),
            P(displayName="NDVI threshold for vegetation bias correction", name="ndvi_thr",
              datatype="GPDouble", parameterType="Optional", direction="Input"),

            # ---- v2.7 marsh = land parameters (19–21) ----
            P(displayName="Fill MNDWI false water inside marsh (needs Red + NIR)", name="marsh_as_land",
              datatype="GPBoolean", parameterType="Optional", direction="Input"),
            P(displayName="NDVI threshold defining marsh vegetation", name="marsh_ndvi_thr",
              datatype="GPDouble", parameterType="Optional", direction="Input"),
            P(displayName="Treat lidar ground at/above MHW as land (sparse marsh / mudflat platforms)",
              name="ground_above_mhw_land", datatype="GPBoolean", parameterType="Optional", direction="Input"),

            # ---- v2.8 marsh exclusion (22–23) ----
            P(displayName="Exclude marsh shorelines (keep beach / non-marsh coast only)", name="exclude_marsh",
              datatype="GPBoolean", parameterType="Optional", direction="Input"),
            P(displayName="Marsh search radius (pixels) — majority of nearby land must be marsh", name="marsh_tol_px",
              datatype="GPDouble", parameterType="Optional", direction="Input"),
        ]

        params[7].value = True  # Default: drop border segments

        for p in params[8:]:
            p.category = "DEM constraint (NOAA lidar)"
        params[9].filter.type = "ValueList"
        params[9].filter.list = list(UNIT_FACTORS.keys())
        params[9].value = "Meters"
        params[11].value = 0.5
        params[12].value = 0.5
        params[13].value = True
        params[14].value = True
        params[17].value = 0.15
        params[18].value = 0.3
        params[19].value = True
        params[20].value = 0.3
        params[21].value = False
        params[22].value = True
        params[23].value = 2.5
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        has_dem = bool(parameters[8].value)
        for p in parameters[9:]:
            p.enabled = has_dem
        return

    def updateMessages(self, parameters):
        if parameters[8].value:
            if parameters[10].value is None:
                parameters[10].setErrorMessage(
                    "Enter MHW elevation in the DEM's vertical datum (metres). "
                    "Get it from the nearest NOAA tide station datums page (MHW minus NAVD88) or NOAA VDatum. "
                    "Enter 0 if the DEM is already referenced to MHW.")
            red, nir = parameters[15].value, parameters[16].value
            if bool(red) != bool(nir):
                (parameters[16] if red else parameters[15]).setWarningMessage(
                    "Both Red and NIR bands are needed for marsh detection and vegetation correction.")
            elif not red and (parameters[19].value or parameters[22].value):
                (parameters[22] if parameters[22].value else parameters[19]).setWarningMessage(
                    "Marsh detection needs Red and NIR bands (e.g., B4 and B8 for Sentinel-2); "
                    "without them marsh shorelines cannot be identified or excluded.")
            if parameters[23].value is not None and parameters[23].value <= 0:
                parameters[23].setErrorMessage("Tolerance must be > 0 pixels.")
            try:   # warn early if the DEM and the image do not overlap
                if parameters[0].value:
                    ie = arcpy.Describe(parameters[0].value).extent
                    dd = arcpy.Describe(parameters[8].value)
                    de = dd.extent.projectAs(ie.spatialReference) if ie.spatialReference else dd.extent
                    if (ie.XMax < de.XMin or ie.XMin > de.XMax or ie.YMax < de.YMin or ie.YMin > de.YMax):
                        parameters[8].setWarningMessage(
                            "This DEM does not overlap the Composite Raster. Use lidar covering this study site; "
                            "otherwise the DEM step is skipped (MNDWI only, no marsh detection).")
            except Exception:
                pass
            for i in (11, 12):
                if parameters[i].value is not None and parameters[i].value < 0:
                    parameters[i].setErrorMessage("Margin must be ≥ 0.")
        return

    def execute(self, parameters, messages):
        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True

        try:
            # Inputs
            composite          = parameters[0].value
            green_band_name    = parameters[1].valueAsText
            swir_band_name     = parameters[2].valueAsText
            study_area         = parameters[3].value
            output_folder      = parameters[4].valueAsText
            prefix             = parameters[5].valueAsText
            output_crs_param   = parameters[6]
            drop_border        = bool(parameters[7].value) if parameters[7].altered else True

            dem_in             = parameters[8].value
            use_dem            = bool(dem_in)
            dem_units          = parameters[9].valueAsText or "Meters"
            mhw                = float(parameters[10].value) if parameters[10].value is not None else None
            land_margin        = float(parameters[11].value) if parameters[11].value is not None else 0.5
            water_margin       = float(parameters[12].value) if parameters[12].value is not None else 0.5
            use_hydroflat      = bool(parameters[13].value) if parameters[13].value is not None else True
            dem_connectivity   = bool(parameters[14].value) if parameters[14].value is not None else True
            red_band_name      = parameters[15].valueAsText
            nir_band_name      = parameters[16].valueAsText
            veg_bias           = float(parameters[17].value) if parameters[17].value is not None else 0.15
            ndvi_thr           = float(parameters[18].value) if parameters[18].value is not None else 0.3
            has_ndvi_bands     = use_dem and bool(red_band_name) and bool(nir_band_name)
            use_veg            = has_ndvi_bands and veg_bias > 0
            marsh_as_land      = bool(parameters[19].value) if parameters[19].value is not None else True
            marsh_ndvi_thr     = float(parameters[20].value) if parameters[20].value is not None else 0.3
            ground_land_rule   = bool(parameters[21].value) if parameters[21].value is not None else False
            use_marsh_veg      = marsh_as_land and has_ndvi_bands
            exclude_marsh      = bool(parameters[22].value) if parameters[22].value is not None else True
            marsh_tol_px       = float(parameters[23].value) if parameters[23].value is not None else 2.5
            do_marsh_excl      = exclude_marsh and has_ndvi_bands
            marsh_excl_mask    = None

            if use_dem and mhw is None:
                raise RuntimeError("A DEM was supplied but MHW elevation is empty.")

            # Describe input & lock envs using dataset path/extent (NOT parameter object)
            try:
                comp_desc   = arcpy.Describe(composite)
                comp_path   = comp_desc.catalogPath
                comp_extent = comp_desc.extent
                comp_sr     = comp_desc.spatialReference
                messages.addMessage(f"Raster: {comp_path} | SR: {comp_sr.name if comp_sr else 'Unknown'}")
            except Exception as e:
                messages.addErrorMessage("Could not describe composite raster: " + str(e))
                raise

            arcpy.env.snapRaster = comp_path
            arcpy.env.cellSize   = comp_path
            arcpy.env.extent     = comp_extent
            arcpy.env.mask       = None

            raster_path   = comp_path
            raster_extent = comp_extent
            raster_sr     = comp_sr

            # Paths
            binary_raster_path   = os.path.join(output_folder, f"{prefix}_binary.tif")
            mndwi_raster_path    = os.path.join(output_folder, f"{prefix}_mndwi.tif")
            mndwi_debug_path     = os.path.join(output_folder, f"{prefix}_mndwi_debug.tif")
            dem_zones_path       = os.path.join(output_folder, f"{prefix}_dem_zones.tif")
            dem_grid_path        = os.path.join(output_folder, f"{prefix}_dem_on_grid.tif")
            dem_work_dir         = os.path.join(output_folder, f"{prefix}_dem_work")
            marsh_mask_path      = os.path.join(output_folder, f"{prefix}_marsh_mask.tif")
            segmentation_path    = os.path.join(output_folder, f"{prefix}_segmentation.tif")
            overrides_path       = os.path.join(output_folder, f"{prefix}_dem_overrides.tif")
            marsh_segments       = os.path.join(output_folder, f"{prefix}_marsh_shoreline_segments.shp")
            marsh_segments_proj  = os.path.join(output_folder, f"{prefix}_marsh_shoreline_segments_projected.shp")
            polygon_path         = os.path.join(output_folder, f"{prefix}_water_polygon.shp")
            valid_polygons       = os.path.join(output_folder, f"{prefix}_valid_polygons.shp")
            dissolved_water      = "in_memory/water_diss"
            shoreline_segments   = os.path.join(output_folder, f"{prefix}_shoreline_segments.shp")
            shoreline_proj       = os.path.join(output_folder, f"{prefix}_shoreline_segments_projected.shp")
            projected_study_area = os.path.join(output_folder, f"{prefix}_projected_study_area.shp")
            mask_boundary        = os.path.join(output_folder, f"{prefix}_mask_boundary.shp")
            extent_fc            = os.path.join(output_folder, f"{prefix}_raster_extent.shp")
            metadata_path        = os.path.join(output_folder, f"{prefix}_metadata.json")

            # Optional clipping (via env; the reference grid below honours extent + mask + snap)
            if study_area:
                try:
                    arcpy.Project_management(study_area, projected_study_area, raster_sr)
                    arcpy.env.mask   = projected_study_area
                    arcpy.env.extent = arcpy.Describe(projected_study_area).extent
                    messages.addMessage("Study area reprojected; env.mask/extent set.")
                except Exception as e:
                    messages.addWarningMessage("Study area project failed, using full raster: " + str(e))
                    arcpy.env.mask   = None
                    arcpy.env.extent = raster_extent
                    study_area = None
            else:
                messages.addMessage("No study area provided. Using full raster.")

            # Read bands on ONE reference grid (v2.6 alignment fix)
            green_ras = _resolve_band(raster_path, green_band_name)
            swir_ras  = _resolve_band(raster_path, swir_band_name)
            messages.addMessage(f"Opening bands: {green_band_name}, {swir_band_name}")

            green_f = Float(green_ras)          # evaluated under env extent/mask/snap
            grid = _make_grid(green_f)
            messages.addMessage(f"Reference grid: {grid['ncols']} x {grid['nrows']} cells, "
                                f"LL=({grid['ll'].X:.3f}, {grid['ll'].Y:.3f})")

            green = _autoscale(_read_on_grid(green_f, grid))
            swir  = _autoscale(_read_on_grid(Float(swir_ras), grid))

            # Cell size & auto-params
            cs = grid["cs"] if grid["cs"] > 0 else 1.0
            messages.addMessage(f"Detected raster cell size (map units): {cs}")

            if cs <= 0.5: gaussian_sigma = 2.0
            elif cs <= 2: gaussian_sigma = 1.5
            elif cs <=10: gaussian_sigma = 1.0
            elif cs <=30: gaussian_sigma = 0.7
            else:         gaussian_sigma = 0.5
            morph_radius = max(1, int(round(gaussian_sigma)))
            messages.addMessage(f"Auto params: sigma={gaussian_sigma}, morph_radius={morph_radius}")

            # MNDWI + smoothing
            eps = 1e-10
            with np.errstate(invalid="ignore", divide="ignore"):
                mndwi = (green - swir) / (green + swir + eps)
            try:
                if gaussian_sigma and gaussian_sigma > 0:
                    nan_mask = ~np.isfinite(mndwi)
                    work = mndwi.copy(); work[nan_mask] = 0.0
                    smoothed = gaussian(work, sigma=gaussian_sigma, preserve_range=True)
                    smoothed[nan_mask] = np.nan
                    mndwi = smoothed
                    messages.addMessage(f"Applied Gaussian smoothing (sigma={gaussian_sigma}).")
            except Exception as e:
                messages.addWarningMessage("Gaussian smoothing failed; proceeding unsmoothed: " + str(e))

            # Otsu threshold (robust)
            finite_mask = np.isfinite(mndwi)
            finite_vals = mndwi[finite_mask]
            if finite_vals.size == 0:
                try:
                    _write_float(mndwi, grid, mndwi_debug_path, raster_sr, messages)
                except Exception:
                    pass
                messages.addErrorMessage("MNDWI array has no valid data to threshold.")
                raise RuntimeError("Empty MNDWI data for thresholding.")

            try:
                clipped_vals, lo_p, hi_p = _clip_percentiles(finite_vals, 0.5, 99.5)
                thresh = float(threshold_otsu(clipped_vals))
                messages.addMessage(f"Otsu threshold (robust): {thresh:.6f} (clip {lo_p:.4f}–{hi_p:.4f})")
            except Exception:
                thresh = None
                messages.addWarningMessage("Otsu thresholding failed; fallback to MNDWI > 0.0")

            with np.errstate(invalid="ignore"):
                binary_bool = (mndwi > thresh) if thresh is not None else (mndwi > 0.0)
            binary_bool &= finite_mask

            # ---------------- DEM constraint (v2.6) ----------------
            dem_meta = {"used": False}
            marsh_meta = {"excluded": False}
            land_sure = water_sure = None
            dem = None
            if use_dem:
                try:
                    dem = _prepare_dem(dem_in, dem_units, grid, raster_sr, comp_path,
                                       dem_work_dir, prefix, messages)

                    ndvi = None
                    if has_ndvi_bands:
                        red = _autoscale(_read_on_grid(Float(_resolve_band(raster_path, red_band_name)), grid))
                        nir = _autoscale(_read_on_grid(Float(_resolve_band(raster_path, nir_band_name)), grid))
                        with np.errstate(invalid="ignore", divide="ignore"):
                            ndvi = (nir - red) / (nir + red + eps)
                        if use_veg:
                            messages.addMessage(f"Vegetation correction: -{veg_bias} m where NDVI > {ndvi_thr}.")
                    elif marsh_as_land:
                        messages.addWarningMessage("Marsh = land skipped: Red and NIR bands not provided.")

                    land_sure, water_sure, zone, veg, dem_z = _dem_zones(
                        dem["mean"], dem["flat"], mhw, land_margin, water_margin,
                        use_flat=use_hydroflat, ndvi=ndvi if use_veg else None, ndvi_thr=ndvi_thr,
                        veg_bias=veg_bias if use_veg else 0.0,
                        flat_ndvi=ndvi, flat_ndvi_max=marsh_ndvi_thr if ndvi is not None else None)
                    land_sure &= finite_mask
                    water_sure &= finite_mask
                    with np.errstate(invalid="ignore"):
                        water_by_elev = water_sure & np.isfinite(dem_z) & (dem_z < mhw - water_margin)
                    water_by_flat = water_sure & ~water_by_elev

                    # v2.8 marsh definition for exclusion: low-lying (uncertain zone) vegetation
                    if do_marsh_excl and ndvi is not None:
                        marsh_excl_mask, _ = _marsh_to_land(zone, ndvi=ndvi, marsh_ndvi_thr=marsh_ndvi_thr,
                                                            use_veg_rule=True)
                        marsh_excl_mask &= finite_mask
                        _write_uint8(np.where(finite_mask, marsh_excl_mask, 255).astype(np.uint8),
                                     grid, marsh_mask_path, raster_sr, messages)
                        _label_classified(marsh_mask_path, MARSH_CLASSES, messages)
                        messages.addMessage(
                            f"Marsh mask: {int(marsh_excl_mask.sum())} px "
                            f"({marsh_excl_mask.sum() * grid['cs'] ** 2 / 1e6:.3f} km²) of low-lying vegetation.")

                    # v2.7 marsh = land (uncertain zone only)
                    marsh_veg, marsh_ground = _marsh_to_land(
                        zone, ndvi=ndvi, marsh_ndvi_thr=marsh_ndvi_thr, use_veg_rule=use_marsh_veg,
                        dem_z=dem_z, flat_frac=dem["flat"] if use_hydroflat else None, mhw=mhw,
                        use_ground_rule=ground_land_rule)
                    marsh_veg &= finite_mask
                    marsh_ground &= finite_mask
                    zone[marsh_veg] = 3
                    zone[marsh_ground] = 4
                    land_force = land_sure | marsh_veg | marsh_ground

                    before = binary_bool.copy()
                    binary_bool = _apply_zones(binary_bool, land_sure, water_sure)
                    to_land = int(np.sum(before & ~binary_bool))
                    to_water = int(np.sum(~before & binary_bool & finite_mask))
                    to_water_elev = int(np.sum(~before & water_by_elev))
                    to_water_flat = int(np.sum(~before & water_by_flat))

                    # Diagnostic: where the DEM changed the MNDWI result
                    ovr = np.zeros(before.shape, dtype=np.uint8)
                    ovr[before & land_sure] = 1                 # MNDWI water -> land (elevation high)
                    ovr[~before & water_by_elev] = 2            # MNDWI land -> water (elevation low)
                    ovr[~before & water_by_flat] = 3            # MNDWI land -> water (flat surface)
                    ovr[before & ~water_sure & (marsh_veg | marsh_ground)] = 4   # MNDWI water -> marsh fill
                    ovr[~finite_mask] = 255
                    try:
                        _write_uint8(ovr, grid, overrides_path, raster_sr, messages)
                        _label_classified(overrides_path, OVERRIDE_CLASSES, messages)
                    except Exception as e:
                        messages.addWarningMessage("Could not save DEM overrides raster: " + str(e))
                    marsh_veg_px = int(np.sum(binary_bool & marsh_veg))      # spectral water -> land (veg rule)
                    marsh_ground_px = int(np.sum(binary_bool & marsh_ground))  # spectral water -> land (ground rule)
                    binary_bool = _apply_zones(binary_bool, land_force, water_sure)
                    n_valid = max(1, int(np.sum(finite_mask)))
                    covered = int(np.sum((zone != 255) & finite_mask))
                    px_area = grid["cs"] * grid["cs"]

                    messages.addMessage(
                        f"DEM zoning (MHW={mhw} m, land >{mhw + land_margin:.2f}, water <{mhw - water_margin:.2f}): "
                        f"land-sure {100.0 * land_sure.sum() / n_valid:.1f}%, "
                        f"water-sure {100.0 * water_sure.sum() / n_valid:.1f}%, "
                        f"DEM coverage {100.0 * covered / n_valid:.1f}%.")
                    messages.addMessage(
                        f"DEM overrides: water→land {to_land} px; land→water {to_water} px "
                        f"(low elevation {to_water_elev}, flat surface {to_water_flat}). "
                        f"Map: {os.path.basename(overrides_path)} (1 →land, 2 →water elev, 3 →water flat, 4 marsh fill).")
                    if use_marsh_veg or ground_land_rule:
                        messages.addMessage(
                            f"Marsh = land: {marsh_veg_px} px ({marsh_veg_px * px_area / 1e6:.3f} km²) of MNDWI 'water' "
                            f"reclassified as vegetated marsh (NDVI > {marsh_ndvi_thr}); "
                            f"{marsh_ground_px} px by lidar ground ≥ MHW.")
                    if covered < 0.5 * n_valid:
                        messages.addWarningMessage("DEM covers less than half of the processing area; "
                                                   "uncovered pixels use MNDWI only.")

                    land_sure = land_force   # re-applied after morphology together with marsh = land
                    zone_out = np.where(finite_mask, zone, 255).astype(np.uint8)
                    _write_uint8(zone_out, grid, dem_zones_path, raster_sr, messages)
                    _label_classified(dem_zones_path, ZONE_CLASSES, messages)
                    _write_float(dem["mean"], grid, dem_grid_path, raster_sr, messages)

                    dem_meta = {
                        "used": True,
                        "dem": str(arcpy.Describe(dem_in).catalogPath),
                        "vertical_units": dem_units,
                        "mhw_m": mhw,
                        "land_threshold_m": mhw + land_margin,
                        "water_threshold_m": mhw - water_margin,
                        "hydroflat_as_water": use_hydroflat,
                        "connectivity": dem_connectivity,
                        "veg_correction_m": veg_bias if use_veg else 0.0,
                        "ndvi_threshold": ndvi_thr if use_veg else None,
                        "aggregation_factor": dem["factor"],
                        "coverage_pct": round(100.0 * covered / n_valid, 2),
                        "land_sure_pct": round(100.0 * land_sure.sum() / n_valid, 2),
                        "water_sure_pct": round(100.0 * water_sure.sum() / n_valid, 2),
                        "override_water_to_land_px": to_land,
                        "override_land_to_water_px": to_water,
                        "marsh_as_land": use_marsh_veg,
                        "marsh_ndvi_threshold": marsh_ndvi_thr if use_marsh_veg else None,
                        "ground_above_mhw_as_land": ground_land_rule,
                        "marsh_veg_water_to_land_px": marsh_veg_px,
                        "marsh_ground_water_to_land_px": marsh_ground_px,
                        "marsh_water_to_land_km2": round((marsh_veg_px + marsh_ground_px) * px_area / 1e6, 4),
                        "zone_codes": {"0": "uncertain (MNDWI decides)", "1": "land-sure (elevation)",
                                       "2": "water-sure (elevation/flattened)", "3": "marsh→land (vegetation)",
                                       "4": "marsh→land (lidar ground ≥ MHW)", "255": "no DEM / outside"},
                        "zones_raster": dem_zones_path,
                    }
                except Exception as e:
                    messages.addWarningMessage("DEM step failed; continuing with MNDWI only: " + str(e))
                    use_dem = False
                    land_sure = water_sure = None
                    dem_meta = {"used": False, "error": str(e)}

            water_count = int(np.sum(np.logical_and(binary_bool, finite_mask)))
            messages.addMessage(f"Water pixel count: {water_count}")
            if water_count == 0:
                try:
                    _write_float(mndwi, grid, mndwi_debug_path, raster_sr, messages)
                    messages.addMessage(f"Saved MNDWI debug raster: {mndwi_debug_path}")
                except Exception:
                    messages.addWarningMessage("Failed to save MNDWI debug raster.")
                raise RuntimeError("No water pixels detected; check band names, clipping, or scene quality.")

            # Morphology
            try:
                if morph_radius and morph_radius > 0:
                    selem = disk(morph_radius)
                    binary_bool = opening(binary_bool, selem)
                    small_pixels = max(8, morph_radius * morph_radius)
                    binary_bool = remove_small_objects(binary_bool.astype(bool), min_size=small_pixels)
                    binary_bool = closing(binary_bool, selem)

                    # Fill holes using PIXELS only (auto)
                    if cs <= 10:
                        hole_px = 200
                    elif cs <= 30:
                        hole_px = 120
                    else:
                        hole_px = 80

                    binary_bool = remove_small_holes(binary_bool.astype(bool), area_threshold=int(hole_px))
                    messages.addMessage(f"Filled holes <= {hole_px}px.")
                    messages.addMessage(f"Morphology: open/remove_small/close (r={morph_radius}px).")
            except Exception as e:
                messages.addWarningMessage("Morphological cleanup failed/partial: " + str(e))

            binary_bool &= finite_mask
            if use_dem and land_sure is not None:
                binary_bool = _apply_zones(binary_bool, land_sure, water_sure)
                messages.addMessage("DEM zones (incl. marsh = land) re-applied after morphology.")

            water_all = binary_bool.copy()   # v2.8: all water before ocean selection (for the segmentation map)

            # DEM-assisted ocean connectivity (raster domain)
            ocean_done = False
            if use_dem and dem_connectivity and dem is not None:
                try:
                    ocean, kept = _ocean_connected(binary_bool, finite_mask, dem["min"], dem["max"],
                                                   mhw, land_margin, water_margin)
                    if ocean.any():
                        binary_bool = ocean
                        ocean_done = True
                        messages.addMessage(f"Ocean-connected water kept via DEM channels: {kept} water bodies.")
                except Exception as e:
                    messages.addWarningMessage("DEM connectivity failed; using largest-polygon rule: " + str(e))

            # Save rasters (on the reference grid)
            try:
                _write_float(mndwi, grid, mndwi_raster_path, raster_sr, messages)
                messages.addMessage(f"MNDWI saved: {mndwi_raster_path}")
            except Exception as e:
                messages.addWarningMessage("Failed to save MNDWI: " + str(e))

            try:
                binary_uint8 = np.full(mndwi.shape, 255, dtype=np.uint8)
                binary_uint8[finite_mask] = binary_bool[finite_mask].astype(np.uint8)
                _write_uint8(binary_uint8, grid, binary_raster_path, raster_sr, messages)
                messages.addMessage(f"Binary saved: {binary_raster_path}")
            except Exception as e:
                messages.addErrorMessage("Failed to save binary raster: " + str(e))
                raise

            # v2.8 segmentation map: 1 open water, 2 inland water, 3 land, 4 marsh
            try:
                if ocean_done:
                    ocean_disp = binary_bool
                else:  # same rule as the polygon step: largest edge-touching water body
                    ocean_disp, _ = _ocean_connected(binary_bool, finite_mask)
                seg = _segmentation(finite_mask, water_all, ocean_disp, marsh_excl_mask)
                _write_uint8(seg, grid, segmentation_path, raster_sr, messages)
                _label_classified(segmentation_path, SEG_CLASSES, messages)
                n_seg = max(1, int(finite_mask.sum()))
                messages.addMessage(
                    "Segmentation: " + ", ".join(
                        f"{SEG_CLASSES[k][0]} {100.0 * np.sum(seg == k) / n_seg:.1f}%" for k in (1, 2, 3, 4)))
            except Exception as e:
                messages.addWarningMessage("Failed to save segmentation map: " + str(e))

            # Raster → Polygon (water)
            try:
                arcpy.RasterToPolygon_conversion(binary_raster_path, polygon_path, "SIMPLIFY", "Value")
                messages.addMessage("RasterToPolygon completed.")
            except Exception as e:
                messages.addErrorMessage("RasterToPolygon failed: " + str(e))
                raise
            try:
                arcpy.RepairGeometry_management(polygon_path)
            except Exception:
                pass

            # Detect value field
            value_field = None
            try:
                fields = list(arcpy.ListFields(polygon_path))
                for cand in ("GRIDCODE", "VALUE", "GRID_CODE", "DN"):
                    for f in fields:
                        if f.name.upper() == cand:
                            value_field = f.name; break
                    if value_field: break
                if value_field is None:
                    for f in fields:
                        if f.type in ("SmallInteger", "Integer", "Double", "Single"):
                            value_field = f.name; break
                messages.addMessage(f"Polygon value field detected: {value_field}")
            except Exception:
                value_field = None

            # Keep water polygons (value == 1)
            try:
                if not value_field:
                    arcpy.AddWarning("Could not detect value field; copying all polygons.")
                    arcpy.CopyFeatures_management(polygon_path, valid_polygons)
                else:
                    arcpy.MakeFeatureLayer_management(polygon_path, "poly_lyr")
                    arcpy.SelectLayerByAttribute_management("poly_lyr", "NEW_SELECTION", f'"{value_field}" = 1')
                    arcpy.CopyFeatures_management("poly_lyr", valid_polygons)
                    arcpy.Delete_management("poly_lyr")
                messages.addMessage("Water polygons isolated.")
            except Exception as e:
                messages.addErrorMessage("Filtering water polygons failed: " + str(e))
                raise

            # Mask boundary (for optional border filtering)
            try:
                if study_area and arcpy.Exists(projected_study_area):
                    if arcpy.Exists(mask_boundary): arcpy.Delete_management(mask_boundary)
                    arcpy.PolygonToLine_management(projected_study_area, mask_boundary)
                    messages.addMessage("Mask boundary from study area created.")
                else:
                    ext = grid["extent"]
                    arr = arcpy.Array([
                        arcpy.Point(ext.XMin, ext.YMin),
                        arcpy.Point(ext.XMin, ext.YMax),
                        arcpy.Point(ext.XMax, ext.YMax),
                        arcpy.Point(ext.XMax, ext.YMin),
                        arcpy.Point(ext.XMin, ext.YMin),
                    ])
                    extent_poly = arcpy.Polygon(arr, raster_sr)
                    if arcpy.Exists(extent_fc): arcpy.Delete_management(extent_fc)
                    arcpy.CopyFeatures_management([extent_poly], extent_fc)
                    if arcpy.Exists(mask_boundary): arcpy.Delete_management(mask_boundary)
                    arcpy.PolygonToLine_management(extent_fc, mask_boundary)
                    messages.addMessage("Mask boundary from raster extent created.")
            except Exception as e:
                messages.addWarningMessage("Failed to create mask boundary: " + str(e))

            # -------------------------------
            # Keep ONLY the ocean-connected water (to avoid inland waterbodies)
            # (skipped when the DEM connectivity step already did this in the raster domain)
            # -------------------------------
            if ocean_done:
                messages.addMessage("Ocean selection done by DEM connectivity; keeping all resulting polygons.")
            else:
                ocean_polygons = "in_memory/ocean_polygons"
                try:
                    if arcpy.Exists(ocean_polygons):
                        arcpy.Delete_management(ocean_polygons)

                    arcpy.MakeFeatureLayer_management(valid_polygons, "wpoly_lyr")

                    if arcpy.Exists(mask_boundary):
                        arcpy.SelectLayerByLocation_management(
                            "wpoly_lyr", "INTERSECT", mask_boundary, selection_type="NEW_SELECTION"
                        )

                    sel_count = int(arcpy.GetCount_management("wpoly_lyr").getOutput(0))
                    if sel_count == 0:
                        messages.addWarningMessage("No boundary-touching water found; falling back to largest water polygon overall.")
                        arcpy.SelectLayerByAttribute_management("wpoly_lyr", "CLEAR_SELECTION")

                    oid_field = arcpy.Describe(valid_polygons).OIDFieldName
                    best_oid, best_area = None, -1.0
                    with arcpy.da.SearchCursor("wpoly_lyr", [oid_field, "SHAPE@AREA"]) as cur:
                        for oid, area in cur:
                            if area is not None and area > best_area:
                                best_oid, best_area = oid, area

                    if best_oid is None:
                        raise RuntimeError("No water polygons available after filtering.")

                    arcpy.SelectLayerByAttribute_management("wpoly_lyr", "NEW_SELECTION", f'"{oid_field}" = {best_oid}')
                    arcpy.CopyFeatures_management("wpoly_lyr", ocean_polygons)
                    arcpy.Delete_management("wpoly_lyr")

                    valid_polygons = ocean_polygons
                    messages.addMessage(f"Ocean water selected (OID={best_oid}, area={best_area:.2f}).")

                except Exception as e:
                    messages.addWarningMessage("Ocean-only filtering failed; continuing with all water polygons: " + str(e))

            # Dissolve water → boundary-only shoreline segments
            try:
                if arcpy.Exists(dissolved_water): arcpy.Delete_management(dissolved_water)
                arcpy.Dissolve_management(valid_polygons, dissolved_water)
                shoreline_lines_tmp = "in_memory/shoreline_lines_tmp"
                if arcpy.Exists(shoreline_lines_tmp): arcpy.Delete_management(shoreline_lines_tmp)
                arcpy.PolygonToLine_management(dissolved_water, shoreline_lines_tmp)

                single_parts = "in_memory/shoreline_singleparts"
                split_segments = "in_memory/shoreline_split_segments"
                for tmp in (single_parts, split_segments):
                    if arcpy.Exists(tmp): arcpy.Delete_management(tmp)

                arcpy.MultipartToSinglepart_management(shoreline_lines_tmp, single_parts)
                arcpy.SplitLine_management(single_parts, split_segments)
                try: arcpy.RepairGeometry_management(split_segments)
                except Exception: pass

                if drop_border and arcpy.Exists(mask_boundary):
                    arcpy.MakeFeatureLayer_management(split_segments, "seg_lyr")
                    arcpy.SelectLayerByLocation_management("seg_lyr", "SHARE_A_LINE_SEGMENT_WITH", mask_boundary, selection_type="NEW_SELECTION")
                    arcpy.SelectLayerByAttribute_management("seg_lyr", "SWITCH_SELECTION")
                    if arcpy.Exists(shoreline_segments): arcpy.Delete_management(shoreline_segments)
                    arcpy.CopyFeatures_management("seg_lyr", shoreline_segments)
                    arcpy.Delete_management("seg_lyr")
                    messages.addMessage("Removed study-area/extent border segments.")
                else:
                    if arcpy.Exists(shoreline_segments): arcpy.Delete_management(shoreline_segments)
                    arcpy.CopyFeatures_management(split_segments, shoreline_segments)

                # v2.8.3: remove segments that trace the edge of the image DATA (NoData margins inside the raster
                # extent), which the extent-rectangle border filter above cannot catch
                if drop_border:
                    try:
                        bx, by, boid = [], [], []
                        with arcpy.da.SearchCursor(shoreline_segments, ["OID@", "SHAPE@"]) as cur:
                            for oid, shp in cur:
                                if shp is None:
                                    continue
                                p = shp.positionAlongLine(0.5, True).firstPoint
                                boid.append(oid); bx.append(p.X); by.append(p.Y)
                        edge_flags = _near_data_edge(np.array(bx), np.array(by), finite_mask, grid, 1.0)
                        edge_oids = {o for o, f in zip(boid, edge_flags) if f}
                        if edge_oids:
                            with arcpy.da.UpdateCursor(shoreline_segments, ["OID@"]) as u:
                                for r in u:
                                    if r[0] in edge_oids:
                                        u.deleteRow()
                        messages.addMessage(f"Removed {len(edge_oids)} segments along the image-data edge (NoData).")
                    except Exception as e:
                        messages.addWarningMessage("Image-edge segment filter failed: " + str(e))

                # v2.8: remove shoreline segments that border marsh
                marsh_meta = {"excluded": False}
                if do_marsh_excl and marsh_excl_mask is not None and marsh_excl_mask.any():
                    try:
                        # v2.8.1 FAST: test each segment's midpoint directly against the marsh pixels
                        # (no raster→polygon / buffer of the marsh mask, which is very slow for large marsh areas)
                        oids, mx, my = [], [], []
                        with arcpy.da.SearchCursor(shoreline_segments, ["OID@", "SHAPE@"]) as cur:
                            for oid, shp in cur:
                                if shp is None:
                                    continue
                                p = shp.positionAlongLine(0.5, True).firstPoint
                                oids.append(oid); mx.append(p.X); my.append(p.Y)
                        land_px = finite_mask & ~water_all
                        flags = _marsh_majority(np.array(mx), np.array(my), marsh_excl_mask, land_px, grid,
                                                marsh_tol_px)
                        marsh_oids = {o for o, f in zip(oids, flags) if f}
                        messages.addMessage(f"Marsh check: {len(marsh_oids)} of {len(oids)} segments border marsh.")

                        total_len = _sum_length(shoreline_segments)
                        oid_name = arcpy.Describe(shoreline_segments).OIDFieldName
                        names = [f.name.upper() for f in arcpy.ListFields(shoreline_segments)]
                        if "IS_MARSH" not in names:
                            arcpy.AddField_management(shoreline_segments, "IS_MARSH", "SHORT")
                        with arcpy.da.UpdateCursor(shoreline_segments, ["OID@", "IS_MARSH"]) as u:
                            for r in u:
                                r[1] = 1 if r[0] in marsh_oids else 0
                                u.updateRow(r)
                        kept_tmp = os.path.join(output_folder, f"{prefix}_nonmarsh_tmp.shp")
                        for tmp in (marsh_segments, kept_tmp):
                            if arcpy.Exists(tmp): arcpy.Delete_management(tmp)
                        arcpy.analysis.Select(shoreline_segments, marsh_segments, '"IS_MARSH" = 1')
                        arcpy.analysis.Select(shoreline_segments, kept_tmp, '"IS_MARSH" = 0')
                        arcpy.Delete_management(shoreline_segments)
                        arcpy.CopyFeatures_management(kept_tmp, shoreline_segments)
                        arcpy.Delete_management(kept_tmp)

                        marsh_len = _sum_length(marsh_segments)
                        kept_len = _sum_length(shoreline_segments)
                        marsh_meta = {
                            "excluded": True,
                            "tolerance_px": marsh_tol_px,
                            "total_length_map_units": round(total_len, 2),
                            "kept_nonmarsh_length_map_units": round(kept_len, 2),
                            "excluded_marsh_length_map_units": round(marsh_len, 2),
                            "excluded_pct": round(100.0 * marsh_len / total_len, 2) if total_len else None,
                            "marsh_segments": marsh_segments,
                        }
                        messages.addMessage(
                            f"Marsh shorelines excluded: {marsh_len:.0f} of {total_len:.0f} map units "
                            f"({marsh_meta['excluded_pct']}%). Kept non-marsh shoreline: {kept_len:.0f}. "
                            f"Excluded segments saved to {os.path.basename(marsh_segments)}.")
                    except Exception as e:
                        messages.addWarningMessage("Marsh exclusion failed; keeping all shoreline segments: " + str(e))
                        marsh_meta = {"excluded": False, "error": str(e)}
                elif exclude_marsh:
                    messages.addWarningMessage("Marsh exclusion skipped (needs DEM + Red + NIR bands, "
                                               "or no marsh was detected).")

                for tmp in (shoreline_lines_tmp, single_parts, split_segments):
                    try:
                        if arcpy.Exists(tmp): arcpy.Delete_management(tmp)
                    except Exception:
                        pass

                seg_count = int(arcpy.GetCount_management(shoreline_segments).getOutput(0)) if arcpy.Exists(shoreline_segments) else 0
                messages.addMessage(f"Final shoreline segments (boundary-only): {seg_count}")
                if seg_count == 0:
                    messages.addWarningMessage("No shoreline segments saved. Inspect intermediate outputs.")
            except Exception as e_seg:
                messages.addErrorMessage("Boundary extraction failed: " + str(e_seg))
                messages.addErrorMessage("GP messages:\n" + arcpy.GetMessages(2))
                raise

            # Define + Project to Output CRS
            proj_sr_name = None
            try:
                if arcpy.Exists(shoreline_segments):
                    arcpy.DefineProjection_management(shoreline_segments, raster_sr)
                    arcpy.Project_management(shoreline_segments, shoreline_proj, output_crs_param.value)
                    messages.addMessage(f"Projected shoreline segments: {shoreline_proj}")
                    if shoreline_proj and arcpy.Exists(shoreline_proj):
                        proj_sr_name = arcpy.Describe(shoreline_proj).spatialReference.name
                else:
                    messages.addWarningMessage("Shoreline segments not found; skipping projection.")
                if marsh_meta.get("excluded") and arcpy.Exists(marsh_segments):
                    arcpy.DefineProjection_management(marsh_segments, raster_sr)
                    arcpy.Project_management(marsh_segments, marsh_segments_proj, output_crs_param.value)
            except Exception as e:
                messages.addWarningMessage("Projecting shoreline segments failed: " + str(e))
                messages.addWarningMessage("GP messages:\n" + arcpy.GetMessages(2))
                shoreline_proj = shoreline_segments if arcpy.Exists(shoreline_segments) else None
                try:
                    if shoreline_proj and arcpy.Exists(shoreline_proj):
                        proj_sr_name = arcpy.Describe(shoreline_proj).spatialReference.name
                except Exception:
                    proj_sr_name = None

            # SEG_ID on both outputs
            try:
                for fc in (shoreline_segments, shoreline_proj):
                    if fc and arcpy.Exists(fc):
                        names = [f.name.upper() for f in arcpy.ListFields(fc)]
                        if "SEG_ID" not in names:
                            arcpy.AddField_management(fc, "SEG_ID", "LONG")
                            with arcpy.da.UpdateCursor(fc, ["SEG_ID"]) as u:
                                i = 1
                                for r in u:
                                    r[0] = i; u.updateRow(r); i += 1
                messages.addMessage("SEG_ID added to shoreline segments (raw & projected).")
            except Exception as e:
                messages.addWarningMessage("Failed to add SEG_ID to outputs: " + str(e))

            # Metadata JSON
            try:
                out_crs_label = proj_sr_name
                if not out_crs_label:
                    try:
                        out_crs_label = output_crs_param.valueAsText
                    except Exception:
                        out_crs_label = None

                meta = {
                    "tool_version": "SE+ v2.0.0",
                    "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z",
                    "composite": str(raster_path),
                    "green_band": green_band_name,
                    "swir_band": swir_band_name,
                    "cellsize_map_units": cs,
                    "grid": {"ncols": grid["ncols"], "nrows": grid["nrows"],
                             "lower_left": [grid["ll"].X, grid["ll"].Y]},
                    "gaussian_sigma": gaussian_sigma,
                    "morph_radius_px": morph_radius,
                    "otsu_threshold": None if thresh is None else float(thresh),
                    "water_pixel_count": water_count,
                    "output_crs": out_crs_label,
                    "drop_border_segments": drop_border,
                    "dem_constraint": dem_meta,
                    "marsh_exclusion": marsh_meta,
                    "ocean_selection": "dem_connectivity" if ocean_done else "largest_boundary_polygon",
                    "paths": {
                        "mndwi": mndwi_raster_path,
                        "binary": binary_raster_path,
                        "water_polygons": valid_polygons,
                        "shoreline_segments": shoreline_segments,
                        "shoreline_segments_projected": shoreline_proj
                    }
                }
                with open(metadata_path, "w") as f:
                    json.dump(meta, f, indent=2)
                messages.addMessage(f"Metadata written: {metadata_path}")
            except Exception as e:
                messages.addWarningMessage("Failed to write metadata JSON: " + str(e))

            # Add outputs to current map
            # Order: diagnostics (hidden) → segmentation map → excluded marsh lines → final shoreline on top
            try:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap

                def _add(path, name, visible=True):
                    lyr = m.addDataFromPath(path)
                    try:
                        lyr.name = name
                        lyr.visible = visible
                    except Exception:
                        pass
                    return lyr

                if arcpy.Exists(mndwi_raster_path):
                    _add(mndwi_raster_path, f"{prefix} MNDWI", visible=False)
                if use_dem and arcpy.Exists(dem_zones_path):
                    _add(dem_zones_path, f"{prefix} DEM zones (lidar decision per pixel)", visible=False)
                if use_dem and arcpy.Exists(overrides_path):
                    _add(overrides_path, f"{prefix} DEM overrides (where lidar changed MNDWI)", visible=False)
                if use_dem and marsh_excl_mask is not None and arcpy.Exists(marsh_mask_path):
                    _add(marsh_mask_path, f"{prefix} Marsh mask", visible=False)

                if arcpy.Exists(segmentation_path):
                    # Drawn with the colormap stored in the file (1 dark blue open water, 2 light blue inland water,
                    # 3 tan land, 4 green marsh). The layer's symbology is deliberately NOT modified here.
                    _add(segmentation_path, f"{prefix} Segmentation")

                if marsh_meta.get("excluded") and arcpy.Exists(marsh_segments_proj):
                    lyr = _add(marsh_segments_proj, f"{prefix} Marsh shoreline (excluded)")
                    try:
                        _style_line(lyr, [56, 168, 0], 1.0)
                    except Exception:
                        pass

                final_fc = shoreline_proj if (shoreline_proj and arcpy.Exists(shoreline_proj)) else (
                    shoreline_segments if arcpy.Exists(shoreline_segments) else None)
                if final_fc:
                    title = "Extracted shoreline (non-marsh)" if marsh_meta.get("excluded") else "Extracted shoreline"
                    lyr = _add(final_fc, f"{prefix} {title}")
                    try:
                        _style_line(lyr, [255, 0, 0], 2.5)
                    except Exception:
                        pass
                messages.addMessage("Added to map: segmentation (open water / inland water / land / marsh) "
                                    "and the extracted shoreline (red, on top). MNDWI and DEM zones added hidden.")
            except Exception as e:
                messages.addWarningMessage("Could not add layers to map: " + str(e))

        finally:
            arcpy.CheckInExtension("Spatial")

    def postExecute(self, parameters):
        "JUST BELIEVE"
        return
