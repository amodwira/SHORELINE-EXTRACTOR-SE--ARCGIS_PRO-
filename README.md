# SE+ Shoreline Extractor (v2.0.0)

<!-- After the first Zenodo release, replace XXXXXXX with your record number: -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Automates shoreline extraction from a multispectral composite using **MNDWI + Otsu**, with adaptive smoothing, morphology cleanup, an optional **lidar DEM land/water constraint** (NOAA lidar) that **excludes marsh shorelines** (keeps beach / non-marsh coast only), polygon→line conversion (keeps boundary), and projection to a target CRS.

## Quick start (ArcGIS Pro)
The toolbox is the file **`Shoreline_Extractor.pyt`**, an ArcGIS Pro *Python toolbox*. No other ArcGIS file is needed.

1. Download this repository (green **Code** button → **Download ZIP**) and unzip it. Keep `Shoreline_Extractor.pyt` and the two `.pyt.xml` files together in one folder.
2. Install **scikit-image** once (see *Install* below).
3. Open any ArcGIS Pro project (or create a new blank one).
4. In the **Catalog** pane, right-click **Toolboxes** → **Add Toolbox**, and select `Shoreline_Extractor.pyt`.
5. Expand the toolbox and double-click **Shoreline Extractor (SE+ v2.0)**.
6. Fill in the inputs and click **Run**. Outputs are added to the map: the segmentation map, the extracted shoreline and, with a DEM, the marsh shoreline.

To keep it in every project, right-click the toolbox in Catalog → **Add To Favorites**.

## Requirements
- **ArcGIS Pro 3.x (Windows)**
- **Spatial Analyst** extension
- Python libs:
  - **Bundled**: `arcpy`, `numpy`
  - **Install**: `scikit-image` (for `threshold_otsu`, `gaussian`, `disk/opening/closing/remove_small_objects`)

### Install (recommended)
```bat
# Open "Python Command Prompt" from ArcGIS Pro
conda create -n shoreline-env --clone arcgispro-py3
conda activate shoreline-env
conda install -c conda-forge scikit-image
```
Quick check:
```bat
python -c "import arcpy, numpy; from skimage.filters import threshold_otsu, gaussian; from skimage.morphology import disk, opening, closing, remove_small_objects; print('OK')"
```

## Inputs (Tool Parameters)
- **Composite Raster** (GPRasterLayer)
- **Green Band Name** (e.g., `B3`, `Band_3`, `Layer_3`)
- **SWIR Band Name** (e.g., `B11`, `Band_11`, `Layer_11`)
- **Study Area** *(optional, polygon)*
- **Output Folder**
- **Output File Prefix**
- **Output Coordinate System** (target CRS)
- **Remove study-area/extent border segments** (default on)

### DEM constraint (NOAA lidar) — optional, v2.6
| Parameter | Default | Meaning |
|---|---|---|
| Lidar DEM | — | Bare-earth DEM (GeoTIFF/IMG). Leave empty to run exactly as v2.5. |
| DEM vertical units | Meters | Meters, US Survey Feet or International Feet. |
| MHW elevation | *(required with DEM)* | MHW in the DEM's vertical datum, in **metres** (e.g., MHW above NAVD88). Enter 0 if the DEM is already referenced to MHW. |
| Land margin | 0.5 m | Elevation > MHW + margin → **land for certain**. |
| Water margin | 0.5 m | Elevation < MHW − margin → **water for certain**. |
| Hydro-flattened surfaces as water | on | Perfectly flat, low DEM areas (flattened lakes/bays) → water. |
| Keep water connected through DEM channels | on | Keeps bays/ponds linked to the Gulf by channels too narrow for the imagery. |
| Red / NIR band | — | Needed for marsh = land and vegetation correction (Sentinel-2: B4 / B8). |
| Vegetation bias | 0.15 m | Subtracted from the DEM where NDVI > threshold (lidar reads high in dense marsh). |
| NDVI threshold | 0.3 | Vegetation cut-off for the bias correction. |
| **Exclude marsh shorelines** *(v2.8)* | on | Removes every shoreline segment that borders marsh (low-lying vegetation near MHW). The final shoreline covers only beach / sandy / hard coast. Removed segments are saved separately. |
| Marsh search radius *(v2.8.2)* | 2.5 px | A segment is marsh shoreline when at least half of the land pixels within this radius are marsh. Looks past the mixed water/grass edge pixel. Beaches ≥ ~3 px (30 m at 10 m) are kept; lower to 1.5 for very narrow beaches. |
| Fill MNDWI false water inside marsh *(v2.7)* | on | Internal clean-up: vegetated pixels near MHW that MNDWI calls water are filled, so they don't create stray lines. |
| NDVI threshold defining marsh vegetation | 0.3 | Raise (0.4) if marsh is over-detected; lower (0.2) for sparse marsh. |
| Treat lidar ground at/above MHW as land *(v2.7)* | off | In the uncertain zone, non-flattened lidar ground at or above MHW is land. For sparse marsh / mudflat platforms NDVI misses. |

**Shoreline definition (v2.8):** the land–water boundary along **non-marsh** coast (beaches, sandy and hard shores). Marsh is identified as vegetation (NDVI above threshold) lying in the DEM-uncertain zone near MHW; shoreline bordering marsh is excluded.

## What the Tool Does (Pipeline)
1. Compute **MNDWI = (Green−SWIR)/(Green+SWIR)**.
2. **Adaptive Gaussian smoothing** (sigma chosen from raster cell size).
3. **Otsu threshold** to binary (fallback: `MNDWI > 0`).
4. *(DEM)* Zone every pixel: **land-sure / water-sure / uncertain**; force land-sure and water-sure. *(v2.7)* Inside the uncertain zone, vegetated marsh (NDVI) — and optionally lidar ground ≥ MHW — is forced to land; the rest is left to MNDWI.
5. Morphology: **opening → small-object removal → closing → hole fill**; *(DEM)* zones re-applied afterwards.
6. *(DEM)* Ocean connectivity in raster domain (4-connected; low DEM cells not crossed by high ground act as channels). Without DEM: largest boundary-touching water polygon, as before.
7. Raster→Polygon, select water class, small-area filter (map-units²).
8. Polygon→**FeatureToLine** (fallback: PolygonToLine + SplitLine).
9. **MultipartToSinglepart** + **SplitLine** → fine segments (keep boundary).
10. Define source projection from raster → **Project** to user CRS.
11. Add outputs to current map (if running inside ArcGIS Pro).

## Outputs (in Output Folder)
- `{prefix}_mndwi.tif`
- `{prefix}_binary.tif`
- `{prefix}_water_polygon.shp`
- `{prefix}_valid_polygons.shp` *(and possibly `{prefix}_filtered_polygons.shp`)*
- `{prefix}_raster_extent.shp` *(when no study area)*
- `{prefix}_mask_boundary.shp` *(reference only; kept)*
- `{prefix}_shoreline_segments.shp` *(source CRS)*
- `{prefix}_shoreline_segments_projected.shp` *(final)*
- `{prefix}_mndwi_debug.tif` *(only if no water pixels were detected)*
- `{prefix}_metadata.json` *(parameters, thresholds, DEM statistics)*
- `{prefix}_dem_zones.tif` *(DEM only: 0 uncertain → MNDWI, 1 land-sure, 2 water-sure, 3 marsh→land by vegetation, 4 marsh→land by lidar ground ≥ MHW, NoData = no DEM)*
- `{prefix}_dem_on_grid.tif` *(DEM only: DEM mean, metres, on the composite grid)*
- `{prefix}_segmentation.tif` *(segmentation map: 1 open water (shoreline-forming), 2 inland / unconnected water, 3 land, 4 marsh (excluded); class names in the CLASS field)*
- `{prefix}_marsh_mask.tif` *(v2.8: 1 = marsh, 0 = not marsh)*
- `{prefix}_dem_overrides.tif` *(v2.8.1: where the DEM changed MNDWI — 1 water→land by elevation, 2 land→water by low elevation, 3 land→water by flat surface, 4 marsh fill)*
- `{prefix}_marsh_shoreline_segments(_projected).shp` *(v2.8: the excluded marsh-edge segments, for reference)*
- `{prefix}_dem_work/` *(DEM only: intermediate clip/project/aggregate rasters; safe to delete)*

## What appears on the map after a run
| Layer | Shown | Style |
|---|---|---|
| `{prefix} Extracted shoreline (non-marsh)` | on, top | red, 2.5 pt |
| `{prefix} Marsh shoreline (excluded)` | on | thin green (only when marsh exclusion ran) |
| `{prefix} Segmentation` | on | dark blue open water · light blue inland water · tan land · green marsh |
| `{prefix} DEM zones`, `{prefix} DEM overrides`, `{prefix} MNDWI` | added but hidden | turn on to diagnose |

If the colours are not applied automatically (depends on ArcGIS Pro version), set *Unique Values* on the `Value` field of the segmentation layer; labels are in the `CLASS` field.

## Class codes in the output rasters
Each classified raster stores these names in its **CLASS** field and these colours in its colormap, so they show in the
legend, the Identify tool and the attribute table.

**`{prefix}_segmentation.tif`**
| Value | Class |
|---|---|
| 1 | Open water (shoreline-forming) |
| 2 | Inland / unconnected water |
| 3 | Land |
| 4 | Marsh (excluded) |

**`{prefix}_dem_zones.tif`** (lidar decision per pixel)
| Value | Class |
|---|---|
| 0 | Near MHW, not classed as marsh – MNDWI decides |
| 1 | High ground – land (lidar) |
| 2 | Low or flat water surface – water (lidar) |
| 3 | Marsh – vegetated, near MHW |
| 4 | Bare ground at/above MHW – land (lidar) |

**`{prefix}_dem_overrides.tif`** (where lidar changed MNDWI)
| Value | Class |
|---|---|
| 0 | No change from MNDWI |
| 1 | MNDWI water → land (high ground) |
| 2 | MNDWI land → water (low elevation) |
| 3 | MNDWI land → water (flat surface) |
| 4 | MNDWI water → land (marsh fill) |

**`{prefix}_marsh_mask.tif`**: 0 not marsh, 1 marsh.

## Preparing NOAA lidar (Louisiana)
1. **Download** from NOAA Digital Coast Data Access Viewer (coast.noaa.gov/dataviewer) → Elevation → Lidar. Request a **bare-earth DEM** (GeoTIFF) rather than points; choose the composite's horizontal system if offered (e.g., NAD83 UTM 15N) — otherwise the tool projects it.
2. **Match dates**: pick a collection as close as possible to the image date (Louisiana subsidence and marsh loss make old DEMs drift).
3. **Vertical datum**: most Louisiana lidar is NAVD88 (GEOID12B/GEOID18). Get **MHW relative to NAVD88** from the nearest NOAA Tides & Currents station *Datums* page (MHW − NAVD88, same epoch/units) or from NOAA VDatum, and enter it in metres.
4. **Check the DEM**: water may be hydro-flattened (flat surface) or NoData. NoData pixels are left to MNDWI.
5. **Inspect** `{prefix}_dem_zones.tif` over the imagery after the first run and adjust the margins.

### Tuning notes
- Louisiana marsh platforms sit within a few tenths of a metre of MHW, so they fall in the **uncertain** zone. In v2.7 the marsh = land rule decides vegetated pixels there; only unvegetated pixels are left to MNDWI. The elevation zones fix ridges/levees/spoil banks (false water) and floating vegetation over flattened/deep water (false land).
- Check `{prefix}_metadata.json` → `marsh_exclusion`: total, kept (non-marsh) and excluded (marsh) shoreline length.
- Vegetated ridges, cheniers and dunes above MHW + land margin are **not** marsh (they are land-sure), so shoreline in front of them is kept.
- Floating vegetation over water that the DEM does **not** mark as water (outside DEM coverage, or unflattened small ponds) will be treated as marsh land. Inspect zone 3 over such areas.
- Leveed areas below sea level are marked water-sure by elevation but are removed by the ocean-connectivity step because the levee blocks them. Keep connectivity on when your area includes polders.
- If too many interior marsh ponds are kept as “ocean-connected”, increase the water margin (channels must be deeper to count).

## Usage (ArcGIS Pro)
1. Add the `.pyt` to your project (Catalog → **Add**).
2. Open **Shoreline Extractor**.
3. Set parameters (see above) and **Run**.

## Tips & Troubleshooting
- **No water found**: Verify band names and AOI; inspect `{prefix}_mndwi_debug.tif`.
- **DEM step failed**: the tool falls back to MNDWI only and says why in the messages; check units, MHW and that the DEM overlaps the image.
- **Band names**: Must match the internal band labels in the composite (e.g., `B3`, `B11`).
- **Performance**: The lidar steps (clip, reproject, 3×3 flatness test at 1 m, aggregation) dominate run time on large scenes; use a Study Area to test quickly. Marsh exclusion is fast since v2.8.1.
- **Licensing**: Requires **Spatial Analyst** (tool checks out at runtime).

## Testing
The core NumPy logic (elevation zones, marsh rules, connectivity, segment labelling, edge filter, class labels) is tested outside ArcGIS with `arcpy` stubbed:
```bash
pip install -r tests/requirements-test.txt
python tests/test_dem.py      # and the other files in tests/
```
The ArcGIS steps themselves (raster I/O, polygon conversion, projection) must be checked in ArcGIS Pro.

## Citing
If you use SE+ in published work, please cite the archived release (DOI on the Zenodo record) using the metadata in `CITATION.cff`. GitHub's **Cite this repository** button produces APA and BibTeX from that file.

## Authors
Developed by Curtis Amo Dwira (Louisiana State University), with supervision from Dr. Ahmmed Abdallah (Louisiana State University).

## License
MIT — see `LICENSE`. Version history: `CHANGELOG.md`.

