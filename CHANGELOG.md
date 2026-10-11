# Changelog

## 2.0.0 (2026-10-11): first public release

Shown in ArcGIS Pro as **Shoreline Extractor (SE+ v2.0)**.

- MNDWI with robust Otsu thresholding, adaptive smoothing and morphology; ocean-connected water; boundary shoreline segments; border and image-edge filtering; projection; metadata JSON.
- Optional lidar DEM constraint referenced to Mean High Water: high ground forced to land, low or flat water surfaces forced to water, the near-MHW zone left to MNDWI.
- NDVI vegetation-bias correction of lidar; flat lidar surfaces treated as water only when not vegetated.
- Marsh identified as vegetation in the near-MHW zone; MNDWI false water inside marsh filled.
- Shoreline segments labelled marsh or non-marsh by a majority vote of nearby land pixels; non-marsh and marsh shorelines delivered as separate layers.
- DEM-assisted connectivity to the open coast through sub-pixel channels, with levees as barriers.
- Segmentation map (open water, inland water, land, marsh) and diagnostic rasters (DEM zones, DEM overrides, marsh mask), each storing class names and colours in the file.
- Checks: warns when the DEM does not overlap the image; validates inputs on the tool form.
- Fix relative to the earlier SE+ v2.5: bands read on one reference grid (Study Area alignment); band NoData carried as NaN.

## Development history (unreleased builds)

These internal builds led to 2.0.0 and were not published.

- **2.8.5**: DEM/image overlap warning.
- **2.8.4**: class names and colours stored in every classified raster.
- **2.8.3**: image-edge segment filter; segmentation colours embedded in the file.
- **2.8.2**: majority-vote marsh test with search radius.
- **2.8.1**: fast marsh exclusion; vegetated-flat guard; DEM overrides raster.
- **2.8**: marsh shorelines excluded and saved separately; segmentation map.
- **2.7**: marsh fill; optional lidar ground ≥ MHW rule.
- **2.6**: lidar DEM constraint, vegetation correction, connectivity; grid-alignment fix.
- **2.5**: original MNDWI–Otsu toolbox (no DEM).
