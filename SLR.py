# =====================================================
# Sea Level Rise + Subsidence Flood Inundation Modelling
# =====================================================

import rasterio
import numpy as np
import geopandas as gpd
from rasterio.mask import mask
from rasterio.features import shapes
from shapely.geometry import shape
import pandas as pd

# -----------------------------
# 1. INPUT FILES
# -----------------------------

dem_path = "lagos_dem.tif"          # DEM in meters
aoi_path = "lagos_boundary.shp"     # Lagos boundary shapefile
output_raster = "lagos_inundation_2050.tif"
output_vector = "lagos_inundation_2050.shp"
output_csv = "inundation_summary.csv"

# -----------------------------
# 2. USER PARAMETERS
# -----------------------------

baseline_year = 2020
target_year = 2050

# Sea-level-rise rate in mm/year
sea_level_rise_rate = 3.6

# Lagos land subsidence rate in mm/year
subsidence_rate = 4.0

# Optional storm surge height in meters
storm_surge_m = 0.5

# -----------------------------
# 3. RELATIVE SEA LEVEL RISE
# -----------------------------

years = target_year - baseline_year

relative_slr_m = ((sea_level_rise_rate + subsidence_rate) * years) / 1000

total_water_level_m = relative_slr_m + storm_surge_m

print(f"Relative SLR by {target_year}: {relative_slr_m:.3f} m")
print(f"Total flood water level: {total_water_level_m:.3f} m")

# -----------------------------
# 4. READ AOI AND DEM
# -----------------------------

aoi = gpd.read_file(aoi_path)

with rasterio.open(dem_path) as src:
    aoi = aoi.to_crs(src.crs)

    dem_clip, transform = mask(
        src,
        aoi.geometry,
        crop=True,
        filled=True,
        nodata=src.nodata
    )

    profile = src.profile.copy()
    profile.update({
        "height": dem_clip.shape[1],
        "width": dem_clip.shape[2],
        "transform": transform,
        "count": 1,
        "dtype": "uint8",
        "nodata": 0
    })

dem = dem_clip[0]

# -----------------------------
# 5. CREATE INUNDATION MODEL
# -----------------------------

# Areas lower than projected water level are exposed
inundation = np.where(
    (dem <= total_water_level_m) & (dem >= 0),
    1,
    0
).astype("uint8")

# -----------------------------
# 6. SAVE INUNDATION RASTER
# -----------------------------

with rasterio.open(output_raster, "w", **profile) as dst:
    dst.write(inundation, 1)

print(f"Saved raster: {output_raster}")

# -----------------------------
# 7. CONVERT RASTER TO VECTOR
# -----------------------------

results = (
    {"properties": {"inundated": value}, "geometry": shape(geom)}
    for geom, value in shapes(
        inundation,
        mask=inundation == 1,
        transform=transform
    )
)

geoms = list(results)

inundation_gdf = gpd.GeoDataFrame.from_features(geoms, crs=aoi.crs)

inundation_gdf.to_file(output_vector)

print(f"Saved vector: {output_vector}")

# -----------------------------
# 8. CALCULATE INUNDATED AREA
# -----------------------------

# Reproject to metric CRS for area calculation
# UTM Zone 31N is suitable for Lagos
inundation_metric = inundation_gdf.to_crs("EPSG:32631")

inundation_area_km2 = inundation_metric.area.sum() / 1e6

summary = pd.DataFrame([{
    "baseline_year": baseline_year,
    "target_year": target_year,
    "sea_level_rise_rate_mm_yr": sea_level_rise_rate,
    "subsidence_rate_mm_yr": subsidence_rate,
    "relative_slr_m": relative_slr_m,
    "storm_surge_m": storm_surge_m,
    "total_water_level_m": total_water_level_m,
    "inundated_area_km2": inundation_area_km2
}])

summary.to_csv(output_csv, index=False)

print(summary)
print(f"Saved summary: {output_csv}")
