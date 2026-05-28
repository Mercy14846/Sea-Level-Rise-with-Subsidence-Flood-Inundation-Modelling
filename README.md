# Sea Level Rise with Subsidence Flood Inundation Modelling
This is a bathtub inundation model. It assumes every low-lying area below the water level floods. That is useful for first-level screening, but it can overestimate flooding because it ignores drainage, barriers, canals, flow connectivity, tide timing, and coastal defenses.

    Bash
    pip install rasterio geopandas shapely pandas numpy


# Lagos Relative Sea-Level Rise Climate Risk Model (Technical Documentation)

## 1. Overview

This script implements a **spatial climate risk assessment model for Lagos, Nigeria** using Google Earth Engine (GEE). The model integrates elevation, water occurrence, built-up extent, population density, and subsidence data to generate a **composite Climate Risk Index (0–1)** and a categorical risk classification (1–4). It also estimates flood exposure, population at risk, and exports geospatial outputs for further analysis.

---

## 2. Study Area Definition

The study area is defined as a rectangular bounding box representing Lagos State using geographic coordinates (longitude and latitude).

```javascript
var lagos = ee.Geometry.Rectangle([2.65, 6.30, 4.35, 6.75]);
Map.centerObject(lagos, 9);
```

* West Boundary: 2.65°E
* South Boundary: 6.30°N
* East Boundary: 4.35°E
* North Boundary: 6.75°N

The map is centered on Lagos at zoom level 9 for visualization.

---

## 3. Digital Elevation Model (DEM)

Elevation data is derived from NASADEM (NASA/NASADEM_HGT/001) and clipped to the study area.

```javascript
var dem = ee.Image('NASA/NASADEM_HGT/001')
  .select('elevation')
  .clip(lagos);
```

Visualization parameters:

* Minimum elevation: 0 m
* Maximum elevation: 20 m

The DEM serves as a primary exposure factor in flood and sea-level rise modelling.

---

## 4. Surface Water Occurrence and Water Mask

Permanent and seasonal water bodies are derived from the JRC Global Surface Water dataset.

```javascript
var waterOccurrence = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
  .select('occurrence')
  .clip(lagos);

var waterMask = waterOccurrence.gt(50);
```

* Pixels with >50% water occurrence are classified as water bodies.
* The resulting binary mask identifies persistent water features.

---

## 5. Distance to Water Bodies

A proximity layer is computed using a Euclidean distance transformation from non-water pixels.

```javascript
var distanceToWater = waterMask.not()
  .fastDistanceTransform(512)
  .sqrt()
  .multiply(30)
  .rename('distance_to_water');
```

This generates a continuous raster representing distance (in meters) to nearest water body.

---

## 6. Built-Up Areas

Built-up intensity is derived from the Global Human Settlement Layer (GHSL P2023A).

```javascript
var builtup = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_S')
  .filterBounds(lagos)
  .mosaic()
  .select(0)
  .clip(lagos);

var builtMask = builtup.gt(0)
  .rename('built_risk');
```

* Binary classification:

  * 1 = built-up area
  * 0 = non-built-up area

---

## 7. Population Density

Population data is sourced from WorldPop (2020 annual dataset).

```javascript
var population = ee.ImageCollection('WorldPop/GP/100m/pop')
  .filterBounds(lagos)
  .filterDate('2020-01-01', '2020-12-31')
  .mean()
  .select(0)
  .clip(lagos);
```

* Represents population per 100m grid cell.
* Used as a socio-economic exposure variable.

---

## 8. Land Subsidence

A constant subsidence rate is assumed for Lagos.

```javascript
var subsidence = ee.Image.constant(4)
  .rename('subsidence_mm_yr')
  .clip(lagos);
```

* Value: 4 mm/year
* Represents vertical land movement contributing to relative sea-level rise.

---

## 9. Sea-Level Rise Scenario

A medium sea-level rise scenario is applied.

```javascript
var years = 2050 - 2020;
var oceanSLR_med = 3.6;
```

* Time horizon: 30 years (2020–2050)
* Sea-level rise rate: 3.6 mm/year (ocean component)

---

## 10. Relative Sea-Level Rise (RSLR)

Relative sea-level rise combines ocean rise and land subsidence.

```javascript
var ocean = ee.Image.constant(oceanSLR_med);
var relativeRate = ocean.add(subsidence);

var rslr = relativeRate
  .multiply(years)
  .divide(1000)
  .rename('rslr_m');
```

This produces cumulative relative sea-level rise in meters.

---

## 11. Flood Exposure Mapping

Flood exposure is derived by comparing elevation with projected relative sea-level rise.

```javascript
var floodExposure = dem.lte(rslr)
  .selfMask()
  .rename('flood_exposure');
```

* Pixels where elevation ≤ RSLR are classified as exposed.
* Output represents potential inundation zones.

---

## 12. Normalized Risk Factors

All risk variables are normalized to a 0–1 scale.

### 12.1 Elevation Risk

```javascript
var elevationRisk = dem.expression(
  '1 - min(elev / 15, 1)', {
    elev: dem
});
```

Lower elevation corresponds to higher risk.

---

### 12.2 Subsidence Risk

```javascript
var subsidenceRisk = subsidence.expression(
  'min(sub / 10, 1)', {
    sub: subsidence
});
```

Higher subsidence increases risk.

---

### 12.3 Water Proximity Risk

```javascript
var waterRisk = distanceToWater.expression(
  '1 - min(dist / 5000, 1)', {
    dist: distanceToWater
});
```

Closer proximity to water increases risk.

---

### 12.4 Population Risk

```javascript
var popRisk = population.expression(
  'min(pop / 1000, 1)', {
    pop: population
});
```

Higher population density increases exposure.

---

## 13. Climate Risk Index

A weighted multi-criteria index is computed as follows:

```javascript
var climateRisk = elevationRisk.multiply(0.30)
  .add(subsidenceRisk.multiply(0.25))
  .add(waterRisk.multiply(0.20))
  .add(popRisk.multiply(0.15))
  .add(builtMask.multiply(0.10));
```

### Weight Distribution:

* Elevation: 30%
* Subsidence: 25%
* Water proximity: 20%
* Population: 15%
* Built-up areas: 10%

The output is a normalized raster (0 = low risk, 1 = high risk).

---

## 14. Risk Classification

The continuous risk index is reclassified into four ordinal categories.

```javascript
var riskClass = climateRisk.expression(
  "(risk <= 0.25) ? 1" +
  ": (risk <= 0.50) ? 2" +
  ": (risk <= 0.75) ? 3" +
  ": 4", {
    risk: climateRisk
});
```

### Classes:

* 1 = Low risk
* 2 = Moderate risk
* 3 = High risk
* 4 = Very high risk

---

## 15. Risk Area Statistics

Area statistics are computed per risk class using zonal aggregation.

```javascript
var areaImage = ee.Image.pixelArea()
  .divide(1e6)
  .rename('area_km2');

var combined = areaImage.addBands(riskClass);
```

A grouped reducer calculates total area per risk category.

---

## 16. Population at Risk

High-risk zones (classes ≥ 3) are used to estimate exposed population.

```javascript
var highRiskMask = riskClass.gte(3);

var populationAtRisk = population
  .updateMask(highRiskMask);
```

Population values are aggregated over high-risk zones.

---

## 17. Data Export

Final outputs are exported to Google Drive for further analysis and visualization.

### 17.1 Climate Risk Index

```javascript
Export.image.toDrive({
  image: climateRisk,
  description: 'Lagos_Climate_Risk_Index'
});
```

### 17.2 Risk Classes

```javascript
Export.image.toDrive({
  image: riskClass,
  description: 'Lagos_Risk_Classes'
});
```

### 17.3 Flood Exposure

```javascript
Export.image.toDrive({
  image: floodExposure,
  description: 'Lagos_Flood_Exposure'
});
```

---

## 18. Summary

This model integrates physical, environmental, and socio-economic datasets to produce a **multi-factor climate risk assessment framework for Lagos**. The outputs support flood risk mapping, urban vulnerability analysis, and climate adaptation planning.
