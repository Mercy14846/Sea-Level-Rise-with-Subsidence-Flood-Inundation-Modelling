// ==========================================================
// LAGOS RELATIVE SEA-LEVEL-RISE CLIMATE RISK MODEL
// ==========================================================


// ----------------------------------------------------------
// 1. STUDY AREA
// ----------------------------------------------------------

var lagos = ee.Geometry.Rectangle([2.65, 6.30, 4.35, 6.75]);

Map.centerObject(lagos, 9);


// ----------------------------------------------------------
// 2. DEM
// ----------------------------------------------------------

var dem = ee.Image('NASA/NASADEM_HGT/001')
  .select('elevation')
  .clip(lagos);

Map.addLayer(
  dem,
  {min: 0, max: 20},
  'Elevation'
);


// ----------------------------------------------------------
// 3. WATER MASK
// ----------------------------------------------------------

var waterOccurrence = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
  .select('occurrence')
  .clip(lagos);

var waterMask = waterOccurrence.gt(50);

Map.addLayer(
  waterMask.selfMask(),
  {palette: ['blue']},
  'Water Bodies'
);


// ----------------------------------------------------------
// 4. DISTANCE TO WATER
// ----------------------------------------------------------

var distanceToWater = waterMask.not()
  .fastDistanceTransform(512)
  .sqrt()
  .multiply(30)
  .rename('distance_to_water');


// ----------------------------------------------------------
// 5. BUILT-UP AREAS
// ----------------------------------------------------------

// IMPORTANT:
// Select only ONE band

var builtup = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_S')
  .filterBounds(lagos)
  .mosaic()
  .select(0)   // FIX
  .clip(lagos);

var builtMask = builtup.gt(0)
  .rename('built_risk');

Map.addLayer(
  builtMask.selfMask(),
  {palette: ['gray']},
  'Built-up Areas'
);


// ----------------------------------------------------------
// 6. POPULATION
// ----------------------------------------------------------

var population = ee.ImageCollection(
  'WorldPop/GP/100m/pop'
)
.filterBounds(lagos)
.filterDate('2020-01-01', '2020-12-31')
.mean()
.select(0) // FIX
.clip(lagos);

Map.addLayer(
  population,
  {min: 0, max: 1000},
  'Population Density'
);


// ----------------------------------------------------------
// 7. SUBSIDENCE
// ----------------------------------------------------------

var subsidence = ee.Image.constant(4)
  .rename('subsidence_mm_yr')
  .clip(lagos);

Map.addLayer(
  subsidence,
  {min: 0, max: 10, palette: ['yellow', 'red']},
  'Subsidence'
);


// ----------------------------------------------------------
// 8. SEA LEVEL RISE SCENARIOS
// ----------------------------------------------------------

var years = 2050 - 2020;

var oceanSLR_med = 3.6;


// ----------------------------------------------------------
// 9. RELATIVE SLR
// ----------------------------------------------------------

var ocean = ee.Image.constant(oceanSLR_med);

var relativeRate = ocean.add(subsidence);

var rslr = relativeRate
  .multiply(years)
  .divide(1000)
  .rename('rslr_m');


// ----------------------------------------------------------
// 10. FLOOD EXPOSURE
// ----------------------------------------------------------

var floodExposure = dem.lte(rslr)
  .selfMask()
  .rename('flood_exposure');

Map.addLayer(
  floodExposure,
  {palette: ['cyan']},
  'Flood Exposure'
);


// ----------------------------------------------------------
// 11. NORMALIZED RISK FACTORS
// ----------------------------------------------------------

// Elevation risk
var elevationRisk = dem.expression(
  '1 - min(elev / 15, 1)', {
    elev: dem
}).rename('elev_risk');


// Subsidence risk
var subsidenceRisk = subsidence.expression(
  'min(sub / 10, 1)', {
    sub: subsidence
}).rename('subsidence_risk');


// Water proximity risk
var waterRisk = distanceToWater.expression(
  '1 - min(dist / 5000, 1)', {
    dist: distanceToWater
}).rename('water_risk');


// Population risk
var popRisk = population.expression(
  'min(pop / 1000, 1)', {
    pop: population
}).rename('pop_risk');


// ----------------------------------------------------------
// 12. CLIMATE RISK INDEX
// ----------------------------------------------------------

var climateRisk = elevationRisk.multiply(0.30)
  .add(subsidenceRisk.multiply(0.25))
  .add(waterRisk.multiply(0.20))
  .add(popRisk.multiply(0.15))
  .add(builtMask.multiply(0.10))
  .rename('climate_risk');

Map.addLayer(
  climateRisk,
  {
    min: 0,
    max: 1,
    palette: ['green', 'yellow', 'orange', 'red']
  },
  'Climate Risk Index'
);


// ----------------------------------------------------------
// 13. RISK CLASSIFICATION
// ----------------------------------------------------------

var riskClass = climateRisk.expression(
  "(risk <= 0.25) ? 1" +
  ": (risk <= 0.50) ? 2" +
  ": (risk <= 0.75) ? 3" +
  ": 4", {
    risk: climateRisk
}).rename('risk_class');

Map.addLayer(
  riskClass,
  {
    min: 1,
    max: 4,
    palette: ['green', 'yellow', 'orange', 'red']
  },
  'Risk Classes'
);


// ----------------------------------------------------------
// 14. RISK AREA STATISTICS
// ----------------------------------------------------------

var areaImage = ee.Image.pixelArea()
  .divide(1e6)
  .rename('area_km2');

var combined = areaImage.addBands(riskClass);

var stats = combined.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'risk_class'
  }),
  geometry: lagos,
  scale: 30,
  maxPixels: 1e13
});

print('Risk Area Statistics', stats);


// ----------------------------------------------------------
// 15. POPULATION AT RISK
// ----------------------------------------------------------

var highRiskMask = riskClass.gte(3);

var populationAtRisk = population
  .updateMask(highRiskMask);

var popStats = populationAtRisk.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: lagos,
  scale: 100,
  maxPixels: 1e13
});

print('Population at High Risk', popStats);


// ----------------------------------------------------------
// 16. EXPORTS
// ----------------------------------------------------------

Export.image.toDrive({
  image: climateRisk,
  description: 'Lagos_Climate_Risk_Index',
  region: lagos,
  scale: 30,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: riskClass,
  description: 'Lagos_Risk_Classes',
  region: lagos,
  scale: 30,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: floodExposure,
  description: 'Lagos_Flood_Exposure',
  region: lagos,
  scale: 30,
  maxPixels: 1e13
});
