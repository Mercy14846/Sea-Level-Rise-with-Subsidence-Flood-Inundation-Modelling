/* ============================================================
   LULC CLASSIFICATION — 1986, 2006, 2025
   Classes: 1=Built-up, 2=Vegetation, 3=Wetland, 4=Water, 5=Bareland
   Sensors : Landsat 5 TM (1986, 2006) | Landsat 9 OLI (2025)
   Method  : Spectral-index threshold classification
             (chosen because no ground-truth/training points exist
             for 1986 — keeps the method identical & comparable
             across all three years, which matters for the Markov
             transition analysis used in the projection step)
   ============================================================ */

// ---------------------------------------------------------
// 1. DEFINE YOUR STUDY AREA  <-- REPLACE THIS
// ---------------------------------------------------------
// Draw a polygon in the GEE code editor (Geometry tools) and
// import it as `studyArea`, OR paste coordinates below.
var studyArea = ee.Geometry.Polygon([
  [[3.80, 7.90], [3.80, 7.60], [4.20, 7.60], [4.20, 7.90]]
]); // <-- placeholder AOI, replace with your actual boundary

Map.centerObject(studyArea, 10);

// ---------------------------------------------------------
// 2. CLOUD MASK FOR LANDSAT COLLECTION 2 LEVEL-2 (SR)
// ---------------------------------------------------------
function maskL2sr(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);

  return image.addBands(opticalBands, null, true)
              .updateMask(qaMask)
              .updateMask(saturationMask);
}

// ---------------------------------------------------------
// 3. BUILD A CLASSIFIED IMAGE FROM A COMPOSITE
//    bandMap = {blue, green, red, nir, swir1, swir2}
// ---------------------------------------------------------
function classifyLULC(image, bandMap) {
  var blue  = image.select(bandMap.blue);
  var green = image.select(bandMap.green);
  var red   = image.select(bandMap.red);
  var nir   = image.select(bandMap.nir);
  var swir1 = image.select(bandMap.swir1);

  var ndvi  = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
  var mndwi = green.subtract(swir1).divide(green.add(swir1)).rename('MNDWI');
  var ndbi  = swir1.subtract(nir).divide(swir1.add(nir)).rename('NDBI');
  var lswi  = nir.subtract(swir1).divide(nir.add(swir1)).rename('LSWI');

  // Start with Bareland (5) as the default class
  var classified = ee.Image.constant(5).rename('LULC').byte();

  // Built-up: NDBI positive, low vegetation
  var builtMask = ndbi.gt(0).and(ndvi.lt(0.3));
  classified = classified.where(builtMask, 1);

  // Vegetation: high NDVI
  var vegMask = ndvi.gt(0.4);
  classified = classified.where(vegMask, 2);

  // Wetland: moderate vegetation + high moisture, not open water
  var wetMask = ndvi.gt(0.1).and(ndvi.lt(0.4))
                 .and(lswi.gt(0.1))
                 .and(mndwi.lt(0.2));
  classified = classified.where(wetMask, 3);

  // Water: high MNDWI (applied LAST so it always wins — water
  // signatures are the most reliable/unambiguous)
  var waterMask = mndwi.gt(0.2);
  classified = classified.where(waterMask, 4);

  return classified.clip(studyArea);
}

var bandMap_L5 = {blue:'SR_B1', green:'SR_B2', red:'SR_B3',
                   nir:'SR_B4', swir1:'SR_B5', swir2:'SR_B7'};
// Landsat 7 ETM+ C02 L2 uses the SAME band names as Landsat 5 TM,
// so bandMap_L5 works for merged L5+L7 collections too.

// ---------------------------------------------------------
// 3b. ROBUST COMPOSITE BUILDER
//     Fixes "Image with no bands" errors caused by an empty
//     ImageCollection (common for Lagos/W.Africa pre-2010,
//     where many Landsat 5 scenes were never downlinked to a
//     ground station and are missing from the archive for a
//     given single calendar year).
//     - Merges Landsat 5 + Landsat 7 for more scene availability
//     - Widens the date window automatically if the first
//       attempt returns zero images
//     - Prints collection size to the Console so you can SEE
//       whether it worked, instead of hitting a silent "no bands"
//       error three steps later
// ---------------------------------------------------------
function maskL2srC2(image) { return maskL2sr(image); }

function getComposite(centerYear, cloudThresh, label) {
  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2');

  function buildFor(startYear, endYear) {
    var start = startYear + '-01-01';
    var end = endYear + '-12-31';
    var merged = l5.merge(l7)
      .filterBounds(studyArea)
      .filterDate(start, end)
      .filter(ee.Filter.lt('CLOUD_COVER', cloudThresh))
      .map(maskL2srC2);
    return merged;
  }

  // Attempt 1: +/- 0 years (single calendar year)
  var col = buildFor(centerYear, centerYear);
  var size = col.size().getInfo();
  print('Scenes found for ' + label + ' (single year, cloud<' + cloudThresh + '):', size);

  // Attempt 2: widen to a 3-year window if empty
  if (size === 0) {
    col = buildFor(centerYear - 1, centerYear + 1);
    size = col.size().getInfo();
    print('Scenes found for ' + label + ' (widened to ' + (centerYear - 1) + '-' + (centerYear + 1) + '):', size);
  }

  // Attempt 3: widen further + relax cloud filter if still empty
  if (size === 0) {
    col = l5.merge(l7)
      .filterBounds(studyArea)
      .filterDate((centerYear - 3) + '-01-01', (centerYear + 3) + '-12-31')
      .filter(ee.Filter.lt('CLOUD_COVER', 80))
      .map(maskL2srC2);
    size = col.size().getInfo();
    print('Scenes found for ' + label + ' (widened +/-3yr, cloud<80):', size);
  }

  if (size === 0) {
    print('WARNING: still zero scenes for ' + label +
      '. Consider sourcing an external product for this period ' +
      '(e.g. your uploaded Lagos_LULC_2005_2006.tif) instead of a live pull.');
  }

  return col.median();
}

// ---------------------------------------------------------
// 4. 1986 — LANDSAT 5 TM (+ L7 fallback, widened window)
// ---------------------------------------------------------
var l5_1986 = getComposite(1986, 30, '1986');
var lulc1986 = classifyLULC(l5_1986, bandMap_L5);

// ---------------------------------------------------------
// 5. 2006 — LANDSAT 5/7 (+ widened window)
//    NOTE: if you already have a trusted Lagos_LULC_2005_2006.tif
//    (e.g. uploaded to Drive/Assets), you can skip this GEE
//    classification entirely and use that file directly as your
//    2006 input for the Markov-CA projection script instead.
// ---------------------------------------------------------
var l5_2006 = getComposite(2006, 30, '2006');
var lulc2006 = classifyLULC(l5_2006, bandMap_L5);

// ---------------------------------------------------------
// 6. 2025 — LANDSAT 9 OLI
//    (reuses maskL2srC2 defined in section 3b — same mask logic
//    works for OLI since QA_PIXEL/QA_RADSAT bands are identical
//    across Collection 2 Level-2 products)
// ---------------------------------------------------------
var l9_2025 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2025-01-01', '2025-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 30))
  .map(maskL2srC2)
  .median();

var bandMap_L9 = {blue:'SR_B2', green:'SR_B3', red:'SR_B4',
                   nir:'SR_B5', swir1:'SR_B6', swir2:'SR_B7'};

var lulc2025 = classifyLULC(l9_2025, bandMap_L9);

// ---------------------------------------------------------
// 7. VISUALIZATION
// ---------------------------------------------------------
var palette = ['d7191c', '1a9641', '80cdc1', '2c7fb8', 'ffffbf'];
// order matches classes 1..5: Built-up, Vegetation, Wetland, Water, Bareland
var visParams = {min: 1, max: 5, palette: palette};

Map.addLayer(lulc1986, visParams, 'LULC 1986');
Map.addLayer(lulc2006, visParams, 'LULC 2006');
Map.addLayer(lulc2025, visParams, 'LULC 2025');

// ---------------------------------------------------------
// 8. AREA STATISTICS PER CLASS (km²) — printed to Console
// ---------------------------------------------------------
function printAreaStats(image, label) {
  var areaImage = ee.Image.pixelArea().divide(1e6).addBands(image);
  var stats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: studyArea,
    scale: 30,
    maxPixels: 1e13
  });
  print('Area (km2) - ' + label, stats);
}
printAreaStats(lulc1986, '1986');
printAreaStats(lulc2006, '2006');
printAreaStats(lulc2025, '2025');

// ---------------------------------------------------------
// 9. EXPORT TO DRIVE AS GEOTIFF
//    Class values are preserved as an integer (byte) raster.
// ---------------------------------------------------------
Export.image.toDrive({
  image: lulc1986,
  description: 'LULC_1986',
  folder: 'LULC_Project',
  fileNamePrefix: 'LULC_1986',
  region: studyArea,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: lulc2006,
  description: 'LULC_2006',
  folder: 'LULC_Project',
  fileNamePrefix: 'LULC_2006',
  region: studyArea,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: lulc2025,
  description: 'LULC_2025',
  folder: 'LULC_Project',
  fileNamePrefix: 'LULC_2025',
  region: studyArea,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

/* ============================================================
   NOTES
   - If cloud-free composites are sparse for a given year,
     widen the date filter (e.g. Jan-Dec of adjacent year) or
     raise the CLOUD_COVER threshold.
   - For higher accuracy than threshold rules, replace
     classifyLULC() with a Random Forest classifier
     (ee.Classifier.smileRandomForest) trained on your own
     labeled points per era — ask if you want that version.
   - Run these three exports, then download the GeoTIFFs from
     Google Drive and use them as inputs to the Markov-CA
     projection script (Python) for 2030/2050/2100.
   ============================================================ */
