//Take Songnen Plain (CBS) in 2021 as an example

//Load the samples file, grid file and shp file for the study area
var agr_region = ee.FeatureCollection('users/xuanhy12345678/shp_file/agriculture_region');
var samples = ee.FeatureCollection('users/xuanhy12345678/samples/2021CBS_samples');
var grid15 = ee.FeatureCollection('users/103b6011/class/grid15');
var region_name = '长白山脉区';
var region = agr_region.filter(ee.Filter.eq('NAME',region_name));

//Load the cropland mask 
var cropland = ee.ImageCollection('users/potapovpeter/Global_cropland_2019').first().clip(region).eq(1);

//Define the date and using bands
var year = 2021;
var bands = ee.List(['B6','B7','NDVI','EVI','LSWI','NDSVI','NDTI','GCVI']); 
var start_day = ee.Date.fromYMD(year,1,1);
var end_day =  ee.Date.fromYMD(year+1,1,1);

function maskL8sr(image) {
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask).divide(10000)
              .set('system:time_start',image.get('system:time_start'))
              .set('system:index',image.get('system:index'));
};
function addIndices(image){
  
  var DOY = image.date().getRelative('day', 'year')
  var year = image.date().get('year')
  
  var ndvi = image.normalizedDifference(['B5', 'B4']).rename(['ndvi']);
  var nir = image.select('B5');
  var red = image.select('B4');
  var blue = image.select('B2');
  var green = image.select('B3');
  var evi = image.expression("2.5 * (B5 - B4) / (B5 + 6*B4 - 7.5*B2 + 1)",
  {
    'B5': image.select('B5'),
    'B4': image.select("B4"),
    'B2': image.select('B2')
  }).rename(['evi']);
  var lswi = image.normalizedDifference(['B3','B6']).rename(['lswi']);
  var ndsvi = image.normalizedDifference(['B6','B4']).rename(['ndsvi']);
  var ndti = image.normalizedDifference(['B7','B6']).rename(['ndti']);
  var gcvi = image.expression("((B5 / B3) - 1)",
  {
    'B5':image.select('B5'),
    'B3':image.select('B3')
  }).rename(['gcvi']);
  return image.addBands(ndvi.toDouble().rename("NDVI")).addBands(evi.toDouble().rename("EVI"))
  .addBands(lswi.toDouble().rename("LSWI")).addBands(gcvi.toDouble().rename("GCVI"))
  .addBands(ndsvi.toDouble().rename("NDSVI")).addBands(ndti.toDouble().rename("NDTI"))
  .addBands(ee.Image(DOY).rename('DOY').toDouble())
  .addBands(ee.Image(year).rename('Year').toDouble())
  .set('DOY',DOY)
};
function mergeBands(image,pre){
  return ee.Image(pre).addBands(image)
};
function addTimeBands(image) {
  var tstamp = ee.Date(image.get('system:time_start'));
  var tdelta = tstamp.difference(start_day, 'year');
  var img_fitting = image.select()
    .addBands(1)
    .addBands(tdelta.multiply(3*Math.PI).sin())
    .addBands(tdelta.multiply(3*Math.PI).cos())
    .addBands(tdelta.multiply(6*Math.PI).sin())
    .addBands(tdelta.multiply(6*Math.PI).cos())
    .addBands(image.select('NDVI','EVI','LSWI'))
    .toDouble();
  return img_fitting;
};

//Select image set
var l8img = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
                  .filterBounds(region)
                  .filterDate(start_day,end_day)
                  .map(maskL8sr)
                  .map(addIndices);
                  
var l8img_col = l8img.select(bands);

//Monthly Composite Image
var monthcomposite = ee.List.sequence(0, 11).map(function(num) {
  var start = ee.Date(start_day).advance(num, 'month');
  var end = start.advance(1, 'month');
  var doy = start.getRelative('day','year');
  var l8_month = l8img_col.filterDate(start, end).median().clip(region);
  var band_length = l8_month.bandNames().length();
  return l8_month.addBands(ee.Image.constant(doy).rename('doy').float())
                .set('system:time_start',ee.Date.fromYMD(year,1,1).advance(doy,'day').millis())
                .set('doy',doy)
                .set('length',band_length);
}).flatten();

//Linear interpolation
var monthcomposite_inter = ee.List.sequence(4,10,1).map(function(month){
  var month = ee.Number(month)
  var back_month = ee.ImageCollection.fromImages(monthcomposite.slice(month.subtract(2),month))
    .filter(ee.Filter.gt('length',0)).mosaic();
  var advance_month = ee.ImageCollection.fromImages(monthcomposite.slice(month.add(1),month.add(3)).reverse())
    .filter(ee.Filter.gt('length',0)).mosaic();
  var back_Y = back_month.select(bands);
  var back_doy = back_month.select('doy');
  var advance_Y = advance_month.select(bands);
  var advance_doy = advance_month.select('doy');
  var inter_img = ee.Image(monthcomposite.get(month));
  var inter_img_doy = ee.Image.constant(inter_img.get('doy')).float();;
  var Y = advance_Y.subtract(back_Y).divide(advance_doy.subtract(back_doy))
      .multiply(inter_img_doy.subtract(back_doy)).add(back_Y);
  var final_img = ee.Image(ee.Algorithms.If({
    condition : ee.Number(inter_img.get('length')).gt(0), 
    trueCase : inter_img.select(bands).unmask(Y),
    falseCase : Y
  }));
  return final_img.clip(region)
    .set('system:time_start',inter_img.get('system:time_start'),'doy',inter_img.get('doy'));
});
monthcomposite_inter = ee.ImageCollection.fromImages(monthcomposite_inter);

//WhittakerSmoothing
var WS = require('users/xuanhy12345678/L1:addon/WS');
var bandList = ['B6','B7','NDVI','EVI','LSWI','NDSVI','NDTI','GCVI'];
var WScollection =ee.ImageCollection([]);
for(var i=0;i<8;i++){
  var band = bandList[i];
  var bandCollection = monthcomposite_inter.select(band);
  var bandName = bandCollection.map(function(img){return img.unmask(bandCollection.mean())});
  var bandName = WS.whittakerSmoothing(bandName)[0];
  var newband = band+'_fitted';
  var collection = bandName.select(newband);
  var WScollection = WScollection.merge(collection);
}
var img_allbands = ee.Image(WScollection.iterate(mergeBands, ee.Image([])));

//Compute the percentile features
var newstart = ee.Date.fromYMD(year,4,1).getRelative('day', 'year');
var l8img_col_new = l8img_col.filter(ee.Filter.dayOfYear(newstart, newstart.add(60*3)));
var percentile = l8img_col_new.reduce(ee.Reducer.percentile([10,25,50,75,90]));
var range = l8img_col_new.reduce(ee.Reducer.percentile([75])).subtract(l8img_col_new.reduce(ee.Reducer.percentile([25])));
var features_combine = percentile.addBands(range);

//Compute the glcms features
var NDVI_h = l8img_col.filter(ee.Filter.dayOfYear(newstart.add(90), newstart.add(180))).select('NDVI').median().rename('NDVI_h');
var EVI_h = l8img_col.filter(ee.Filter.dayOfYear(newstart.add(90), newstart.add(180))).select('EVI').median().rename('EVI_h');
var glcms_NDVI = NDVI_h.multiply(10000).uint16().glcmTexture();
var glcms_EVI = EVI_h.multiply(10000).uint16().glcmTexture();
var glcms = glcms_NDVI.addBands(glcms_EVI);

//Harmonic regression
var l8 = l8img.map(addTimeBands);
var indexs = ee.List(['NDVI','EVI','LSWI']);
var component = ee.List(['constant', 'cos3', 'sin3', 'cos6' , 'sin6']);
var harmonic = l8.reduce(ee.Reducer.linearRegression(5,3));
var coefficients = harmonic.select('coefficients').matrixTranspose().arrayFlatten([indexs,component]).clip(region);  

//Final image composite and cropland mask
var finalimage = img_allbands.addBands(features_combine).addBands(glcms).addBands(coefficients).clip(region);
var finalImage = finalimage.updateMask(cropland);

var regionGrid = grid15.filterBounds(region);            
var gridList = regionGrid.toList(regionGrid.size());

//80% samples for training and 20% samples for testing
var sampleData = samples.filterBounds(region).randomColumn('random');
var training_samples_region = sampleData.filter(ee.Filter.lte("random", 0.8));
var testinf_samples_region = sampleData.filter(ee.Filter.gt("random", 0.8));

//Tile-based classification
var ls = [];
for(var i=0;i<gridList.size().getInfo();i++){
  //Compute neighboring Grids 
  var coordinates = ee.Feature(gridList.get(i)).geometry().centroid().coordinates();
  var x = ee.Number(coordinates.get(0));
  var y = ee.Number(coordinates.get(1));
  var roi = ee.FeatureCollection(ee.Geometry.Polygon(
        [[[x.subtract(1.49), y.add(1.49)],
          [x.add(1.49), y.add(1.49)],
          [x.add(1.49), y.subtract(1.49)],
          [x.subtract(1.49), y.subtract(1.49)],
          [x.subtract(1.49), y.add(1.49)]]]));
  var gridNear = regionGrid.filterBounds(roi).union();
  var clipImage = finalImage.clip(gridNear);
  var training_samples_gridNear = training_samples_region.filterBounds(gridNear);
  var training = clipImage.sampleRegions({
    collection: training_samples_gridNear, 
    properties: ["class"], 
    scale: 30,
    tileScale:16
  });
  //Classifier training
  var classifier = ee.Classifier.smileRandomForest(100)
    .train({
      features: training, 
      classProperty: 'class', 
      inputProperties: clipImage.bandNames()});
  var Classified_RF = finalImage.clip(ee.Feature(gridList.get(i))).classify(classifier);
  ls.push(ee.Image(Classified_RF));
};

var class_img = ee.ImageCollection(ls).mosaic();

//Export the classification maps
Export.image.toAsset({
  image: class_img,
  description: class_img+year+region_name,
  scale: 30,
  region: region,  
  maxPixels : 1e13
});

Export.table.toAsset({
  collection: testinf_samples_region,
  description: 'testinf_samples_region'+year+region_name,
  maxVertices:1e13
});
