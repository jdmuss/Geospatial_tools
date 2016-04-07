#!/usr/bin/env python
"""
rasterClass.py
Description:
   These classes and routines use GDAL to get file information ('header information')
   for a raster named <file_name>. They also load raster data into numpy arrays
   for manipulation and merging.

Dependencies: 
         Public:     datetime, gdal, gdalnumeric, numpy, osr, re, sys, __future__
         Private:    Py_utils
		    Raster_tools requires scikit-learn (install with 'pip install -U scikit-learn')
   
Created by: Jordan Muss
            Forest Landscape Ecology Lab
            Department of Forestry & Wildlife Ecology
            UW Madison
Creation Date:  7-13-2010
Updated:       12-20-2012   Added a clip method so that a raster can be clipped
                                using another raster as the clip mask.
		2-01-2013   Changed script name from "clipRaster.Py"
                                to "rasterClass.py"
                7-15-2013   Added a reprojection method
                3-04-2014   Added:
		              'write_raster' - procedure to write a raster to a file
			      'clipClass' - a class that sets up an environment in
			          order to reproject and clip one raster using
				  another raster.
			      'read_and_clip_raster' - procedure that reads a band
			          from an hdf file and clips it using a clipping
				  raster. TODO: make this more robust
			      'fill_band' - procedure that uses data and a mask from
				   another scene to fill in pixels that are flagged
				   by a source mask
			      'fill_all_bands' - procedure, that uses data from other
			           scenes to fill in pixels obscured by snow, cloud,
				   and/or shadow
			      'read_band' methods to the 'Raster' and 'Raster_hdf' classes.
				   These methods return the data array for the specified
				   band number
                7-24-2014   Updated: clipClass to calculate the number of cols correctly;
		                   changed the use of 'np.floor' to 'np.ceil'.
                8-19-2014   Added: 'ReprojectClass' to calculate the number of cols correctly;
			      'ReprojectClass' - a class that sets up the environment in
			          order to reproject a raster to a user specified projection.
			     Updated 'RasterClass', making the 'Reproject' method functional
			          using the new 'ReprojectClass'.
                9-15-2015   Added: 'GetSingleNoData' method to RasterClass to
		                 return a single no_data for the entire raster.
                3-27-2016   Added: 'proximity' function, which calculates a
		                 proximity matrix for a raster or band using
				 GDAL's ComputeProximity method. Distances are
				 in number of pixels. Multiply by pixel size to
				 get a real distance ***This could be folded
				 into RasterClass as method ***
                4-01-2016   Deleted: 'read_and_clip_raster' function, because
		                 it's functionality is accomplished using class
				 'clipClass'.
			     Deleted: 'clip_and_rescale_raster' method in class
			         'Raster_hdf'. It was redundant.
			     Updated: 'clip' & 'clip_raster' methods in 'RasterClass'
			         to utilize clip_class (removed redundant code).
			     Updated: 'clip_raster' method in 'Raster_hdf' to
			         utilize clip_class (removed redundant code).
			     Updated: 'clipClass' to return True/False on success/failure,
			         and added 'get_clipped_raster' method to return
				 the clipped raster as a numpy structure. Now allow a file
				 name to be passed for the raster that will be clipped.
                4-07-2016   Merged this version with one that storesreferences
				 dicts, etc. in gdal_defs.py

ToDo:
      1) replace calls to 'read_and_clip_raster' & 'clip_and_rescale_raster' in
             other scripts with calls to 'clipClass' or 'Raster_hdf'
      2) Combine/reconcile 'Reproject_Class' & 'clipClass' & calls by the 'Repreject'
             method in 'RasterClass'
      3) Look into 'LoadRaster' class & aggFunc (dict is located in gdal_defs)
"""
from __future__ import print_function
import sys
import datetime, re
import numpy as np
from PIL import Image, ImageDraw
from gdal_defs import *
from Raster_tools import mask_type
try:
    from osgeo import gdal, gdalnumeric, osr
except:
    import gdal, gdalnumeric, osr

gdal.UseExceptions()

def transformPoint(pt, from_cs, to_cs):
    sourceSR = osr.SpatialReference(wkt=from_cs)
    destSR = osr.SpatialReference(wkt=to_cs)
    transform = osr.CoordinateTransformation(sourceSR, destSR)
    return transform.TransformPoint(pt['X'], pt['Y'], pt['Z']) 

def fill_band(source_band, f_band, source_mask, fill_mask):
    ''' Fill in data that satisfy a masking criteria with data from another scene
        if the data in the other scene are valid (not masked by their own criteria). '''
    f_mask = source_mask & np.invert(fill_mask)
    remaining_mask = source_mask & fill_mask
    source_band[f_mask] = f_band[f_mask]
    return source_band

def GetEnd(sizeOffset):
    return (sizeOffset['size'] - sizeOffset['offset'])

''' Deprecated: def read_and_clip_raster(hdf_raster, band_num, clip_raster) '''

def write_raster(raster_bands, dest_geotransform, dest_projection, outRast, no_data=-9999, rasterFormat='GTiff'): 
    '''Create and write to a new multi-band raster. '''
    rastDriver = gdal.GetDriverByName(rasterFormat) 
    metadata = rastDriver.GetMetadata()
    if metadata.has_key(gdal.DCAP_CREATE) and metadata[gdal.DCAP_CREATE] == 'YES':
        '''Create the new multi-band raster: '''
        nBands = len(raster_bands)
        destRaster = rastDriver.Create(outRast, raster_bands[0].shape[1], \
                                       raster_bands[0].shape[0], nBands, \
                                       gdalDataTypes[str(raster_bands[0].dtype)]) 
        destRaster.SetGeoTransform(dest_geotransform) 
        destRaster.SetProjection(dest_projection) 
        '''Write each band to the new raster: '''
        for b in range(nBands):
            destBand = destRaster.GetRasterBand(b+1)
            destBand.WriteArray(raster_bands[b])
            if (no_data is not None) and (no_data is not np.nan):
		destBand.SetNoDataValue(no_data)
        destRaster = None
        ret_status = True
    else:
        print("The %s driver does not support a Create() method." % rasterFormat)
        ret_status = False  
    return ret_status

def proximity(rast):
    ''' Get the proximity array for a raster using GDAL's ComputeProximity
        method. Distances are in number of pixels. Multiply by pixel size to
	get a real distance: '''
    mem_driver = gdal.GetDriverByName('MEM')
    if isinstance(rast, RasterClass):
	inBand = rast.raster.GetRasterBand(1)
	geotransform = rast.geotransform
	nRow = rast.rows
	nCol = rast.cols
    elif isinstance(rast, np.ndarray):
	nRow, nCol = rast.shape
	d_Type = rast.dtype.name
	if d_Type in gdalDataTypes:
	    gdal_Type = gdalDataTypes[d_Type]
	    inRaster = mem_driver.Create('', nCol, nRow, 1, gdal_Type)
	    geotransform = (0, 1, 0, 0, 0, 1)
	    inRaster.SetGeoTransform(geotransform)
	    inBand = inRaster.GetRasterBand(1)
	    inBand.WriteArray(rast)
    else:
	print("An invalid data type was passed (raster class or numpy array is required).")
	return None
    outRaster = mem_driver.Create('', nCol, nRow, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(geotransform)
    dst = np.zeros((nRow, nCol))
    outBand = outRaster.GetRasterBand(1)
    outBand.WriteArray(dst)
    
    gdal.ComputeProximity(inBand, outBand)
    prox_result = outBand.ReadAsArray()
    del inBand, outRaster
    return prox_result

class RasterClass:
    def __init__(self, file_name, access=gdal.GA_ReadOnly, noDataVal=np.NaN):
        try:
            self.raster = gdal.Open( file_name, access )
	    self.rasterPath = file_name
	    self.isRaster = True
	    self.driver = self.raster.GetDriver().ShortName
	    self.numBands = self.raster.RasterCount
	    self.cols = self.raster.RasterXSize
	    self.rows = self.raster.RasterYSize
	    self.band_type = self.raster.GetRasterBand(1).DataType
	    self.dataType = gdal.GetDataTypeName(self.band_type)
	    self.projection = self.raster.GetProjection()
	    srs = osr.SpatialReference()
	    srs.ImportFromWkt(self.projection)
	    srs.AutoIdentifyEPSG()
	    self.EPSG = srs.GetAuthorityCode(None)
	    if self.EPSG is not None : self.EPSG = int(self.EPSG)
	    self.geotransform = self.raster.GetGeoTransform()
	    self.pixelSize = self.geotransform[1]
	    self.AltPixelSize = self.geotransform[5]
	    self.ulX = self.geotransform[0]
	    self.ulY = self.geotransform[3]
	    self.lrX = self.calc_lrX()
	    self.lrY = self.calc_lrY()
	    self.MD = self.raster.GetMetadata()
	    self.bands = {}
	    for band in self.MD:
		if band.partition('_')[2].isdigit():
		    self.bands[self.MD[band].partition('(')[2].partition(':')[0]] = \
			int(band.partition('_')[2])
    
	    self.NoData = []
	    for band_num in range(1, self.numBands+1):
		b = self.raster.GetRasterBand(band_num)
		bNoData = b.GetNoDataValue()
		if (bNoData is None) | (bNoData is np.nan):
		    self.NoData.append(noDataVal)
		else:
		    self.NoData.append(bNoData)
		b = None
	    if len(self.NoData) > 0:
		self.nodata_value = self.NoData[0]
	    else:
		self.nodata_value = None
	    CT = self.raster.GetRasterBand(1).GetRasterColorTable()
	    if CT is not None:
		self.CT = CT.Clone()
	    else:
		self.CT = None
        except:
            self.EmptyInfoFile()
            self.isRaster = False
	    print(("WARNING: either %s does not exist or there was a problem " \
	          + "opening/reading the raster.") % file_name)
            return
    def EmptyInfoFile(self):
        self.raster = None
        self.isRaster = True
        self.rasterPath = None
        self.nodata_value = None
        self.driver = None
        self.numBands = None
        self.cols = None
        self.rows = None
        self.band_type = None
        self.dataType = None
        self.projection = None
        self.geotransform = None
        self.pixelSize = None
        self.AltPixelSize = None
        self.ulX = None
        self.ulY = None
        self.lrX = None
        self.lrY = None
        self.CT = None
    def copy(self):
        copy = RasterClass('', 0)
        for key in self.__dict__.keys():
            exec('copy.' + key + '= self.' + key)
        return copy
    def calc_lrX(self):
        return self.ulX + self.pixelSize * self.cols
    def calc_lrY(self):
        return self.ulY + self.AltPixelSize * self.rows
    def clip(self, clipper, outRasterName=None):
	''' Clip a raster using the extent of a shapefile or another raster.
	    This assumes that clipping shapefile/raster & original raster have
	    identical projections: '''
        pixel = {'width':self.geotransform[1], 'height':self.geotransform[5]}
        origin = {'x':self.ulX, 'y':self.ulY}
        cols = self.cols
        rows = self.rows
	if type(clipper) is type(()):
	    ''' A tuple of (ulX, lrX, lrY, ulY) format was passed'''
	    intersect = self.GetShapeIntersect(clipper)
	    if intersect != None:
		clipULoffset = {'x':int(np.floor((intersect['ulX'] - origin['x']) / pixel['width'])),
		                'y':int(np.floor((intersect['ulY'] - origin['y']) / pixel['height']))}
		clipLRoffset = {'x':int(np.floor((intersect['lrX'] - origin['x']) / pixel['width'])),
		               'y':int(np.floor((intersect['lrY'] - origin['y']) / pixel['height']))}
	    new_geotransform = (intersect['ulX'], self.geotransform[1], self.geotransform[2], \
		                intersect['ulY'], self.geotransform[4], self.geotransform[5])
	    new_height = clipLRoffset['y'] - clipULoffset['y'] + 1
	    new_width = clipLRoffset['x'] - clipULoffset['x'] + 1
	    clipped_rast = np.full(shape=(self.numBands, new_height, new_width), \
		            fill_value=self.nodata_value, dtype=gdalNumpyTypes[self.dataType])
	    
	    xStart = clipULoffset['x'] + 1
	    yStart = clipULoffset['y'] + 1
	    ''' Get each band, clip it, and load it into the out raster: '''
	    for band in range(self.numBands):
		srcBand = self.read_band(band+1)
		clipped_rast[band] = srcBand[clipULoffset['y']:(clipLRoffset['y']+1), \
		                        clipULoffset['x']:(clipLRoffset['x']+1)]
	else:
	    if isinstance(clipper, str) and os.path.isfile(clipper):
		clipper = RasterClass(clipper)
	    if isinstance(clipper, RasterClass):
		'''clipper is an instance of RasterClass'''
		clipped_rast = self.clip_raster(clipper)
		xStart = clipper.ulX
		yStart = clipper.ulY
	    else:
		print("Warning (RasterClass.clip): the clipper provided was not valid (%s)" % self.rasterPath)
		return None
	    if clipped_rast is not None:
		if outRasterName is not None:
		    if len(clipped_rast.shape) == 2 : clipped_rast = [clipped_rast]
		    write_raster(clipped_rast, new_geotransform, self.projection, \
			         outRasterName, no_data=self.nodata_value, \
		                 rasterFormat=self.driver)
		return {'xStart':xStart, 'yStart':yStart, 'rast':clipped_rast}
	    else:
		return None
    def clip_raster(self, clip_rast, trans='NN'):
	''' Reproject, and clip a raster using a template raster: '''
	clipRaster = clipClass(clip_rast, self.pixelSize, num_bands=self.numBands, rast_type=self.dataType)
	clip = clipRaster.clip_raster(self.raster, trans=trans)
	if not clip:
	    ''' There was an error clipping/reprojecting the raster. The error
	        message was printed in 'clipcClass':'''
	    return None
	else:
	    return clip.get_clipped_raster()
    def read_band(self, band_num):
	return self.raster.GetRasterBand(band_num).ReadAsArray()
    def read_all_bands(self):
	bands = np.empty([self.numBands, self.rows, self.cols], dtype=gdalNumpyTypes[gdal.GetDataTypeName(self.band_type)])
	for b in range(1, (self.numBands+1) ):
	    if len(self.NoData) == 1:
		bands[b-1].fill(self.NoData)
	    else:
		bands[b-1].fill(self.NoData[b-1])
	    bands[b-1] = self.raster.GetRasterBand(b).ReadAsArray()
	return bands
    def Reproject(self, newProj, transformation):
	''' Reproject a raster using the specified transformations: '''
	reproj = ReprojectClass(self, newProj, rast_type = 'MEM', trans = transformation)
	r = reproj.reproject(self)
	return r
    def Rescale(self, newPixelSize):
	'''ToDo: this looks like it needs fixing. repeat/aggregate pixel values?'''
        # Reset cols, rows, & lower right coordinates for a new raster resolution
        self.pixelSize = newPixelSize
        newDims = self.GetAltDims(newPixelSize)
        self.cols = newDims['cols']
        self.rows = newDims['rows']
        self.lrX = self.calc_lrX()
        self.lrY = self.calc_lrY()
    def GetIntersect(self, Rast2_Info):
        # This assumes pixels are square
        int_ulX = max(self.ulX, Rast2_Info.ulX)
        int_ulY = min(self.ulY, Rast2_Info.ulY)
        int_lrX = min(self.lrX, Rast2_Info.lrX)
        int_lrY = max(self.lrY, Rast2_Info.lrY)
        if (int_ulX > self.lrX) | (int_lrY < self.lrY):
            return None
        return {'ulX':int_ulX, 'ulY':int_ulY, 'lrX':int_lrX, 'lrY':int_lrY}
    def GetShapeIntersect(self, shape_extent):
	''' Get the intersection of the raster with the extent of a shapefile.
	    This assumes that the raster and shapefile have identical projections:'''
	int_ulX = max(self.ulX, shape_extent[0])
	int_ulY = min(self.ulY, shape_extent[3])
	int_lrX = min(self.lrX, shape_extent[1])
	int_lrY = max(self.lrY, shape_extent[2])
	if (int_ulX > self.lrX) | (int_lrY < self.lrY):
	    print("Warning: The shapefile and raster do not intersect")
	    return None
        return {'ulX':int_ulX, 'ulY':int_ulY, 'lrX':int_lrX, 'lrY':int_lrY}
    def GetNewDims(self, winDims):
        int_xSize = int((winDims['lrX'] - winDims['ulX'])/self.pixelSize)
        int_ySize = int((winDims['ulY'] - winDims['lrY'])/self.pixelSize)
        return {'cols':int_xSize, 'rows':int_ySize}
    def GetAltDims(self, pixelSize):
        alt_xSize = int(round((self.lrX - self.ulX)/pixelSize))
        alt_ySize = int(round((self.ulY - self.lrY)/pixelSize))
        return {'cols':alt_xSize, 'rows':alt_ySize}
    def GetSingleNoData(self, default_no_data=-9999):
	no_data = np.unique(self.NoData)
	if len(no_data) > 0:
	    if len(no_data) != 1:
		print("Warning: %s has more than one no_data value! The first value (%s) will be used" % (self.rasterPath, no_data[0]), file=sys.stderr)
	    no_data = no_data[0]
	    if no_data != default_no_data:
		print("Error: the no_data value for %s does not match the default no_data value (%s)!" % (self.rasterPath, default_no_data), file=sys.stderr)
		return(False, -2, no_data)
	    else:
		return(True, no_data)
	else:
	    print("Warning: No no_data values were found for %s! The default value (%s) will be used" % (self.rasterPath, default_no_data), file=sys.stderr)
	    return(False, -1)
    def GetOffsets(self, point):
	''' Return the x & y indices (offsets) for a point in the raster's
	    coordinate system: '''
        x = int((point['X'] - self.ulX)/self.pixelSize)
        y = int((point['Y'] - self.ulY)/self.AltPixelSize)
        return {'X':x, 'Y':y}
    def GetCoordinate(self, point):
	''' Return the coordinates (in the raster's coordinate system) of a
	    point's offset (the array's x & y indices): '''
	proj_Y = self.ulY + self.AltPixelSize*(2*point['Y'] + 1)/2
	proj_X = self.ulX + self.pixelSize*(2*point['X'] + 1)/2
        return {'X':proj_X, 'Y':proj_Y}
    def Update_Band(self, band_num, datta_array):
	band = self.raster.GetRasterBand(band_num)
        result = band.WriteArray(datta_array)
	band.FlushCache()
        return result
    def __str__(self):
        return ('Filename: %s\n\
Format: %s\n\
File Size: %d x %d x %d\n\
Pixel Size: (%f x %f)\n\
UL:(%f, %f)\n\
LR:(%f, %f)' % (self.rasterPath, self.driver, self.cols, self.rows, self.numBands, \
                self.pixelSize, self.AltPixelSize, self.ulX, self.ulY, \
                self.lrX, self.lrY))
    def report(self):
        print('Filename: %s' % self.rasterPath)
        print('Driver: %s' % self.driver)
        print('File Size: %d x %d x %d' % (self.cols, self.rows, self.numBands))
        print('Pixel Size: (%f x %f)' % (self.pixelSize, self.AltPixelSize))
        print('UL:(%f, %f)' % (self.ulX, self.ulY))
        print('LR:(%f, %f)' % (self.lrX, self.lrY))
    def close(self):
        self.raster = None
    def __del__(self):
        self.raster = None

class Raster_hdf:
    def __init__(self, file_name):
	self.hdf_raster = gdal.Open(file_name)
	''' Each band has two entries [name, description] '''
	self.bands = self.hdf_raster.GetSubDatasets()
	self.num_bands = len(self.bands)
	self.metadata = self.hdf_raster.GetMetadata()
	if 'AcquisitionDate' in self.metadata:
	    self.dt_time = datetime.datetime.strptime(self.metadata['AcquisitionDate'], '%Y-%m-%dT%H:%M:%S.%fZ')
	    self.long_date = self.dt_time.strftime('%Y-%m-%d')
	    self.date = self.dt_time.strftime('%Y%j')
	    self.year = self.dt_time.strftime('%Y')
	    self.doy = self.dt_time.strftime('%j')
	else:
	    self.long_date = '' # No Acquisition date information available
    def clip_raster(self, clip_rast, band_num, trans='NN', pix_size=None):
	''' Reproject, and clip a raster using a template raster: '''
	band = RasterClass(self.bands[band_num][0], gdal.GA_ReadOnly)
	if pix_size is not None:
	    self.pixelSize = pix_size
	else:
	    self.pixelSize = band.pixelSize
	clipRaster = clipClass(clip_rast, band.pixelSize, num_bands=1, rast_type=band.dataType)
	clip = clipRaster.clip_raster(band.raster, trans=trans)
	band.close()
	if not clip:
	    ''' There was an error clipping/reprojecting the raster:'''
	    m = re.search(r'(?<=band)\d+', self.bands[band_num][1], flags=re.IGNORECASE)
	    if m:
		band = m.group()
	    else:
		band = self.bands[band_num][1].split(' ')[1]
	    print("Error reprojecting Band %s" % band)
	    return None
	else:
	    return clip.get_clipped_raster()
    '''Deprecated: def clip_and_rescale_raster(self, clip_rast, band_num, pix_size, trans='NN'):'''
    def print_metadata(self):
	for k, v in self.metadata.iteritems(): print("%s:\t%s" % (k, v))
    def print_band_names(self):
	for v in self.bands: print(v[0])
    def print_band_description(self):
	for v in self.bands: print(v[1])
    def close(self):
        self.hdf_raster = None
    def __del__(self):
        self.hdf_raster = None

class Band:
    def __init__(self, fileInfo):
        self.Band = np.zeros((fileInfo.ySize,fileInfo.xSize),
                             dtype = gdalNumpyTypes[fileInfo.dataType])
        self.FillBand(fileInfo.nodata_value)
    def FillBand(self, fillValue):
        self.Band.fill(fillValue)
    def LoadData(self, file_name):
        try:
            self.Band = gdalnumeric.LoadFile(file_name)
        except RuntimeError:
            return False
        return True
    def FillGreater(self, ceiling, fillValue):
        self.Band = np.choose(np.greater(self.Band, ceiling), (self.Band, fillValue))
    def CopyFrom(self, fromBand):
        """ This is a simple array copy. It first tests that the arrays have the
            same dimensions.
        """
        if self.Bandshape == fromBand.shape:
            self.Band =  fromBand.copy()
            return True
        else:
            print("Copy Failed. The Band dimensions do not match.", file=sys.stderr)
            return False

class LoadRaster:
    def __init__(self, file_name, ndVal):
        if file_name == None:
            self.Info = None
            self.Rast = None
        else:
            self.Info = RasterClass(file_name, ndVal)
            self.Rast = Band(self.Info)
            self.Rast.LoadData(self.Info.rasterPath)
    def FillGreater(self, ceiling, fillValue):
        self.Rast.Band = np.choose(np.greater(self.Rast.Band, ceiling), \
                                   (self.Rast.Band, fillValue))
    def CopyFrom(self, raster):
        intersect = self.Info.GetIntersect(raster.Info)
        intersect = {'X':intersect['ulX'], 'Y':intersect['ulY']}
        toOffset = self.Info.GetOffsets(intersect)
        fromOffset = raster.Info.GetOffsets(intersect)
        for row in range(min(GetEnd({'size':self.Info.rows, 'offset':toOffset['Y']}),
                        GetEnd({'size':raster.Info.rows, 'offset':fromOffset['Y']}))):
            for col in range(min(GetEnd({'size':self.Info.cols, 'offset':toOffset['X']}),
                            GetEnd({'size':raster.Info.cols, 'offset':fromOffset['X']}))):
                self.Rast.Band[row + toOffset['Y'], col+toOffset['X']] = \
                    raster.Rast.Band[row + fromOffset['Y'], col + fromOffset['X']]
    def Align(self, toRast):
        Aligned = LoadRaster(None, 0)
        Aligned.Info = toRast.Info.copy()
        Aligned.Info.rasterPath = None
        Aligned.Info.Rescale(self.Info.pixelSize)
        Aligned.Rast = Band(Aligned.Info)
        Aligned.CopyFrom(self)
        return Aligned
    def Rescale(self, newPixelSize, NDval, Aggregation = 'mean'):
        CF = int(newPixelSize/self.Info.pixelSize)
        newRast = LoadRaster(None, NDval)
        newRast.Info = self.Info.copy()
        if Aggregation == 'mean':
            newRast.Info.dataType = 'Float32'
        newRast.Info.Rescale(newPixelSize)
        newRast.Rast = Band(newRast.Info)
        for row in range(newRast.Info.rows):
            for col in range(newRast.Info.cols):
                tmp=np.array(self.Rast.Band[row*CF:(row+1)*CF, col*CF:(col+1)*CF])
                newRast.Rast.Band[row,col] = eval(aggFunc['mean'] + \
                                    '( tmp[~np.equal(tmp, NDval)] )')
        return newRast

class ReprojectClass:
    def __init__(self, rasterCl, new_projection, rast_type = 'MEM', trans = 'cubicspline'):
	self.rast_driver = gdal.GetDriverByName( rast_type )
	self.from_projection = rasterCl.projection
	self.new_projection = new_projection
	self.trans = trans

        self.ulX, self.ulY, self.ulZ = transformPoint({'X':rasterCl.ulX, \
	        'Y':rasterCl.ulY, 'Z':0}, self.from_projection, self.new_projection)
	self.lrX, self.lrY, self.lrZ = transformPoint({'X':rasterCl.lrX, \
	        'Y':rasterCl.lrY, 'Z':0}, self.from_projection, self.new_projection)
	
	if (osr.SpatialReference(self.from_projection).GetAttrValue("UNIT") == 'degree'):
		spheroid_radius = float(osr.SpatialReference(self.from_projection).GetAttrValue("SPHEROID",1))
		self.npx = rasterCl.pixelSize/(np.radians(1)*spheroid_radius)
	else:
		self.npx = rasterCl.pixelSize

	self.nCols = int(np.ceil(np.round((self.lrX - self.ulX)/self.npx, 1)))
	self.nRows = int(np.ceil(np.round((self.ulY - self.lrY)/self.npx, 1)))
	self.lrX = self.ulX + self.npx * self.nCols
	self.lrY = self.ulY - self.npx * self.nRows
	''' Calculate the new geotransform '''
	self.new_geotransform = ( self.ulX, self.npx, rasterCl.geotransform[2], self.ulY, rasterCl.geotransform[4], -self.npx )
    def reproject(self, rasterCl, out_file = '', num_bands=1):
	from_rast = rasterCl.raster
	
	''' Create a raster (default=in memory) with the pixel-size & the boundary set for the new projection: '''
	self.dest = self.rast_driver.Create(out_file, self.nCols, self.nRows, num_bands, rasterCl.band_type)
	''' Set the geotransform '''
	self.dest.SetGeoTransform(self.new_geotransform)
	self.dest.SetProjection(self.new_projection)
	res = gdal.ReprojectImage( from_rast, self.dest, rasterCl.projection, self.new_projection, transforms[self.trans] )
	if (res != 0):
	    print("Error reprojecting raster %s" % rasterCl.rasterPath)
	    sys.exit(2)
	else:
	    return self.dest

class clipClass:
    def __init__(self, clip_raster_name, nom_pix_sz, num_bands=1, rast_type='Byte'):
	memory = gdal.GetDriverByName( 'MEM' )
	clipRaster = RasterClass(clip_raster_name, gdal.GA_ReadOnly)
	self.ulX, self.ulY, ulZ = clipRaster.ulX, clipRaster.ulY, 0
	lrX, lrY, lrZ = clipRaster.lrX, clipRaster.lrY, 0

	if (osr.SpatialReference(clipRaster.projection).GetAttrValue("UNIT") == 'degree'):
	    spheroid_radius = float(osr.SpatialReference(clipRaster.projection).GetAttrValue("SPHEROID",1))
	    self.npx = nom_pix_sz/(np.radians(1)*spheroid_radius)
	else:
	    self.npx = nom_pix_sz

	self.nCols = int(np.ceil(np.round((lrX-self.ulX)/self.npx, 1)))
	self.nRows = int(np.ceil(np.round((self.ulY-lrY)/self.npx, 1)))
	self.lrX = self.ulX + self.npx*self.nCols
	self.lrY = self.ulY - self.npx*self.nRows

	''' Calculate the new geotransform '''
	self.geotransform = ( self.ulX, self.npx, clipRaster.geotransform[2], \
	                      self.ulY, clipRaster.geotransform[4], -self.npx )
	self.projection = clipRaster.projection
	''' Create a raster in memory with new pixel-size & the boundary of the clip raster: '''
	self.dest = memory.Create('', self.nCols, self.nRows, num_bands, gdalDataTypes[rast_type])
	''' Set the geotransform '''
	self.dest.SetGeoTransform(self.geotransform)
	self.dest.SetProjection(self.projection)
	clipRaster = None
    def clip_raster(self, rast_to_clip, trans='NN'):
	''' Reproject, and clip a raster using a template raster: '''
	#clipRaster = clipClass(clip_rast, self.pixelSize, num_bands=self.numBands, rast_type=self.dataType)
	if isinstance(rast_to_clip, str) and os.path.isfile(rast_to_clip):
	    rast_to_clip = RasterClass(rast_to_clip)
	if isinstance(rast_to_clip, RasterClass):
	    inRast = rast_to_clip
	    rastPath = inRast.rasterPath
	    res = gdal.ReprojectImage( inRast.raster, self.dest, inRast.projection, \
		                             self.projection, transforms[trans] )
	    if res != 0:
		print("Error clipping/reprojecting raster %s" % rastPathh)
		return False
	    else:
		return True
	else:
	    print("Warning: ", rast_to_clip, " is an invalid raster/raster file name and can not be clipped")
	    return False
    def get_clipped_raster(self):
	return self.dest.ReadAsArray()

class rasterizeClass:
    def __init__(self, to_raster_class):
	''' Initialize the rasterization class. Set the to_raster class and
	    projection. Set all other fields to None: '''
	self.base_raster = to_raster_class
	self.projection = to_raster_class.projection
	self.geotransform = None
	self.envelope = None
	self.pixel_width = to_raster_class.geotransform[1]
	self.pixel_height = to_raster_class.geotransform[5]
	self.mask = None
	self.mask_shape = None
    def rasterize_feature(self, feature, to_raster_name, write=True):
	''' Rasterize the selected feature into the projection and geotransform of
	    the initialized raster class. Do not create a GeoTiff if write=False:'''
	cs = feature.GetGeometryRef().GetSpatialReference().ExportToWkt()
	''' Map points to pixels for drawing the boundary on a blank 8-bit,
	    black and white image mask:                                   '''
	points = []
	pixels = []
	''' Get the feature envelope from the nodes, because "Envelope" does 
	    not seem to capture all of the points:                         '''
	ul_x, lr_x, lr_y, ul_y = (None, None, None, None)
	
	''' Get the coordinates of each node of the feature and convert them
	    in to the coordinate system of the reference (to-) raster:     '''
	g = feature.GetGeometryRef()
	pts = g.GetGeometryRef(0)
	for p in range(pts.GetPointCount()):
	    pt = (pts.GetX(p), pts.GetY(p))
	    (x, y, z) = transformPoint({'X':pt[0], 'Y':pt[1], 'Z':0}, cs, self.projection)
	    if ul_x is None:
		ul_x, lr_x, lr_y, ul_y = x, x, y, y
	    else:
		if x < ul_x : ul_x = x
		if x > lr_x : lr_x = x
		if y < lr_y : lr_y = y
		if y > ul_y : ul_y = y
	    points.append((x, y))
	''' Create the geotransform of the new raster mask for the rasterized
	    feature, determine its shape, & get the relative indices for the 
	    nodes of the feature:                                          '''	
	self.geotransform = (ul_x, self.pixel_width, 0, ul_y, 0, self.pixel_height)
	UL_idx = self.base_raster.GetOffsets({'X':ul_x, 'Y':ul_y})
	LR_idx = self.base_raster.GetOffsets({'X':lr_x, 'Y':lr_y})
	self.envelope = {'ul':UL_idx, 'lr':LR_idx}
	nRows = LR_idx['Y'] - UL_idx['Y']
	nCols = LR_idx['X'] - UL_idx['X']
	self.mask_shape = (nRows, nCols)

	for p in points:
	    pt = self.base_raster.GetOffsets({'X':p[0], 'Y':p[1]})
	    pixels.append((pt['X']-UL_idx['X'], pt['Y']-UL_idx['Y']))
	''' Create a new 8-bit black & white PIL image (mode='L') that is the
	    shape of the feature's 'envelope', with color=0. Load the feature's
	    nodes into the image, and determine the polygon. Save the new
	    polygon feature as an array mask (1 = polygon, 0 = empty space:'''
	rasterPoly = Image.new("L", (nCols, nRows), 0)
	ImageDraw.Draw(rasterPoly).polygon(pixels, outline=1, fill=1)
	self.mask = np.array(rasterPoly)
	if write:
	    write_raster([self.mask], self.geotransform, self.projection, outRast=to_raster_name, no_data=0)
    def reset(self):
	''' Reset the mask variables of the class instance:                '''
	self.geotransform = None
	self.envelope = None
	self.mask = None
	self.mask_shape = None
'''------------------------------------------------------------------------------------
   The following procedures are specific to landsat hdf files with a cfMask/fmask band:
   ------------------------------------------------------------------------------------'''
def fill_all_bands(source_file, fill_files, output_raster, out_type = 'GTiff', 
                   fill_val = 1, fill_crit = 'gt', no_data = -9999):
    ''' Fill in data that satisfy a masking criteria with data from another scene
        if the data in the other scene are valid (not masked by their own criteria).
        This expects Landsat hdf files with fmask bands as the source and fill files. '''
    source_raster = Raster_hdf(source_file)
    ''' Get the source bands; blue (band 1) = 0, green  (band 2) = 1, red (band 3) = 2,
        NIR (band 4) = 3, SWIR (band 5) = 4, SWIR2 (band 7) = 5, & cfmask (cloud,
        shadow, snow, & water) = 16. : '''
    bands = {0:'', 1:'', 2:'', 3:'', 4:'', 5:''}
    b_idx = bands.keys()
    b_idx.sort()
    for b in b_idx:
        band = RasterClass(source_raster.bands[b][0], gdal.GA_ReadOnly)
        bands[b] = band.read_band(1)
        band.close()
    ''' Get the source cloud mask: '''
    source_mask_band = RasterClass(source_raster.bands[16][0], gdal.GA_ReadOnly).read_band(1)
    source_mask = mask_type(fill_crit, source_mask_band, [fill_val]).mask
    del(source_mask_band)
    ''' Set the clip raster: '''
    band = RasterClass(source_raster.bands[0][0], gdal.GA_ReadOnly)
    clipRaster = clipClass(source_raster.bands[0][0], band.pixelSize, num_bands=1, rast_type=band.dataType)
    cloudClipRaster = clipClass(source_raster.bands[16][0], band.pixelSize, num_bands=1, rast_type='Byte')
    band.close()
    for ff in fill_files:
        fill_raster = Raster_hdf(ff)
        ''' Get the cloud mask of the infill raster (band 16): '''
        fill_mask_band = read_and_clip_raster(fill_raster, 16, cloudClipRaster)
        fill_mask = mask_type(fill_crit, fill_mask_band, [fill_val]).mask
        del(fill_mask_band)
        f_mask = source_mask & np.invert(fill_mask)
        for b in b_idx:
            infill_band = read_and_clip_raster(fill_raster, b, clipRaster)
            bands[b][f_mask] = infill_band[f_mask]
        ''' Determine which pixels remain to be filled/replaced: '''
        source_mask = source_mask & fill_mask
    status = write_raster([bands[b] for b in b_idx], clipRaster.geotransform, \
                          clipRaster.projection, output_raster, no_data = no_data, \
                          rasterFormat = out_type)
    return status
