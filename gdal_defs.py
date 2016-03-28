#!/usr/bin/env python
"""
gdal_defs.py
Description:
   These are common dicts, variables, etc. used by my gdal helper libraries.
   They have been pulled out of the rasterClass and shapeFileClass libraries
   in order to reduce errors & eliminate duplication.

Dependencies: 
         Public:     datetime, gdal, gdalnumeric, numpy, osr, re, sys, __future__
         Private:    
   
Created by: Jordan Muss
            
Creation Date:  3-25-2016
Updated:       xx-xx-20xx   
"""
import numpy as np
from Raster_tools import mask_type
try:
    from osgeo import gdal, gdalnumeric, osr
except:
    import gdal, gdalnumeric, osr
    
'''--------------------------------------------------------------------------
      Define variables/dicts.
      Create numpy/gdal type mapping table: 
   --------------------------------------------------------------------------'''
if gdal.__version__[:gdal.__version__.rfind('.')] == "1.9":
    transforms = { 'NearestNeighbor'  : gdal.GRA_NearestNeighbour,
                   'near'             : gdal.GRA_NearestNeighbour,
                   'NN'               : gdal.GRA_NearestNeighbour,
                   'bilinear'         : gdal.GRA_Bilinear,
                   'BL'               : gdal.GRA_Bilinear,
                   'cubic'            : gdal.GRA_Cubic,
                   'C'                : gdal.GRA_Cubic,
                   'cubicspline'      : gdal.GRA_CubicSpline,
                   'CS'               : gdal.GRA_CubicSpline,
                   'L'                : gdal.GRA_Lanczos,
                   'lanczos'          : gdal.GRA_Lanczos}
else:
    transforms = { 'NearestNeighbor'  : gdal.GRA_NearestNeighbour,
                   'near'             : gdal.GRA_NearestNeighbour,
                   'NN'               : gdal.GRA_NearestNeighbour,
                   'bilinear'         : gdal.GRA_Bilinear,
                   'BL'               : gdal.GRA_Bilinear,
                   'cubic'            : gdal.GRA_Cubic,
                   'C'                : gdal.GRA_Cubic,
                   'cubicspline'      : gdal.GRA_CubicSpline,
                   'CS'               : gdal.GRA_CubicSpline,
                   'L'                : gdal.GRA_Lanczos,
                   'lanczos'          : gdal.GRA_Lanczos,
                   'average'          : gdal.GRA_Average,
                   'avg'              : gdal.GRA_Average,
                   'mode'             : gdal.GRA_Mode,
                   'M'                : gdal.GRA_Mode}
# also check 'average' and 'mode' resampling

gdalDataTypes = { 'Unknown'  : gdal.GDT_Unknown,
                  'Byte'    : gdal.GDT_Byte,
                  'byte'    : gdal.GDT_Byte,
                  'Int8'    : gdal.GDT_Byte,
                  'int8'    : gdal.GDT_Byte,
                  'UInt16'   : gdal.GDT_UInt16,
                  'uint16'   : gdal.GDT_UInt16,
                  'Int16'    : gdal.GDT_Int16,
                  'int16'    : gdal.GDT_Int16,
                  'UInt32'   : gdal.GDT_UInt32,
                  'uint32'   : gdal.GDT_UInt32,
                  'Int32'    : gdal.GDT_Int32,
                  'int32'    : gdal.GDT_Int32,
                  'Float32'  : gdal.GDT_Float32,
                  'float32'  : gdal.GDT_Float32,
                  'Float64'  : gdal.GDT_Float64,
                  'float64'  : gdal.GDT_Float64,
                  'CInt16'   : gdal.GDT_CInt16,
                  'cint16'   : gdal.GDT_CInt16,
                  'CInt32'   : gdal.GDT_CInt32,
                  'cint32'   : gdal.GDT_CInt32,
                  'CFloat32' : gdal.GDT_CFloat32,
                  'cfloat32' : gdal.GDT_CFloat32,
                  'CFloat64' : gdal.GDT_CFloat64,
                  'cfloat64' : gdal.GDT_CFloat64}
  
gdalNumpyTypes = {'Unknown'  : 'uint8',
                  'Byte'     : 'int8',
                  'UInt16'   : 'uint16',
                  'Int16'    : 'int16',
                  'UInt32'   : 'uint32',
                  'Int32'    : 'int32',
                  'Float32'  : 'float32',
                  'Float64'  : 'float64',
                  'CInt16'   : np.dtype('int8, int8'),
                  'CInt32'   : np.dtype('int16, int16'),
                  'CFloat32' : 'complex',
                  'CFloat64' : 'complex'}

aggFunc = {'mean': 'np.mean',
           'max' : 'np.max',
           'min' : 'np.min'}
'''--------------------------------------------------------------------------
      Create special projections:
   --------------------------------------------------------------------------'''
'''------------------MODIS sinusoidal projection (6974):---------------------'''
modis_srs = osr.SpatialReference()
modis_srs.SetWellKnownGeogCS( "WGS84" )
modis_srs.SetProjCS( "MODIS Sinusoidal" )
modis_srs.SetWellKnownGeogCS( "EPSG:4326" )
modis_srs.SetProjection("Sinusoidal")
modis_srs.SetProjParm("false_easting", 0.0)
modis_srs.SetProjParm("false_northing", 0.0)
modis_srs.SetProjParm("central_meridian", 0.0)
modis_srs.SetProjParm("Semi_Major", 6371007.181)
modis_srs.SetProjParm("Semi_Minor", 6371007.181)
modis_srs.SetLinearUnits("Meter",1.0)
