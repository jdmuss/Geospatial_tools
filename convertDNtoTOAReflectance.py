#!/usr/bin/env python
"""
convertDNtoTOAReflectance.py

Description:
   This python program uses meta file information stored in a class (provided by
   "readLandsatMTL.py") to convert DN values to top of atmosphere (TOA) reflectance
   values. It currently only works with Landsat 7 data. Data are processed and 
   written to the raster on the fly due to memory constraints within Python.

Dependencies: 
   Dependencies:
         Public:     datetime, gdal, gdalnumeric, numpy
         Private:    readLandsatMTL, rasterHeader, rasterClass
    
Created by: Jordan Muss

Creation Date: 2-1-2013
Updated:       ??-??-201?  

ToDo:
          Modify the write method to reconstruct the metafile
"""
import datetime
from numpy import cos, pi, radians, sin  # can use math instead
import numpy as np

try:
    from osgeo import gdal, gdalnumeric
except:
    import gdal, gdalnumeric

from readLandsatMTL import *
from rasterHeader import *
from rasterClass import *

esun = {
    1 : 1996.6, # 1997
    2 : 1812, # 1812
    3 : 1533, # 1533
    4 : 1039, # 1039
    5 : 230.8, # 230.8
    7 : 84.90, # 84.90
    8 : 1362 # 1362
}

L7bandInfo = {
    1  : ['Band_1','Band 1 (Blue: 450-520 nm)'],
    2  : ['Band_2','Band 2 (Green: 520-600 nm)'],
    3  : ['Band_3','Band 3 (Red: 630-690 nm)'],
    4  : ['Band_4','Band 4 (NIR: 770-900 nm)'],
    5  : ['Band_5','Band 5 (SWIR1: 1550-1750 nm)'],
    6  : ['Band_6','Band 6 (Thermal : 10400-125000 nm)'],
    61 : ['Band_61','Band 6L (Thermal low gain: 10400-125000 nm)'],
    62 : ['Band_62','Band 6H (Thermal low gain: 10400-125000 nm)'],
    7  : ['Band_7','Band 7 (SWIR2: 2060-2350 nm)'],
    8  : ['Band_8','Band 8 (Panchromatic: 520-900 nm)}']
}

def stripNonAlpha(iStr):
    return ''.join(ch for ch in iStr if ch.isalnum())

def getDOY(dateTuple):
    return dateTuple.timetuple().tm_yday

def earthSunDist(DOY, alt=False):
    if alt:
        d = (1+0.01672*sin(2*pi*(DOY-93.5)/365))
    else:
        d = (1-0.01672*cos(radians(0.9856*(DOY-4))))
    return d

def getTOArefl(DN, radMax, radMin, dnMax, dnMin, d, esun, cosZenAngle):
    radiance = ((radMax - radMin)/(dnMax - dnMin))*(DN - dnMin) + radMin
    return (pi*radiance*d*d)/(esun*cosZenAngle)

def DNtoTOAReflectance(outFileName, metaFileName, controlFileName,
                       bandsToProc, NoData = 0, outRasterFormat = 'ENVI'):
    meta=L7metaFileClass(metaFileName, controlFileName)
    sensor = stripNonAlpha(meta.spacecraftID).lower()
    
    if sensor == 'landsat7':
        bandFiles = {
            1 : meta.band1FileName, 2 : meta.band2FileName, 3 : meta.band3FileName,
            4 : meta.band4FileName, 5 : meta.band5FileName, 61 : meta.band61FileName,
            62 : meta.band62FileName, 7 : meta.band7FileName, 8 : meta.band8FileName
        }
    
        bands = {
            1 : {'qCalMin':meta.quantizeCalMinBand1, 'qCalMax':meta.quantizeCalMaxBand1,
                 'radMin':meta.radianceMinBand1, 'radMax':meta.radianceMaxBand1},
            2 : {'qCalMin':meta.quantizeCalMinBand2, 'qCalMax':meta.quantizeCalMaxBand2,
                 'radMin':meta.radianceMinBand2, 'radMax':meta.radianceMaxBand2},
            3 : {'qCalMin':meta.quantizeCalMinBand3, 'qCalMax':meta.quantizeCalMaxBand3,
                 'radMin':meta.radianceMinBand3, 'radMax':meta.radianceMaxBand3},
            4 : {'qCalMin':meta.quantizeCalMinBand4, 'qCalMax':meta.quantizeCalMaxBand4,
                 'radMin':meta.radianceMinBand4, 'radMax':meta.radianceMaxBand4},
            5 : {'qCalMin':meta.quantizeCalMinBand5, 'qCalMax':meta.quantizeCalMaxBand5,
                 'radMin':meta.radianceMinBand5, 'radMax':meta.radianceMaxBand5},
            61 : {'qCalMin':meta.quantizeCalMinBand61, 'qCalMax':meta.quantizeCalMaxBand61,
                 'radMin':meta.radianceMinBand61, 'radMax':meta.radianceMaxBand61},
            62 : {'qCalMin':meta.quantizeCalMinBand62, 'qCalMax':meta.quantizeCalMaxBand62,
                 'radMin':meta.radianceMinBand62, 'radMax':meta.radianceMaxBand62},
            7 : {'qCalMin':meta.quantizeCalMinBand7, 'qCalMax':meta.quantizeCalMaxBand7,
                 'radMin':meta.radianceMinBand7, 'radMax':meta.radianceMaxBand7},
            8 : {'qCalMin':meta.quantizeCalMinBand8, 'qCalMax':meta.quantizeCalMaxBand8,
                 'radMin':meta.radianceMinBand8, 'radMax':meta.radianceMaxBand8}
        }

    if sensor == 'landsat5':
        bandFiles = {
            1 : meta.band1FileName, 2 : meta.band2FileName, 3 : meta.band3FileName,
            4 : meta.band4FileName, 5 : meta.band5FileName, 6 : meta.band6FileName,
            7 : meta.band7FileName }
    
        bands = {
            1 : {'qCalMin':meta.quantizeCalMinBand1, 'qCalMax':meta.quantizeCalMaxBand1,
                 'radMin':meta.radianceMinBand1, 'radMax':meta.radianceMaxBand1},
            2 : {'qCalMin':meta.quantizeCalMinBand2, 'qCalMax':meta.quantizeCalMaxBand2,
                 'radMin':meta.radianceMinBand2, 'radMax':meta.radianceMaxBand2},
            3 : {'qCalMin':meta.quantizeCalMinBand3, 'qCalMax':meta.quantizeCalMaxBand3,
                 'radMin':meta.radianceMinBand3, 'radMax':meta.radianceMaxBand3},
            4 : {'qCalMin':meta.quantizeCalMinBand4, 'qCalMax':meta.quantizeCalMaxBand4,
                 'radMin':meta.radianceMinBand4, 'radMax':meta.radianceMaxBand4},
            5 : {'qCalMin':meta.quantizeCalMinBand5, 'qCalMax':meta.quantizeCalMaxBand5,
                 'radMin':meta.radianceMinBand5, 'radMax':meta.radianceMaxBand5},
            6 : {'qCalMin':meta.quantizeCalMinBand6, 'qCalMax':meta.quantizeCalMaxBand6,
                 'radMin':meta.radianceMinBand6, 'radMax':meta.radianceMaxBand6},
            7 : {'qCalMin':meta.quantizeCalMinBand7, 'qCalMax':meta.quantizeCalMaxBand7,
                 'radMin':meta.radianceMinBand7, 'radMax':meta.radianceMaxBand7},
        }
    
    DOY = getDOY(meta.acquisitionDate)
    d = earthSunDist(DOY)
    cosZenAngle = cos(radians(90 - meta.sunElevation))
    
    # Get raster projection and geotransformation info:
    rasterInfo = RasterClass(bandFiles[1], gdal.GA_ReadOnly)
    # Create, write, and close the destination rasters:
    rastDriver = gdal.GetDriverByName(outRasterFormat)
    destRaster = rastDriver.Create(outFileName, rasterInfo.cols, rasterInfo.rows,
                                   len(bandsToProc), gdal.GDT_Float32)
    destRaster.SetGeoTransform(rasterInfo.geotransform)
    destRaster.SetProjection(rasterInfo.projection)
    
    # Process data each band:
    reflTable = {}
    for band in bandsToProc:
        for DN in range(int(bands[band]['qCalMin']), int(bands[band]['qCalMax'] + 1)):
            reflTable[DN] = getTOArefl(DN, bands[band]['radMax'], bands[band]['radMin'],
                                bands[band]['qCalMax'], bands[band]['qCalMin'], d,
                                esun[band] ,cosZenAngle)
    
        # Read the band into memory:
        srcRaster = gdal.Open(bandFiles[band], gdal.GA_ReadOnly)
        data = srcRaster.GetRasterBand(1).ReadAsArray()
        srcRaster = None
    
        # Change all DN values TOA reflectance:
        TOAraster = np.ndarray(shape=data.shape, dtype=np.float32)
        TOAraster.fill(NoData)
        for DN, refVal in reflTable.iteritems(): TOAraster[data==DN] = refVal
    
        destRaster.GetRasterBand(bandsToProc.index(band)+1).WriteArray(TOAraster)
        del(TOAraster)
    
    MD = destRaster.GetMetadata()
    bnames = {}
    for band in bandsToProc:
        bnames[L7bandInfo[band][0]] = L7bandInfo[band][1]
        MD['Band_'+str(bandsToProc.index(band)+1)] = L7bandInfo[band][1]

    destRaster.SetMetadata(MD)
    destRaster = None
    
    outDesc = 'Landsat TOA Reflectance Calibration [R' + str(meta.row) + 'P' \
        + str(meta.path) + ' ' + meta.acquisitionDate.strftime('%d%b%Y') + ']'
    sensor = meta.spacecraftID + ' ' + meta.sensorID
    mapInfo = meta.mapProjection + ', 1, 1, ' + str(rasterInfo.geotransform[0]) + \
       ', ' + str(rasterInfo.geotransform[3]) + ', ' + str(rasterInfo.geotransform[1]) \
       + ', ' + str(rasterInfo.geotransform[1]) + ', ' + str(meta.UTMZone) + ', ' + \
       orientation[meta.orientation] + ', ' + meta.datum + ', units=Meters'
    
    enviHdr = rasterHeader()
    enviHdr.load(outRasterFormat, outDesc, meta.samplesRef, meta.linesRef,
            len(bandsToProc), 'float32', sensor, 0, mapInfo, rasterInfo.projection,
            bnames, 'nanometers')
    hdrFile = outFileName + '.hdr'
    enviHdr.write(hdrFile)
