#!/usr/bin/env python
"""
ASCII_Grid_class.py

Description:
   This python script will copy the header of one ascii grid file to another,
   creating a new ascii grid. It assumes that the two files have the same
   dimensions, just that there is a problem with the header in the original file.
 
Dependencies: 
         Public:     numpy, os, sys, __future__
         Private:    Py_utils
     
Created by: Jordan Muss
 
Creation Date: 11-3-2015
Updated:       
   Author:     Jordan Muss 
   Mod. Date:  ??-??-20?? 
 
ToDo: 1) This is a work in progress. It is a one-off, but could be expanded
         into a general utility if a need arises.
"""
from __future__ import print_function
import os, sys
import numpy as np
from rasterClass import gdalDataTypes
from Py_utils import pathExpander
'''--------------------------------------------------------------------------'''
'''-------------------End functions-------------------------------------------'''
class GridClass:
    def __init__(self, file_name=None, **kwargs):
	self.EmptyInfoFile()
        if file_name is not None:
	    if os.path.exists(file_name) and os.path.isfile(file_name):
		self.file_name = file_name
		if self.GetHeader(**kwargs):
		    self.GetData()
		    self.isGrid = True
	    else:
		print("WARNING: Either %s does not exist or is not an ascii grid." % file_name)
		self.EmptyInfoFile()

    def EmptyInfoFile(self):
	self.empty = True
        self.file_name = None
        self.isGrid = False
	self.cols = None
	self.rows = None
	self.cell_size = None
	self.xllcorner = None
	self.yllcorner = None
	self.xllcenter = None
	self.yllcenter = None
        self.nodata_value = None
        self.data_type = None
	self.lines = None
	self.projection = None
	self.geotransform = None
	self.raster_band = None

    def GetHeader(self, **kwargs):
	projection_file = kwargs.pop("proj_file", format("%s.prj" % self.file_name.split('.asc')[0]))
	success = False
	with open(self.file_name, 'r') as infile:
	    ''' Read the header: '''
	    ''' ToDo: Test that value types are correct. Check that both 
		corners/or centers are provided.  '''
	    l = infile.readline().strip().split()
	    if len(l) == 2:
		if l[0].lower() == 'ncols':
		    self.cols = int(l[1])
		else:
		    self.header_error('ncols')
		    return False
	    l = infile.readline().strip().split()
	    if len(l) == 2:
		if l[0].lower() == 'nrows':
		    self.rows = int(l[1])
		else:
		    self.header_error('nrows')
		    return False
	    l = infile.readline().strip().split()
	    if len(l) == 2:
		if (l[0].lower() == 'xllcenter') or (l[0].lower() == 'xllcorner'):
		    if l[0].lower() == 'xllcenter':
			self.xllcenter = float(l[1])
		    else:
			self.xllcorner = float(l[1])
		else:
		    self.header_error('xllcenter/xllcorner')
		    return False
	    l = infile.readline().strip().split()
	    if len(l) == 2:
		if (l[0].lower() == 'yllcenter') or (l[0].lower() == 'yllcorner'):
		    if l[0].lower() == 'yllcenter':
			self.yllcenter = float(l[1])
		    else:
			self.yllcorner = float(l[1])
		else:
		    self.header_error('yllcenter/yllcorner')
		    return False
	    l = infile.readline().strip().split()
	    if len(l) == 2:
		if (l[0].lower() == 'cellsize') and (float(l[1]) > 0):
		    self.cell_size = float(l[1])
		else:
		    self.header_error('cell_size')
		    return False
	    l = infile.readline().strip().split()
	    if len(l) == 2:
		if (l[0].lower() == 'nodata_value') or (l[0].lower() == 'nodatavalue'):
		    first_value = infile.readline().strip().split()[0]
		    try:
			self.nodata_value = int(l[1])
			fv = int(first_value)
		    except:
			self.nodata_value = float(l[1])
			fv = float(first_value)
		    self.data_type = type(self.nodata_value).__name__
		else:
		    self.header_error('nodata_value')
		    return False
	    success = True
	    ''' Get projection information if it exists:'''
	    self.readProjectionFile(projection_file)
	return success

    def GetData(self):
	status = False
	with open(self.file_name, 'r') as infile:
	    status = True
	    ''' ToDo: Test that each line has the correct number of columns. '''
	    self.lines = infile.readlines()
	    ''' Remove all blank lines (=='\n') from the list:'''
	    self.lines = filter(lambda x:x!='\n', self.lines)
	    ''' Skip the header: '''
	    self.lines = self.lines[6:]
	    if self.rows != len(self.lines):
		print("Error: the number of rows in the grid does not match the header.")
		#self.EmptyInfoFile()
		return False
	    for l in self.lines:
		if self.cols != len(l.strip().split()):
		    print("Error: the number of columns in the grid does not match the header.")
		    status = False
		    self.EmptyInfoFile()
		    break
	    self.empty = False
	return status

    def GetValue(self, x=None, y=None):
	if (x is None) or (y is None):
	    print("ASCII_GRID_CLass.GetValue requires  x & y index values")
	    return None
	else:
	    if self.raster_band is None : self.GridToRaster()
	    return self.raster_band[y, x]

    def RasterToGrid(self, **kwargs):
	''' Load a raster into an ASCII GRID: '''
	raster_band = kwargs.pop("raster_band", None)
	geotransform = kwargs.pop("geotransform", None)
	llX = kwargs.pop("llX", None)
	llY = kwargs.pop("llY", None)
	corner = kwargs.pop("corner", None)
	noData = kwargs.pop("noData", -9999)
	data_type = kwargs.pop("data_type", None)
	projection = kwargs.pop("projection", None)
	if raster_band is not None:
	    self.empty = False
	    self.isGrid = True
	    self.cols = raster_band.shape[1]
	    self.rows = raster_band.shape[0]
	    if data_type is not None:
		self.data_type = data_type
	    else:
		self.data_type = gdalDataTypes[raster_band.dtype.name]
	    self.arrayToGridStruct(raster_band)
	    self.nodata_value = noData
	    if geotransform is not None:
		self.cell_size = geotransform[1]
		llX = geotransform[0]
		llY = geotransform[3] + self.rows*geotransform[5]
	    else:
		self.cell_size = 1
		if llX is None:
		    print("Warning: Neither a geotransform or lower left x coordinate were provided. 0 will be assigned.")
		    llX = 0
		if llY is None:
		    print("Warning: Neither a geotransform or lower left y coordinate were provided. 0 will be assigned.")
		    llY = 0
	    if corner:
		self.xllcorner = llX
		self.yllcorner = llY
	    else:
		self.xllcenter = llX
		self.yllcenter = llY
	    if projection is not None:
		self.projection = projection
	    return True
	else:
	    print("Error: RasterToGrid requires a raster band.")
	    self.empty = True
	    return False

    def arrayToGridStruct(self, arr):
	if arr.dtype.name.lower().find('int') > -1:
	    format_string = '{0}'
	else:
	    format_string = '{:f}'
	self.lines = []
	for row in arr:
	    self.lines.append('\t'.join(format_string.format(item) for item in row))

    def GridToRaster(self, **kwargs):
	''' Load the data from an ASCII GRID into a 2D numpy area & prepare
	    geotransform and projection files for future GDAL raster operations: '''
	projection_file = kwargs.pop("projection_file", None)
	if (((self.xllcorner is not None) and (self.yllcorner is not None)) or \
	    ((self.xllcenter is not None) and (self.yllcenter is not None))) and \
	   (self.cols is not None) and (self.rows is not None) and(self.cell_size is not None):
	    if (self.xllcorner != None) and (self.yllcorner != None):
		ulX = self.xllcorner
		ulY = self.yllcorner
	    else:
		ulX = self.xllcenter
		ulY = self.yllcenter
	    ulY += self.rows*self.cell_size
	    self.geotransform = (ulX, self.cell_size, 0, ulY, 0, -self.cell_size)
	    self.raster_band = np.full((self.rows, self.cols), self.nodata_value, dtype=self.data_type)
	    if self.data_type == 'int':
		toNumeric = lambda x:int(float(x))
	    else:
		toNumeric = lambda x:float(x)
	    for r in range(self.rows):
		self.raster_band[r] = map(toNumeric, self.lines[r].split())
	    ''' Get the projection from an external file:'''
	    if projection_file is not None : self.readProjectionFile(projection_file)
	    return True
	else:
	    print("Error: GridToRaster failed, bad geotransform data.")
	    return False

    def header_error(self, line):
	print("There was an error reading the header (%s) of grid (%s)" % (line, self.file_name))
	self.EmptyInfoFile()
	return

    def WriteGrid(self, out_file_name):
	self.file_name = out_file_name
	if not self.empty:
	    with open(out_file_name, 'w') as outfile:
		''' Write the header: '''
		print("ncols %d" % self.cols, file=outfile)
		print("nrows %d" % self.rows, file=outfile)
		if self.xllcorner != None:
		    print("xllcorner %s" % repr(self.xllcorner), file=outfile)
		else:
		    ''' This assumes the header is correct: '''
		    print("xllcenter %s" % repr(self.xllcenter), file=outfile)
		if self.yllcorner != None:
		    print("yllcorner %s" % repr(self.yllcorner), file=outfile)
		else:
		    ''' This assumes the header is correct: '''
		    print("yllcenter %s" % repr(self.yllcenter), file=outfile)
		print("cellsize %s" % repr(self.cell_size), file=outfile)
		print("nodata_value %s" % repr(self.nodata_value), file=outfile)
		for l in self.lines:
		    print(l, file=outfile)
		if self.projection is not None:
		    self.writeProjectionFile(out_file_name=format("%s.prj" % self.file_name.split('.asc')[0]))
		return True
	return False

    def writeProjectionFile(self, **kwargs):
	out_file_name = kwargs.pop("out_file_name", None)
	if out_file_name is None : out_file_name = self.file_name.split('.asc')[0]
	if self.projection is not None:
	    with open(out_file_name, 'w') as outfile:
		print(self.projection, file=outfile)

    def readProjectionFile(self, proj_file, **kwargs):
	error_msg = kwargs.pop("error_msg", False)
	if os.path.exists(proj_file):
	    with open(proj_file, 'r') as prjFile:
		''' Read the projection file: '''
		self.projection = prjFile.readline().strip()
		return True
	else:
	    if error_msg:
		print("GRID Error: projection file does not exist.")
	    return False
