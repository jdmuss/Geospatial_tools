#!/usr/bin/env python
"""
Summarize_rescale_raster.py
 
Description:
   This python program will open an envi raster and extract summary statistics.
 
Dependencies: 
         Public:     getopt, numpy, os, sys, time, __future__
         Private:    Py_utils, rasterClass
     
Created by: Jordan Muss
            
 
Creation Date: 8-1-2013
Updated:       3-01-2015   Generalized the routine to be callable from another
                              script/procedure
 
ToDo: This is for single band rasters, and should be generalized for multiband
         rasters.          
"""
from __future__ import print_function
import getopt, os, sys
import numpy as np
from Py_utils import pathExpander
from rasterClass import *
from time import time, clock

def summarize_raster(src_raster_name, dest_raster_name, rescale_raster_template, binary_flag=False):
    '''----------------------------------------------------------------------------'''
    ''' Open and read the source raster to be summarized: '''
    src_rast = RasterClass(src_raster_name, gdal.GA_ReadOnly)
    src_geotransform = src_rast.geotransform
    src_projection = src_rast.projection
    src_data = src_rast.read_band(1)

    rescale_rast = RasterClass(rescale_raster_template, gdal.GA_ReadOnly)
    rescale_geotransform = rescale_rast.geotransform
    rescale_projection = rescale_rast.projection
    rescale_rast.close()

    if rescale_projection == src_projection:
        x_dim = int(np.abs(np.ceil(src_data.shape[1]/rescale_geotransform[1]*src_geotransform[1])))
        y_dim = int(np.abs(np.ceil(src_data.shape[0]/rescale_geotransform[5]*src_geotransform[5])))
        x_rescale_ratio = np.abs(rescale_geotransform[1]/src_geotransform[1])
        y_rescale_ratio = np.abs(rescale_geotransform[5]/src_geotransform[5])
        summary_rast = np.zeros(shape=(y_dim, x_dim))
        for y in range(y_dim):
            y_start = y*y_rescale_ratio
            y_end = int(np.floor(np.min([( y_start + y_rescale_ratio), src_data.shape[0] ])))
            y_start = int(np.floor(y_start))
            for x in range(x_dim):
                x_start = x*x_rescale_ratio
                x_end = int(np.floor(np.min([( x_start + x_rescale_ratio), src_data.shape[1] ])))
                x_start = int(np.floor(x_start))
                block = src_data[y_start:y_end, x_start:x_end]
                summary_rast[y,x] = float(np.sum(block))/np.product(block.shape)
        print(write_raster([summary_rast], rescale_geotransform, rescale_projection, outRast=dest_raster_name, no_data=0))
    else:
        print("Error: the projections of the source and rescale rasters do not match.")
    src_rast.close()
    if binary_flag:
        bin_dest_raster_name = '_bin.'.join(dest_raster_name.rsplit('.'))
        cutoff_mask = np.ma.masked_greater_equal(summary_rast, 0.5).mask
        binary_rast = np.zeros(shape = summary_rast.shape, dtype='int8')
        binary_rast[cutoff_mask] = 1
        print(write_raster([binary_rast], rescale_geotransform, rescale_projection, outRast=bin_dest_raster_name, no_data=0))

if __name__ == "__main__":
    keep_time = False
    convert_to_binary = False
    usage = '\n\tSummarize_rescale_raster.py  -s <source-raster> -d <destination-raster> -c <clip-raster> [-b, -t, --Time]\n'
    description = "Summarize and rescale a raster."

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hs:d:c:bt', ["source=", "destination=", "clip=", "AsBinary", "Time"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if not opts:
        print(usage)
        sys.exit(2)
    for opt in opts:
	if opt[0] == '-h':
	    print(usage)
	    sys.exit()
	elif opt[0] in ("-s"):
	    src_name = pathExpander(opt[1])
	elif opt[0] in ("-d"):
	    dest_name = pathExpander(opt[1])
	elif opt[0] in ("-c"):
	    template_name = pathExpander(opt[1])
	elif opt[0] in ("-b", "--AsBinary"):
	    convert_to_binary = True
	elif opt[0] in ("-t", "--Time"):
	    keep_time = True
    if keep_time:
	w_t0 = time.time()
	p_t0 = time.clock()
    try:
        src_name, opt, template_name
    except:
        print("Error: required parameters are missing.")
	print(usage)
        sys.exit(2)
    result = summarize_raster(src_name, dest_name, template_name, binary_flag=convert_to_binary)
    if keep_time:
	print ("\n %f seconds of processor time" % (time.clock() - p_t0), file=sys.stderr)
	print ("\n %f seconds of real time" % (time.time() - w_t0), file=sys.stderr)    