# Geospatial_tools
Python-based geospatial tools that rely heavily on GDAL libraries
This is a work in progress. These libraries began as helper scripts to open and process rasters and shapefiles
using GDAL, and have slowly evolved int more robust scripts. The primary caveat is that there may be legacy
sections, written for specific purposes, which have not yet been updated or tested for more general purposes.

It should be noted that these scripts require GDAL and Numpy at a minimum. Some classes or functions
may also require SciPy, Scikit-learn, and PIL

The main files are:

1) gdal_defs.py: translation dictionaries and other variables that are referenced by classes in each library
2) rasterClass.py: general tools to open and process raster files
3) Raster_tools.py: tools to mask rasters, perform stretches, and the beginning of some machine-learning
      classification techniques (this is a work in progress and requires some re-orginization)
4) veg_indices.py: a class to calculate some useful VIs
5) shapeFileClass.py: Classes and functions to open and process shapefiles
6) convertDNtoTOAReflectance.py: Class to convert raw landsat DN values to top of atmosphere reflectance values
7) ASCII_Grid_class.py: Class to read and convert between ESRI ASCII Grids & GeoTiffs

Summarize_rescale_Raster.py has been deprecated

# Updated 4-01-2015
