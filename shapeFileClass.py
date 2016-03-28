#!/usr/bin/env python
"""
shapeFileClass.py
Description:
   These classes use OGR to create, or get file information from vectors.
   It uses the vector named <file_name>. Currently it has only been tested with
   Shapefiles.

Dependencies: 
         Public:     gdal, ogr, os, osr, __future__
         Private:    
   
Created by: Jordan Muss
Creation Date:  02-05-2014
Updated:        03-24-2016   Updated the feature class to get feature values and keys.
                             Added a buffer method to the ShapeFileClass to
			        create a new shapefile of buffered features
				based on the original shapefile in the class.
				Added 'rasterize' function.

ToDo:
      1) ToDo: add an option/method to merge the buffered features.
      2) Merge features within a specified distance tolerance.
      3) Create an OGR_geometry-type dictionary. Types are:
          ogr.wkbGeometryCollection, ogr.wkbLineString, ogr.wkbLinearRing, ogr.wkbMultiLineString,
	  ogr.wkbMultiPoint, ogr.wkbMultiPolygon, ogr.wkbNDR, ogr.wkbNone, ogr.wkbPoint,
	  ogr.wkbPolygon, ogr.wkbUnknown, ogr.wkbXDR
      4) Merge/make universal & test calls to feature class
      5) Push rasterize into the class? Merge features from earlier version of
           function

"""
from __future__ import print_function
import os, sys
from gdal_defs import *
try:
    from osgeo import gdal, ogr, osr
except:
    import gdal, ogr, osr

ESRI = 'ESRI Shapefile'
e00 = 'AVCE00'

vec_open_methods = { 'read'        : gdal.GA_ReadOnly,
                     'r'           : gdal.GA_ReadOnly,
                     'write'       : gdal.GA_Update,
                     'w'           : gdal.GA_Update,
                     'read/write'  : gdal.GA_Update,
                     'rw'          : gdal.GA_Update,
                     'update'      : gdal.GA_Update,
                     'u'           : gdal.GA_Update
                   }

ogrDataTypes = { 'integer'        : ogr.OFTInteger,       # Simple 32-bit integer
                 'int'            : ogr.OFTInteger,       # Simple 32-bit integer
                 'integerlist'    : ogr.OFTIntegerList,   # List of 32-bit integers
                 'intlist'        : ogr.OFTIntegerList,   # List of 32-bit integers
                 'real'           : ogr.OFTReal,          # Double precision floating point numbers
                 'reallist'       : ogr.OFTRealList,      # List of double precision floating point numbers
                 'string'         : ogr.OFTString,        # String of ASCII chars
                 'stringlist'     : ogr.OFTStringList,    # Array of strings
                 'widestring'     : ogr.OFTString,        # String of ASCII chars
                 'widestringlist' : ogr.OFTStringList,    # Deprecated
                 'binary'         : ogr.OFTBinary,        # Raw binary data
                 'date'           : ogr.OFTDate,          # Date
                 'time'           : ogr.OFTTime,          # Time
                 'datetime'       : ogr.OFTDateTime       # SDate & time
               }
 
def dig(i, master, sub):
    ''' This routine will recursively search through a dictionary of shape FIDs.
        Each entry (key) in the dictionary contains a list of overlapping/touching
	FIDs for that key (feature FID). The recursion returns a list of all FIDs
        that touch or overlap the initial FID (KEY). Be sure to remove the source
	key from the resulting list, and then use set and sort to remove duplicate
	values and sort the list'''

    for v in master[i]:
	if v not in sub:
	    sub.append(v)
	    sub = dig(v, master, sub)
    return sub

class ShapeFileClass:
    def __init__(self, file_name, shpFormat=ESRI, open_method='r', \
                 spatial_ref=4326, geometry_type=ogr.wkbPolygon):
	''' The default spatial reference (4326) is WGS84 '''
        outDriver = ogr.GetDriverByName(shpFormat)
	self.file_name = file_name
	self.name = os.path.basename(file_name)
	self.shapeDriverName = shpFormat
        if os.path.exists(file_name):
	    self.vector = outDriver.Open(file_name, vec_open_methods[open_method.lower()])
	    self.num_layers = self.vector.GetLayerCount()
	    self.layer = self.vector.GetLayer()
	    self.layer_name = self.layer.GetName()
	    self.featureDefn = self.layer.GetLayerDefn()
	    self.extent = self.layer.GetExtent()
	    self.num_features = self.layer.GetFeatureCount()
	    self.num_fields = self.layer.GetLayerDefn().GetFieldCount()
	    self.srs = self.layer.GetSpatialRef()
	    self.wkt = self.srs.ExportToPrettyWkt()
	    if self.srs is None:
		print("Warning: %s does not have an assigned spatial reference. %d will be assigned." % (self.name, 4326))
		self.create_projection(spatial_ref)
		self.EPSG = spatial_ref
		self.srs= osr.SpatialReference()
		self.srs.ImportFromEPSG(spatial_ref)		
	    else:
		if self.layer.GetSpatialRef().AutoIdentifyEPSG() == 0:
		    self.EPSG = int(self.layer.GetSpatialRef().GetAttrValue("Authority",1))
		else:
		    self.EPSG = 0
	else:
	    print("The vector file (%s) does not exist. It will be created" % file_name)
	    self.vector = outDriver.CreateDataSource(file_name) 
	    srs_type = type(spatial_ref)
	    if srs_type is int:
		''' srs is in EPSG format:                 '''
		srs = osr.SpatialReference()
		srs.ImportFromEPSG(spatial_ref)
	    elif srs_type is str:
		''' srs is in Wkt format:                 '''
		srs = osr.SpatialReference()
		srs.ImportFromWkt(spatial_ref)
	    elif srs_type is type(osr.SpatialReference()):
		srs = spatial_ref
	    else:
		print("%s is an unrecognized SRS type. WGS84 will be assigned as a default." % 4326)
		srs = osr.SpatialReference()
		srs.ImportFromEPSG(4326)
	    self.layer = self.vector.CreateLayer(file_name, srs, geom_type = geometry_type) 
	    self.featureDefn = self.layer.GetLayerDefn()
	return
    def add_polygon(self, vertices, field_name_lst=[], field_value_lst=[]):
	''' Add vertices in (X, Y) order '''
	if self.layer.TestCapability( ogr.OLCCreateField ):
	    ring = ogr.Geometry(ogr.wkbLinearRing)
	    for v in vertices:
		ring.AddPoint(v[0], v[1])
	    ''' Close the polygon'''
	    ring.AddPoint(vertices[0][0], vertices[0][1])
	    poly = ogr.Geometry(ogr.wkbPolygon)
	    poly.AddGeometry(ring)
	else:
	    print("The vector file (%s) was not opened with create/write priveleges" % self.name)
	''' add new geometry to the layer: '''
	outFeature = ogr.Feature(self.featureDefn)
	outFeature.SetGeometry(poly)
	if field_name_lst and field_value_lst:
	    if len(field_name_lst) == len(field_value_lst):
		for i in range(len(field_name_lst)):
		    fld_nm = field_name_lst[i].lower()
		    if (fld_nm == 'area') & (field_value_lst[i] == 0):
			val = outFeature.GetGeometryRef().Area()
		    elif (fld_nm == 'perimeter') & (field_value_lst[i] == 0):
			    val = outFeature.GetGeometryRef().Boundary().Length()
		    else:
			val = field_value_lst[i]
		    outFeature.SetField(field_name_lst[i], val)
	    else:
		print("Warning: Field name and value list lengths differ; they will not be assigned." % self.name)
	self.layer.CreateFeature(outFeature)
	outFeature.Destroy
    def add_feature(self, feature_geom, field_name_lst=[], field_value_lst=[]):
	''' Add a feature and its attributes to the layer '''
	if self.layer.TestCapability( ogr.OLCCreateField ):
	    feature = ogr.Feature(self.featureDefn)
	    if len(field_name_lst) == len(field_value_lst):
		for i in range(len(field_name_lst)):
		    fld_nm = field_name_lst[i].lower()
		    if (fld_nm == 'area') & (field_name_lst[i] == 0):
			val = feature_geom.Area()
		    elif (fld_nm == 'perimeter') & (field_name_lst[i] == 0):
			    val = feature_geom.Boundary().Length()
		    else:
			val = field_value_lst[i]
		    feature.SetField(field_name_lst[i], val)
	    else:
		print("Warning: Field name and value list lengths differ; they will not be assigned." % self.name)
	    feature.SetGeometry(feature_geom)
	    self.layer.CreateFeature(feature)
	    ''' Destroy the feature to free resources: '''
	    feature.Destroy()
	else:
	    print("The vector file (%s) was not opened with create/write priveleges" % self.name)
    def buffer(self, bufferDist, out_file):
	''' ToDo: add an option/method to merge the buffered features.'''
	new_layer = ShapeFileClass(out_file, shpFormat=self.shapeDriverName, \
	                    spatial_ref=self.srs, geometry_type=self.layer.GetGeomType())
	for fid in range(self.num_features):
	    f = FeatureClass(self.layer, fid)
	    buffered_geom = f.geometry.Buffer(bufferDist)
	    new_layer.add_feature(buffered_geom, f.field_list, f.val_list)

    def create_field(self, fld_name, fld_type=None, width=0, precision=0, justify=0):
	''' Add a field to the layer. Types are: 'int', 'integer', 'integerlist',
	     'intlist','real', 'reallist', 'string', 'stringlist', 'widestring',
	     'widestringlist', 'binary', 'date', 'time', 'datetime         '''
	fields = LayerClass(self.vector).field_dict.keys()
	if self.layer.TestCapability( ogr.OLCCreateField ):
	    if fld_name in fields : 
		print("%s already has a field named '%s', which will be used in the count operations." % (self.name, fld_name))
	    else:
		if type(fld_name) is ogr.FieldDefn:
		    field_name = fld_name
		else:
		    field_name = ogr.FieldDefn(fld_name, ogrDataTypes[fld_type.lower()])
		    field_name.SetWidth(width)
		    field_name.SetPrecision(precision)
		    field_name.SetJustify(justify)
		self.layer.CreateField(field_name)
	else:
	    print("The vector file (%s) was not opened with create/write priveleges" % self.name)
    def create_projection(self, new_spatial_ref):
	spatialRef = osr.SpatialReference()
	spatialRef.ImportFromEPSG(new_spatial_ref)
	
	spatialRef.MorphToESRI()
	file = open('.'.join([self.file_name.split('.')[0], 'prj']), 'w')
	file.write(spatialRef.ExportToWkt())
	file.close()	
    def dissolve(self, new_layer_name):
	''' Dissolve overlapping and touching features into single features in
	    a layer. The results will be written into a new layer with only
	    'Area' and 'Perimeter' fields. '''
	self.layer.GetSpatialRef().AutoIdentifyEPSG()
	new_layer = ShapeFileClass(new_layer_name, shpFormat=self.shapeDriverName, open_method='u', \
	                           spatial_ref=self.srs, geometry_type=self.layer.GetGeomType())
	new_layer.create_field('Area', 'real', width=15, precision=4)
	new_layer.create_field('Perimeter', 'real', width=15, precision=4)

	FID_list = {i:[] for i in range(self.layer.GetFeatureCount())}
	FID_delete_list = []

	for FID in range(self.layer.GetFeatureCount()):
	    f = self.layer.GetFeature(FID)
	    g = f.GetGeometryRef()
	    for FID_2 in FID_list:
		if FID != FID_2:
		    f2 = self.layer.GetFeature(FID_2)
		    g2 = f2.GetGeometryRef()
		    if g.Equals(g2) or g.Contains(g2):
			FID_delete_list.append(FID_2)
		    elif g.Intersects(g2) or g.Touches(g2) or g.Contains(g2):
			FID_list[FID].append(FID_2)
		    del(g2)
		    f2.Destroy()
		    del(f2)
	    f.Destroy()
	''' Remove items in 'FID_delete_list' (the dissolved FIDs) from the layer. '''
	FID_delete_list=list(set(FID_delete_list))
	FID_delete_list.sort()
	for dFID in FID_delete_list:
	    for dv in FID_list[dFID]:
		if dFID in FID_list[dv] : FID_list[dv].remove(dFID)
	    del(FID_list[dFID])
	
	for k in range(len(FID_list)):
	    if k in FID_list:
		if FID_list[k]:
		    s = dig(k, FID_list, [])
		    s.remove(k)
		    for sub_v in s:
			del(FID_list[sub_v])
		    FID_list[k] = list(set(s))
		    FID_list[k].sort()
	self.FID_delete_list = FID_delete_list
	self.FID_list = FID_list
	''' Add dissolved features with updated perimeter and area to 'new_layer. '''
	for FID, merge_list in FID_list.iteritems():
	    f = self.layer.GetFeature(FID)
	    g = f.GetGeometryRef()
	    if merge_list:
		for FID_2 in merge_list:
		    f2 = self.layer.GetFeature(FID_2)
		    g2 = f2.GetGeometryRef()
		    g = g.Union(g2)
		    f2.Destroy()
	    new_layer.add_feature(g, ['Area', 'Perimeter'], [0, 0])
	    f.Destroy()
	return
    def merge_shapefiles(self, file_list=[], out_file=None):
	''' Merge layers into one. A copy is made if no additioanllayers are
	    passed in file_list, and the name wile be <self.file_name>_copy.  '''
	field_list = []
	field_dict = {}
	if not out_file:
	    out_file = self.file_name.rsplit('.')
	    out_file = '.'.join(['_'.join([out_file[0],'copy']), out_file[1]])
	new_layer = ShapeFileClass(out_file, shpFormat=self.shapeDriverName, open_method='u', \
	                    spatial_ref=self.srs, geometry_type=self.layer.GetGeomType())
	''' Create fields in the new layer:  '''
	field_dict[-1] = {}
	ld = self.layer.GetLayerDefn()	
	for fld_idx in range(ld.GetFieldCount()):
	    field_defn = ld.GetFieldDefn(fld_idx)
	    field_name = field_defn.GetName()
	    field_dict[-1][field_name] = field_name
	    field_list.append(field_name)
	    new_layer.create_field(field_defn)
	''' Open the layers to join into the new layer: '''
	j_files = {-1:None}
	for file_idx in range(len(file_list)):
	    field_dict[file_idx] = {}
	    j_files[file_idx] = ShapeFileClass(file_list[file_idx], shpFormat=self.shapeDriverName, \
	                            spatial_ref=self.srs, geometry_type=self.layer.GetGeomType())
	    ''' Create fields in the new layer:  '''
	    ld = j_files[file_idx].layer.GetLayerDefn()	
	    for fld_idx in range(ld.GetFieldCount()):
		field_defn = ld.GetFieldDefn(fld_idx)
		field_name = field_defn.GetName()
		old_name = field_defn.GetName()
		if field_name in [f.lower() for f in field_list]:
		    i_str = str(fld_idx)
		    istr_len = len(istr) + 1
		    if len(field_name) > (10 - istr_len):
			field_name = '_'.join([field_name[:-istr_len], i_str])
		    field_defn.SetName(field_name)		
		field_dict[file_idx][old_name] = field_name
		field_list.append(field_name)
		new_layer.create_field(field_defn)
	''' Loop through each layer and add its features and feature values to the new layer: '''
	for file_idx in j_files.keys():
	    if file_idx == -1:
		''' Copy the main layer:'''
		layer = self.layer
	    else:
		layer = j_files[file_idx].layer
	    feature = layer.GetNextFeature()
	    while feature:
		g = feature.GetGeometryRef()
		feature_field_list = []
		feature_val_list = []
		for k, v in field_dict[file_idx].iteritems():
		    feature_field_list.append(v)
		    feature_val_list.append(feature.GetField(k))
		new_layer.add_feature(g, feature_field_list, feature_val_list)
		feature.Destroy()
		feature = layer.GetNextFeature()
	    if file_idx != -1:
		j_files[file_idx].close()
	new_layer.close()	    
    def close(self):
	self.vector.Destroy()
	self.vector = None
    def __del__(self):
	#self.vector.Destroy()
	self.vector = None

class LayerClass:
    def __init__(self, vector, layer_num = 0):
	self.layer = vector.GetLayer(layer_num)
	self.layer_type = self.layer.GetGeomType()
	self.layer_type_name = ogr.GeometryTypeToName(self.layer_type)
	self.name = self.layer.GetName()
	self.num_features = self.layer.GetFeatureCount()
	self.srs = self.layer.GetSpatialRef()
	self.prettyWkt = self.srs.ExportToPrettyWkt()
	self.EPSG = self.srs.GetAuthorityCode(None)
	if self.EPSG is not None : self.EPSG = int(self.EPSG)
	self.num_fields = self.layer.GetLayerDefn().GetFieldCount()
	self.field_dict = { self.layer.GetLayerDefn().GetFieldDefn(i).GetName():
	                   {'id':i,
	                    'type':self.layer.GetLayerDefn().GetFieldDefn(i).GetTypeName(),
	                    'typenum':self.layer.GetLayerDefn().GetFieldDefn(i).GetType(),
	                    'width':self.layer.GetLayerDefn().GetFieldDefn(i).GetWidth(),
	                    'precision':self.layer.GetLayerDefn().GetFieldDefn(i).GetPrecision()
	                   } for i in range(self.num_fields) }

class FeatureClass:
    def __init__(self, layer, feature_num = 0):
	self.feature = layer.GetFeature(feature_num)
	self.FID = self.feature.GetFID()
	self.num_fields = self.feature.GetFieldCount()
	self.geometry = self.feature.GetGeometryRef()
	self.field_list = []
	self.val_list = []
	for k, v in self.feature.items().iteritems():
	    self.field_list.append(k)
	    self.val_list.append(v)
	geom_type = ogr.GeometryTypeToName(layer.GetGeomType())
	if geom_type.lower() == 'point':
	    self.X = self.geometry.GetX()
	    self.Y = self.geometry.GetY()
    def __del__(self):
        self.feature = None
	self.FID = None
	self.num_fields = None
	self.geometry = None
	geom_type = None
	self.X = None
	self.Y = None

def nodes_to_ring(nodes, calc_geom_measures=False):
    perim = area = 0
    pt1_lon = nodes[0][0]
    pt1_lat = nodes[0][1]
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for pt in nodes:
	if len(pt) == 3:
	    z = pt[2]
	else:
	    z = 0
	ring.AddPoint(pt[0], pt[1], z)
	if calc_geom_measures:
	    ''' TODO: Calculate the area in some unit other than degrees.  '''    
	    perim += haversine_dist(pt1_lon, pt1_lat, pt[0], pt[1])
	    pt1_lon = pt[0]
	    pt1_lat = pt[1]	    
    ''' Close ring if it is still open:                                '''
    if nodes[0] != nodes[-1]:
	print("Warning (nodes_to_ring): the the supplied nodes do not close the ring, it will be forced closed.", file=sys.stderr)
	ring.CloseRings()
    ''' Create the polygon:                                            '''
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    return {'geom':poly, 'perim':perim, 'area':area}

def rasterize_layer(shape_fl_nm, to_rast, shpFormat=ESRI, inMem=True):
    ''' ToDo: import raster class & do it using the class.'''
    rast_driver = gdal.GetDriverByName('MEM')    
    dest_rast = rast_driver.Create('', to_rast.cols, to_rast.rows, 1, gdalDataTypes['int8'])
    dest_rast.SetGeoTransform(to_rast.geotransform)
    dest_rast.SetProjection(to_rast.projection)
    srcVect = ShapeFileClass(shape_fl_nm, shpFormat=shpFormat, open_method='r')
    gdal.RasterizeLayer(dest_rast, [1], srcVect.layer, burn_values=[1])
    band = dest_rast.GetRasterBand(1).ReadAsArray().astype('int8')
    if not inMem:
	write_raster([band], to_rast.geotransform, to_rast.projection, outRast = shape_fl_nm.replace('.shp','.tif'), no_data = -9999)
    srcVect.close()
    del srcVect, dest_rast
    return band

# Get First feature:
#shapeFile = ShapeFileClass(shpFile_name, open_method=gdal.GA_Update)
#shapeFile.layer.DeleteFeature(27)
#f = shapeFile.layer.GetFeature(12)
#g = f.GetGeometryRef()

#print(g.Distance(g2))

#geom.Area()
#geom.Boundary().Length()

#f.GetFieldType(0)
#f.GetDefnRef()
#f.GetFID()
#f.GetField('Area') # f.GetField(0)
#f.GetFieldCount()
#f.GetFieldDefnRef(2)
#f.GetFieldIndex(0)


#feature = self.layer.GetNextFeature()
#f.SetGeometry(g)
#f.SetField('Area', g.GetArea())
#f.SetField('Perimeter', g.Boundary().Length())
#new_layer.layer.SetFeature(f)
'''--------------------------------------------------------------------------
   This is for future reference. Especially when dealing with angular projections:
   
           ref_geometry = ref_feature.GetGeometryRef()
        pts = ref_geometry.GetGeometryRef(0)
        points = []
        for p in xrange(pts.GetPointCount()):
            points.append((pts.GetX(p), pts.GetY(p)))

def edges_index(points):
    """
    compute edges index for a given 2D point set

    1- The number of edges which form the polygon
    2- Perimeter
    3- The length of the longest edge in a polygon
    4- The length of the shortest edge in a polygon
    5- The average length of all of edges in a polygon
    6- The lengths of edges deviate from their mean value
    """
    Nedges = len(points)-1
    length = []
    for i in xrange(Nedges):
        ax, ay = points[i]
        bx, by = points[i+1]
        length.append(math.hypot(bx-ax, by-ay))
    edges_perimeter = numpy.sum(length)
    edges_max = numpy.amax(length)
    edges_min = numpy.amin(length)
    edges_average = numpy.average(length)
    edges_std = numpy.std(length)
    return (Nedges,edges_perimeter,edges_max,edges_min,edges_average,edges_std)
    
--------------------------------------------------------------------------'''
