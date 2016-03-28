#!/usr/bin/env python
"""
Raster_tools.py
 
Description:
   This python script contains a variety of raster processing tools, including
   y-by-x despeckling filters and masking operations.
 
Dependencies: 
         Public:     numpy, __future__
         Private:    Py_utils
     
Created by: Jordan Muss

Creation Date: 5-8-2014
Updated:       11-20-2014 Renamed the script to 'Raster_tools'. Moved the masking
                            classes and functions from 'Py_utils.py' and
			    'Display_Raster.py' into this script.

	        03-28-2016 Moved stretch & display classes and functions from
		            'Display_Raster.py'
ToDo: 
          1) Add a cummulative linear (or percentile) stretch
	  2) Add more clustering/machine learning methods
	  3) Allow the use of the full box (use_diags=True)
"""
from __future__ import print_function
#import sys
import numpy as np
#from Py_utils import *

'''-------------------------------------------------------------------------
    Masking tools:
-------------------------------------------------------------------------'''
class mask_type:
    def __init__(self, condition, band_array, mval):
	# Mask options:
	mask_opts = {'equal'      : self.equal,
	             'eq'         : self.equal,
	             'not_eq'     : self.not_equal,
	             'neq'        : self.not_equal,
	             'less'       : self.less_than,
	             'lt'         : self.less_than,
	             'less_eq'    : self.less_than_equal,
	             'leq'        : self.less_than_equal,
	             'greater'    : self.greater_than,
	             'gt'         : self.greater_than,
	             'greater_eq' : self.greater_than_equal,
	             'geq'        : self.greater_than_equal,
	             'invalid'    : self.invalid,
	             'inv'        : self.invalid
                }
	self.mval = mval
	self.mask = mask_opts[condition](band_array)
    def equal(self, b_array):
	return np.ma.masked_equal(b_array, self.mval).mask
    def not_equal(self, b_array):
	return np.ma.masked_not_equal(b_array, self.mval).mask
    def less_than(self, b_array):
	return np.ma.masked_less(b_array, self.mval).mask
    def less_than_equal(self, b_array):
	return np.ma.masked_less_equal(b_array, self.mval).mask
    def greater_than(self, b_array):
	return np.ma.masked_greater(b_array, self.mval).mask
    def greater_than_equal(self, b_array):
	return np.ma.masked_greater_equal(b_array, self.mval).mask
    def invalid(self, b_array):
	return np.ma.masked_invalid(b_array).mask

def OR_mask(src_rast, condition_list, value_list):
    out_mask = np.zeros_like(src_rast, dtype='bool')
    if len(condition_list) == len(value_list):
	for cond, val in zip(condition_list, value_list):
	    out_mask = out_mask | mask_type(cond, src_rast, val).mask
	return(out_mask)
    else:
	print("Error: The condition list length (%d) does not equal the value list length(%d)" % (en(condition_list), len(value_list)))
	return None

def AND_mask(src_rast, condition_list, value_list):
    out_mask = np.ones_like(src_rast, dtype='bool')
    if len(condition_list) == len(value_list):
	for cond, val in zip(condition_list, value_list):
	    out_mask = out_mask &  mask_type(cond, src_rast, val).mask
	return(out_mask)
    else:
	print("Error: The condition list length (%d) does not equal the value list length(%d)" % (en(condition_list), len(value_list)))
	return None
'''-------------------------------------------------------------------------
    Stretch/display tools:
-------------------------------------------------------------------------'''
def create_rgba_array(Red, Green, Blue, alpha = 0, mask_values = [0, 0, 0], \
                      m_type = 'equal', data_type = 'uint8'):
    ''' Create a 4-band RGBA array:                                        '''
    # Test whether the red, green, & blue channels are all the same size:
    if Red.shape == Green.shape == Blue.shape:
	rgb_arr = np.zeros((Red.shape[0], Red.shape[1], 4), data_type)
	# Create a mask for all pixels that are missing data in at least one band:
	red_mask = mask_type(m_type, Red, mask_values[0]).mask
	green_mask = mask_type(m_type, Green, mask_values[1]).mask
	blue_mask = mask_type(m_type, Blue, mask_values[2]).mask
	
	# Load the RGB channels into the RGB array:
	rgb_arr[:,:,0] = (255. * np.where(red_mask, 0, Red))/np.max(Red)
	rgb_arr[:,:,1] = (255. * np.where(green_mask, 0, Green))/np.max(Green)
	rgb_arr[:,:,2] = (255. * np.where(blue_mask, 0, Blue))/np.max(Blue)
	# Set the alpha channel (transparency) for the masked values. All others are opaque
	if type(alpha) == np.ndarray:
	    if alpha.shape == Red.shape:
		rgb_arr[:,:,3] = alpha
	    else:
		print("Warning: the alpha channel is not the same size as the RGB channels. Using alpha = 0")
		rgb_arr[:,:,3] = np.where(red_mask | green_mask | blue_mask, 0, 255)
	else:
	    rgb_arr[:,:,3] = np.where(red_mask | green_mask | blue_mask, alpha, 255)
    else:
	print("Error: the red, green, & blue channels are not the same size")
	rgb_arr = None
    return rgb_arr

def linear_stretch_perc(band, percentile, **kwargs):
    out_min = kwargs.pop("out_min", 0.)
    out_max = kwargs.pop("out_max", 255.)
    in_lo = np.percentile(band, percentile)
    in_hi = np.percentile(band, 100 - percentile)
    out = np.where(band < in_lo, in_lo, band)
    out = np.where(out > in_hi, in_hi, out)
    if (in_hi - in_lo) != 0:
	out = np.float_(out - in_lo)/(in_hi - in_lo) * out_max
    return out

def histogram_stretch(band, **kwargs):
    num_bins = kwargs.pop("num_bins", 10)
    invert = kwargs.pop("invert", False)
    normalize = kwargs.pop("normalize", False)
    out_min = kwargs.pop("out_min", 0.)
    out_max = kwargs.pop("out_max", 255.)
    clip_range = kwargs.pop("range", (band.min(), band.max()))
    stretch_bins = np.histogram(range(int(out_min), int(out_max)), bins=num_bins)[1].astype('int')
    if out_max not in stretch_bins:
	max_diff = out_max - stretch_bins[-1]
	stretch_bins = [v+max_diff if v>out_min else v for v in stretch_bins]
    counts, breaks = np.histogram(band, bins=num_bins, range=clip_range)
    if normalize:
	cdf = counts.cumsum()
	stretch_bins = out_max*cdf/cdf[-1]
	breaks = breaks[:-1]
    #out = np.interp(band.flatten(), breaks, stretch_bins).reshape(band.shape)
    if invert : stretch_bins = stretch_bins[::-1]
    out = np.interp(band, breaks, stretch_bins)
    return out

def normalize_raster(src_raster, axis=None):
    ''' This assumes that multi-band rasters are dimesioned (bands, rows, cols) '''
    fail = False
    message = {'a':"Normalization will occur using the entire raster",
               'b':"Normalization will occur by band",
               'p':"Normalization will occur across pixels",
               'r':"Normalization will occur across rows",
               'c':"Normalization will occur across columns"}
    rast_dims = src_raster.shape
    num_dims = len(rast_dims)
    if num_dims > 3:
        print("Warning:function can only normalize two and three dimension rasters")
        return False
    else:
        if axis is None:
            axis = 'a'
            print(message[axis])
        elif axis is int():
            if axis >= num_dims:
                alt_msg = " (axis value exceeds the raster's dimensions)"
                fail = True
            elif axis < 0:
                alt_msg = ''
                fail = True
            elif num_dims == 2:
                axis = {0:'r', 1:'c'}[axis]
            else:
                axis = {0:'p', 1:'r', 2:'c'}[axis]
        if fail or axis not in message:
            print("Warning:invalid axis value%s." % alt_msg)
            print("\tChoose a value from the following list:")
            for k, v in  message.iteritems() : print("\t\t'%s' for %s" % (k, v.replace('will occur ', '').lower()))
            return False
        axis = axis.lower()[0]
        if num_dims == 3:
            axes= {'a':None, 'p':0, 'c':1, 'r':2}
            shape = {'a':rast_dims, 'p':rast_dims[0], 'c':rast_dims[1], 'r':rast_dims[2]}
        else:
            axes= {'a':None, 'b':None, 'p':None, 'c':0, 'r':1}
            shape = {'a':rast_dims, 'p':rast_dims, 'c':rast_dims[0], 'r':rast_dims[1]}
            #d, nRows, nCols = (0,) + rast_dims
            if axis == 'p':
                print("%s %s" %
                      ("Warning:two-dimensional rasters can not be normalized at the pixel-level.", message['a']))
    print("Normalizing %s" % message[axis][25:])
    if (num_dims == 3) & (axis == 'b'):
	mean_=np.array([np.tile(np.nanmean(src_raster[b]), rast_dims[1:]) for b in range(rast_dims[0])])
	std_=np.array([np.tile(np.nanstd(src_raster[b]), rast_dims[1:]) for b in range(rast_dims[0])])
    else:
	mean_ = np.tile(np.nanmean(src_raster, axis=axes[axis]), shape[axis]).reshape(rast_dims)
	std_ = np.tile(np.nanstd(src_raster, axis=axes[axis]), shape[axis]).reshape(rast_dims)
    if (std_ == 0).any():
        print("Warning:Divide by zero error (there is a zero value in the standard deviation)")
        return False
    else:
        norm_img = (src_raster - mean_)/std_    
    return norm_img
'''-------------------------------------------------------------------------
    Moving window, despeckling,... tools:
-------------------------------------------------------------------------'''
def despeckle_3x3(src_array, fg=None, bg=0, num_passes=1, use_diags=False):
    ''' A 'fast' a 3x3 despekling filter that does not rely on looping.
        TODO: offer the option of using a complete filtering box that includes
	the diagonals. '''
    ''' Get the foreground value if it is not specified. This assumes that the
        array values are only binary. Unexpected results will occur if this is
        not the case.                                                      '''
    if fg is None : fg = [v for v in np.unique(src_array) if v != bg][0]
    if use_diags : print("Warning:The 'use diagonals' feature is not yet implemented.")
    t = np.copy(src_array)
    for p in range(num_passes):
	up = np.roll(t, -1, axis = 0)
	up[-1] = bg
    
	down= np.roll(t, 1, axis = 0)
	down[0] = bg
    
	left = np.roll(t, -1, axis = 1)
	left[:,-1] = bg
    
	right = np.roll(t, 1, axis = 1)
	right[:,0] =bg
    
	toggle_off_mask = np.ma.masked_equal(up, bg).mask & np.ma.masked_equal(down, bg).mask & \
	           np.ma.masked_equal(right, bg).mask & np.ma.masked_equal(left, bg).mask
	t[toggle_off_mask] = bg
	toggle_on_mask = np.ma.masked_equal(up, fg).mask & np.ma.masked_equal(down, fg).mask & \
	           np.ma.masked_equal(right, fg).mask & np.ma.masked_equal(left, fg).mask
	t[toggle_on_mask] = fg
    return t

def despeckle(src_array, bg=0, size=[3,3], num_passes=1, use_diags=False):
    ''' Pass a 3x3 filter across the scene first, and then grow the filtering
        box until it is the size of the target filter. This will ensure that
	smaller	elements are eliminated. TODO: offer the option of using a complete
	filtering box that includes the diagonals. '''
    if use_diags : print("Warning:The 'use diagonals' feature is not yet implemented.")
    q, r = divmod(size[0], 2)
    yStart = -q
    expanded_array = np.insert(src_array, src_array.shape[0], np.tile(np.repeat(bg, src_array.shape[1]), (q,1)), axis=0)
    if r == 0 : q -= 1
    expanded_array = np.insert(expanded_array, 0, np.tile(np.repeat(bg, expanded_array.shape[1]), (q,1)), axis=0)
    yStart += expanded_array.shape[0] - src_array.shape[0]

    q, r = divmod(size[1], 2)
    xStart = -q
    expanded_array = np.insert(expanded_array, expanded_array.shape[1], np.tile(np.repeat(bg, expanded_array.shape[0]), (q,1)), axis=1)
    if r == 0 : q -= 1
    expanded_array = np.insert(expanded_array, 0, np.tile(np.repeat(bg, expanded_array.shape[0]), (q,1)), axis=1)
    xStart += expanded_array.shape[1] - src_array.shape[1]

    for p in range(num_passes):
	for ySize in range(3, size[0]+1):
	    for xSize in range(3, size[1]+1):
		if [ySize, xSize] == [3, 3]:
		    expanded_array = despeckle_3x3(expanded_array, bg=bg, num_passes=1, use_diags=use_diags)
		else:
		    box_mask = np.ones([ySize, xSize], dtype=bool)
		    box_mask[[0,0,-1,-1],[0,-1,0,-1]] = False
		    test_mask = np.copy(box_mask)
		    test_mask[1:-1,1:-1] = False
		    test_perim = float(len(np.where(test_mask)[0]))
    
		    q, r = divmod(xSize, 2)
		    xEnd = expanded_array.shape[1] - q
		    box_cols = range(-q, q+1)
		    if r == 0 : box_cols = box_cols[1:]
		    box_cols = np.array(box_cols*ySize)
		    
		    q, r = divmod(ySize, 2)
		    yEnd = expanded_array.shape[0] - q
		    box_rows = range(-q, +q+1)
		    if r == 0 : box_rows = box_rows[1:]
		    box_rows = np.repeat(box_rows, xSize)
    
		    non_bg_idx = np.where(np.ma.masked_not_equal(expanded_array, bg).mask)
		    for i in range(len(non_bg_idx[0])):
			box = expanded_array[box_rows+non_bg_idx[0][i], box_cols+non_bg_idx[1][i]].reshape([ySize, xSize])
			if np.sum(box[test_mask])/test_perim == bg:
			    box[box_mask] = bg
			    expanded_array[box_rows+non_bg_idx[0][i], box_cols+non_bg_idx[1][i]] = box.ravel()
    return expanded_array[yStart:yStart+src_array.shape[0], xStart:xStart+src_array.shape[1]]

def reravel_image(raveled_array, w, h, mask=None, d=None):
    '''Recreate the image from an unraveled array:'''
    if mask is None:
        valid_indices = np.where(np.ones(shape=(h, w), dtype='bool'))
    else:
        valid_indices = np.where(mask)
    valid_indices = zip(valid_indices[0], valid_indices[1])
    if d is None:
        image = np.zeros(shape=(h, w))
    else:
        image = np.zeros(shape(h, w, d))
    label_idx = 0
    '''TODO: fix this for multi-band images!                              '''
    for idx in valid_indices:
        image[idx] = raveled_array[label_idx]
        label_idx += 1
    return image

def neighbor_pixel(rast_mat, flag):
    perm_mask = np.ma.masked_equal(rast_mat, flag).mask | np.ma.masked_equal(rast_mat, 0).mask
    
    up = np.roll(rast_mat, -1, axis = 0)
    up[-1] = 0
    up[perm_mask] = rast_mat[perm_mask]
    
    down= np.roll(rast_mat, 1, axis = 0)
    down[0] = 0
    down[perm_mask] = rast_mat[perm_mask]
    
    left = np.roll(rast_mat, -1, axis = 1)
    left[:,-1] = 0
    left[perm_mask] = rast_mat[perm_mask]
    
    right = np.roll(rast_mat, 1, axis = 1)
    right[:,0] = 0
    right[perm_mask] = rast_mat[perm_mask]
    
    adj_rast = np.ma.masked_equal(up, flag).mask | np.ma.masked_equal(down, flag).mask | \
               np.ma.masked_equal(right, flag).mask | np.ma.masked_equal(left, flag).mask
    rast_mat[adj_rast] = flag
    return(rast_mat)
    #return(adj_rast)

def grow_disturbance(dist_mat, flag):
    flag_pixel_ct = 0
    i = 0
    while len(np.where(dist_mat==flag)[0]) > flag_pixel_ct:
	i += 1
	flag_pixel_ct = len(np.where(dist_mat==flag)[0])
	dist_mat = neighbor_pixel(dist_mat, flag)	
    print("Convergence occured in %d passes." %  i)
    return(dist_mat)

'''-------------------------------------------------------------------------
    Cluster analysis tools:
-------------------------------------------------------------------------'''
'''------------------------------------------------------------------------'''
''' k-means clustering:   '''
'''------------------------------------------------------------------------'''
from sklearn.cluster import KMeans
#from sklearn.metrics import pairwise_distances_argmin
from sklearn.utils import shuffle
from time import time

class clustering:
    def __init__(self, src_image, mask_val=None, time_it=False):
	'''Initialize the clustering class'''
	self.messages = {'fim':"Clustering will be performed using the full data set.",
                'ssp_m':"Sub-sampling using a percentage of the data set.",
                'ssp_em':"Sub-sampling percentage must be between 0 and 1.",
                'sss_m':"Sub-sampling using a fixed number of data elements.",
                'sss_em':"Sub-sampling size must be between the number of clusters and the size of the valid data set."
               }
	self.image = src_image
	self.mask_val = mask_val
	self.time_it = time_it
    def flatten(self):
	'''--------------------------------------------------------------------'''
	''' Flatten the image into a 1-dimensional array by d-bands:           '''
	'''--------------------------------------------------------------------'''
	original_shape = self.image.shape
	valid_mask = np.ma.masked_not_equal(self.image, self.mask_val).mask & \
	             ~np.ma.masked_invalid(self.image).mask
	if len(original_shape) == 2:
	    d = 1
	    h, w = original_shape
	    image_array = np.zeros(shape=(valid_mask.sum(), d))
	    image_array[:, 0] = np.ravel(self.image[valid_mask])
	else:
	    valid_mask = np.any(valid_mask, axis = 0)
	    d, h, w = original_shape
	    image_array = np.zeros(shape=(valid_mask.sum(), d))
	    for i in range(d):
		image_array[:, i] = np.ravel(self.image[i, valid_mask])
	self.valid_data_size = valid_mask.sum()
	self.d, self.h, self.w = d, h, w
	self.valid_mask = valid_mask
	self.image_array = image_array
    def kMeans(self, n_clusters=3, sample_perc=None, sample_size=None, ordinal=False):
	''' Perform k-means clustering of an HxW or BxHxW image. This assumes the
	    image has been standardized. This may work for non-raster (image) data,
	    but has not been tested.
	-----------------------------------------------------------------------'''
	self.flatten()
	'''--------------------------------------------------------------------'''
	''' Sub-sample the data:                                               '''
	'''--------------------------------------------------------------------'''
	t0 = time()
	if (sample_size is None) and (sample_perc is None):
	    print("Warning:a sub-sampling parameter was not provided. %s" % messages['fim'])
	    image_sample = self.image_array
	else:
	    if (sample_size is not None):
		if (sample_perc is not None):
		    print("Warning:both sub-sampling size and percentage parameters were provided. Only sample-size will be used.")
		    sample_perc = None
		if (sample_size < n_clusters) or (sample_size > self.valid_data_size):
		    print("Warning:%s %s" % (messages['sss_em'], messages['fim']))
		    image_sample = self.image_array
		    sample_size = None
	    else:
		if (sample_perc < 0) or (sample_perc > 1):
		    print("Warning:an invalid sub-sampling percentage was entered. %s %s" % (messages['ssp_em'], messages['fim']))
		    image_sample = self.image_array
		    sample_size = sample_perc = None
		else:
		    sample_size = int(self.valid_mask.sum()*sample_perc/100)*100
	if (sample_size is not None):
	    print("\nFitting k-means model on a sub-sample of the image (~%d samples)" % sample_size)
	    image_sample = shuffle(self.image_array, random_state=0)[:sample_size]
	k_means = KMeans(n_clusters=n_clusters, random_state=0).fit(image_sample)
	if self.time_it : print("Fitting took %0.3fs." % (time() - t0))
	'''--------------------------------------------------------------------'''
	''' Get labels for all points:                                         '''
	'''--------------------------------------------------------------------'''
	print("Predicting clusters on the full image using the sub-sampled results")
	t0 = time()
	labels = k_means.predict(self.image_array)
	if self.time_it : print("Prediction took %0.3fs." % (time() - t0))
	'''--------------------------------------------------------------------'''
	''' Build the new k-means image with 'num_clusters':                   '''
	'''--------------------------------------------------------------------'''
	print("Constructing the new clustered image:")
	t0 = time()
	km = self.recreate_image(k_means.cluster_centers_, labels, ordinal=ordinal)
	if self.time_it : print("Construction of the clustered image took %0.3fs." % (time() - t0))
    
	return km
    def recreate_image(self,codebook, labels, ordinal=True):
	'''Recreate the image from the code-book & labels'''
	d = codebook.shape[1]
	label_idx = 0
	if self.valid_mask is None:
	    valid_indices = np.where(np.ones(shape=(self.h, self.w), dtype='bool'))
	else:
	    valid_indices = np.where(self.valid_mask)
	valid_indices = zip(valid_indices[0], valid_indices[1])
	if (d == 1) or ordinal:
	    image = np.zeros(shape=(self.h, self.w))
	else:
	    image = np.zeros((d, self.h, self.w))
	'''TODO: fix this for single-/multi-band images!
		use swapaxes?                                  '''
	if ordinal:
	    image = image.astype('int')
	    for idx in valid_indices:
		image[idx] = labels[label_idx] + 1
		label_idx += 1
	else:
	    if d == 3 : image = np.rollaxis(image, 0, 3)
	    for idx in valid_indices:
		image[idx] = codebook[labels[label_idx]]
		label_idx += 1        
	return image

