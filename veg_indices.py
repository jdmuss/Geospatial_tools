#!/usr/bin/env python
"""
veg_indices.py

Description:
   This python program calculate vegetation indices (EVI, NBR, NDII, NDVI, SAVI, 
   SAVI2, Tasseled Cap(brightness, greenness, & wetness), and WDRVI).

Dependencies: 
   Dependencies:
         Public:     numpy
         Private:    
    
Created by: Jordan Muss

Creation Date: 1-21-2014
Updated:       ??-??-201?  

ToDo:
          Modify the write method to reconstruct the metafile
"""
import numpy as np
""" Tell numpy to ignore errors (specifically divide by zero errors. These will 
    be caught and replaced by 'np.isnan'.
"""
np.seterr(all='ignore')

class vegIndex:
    def __init__(self, noDataVal = np.NaN):
        #self.raster = np.array([])
        self.NoData = noDataVal
    def EVI(self, NIR, Red, Blue, missing_val=None):
        if not missing_val : missing_val = self.NoData
        # Enhanced Vegetation Ratio
        c1 = 6.0
        c2 = 7.5
        L = 1
        mask = np.where(np.ma.masked_equal(Blue, self.NoData).mask | 
                        np.ma.masked_equal(Red, self.NoData).mask |
                        np.ma.masked_equal(NIR, self.NoData).mask)
        EVI = 2.5*(NIR-Red)/(c1*NIR - c2*Blue + L)
        if np.isfinite(self.NoData) : EVI[np.isnan(EVI)] = missing_val
        EVI[mask] = missing_val
        return EVI
    def NBR(self, NIR, SWIR2, rescale = True, missing_val=None):
        if not missing_val : missing_val = self.NoData
        # Normalized Burn Ratio; Roy et al (2006) GRSL
        mask = np.where(np.ma.masked_equal(NIR, self.NoData).mask |
                        np.ma.masked_equal(SWIR2, self.NoData).mask)
        NBR = 1.*(NIR-SWIR2)/(NIR+SWIR2)
        if rescale:
            NBR[NBR < -1] = -1
            NBR[NBR > 1] = 1
        if np.isfinite(self.NoData) : NBR[np.isnan(NBR)] = missing_val
        NBR[mask] = missing_val
        return NBR
    def NDII(self, SWIR, NIR, rescale = True, missing_val=None):
        if not missing_val : missing_val = self.NoData
        # Normalized difference infrared index also called Normalized difference water
        #  index (Gao, B. (1996) RSE); sensitive to leaf water content
        mask = np.where(np.ma.masked_equal(NIR, self.NoData).mask |
                        np.ma.masked_equal(SWIR, self.NoData).mask)
        NDII = 1.*(NIR-SWIR)/(NIR+SWIR)
        if rescale:
            NDII[NDII < -1] = -1
            NDII[NDII > 1] = 1
        if np.isfinite(self.NoData) : NDII[np.isnan(NDII)] = missing_val
        NDII[mask] = missing_val
        return NDII
    def NDVI(self, NIR, Red, rescale = True, missing_val=None):
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(NIR, self.NoData).mask |
                        np.ma.masked_equal(Red, self.NoData).mask)
        NDVI = 1.*(NIR-Red)/(NIR+Red)
        if rescale:
            NDVI[NDVI < -1] = -1
            NDVI[NDVI > 1] = 1
        if np.isfinite(self.NoData) : NDVI[np.isnan(NDVI)] = missing_val
        NDVI[mask] = missing_val
        return NDVI
    def SAVI(self, NIR, Red, L=0.5, missing_val=None):
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(NIR, self.NoData).mask |
                        np.ma.masked_equal(Red, self.NoData).mask)
        SAVI = (1.+L) * (NIR-Red)/(NIR + Red + L)
        if np.isfinite(self.NoData) : SAVI[np.isnan(SAVI)] = missing_val
        SAVI[mask] = missing_val
        return SAVI
    def SAVI2(self, NIR, Red, missing_val=None):
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(NIR, self.NoData).mask |
                        np.ma.masked_equal(Red, self.NoData).mask)
        SAVI2 = (2*NIR + 1 - np.sqrt((2*NIR + 1)**2 - 8*(NIR - Red)))/2.0
        if np.isfinite(self.NoData) : SAVI2[np.isnan(SAVI2)] = missing_val
        SAVI2[mask] = missing_val
        return SAVI2
    def TC_Brightness(self, b1, b2, b3, b4, b5, b7, missing_val=None):
        '''Tasseled Cap Brightness. Coefficients are from Crist & Cicone (TGRS 1984). These
           should be for atmospherically corrected landsat TM & ETM+ data where:
           b1 = blue (450 - 520 nm), b2 = green (520 - 600 nm), b3 = red (630 - 690 nm), 
           b4 = NIR (760 - 900 nm), b5 = SWIR1 (1550 - 1750 nm), & b7 = SWIR2 (2080 - 2350 nm)
        '''
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(b1, self.NoData).mask | 
                        np.ma.masked_equal(b2, self.NoData).mask |
                        np.ma.masked_equal(b3, self.NoData).mask |
                        np.ma.masked_equal(b4, self.NoData).mask |
                        np.ma.masked_equal(b5, self.NoData).mask |
                        np.ma.masked_equal(b7, self.NoData).mask)

        c1 = 0.2043 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.3037
        c2 = 0.4158 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.2793
        c3 = 0.5524 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.4743
        c4 = 0.5741 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.5585
        c5 = 0.3124 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.5082
        c6 = 0.2303 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.1863

        Brightness = c1*b1 + c2*b2 + c3*b3 + c4*b4 + c5*b5 + c6*b7
        if np.isfinite(self.NoData) : Brightness[np.isnan(Brightness)] = missing_val
        Brightness[mask] = missing_val
        return Brightness
    def TC_Greenness(self, b1, b2, b3, b4, b5, b7, missing_val=None):
        '''Tasseled Cap Greenness. Coefficients are from Crist & Cicone (TGRS 1984). These
           should be for atmospherically corrected landsat TM & ETM+ data where:
           b1 = blue (450 - 520 nm), b2 = green (520 - 600 nm), b3 = red (630 - 690 nm), 
           b4 = NIR (760 - 900 nm), b5 = SWIR1 (1550 - 1750 nm), & b7 = SWIR2 (2080 - 2350 nm)
        '''
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(b1, self.NoData).mask | 
                        np.ma.masked_equal(b2, self.NoData).mask |
                        np.ma.masked_equal(b3, self.NoData).mask |
                        np.ma.masked_equal(b4, self.NoData).mask |
                        np.ma.masked_equal(b5, self.NoData).mask |
                        np.ma.masked_equal(b7, self.NoData).mask)

        c1 = -0.1603 # from Patrick Griffiths/Tobias Kuemmerle instead of -0.2848
        c2 = -0.2819 # from Patrick Griffiths/Tobias Kuemmerle instead of -0.2435
        c3 = -0.4934 # from Patrick Griffiths/Tobias Kuemmerle instead of -0.5436
        c4 = 0.7940  # from Patrick Griffiths/Tobias Kuemmerle instead of 0.7243
        c5 = -0.0002 # from Patrick Griffiths/Tobias Kuemmerle instead of 0.0840
        c6 = -0.1446 # from Patrick Griffiths/Tobias Kuemmerle instead of -0.1800

        Greenness = c1*b1 + c2*b2 + c3*b3 + c4*b4 + c5*b5 + c6*b7
        if np.isfinite(self.NoData) : Greenness[np.isnan(Greenness)] = missing_val
        Greenness[mask] = missing_val
        return Greenness
    def TC_Wetness(self, b1, b2, b3, b4, b5, b7, missing_val=None):
        '''Tasseled Cap Wetness. Coefficients are from Crist & Cicone (TGRS 1984). These
           should be for atmospherically corrected landsat TM & ETM+ data where:
           b1 = blue (450 - 520 nm), b2 = green (520 - 600 nm), b3 = red (630 - 690 nm), 
           b4 = NIR (760 - 900 nm), b5 = SWIR1 (1550 - 1750 nm), & b7 = SWIR2 (2080 - 2350 nm)
        '''
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(b1, self.NoData).mask | 
                        np.ma.masked_equal(b2, self.NoData).mask |
                        np.ma.masked_equal(b3, self.NoData).mask |
                        np.ma.masked_equal(b4, self.NoData).mask |
                        np.ma.masked_equal(b5, self.NoData).mask |
                        np.ma.masked_equal(b7, self.NoData).mask)

        c1 = 0.0315  # from Patrick Griffiths/Tobias Kuemmerle instead of 0.1509
        c2 = 0.2021  # from Patrick Griffiths/Tobias Kuemmerle instead of 0.1973
        c3 = 0.3102  # from Patrick Griffiths/Tobias Kuemmerle instead of 0.3279
        c4 = 0.1594  # from Patrick Griffiths/Tobias Kuemmerle instead of 0.3406
        c5 = -0.6806 # from Patrick Griffiths/Tobias Kuemmerle instead of -0.7112
        c6 = -0.6109 # from Patrick Griffiths/Tobias Kuemmerle instead of -0.4572

        Wetness = c1*b1 + c2*b2 + c3*b3 + c4*b4 + c5*b5 + c6*b7
        if np.isfinite(self.NoData) : Wetness[np.isnan(Wetness)] = missing_val
        Wetness[mask] = missing_val
        return Wetness
    '''--------------------------------------------------------------------------'''
    def MODIS_TC(self, b1, b2, b3, b4, b5, b6, b7, transform, missing_val=None):
        '''Tasseled Cap Brightness. MODIS specific oefficients are from 
           Lobsler & Ciohen (IJRS 2007). These seem to be the most reliable
           coefficients, but others are possible (e.g. Zhang et al. (GRSS 2002)).
           b1 = red (620 - 670nm), b2 = NIR (841 - 876 nm), b3 = blue (459 - 479 nm), 
           b4 = Green (545 - 565 nm), b5 = MIR (1230 - 1240 nm),
           b6 = SWIR1 (1628 - 1652 nm), & b7 = SWIR2 (2105 - 2155 nm)
        '''
        Coefficients = {'brightness':{'c1':0.4395, 'c2':0.5945, 'c3':0.2460,
                                     'c4':0.3918, 'c5':0.3506, 'c6':0.2136,
                                     'c7':0.2678},
                       'greenness':{'c1':-0.4064, 'c2':0.5129, 'c3':-0.2744,
                                    'c4':-0.2893, 'c5':0.4822, 'c6':-0.0036,
                                    'c7':-0.4169},
                       'wetness':{'c1':0.1147, 'c2':0.2489, 'c3':0.2408,
                                  'c4':0.3132, 'c5':-0.3122, 'c6':-0.6416,
                                  'c7':-0.5087}
                      }
        t = transform.lower()
        if transform.lower() in Coefficients:
            if not missing_val : missing_val = self.NoData
            mask = np.where(np.ma.masked_equal(b1, self.NoData).mask | 
                            np.ma.masked_equal(b2, self.NoData).mask |
                            np.ma.masked_equal(b3, self.NoData).mask |
                            np.ma.masked_equal(b4, self.NoData).mask |
                            np.ma.masked_equal(b5, self.NoData).mask |
                            np.ma.masked_equal(b6, self.NoData).mask |
                            np.ma.masked_equal(b7, self.NoData).mask)
            TC_band = Coefficients[t]['c1']*b1 + Coefficients[t]['c2']*b2 + \
                      Coefficients[t]['c3']*b3 + Coefficients[t]['c4']*b4 + \
                      Coefficients[t]['c5']*b5 + Coefficients[t]['c6']*b6 + \
                      Coefficients[t]['c7']*b7
            if np.isfinite(self.NoData) : TC_band[np.isnan(TC_band)] = missing_val
            TC_band[mask] = missing_val
        else:
            print("Error: %s is not a valid MODIS tasseled cap transformation" % transform)
            TC_band = False
        return TC_band
    '''--------------------------------------------------------------------------'''
    def TC_BG_distance(self, greenness, brightness, missing_val=None):
        '''Tasseled Cap brightness-greenness distance (from Duane et al, Forest
           Science 2010).
        '''
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(brightness, self.NoData).mask | 
                        np.ma.masked_equal(greenness, self.NoData).mask)

        distance = np.sqrt(brightness**2 + greenness**2)
        if np.isfinite(self.NoData) : distance[np.isnan(distance)] = missing_val
        distance[mask] = missing_val
        return distance
    def WDRVI(self, NIR, Red, alpha=1, missing_val=None):
        # Calculate NDVI with a wide dynamic range:
        if not missing_val : missing_val = self.NoData
        mask = np.where(np.ma.masked_equal(NIR, self.NoData).mask |
                        np.ma.masked_equal(Red, self.NoData).mask)
        WDRVI = (NIR-Red)/(NIR+Red)
        if np.isfinite(self.NoData) : WDRVI[np.isnan(WDRVI)] = missing_val
        WDRVI[mask] = missing_val
        return WDRVI
    def close(self):
        self.raster.close()
