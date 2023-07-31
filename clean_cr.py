"""
Cosmic Ray Filter
"""

# IMPORTS
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from IPython import embed
import ipdb
import time


# FUNCTIONS
def difference(data):
    """
    Take the difference of the successive points in each array
    """
    return np.diff(data)


def robust_sigma(diff, zero=0):
    """
    Calculate a resistant estimate of the dispersion of
    a distribution. For an uncontaminated distribution,
    this is identical to the standard deviation.
 
    Use the median absolute deviation as the initial
    estimate, then weight points using Tukey Biweight.
    See, for example, Understanding Robust and
    Exploratory Data Analysis, by Hoaglin, Mosteller
    and Tukey, John Wiley and Sons, 1983.
 
    .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.
 
    :History:
        * H Freudenreich, STX, 8/90
        * Replace MED call with MEDIAN(/EVEN), W. Landsman, December 2001
        * Converted to Python by P. L. Lim, 11/2009
 
    Examples
    --------
    >>> result = robust_sigma(diff, zero=1)
    ----------
    Parameters
    ----------
    diff: array_like
        Vector of quantity for which the dispersion is
        to be calculated
 
    zero: int
        If set, the dispersion is calculated w.r.t. 0.0
        rather than the central value of the vector. If
        Y is a vector of residuals, this should be set.
    -------
    Returns
    -------
    out_val: float
        Dispersion value. If failed, returns -1.
 
    """
    # Flatten array
    y = diff.reshape(diff.size, )
 
    eps = 1.0E-20
    c1 = 0.6745
    c2 = 0.80
    c3 = 6.0
    c4 = 5.0
    c_err = -1.0
    min_points = 3
 
    if zero:
        y0 = 0.0
    else:
        y0 = np.median(y)
 
    dy    = y - y0
    del_y = abs( dy )
 
    # First, the median absolute deviation MAD about the median:
 
    mad = np.median( del_y ) / c1
 
    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps:
        mad = np.mean( del_y ) / c2
    if mad < eps:
        return 0.0
 
    # Now the biweighted value:
    u  = dy / (c3 * mad)
    uu = u*u
    q  = np.where(uu <= 1.0)
    count = len(q[0])
    if count < min_points:
        #print 'ROBUST_SIGMA: This distribution is TOO WEIRD! Returning', c_err
        return c_err
 
    numerator = np.sum( (y[q]-y0)**2.0 * (1.0-uu[q])**4.0 )
    n    = y.size
    den1 = np.sum( (1.0-uu[q]) * (1.0-c4*uu[q]) )
    siggma = n * numerator / ( den1 * (den1 - 1.0) )
 
    if siggma > 0:
        out_val = np.sqrt( siggma )
    else:
        out_val = 0.0
 
    return out_val


def sig_clip(diff, sigma):
    """
    Perform sigma-clipping of (sigma calculated in robust_std_deviation)
    to the data using the differences between each image flux count
    
    If difference value is well above 5sigma,
    it will get flagged and put into a flagged list
    """
    
    cr_jump_detected = []
    
    for val in range(len(diff)):
        if diff[val] >= sigma*5:
            cr_jump_detected.append(val+1)
            
    return cr_jump_detected


def line_fitting_full(data, flag_list, frames):
    """
    Full line fitting function given a certain circumstance of where the cosmic ray jump is
    """
    crJumpIdx = 0
    xsBefore = 0
    xsAfter = 0
    ysBefore = 0
    ysAfter = 0
    slopeBefore = 0
    slopeAfter = 0

    flag = np.array(flag_list)

    crJumpIdx, xsBeforeCR, ysBeforeCR, xsAfterCR, ysAfterCR = params(data, flag)

    # if cr jump is at the very beginning
    if len(ysBeforeCR) < 20:
        xsBefore, xsAfter, slopeBefore, slopeAfter = lenBefore_less(data, crJumpIdx, xsBeforeCR, ysBeforeCR, xsAfterCR, ysAfterCR)
    # if cr jump is at the very end
    elif len(ysAfterCR) < 20:
        xsBefore, xsAfter, slopeBefore, slopeAfter = lenAfter_less(data, crJumpIdx, xsBeforeCR, ysBeforeCR, xsAfterCR, ysAfterCR)
    # cr jump doesnt apply to any of the above
    else:
        xsBefore, xsAfter, slopeBefore, slopeAfter = line_fitting(xsBeforeCR, ysBeforeCR, xsAfterCR, ysAfterCR)

    xsBefore = np.array(xsBefore)
    xsAfter = np.array(xsAfter)
    slopeBefore = np.array(slopeBefore)
    slopeAfter = np.array(slopeAfter)
        
    return xsBefore, xsAfter, ysBefore, ysAfter, slopeBefore, slopeAfter


def params(data, flag):
    """ Setting the x and y values before and after the jump """
    
    crJumpIdx = int(flag)
    
    beforeCR = int(crJumpIdx-20)
    afterCR = int(crJumpIdx+20)

    xsBeforeCR = np.arange(beforeCR, crJumpIdx, 1)
    ysBeforeCR = data[beforeCR:crJumpIdx]
    
    xsAfterCR = np.arange(crJumpIdx, afterCR, 1)
    ysAfterCR = data[crJumpIdx:afterCR]
    
    return crJumpIdx, xsBeforeCR, ysBeforeCR, xsAfterCR, ysAfterCR
    
    
def line_fitting(xs_before_jump, ys_before_jump ,xs_after_jump, ys_after_jump):
    """ Line fitting of 20 points before and after the CR jump """
    
    linregBefore = stats.linregress(xs_before_jump, ys_before_jump)
    linregAfter = stats.linregress(xs_after_jump, ys_after_jump)
    
    slopeFitBefore = linregBefore.intercept + linregBefore.slope * xs_before_jump
    slopeFitAfter = linregAfter.intercept + linregAfter.slope * xs_after_jump
    
    return xs_before_jump, xs_after_jump, slopeFitBefore, slopeFitAfter
    
    
def lenBefore_less(data, jump_idx, xs_before_jump, ys_before_jump, xs_after_jump, ys_after_jump):
    """
    If the length of values before the CR jump is less than 20,
    then this function takes what the length of the values before the CR and takes the slope from that set of points
    """
    xsBefore = np.arange(0, jump_idx, 1)
    ysBefore = data[:jump_idx]
    
    xsAfter = xs_after_jump
    ysAfter = ys_after_jump
    
    slopes = line_fitting(xsBefore, ysBefore, xsAfter, ysAfter)
    
    return slopes


def lenAfter_less(data, jump_idx, xs_before_jump, ys_before_jump ,xs_after_jump, ys_after_jump):
    """
    If the length of values after the CR jump is less than 20,
    then this function takes the length of the values after the CR and takes the slope from that set of points
    """
    
    afterCR = int(jump_idx + len(ys_after_jump))
    
    xsAfter = np.arange(jump_idx, afterCR, 1)
    ysAfter = data[jump_idx:afterCR]
    
    xsBefore = xs_before_jump
    ysBefore = ys_before_jump
    
    slopes = line_fitting(xsBefore, ysBefore, xsAfter, ysAfter)
    
    return slopes


def jump_frame1(frames, data, jump_idx, xs_before_jump, ys_before_jump ,xs_after_jump, ys_after_jump):
    """
    If a jump is detected at frame 1
    then it finds the offset btwn frames 1 and 0, and brings every frame after down then finds actual offset for an initial correction
    """
    
    # value at f1 - value at f0 = offset from frame from rest of data
    offset = data[1] - data[0]
    # subtract offset from the rest of the data
    data[1:] = data[1:] - offset
    
    return data

    
def offset_correct(data, crJump, slopeBeforeCR, slopeAfterCR):
    """
    Calculating the offset between the fitted slopes (first point of slope after CR minus the last point of slope before CR),
    and applying that offset to the points after the CR jump
    """
    
    crJump = int(crJump)
    
    offset =  slopeAfterCR[0] - slopeBeforeCR[-1]
    
    ysAfterCR_corr = data[crJump:] - offset
    
    corrected_pix = np.concatenate((data[0:crJump], ysAfterCR_corr))
    
    return corrected_pix



# MAIN
if __name__ == "__main__":
 
    # Open fits file
    infile = fits.open('in_filename')   # INPUT FILE NAME

    # Name of corrected outfile name 
    outfile = 'out_filename'
    
    img_stack = infile[0].data
    frames = np.arange(0, img_stack.shape[0], 1)
    
    # Keeps pixels that do not need to be corrected
    corrected_det = img_stack

    print('\n')
    print('---------------------------------------')
    print('CORRECTING...')
    print('---------------------------------------')

    # Begin timing code
    t1 = time.time()
    
    # Iterating through each pixel in the detector
    for i in range(img_stack.shape[2]):
        for j in range(img_stack.shape[1]):
            
            singlePixel = img_stack[:, j, i]
            
            # Gets difference between each pixel
            diff = difference(singlePixel)
            
            # Computing the robust standard deviation
            sigma = robust_sigma(diff, zero=0)
            
            # Flagging pixels that have differences above 5sigma
            flag = sig_clip(diff, sigma)
            
            if len(flag) > 0:
                print('\n')
                print('PIXEL:', i, j)
                print('Cosmic rays detected at: ', flag)
                
                for jump in flag:
                    jump_amp = diff[jump-1]
                    print('JUMP AMPLITUDE:', jump_amp)
                    # Corrects for each flag found
                    xsBefore, xsAfter, ysBefore, ysAfter, slopeBefore, slopeAfter = line_fitting_full(singlePixel, jump, frames)
                    
                    corr_pix = offset_correct(singlePixel, jump, slopeBefore, slopeAfter)
                    
                    # Jump at frame 1 -- initial correction
                    if jump == 1:
                        corr_pix = jump_frame1(frames, singlePixel, jump, xsBefore, ysBefore, xsAfter, ysAfter)
                    else:
                        continue
                        
                    corrected_det[:, j, i] = corr_pix

    
    fits.writeto(outfile, corrected_det)
    
    
    # Total time to correct for cosmic rays
    t2 = time.time()
    total = t2 - t1

    print('\n')
    print('---------------------------------------------------')
    print('CORRECTIONS COMPLETE')
    print('Corrections took:', total/60, 'minutes')
    print('---------------------------------------------------')
    

