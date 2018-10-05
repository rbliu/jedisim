# -*- coding: utf-8 -*-
#!/usr/bin/env python
## last modified by Robert Liu at 7/3/2018
## This script is used to scale UP images (eg. from 0.03 to 0.004 pixel scale).

import cv2, re, sys
import numpy as np
from astropy.io import fits
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
    print >>sys.stderr, "Usage: python  scaleUp.py  old_image  old_res  new_res  new_image"
    exit(1);
oldFits = sys.argv[1]
oldRes  = sys.argv[2]
newRes  = sys.argv[3]
newFits = sys.argv[4]

sclaeFactor = float(oldRes) / float(newRes)

# Read in the original image
originFits = fits.open(oldFits)

height, width = originFits[0].shape[:2]
print "The original image has the size of " + str(originFits[0].data.shape) + " pixel."

# Extract the central 1168x1168 pixel area
centralArea = originFits[0].data[1464:2632,1464:2632]

# Scale by a fixed factor in x and y axis
outputFits = cv2.resize(centralArea*6000.0, None, fx=sclaeFactor, fy=sclaeFactor,interpolation=cv2.INTER_CUBIC)

# Save the FITS image
hdu = fits.PrimaryHDU(outputFits/6000.0)
hdu.writeto(newFits)

print "The output image has the size of " + str(outputFits.shape) + " pixel."
