#!/usr/bin/env python
import re, os, sys
import numpy as np
from astropy import wcs
from astropy.io import fits

if len(sys.argv) != 2:
    print >>sys.stderr, "Usage: python jediWCS.py {/path/to/LSST/simulation} "
    exit(1);
sim = sys.argv[1]

w = wcs.WCS(naxis=2)
# what is the center pixel of the X-Y grid
w.wcs.crpix = [1800.0, 1800.0]
# what is the coordinate of that pixel
w.wcs.crval = [0.1, 0.1]
# what is the pixel scale (deg/pix) in longitude and latitude -- 0.2arcsec/pix
w.wcs.cdelt = np.array([-5.55555555555556E-05,5.55555555555556E-05])
# it is a tangential projection
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

wcsheader = w.to_header()

# write wcs to header
d = fits.open('%s' %(sim), mode='update')
h = d[0].header
h += wcsheader

d.flush()
d.close()
print("Fake WCS infomation added to %s!" %(sim))
