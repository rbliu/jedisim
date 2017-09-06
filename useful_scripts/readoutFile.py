#!/usr/bin/env python
## last modified by Robert Liu at 9/30/2016
import re
import sys
import numpy as np
import pyfits

if len(sys.argv) != 2:
    print >>sys.stderr, "Usage: readoutFile {source_fits_file} > {catalog_file}"
    exit(1);
srcfits = sys.argv[1]

# Load sources and print interesting columns
if srcfits=='':
	srcfits = 'src.fits'

data_table, header_table = pyfits.getdata('%s' %(srcfits), 1, header=True)

cols = ("id",
        "coord_ra",
        "coord_dec",
        "ext_shapeHSM_HsmSourceMoments_x",
        "ext_shapeHSM_HsmSourceMoments_y",
        "ext_shapeHSM_HsmSourceMoments_xx",
        "ext_shapeHSM_HsmSourceMoments_yy",
        "ext_shapeHSM_HsmSourceMoments_xy",
        "ext_shapeHSM_HsmShapeRegauss_e1",
        "ext_shapeHSM_HsmShapeRegauss_e2",
        "ext_shapeHSM_HsmShapeRegauss_sigma",
        "base_ClassificationExtendedness_value",
        )

print '#' + ' '.join(cols)

vecs = []
for col in cols:
    if col.endswith(".ra") or col.endswith(".dec") or col.endswith("_ra") or col.endswith("_dec"):
        v = np.rad2deg(data_table.field(col))

    else:
        v = data_table.field(col)

    vecs.append(v)

for vals in zip(*vecs):
    print ' '.join([str(el) for el in vals])
