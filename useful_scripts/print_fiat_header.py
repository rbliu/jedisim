#!/usr/bin/env python
## last modified by Robert Liu at 7/21/2017
import re
import sys
import numpy as np
import pyfits

if len(sys.argv) != 2:
    print >>sys.stderr, "Usage: python print_fiat_header.py {source_fits_file} > {fiat_header_file}"
    exit(1);
srcfits = sys.argv[1]

# Load sources and print column names
if srcfits=='':
	srcfits = 'src.fits'

data_table, header_table = pyfits.getdata('%s' %(srcfits), 1, header=True)

print '# fiat 1.0'

# Print all the flag names
for i in xrange(len(header_table['TFLAG*'])):
    print '# TTYPE' + str(i+1) + ' = ' + header_table['TFLAG*'][i]

# Print all the column names except the 1st column "flags"
for i in xrange(len(header_table['TTYPE*'])-1):
    print '# TTYPE' + str(len(header_table['TFLAG*'])+i+1) + ' = ' + header_table['TTYPE*'][i+1]
