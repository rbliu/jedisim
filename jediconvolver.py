#! /usr/bin/env python
#last modified 8 July 2013
import os, sys, subprocess, math, re, shutil,copy

#parse command line inputs
if(len(sys.argv) != 5):
    print("Usage: python jediconvolver.py HST_image_to_convolve convolvedlist_file output_folder psf_file")
    sys.exit(1);

"""
REQUIREMENTS:
The original folder should have a file in it called "*HST.fits" which is to be convolved. The original folder should have a sub-folder called "convolved/", which should be EMPTY, and a file "convolvedlist.txt" which gives the names of the bands for convolution.

"""



config = {}

#name for this file, perhaps trial7, or something. extract from the non-"HST.fits" part of the HST scale image
config['original_folder'] = sys.argv[1]
config['HST_image'] = config['original_folder']+sys.argv[2]
config['output_folder'] = sys.argv[3]
config['psf_file'] = sys.argv[4]


##user-set values
config['nx'] = "12288"
config['ny'] = "12288"
config['pix_scale'] = "0.03"
config['final_pix_scale'] = "0.27"
config['x_trim'] = "0"
config['y_trim'] = "0"
config['exp_time'] = "6000"
config['noise_mean'] = "21.63"
#######end user set values


config['prefix'] = sys.argv[2].split("HST")[0]
config['HST_convolved_image']="HST_convolved.fits"
config['LSST_convolved_image']="DECam_convolved.fits"
config['LSST_convolved_noise_image']="DECam_convolved_noise.fits"


config['convolved_path'] = config['original_folder'] + "convolved/"
config['convolvedlist_file'] = config['original_folder']+config['prefix']+"convolvedlist.txt"


config['HST_convolved_image'] = config['output_folder'] +config['prefix']+config['HST_convolved_image']
config['LSST_convolved_image'] = config['output_folder'] +config['prefix']+config['LSST_convolved_image']
config['LSST_convolved_noise_image'] = config['output_folder'] +config['prefix']+config['LSST_convolved_noise_image']


print config


#function to run a program and write output to the shell
def run_process(name, args,):
    print "------------------------------------------------------------------------"
    print "Running: %s\nCommand:"%name
    for arg in args:
        print arg,
    print ""
    print "------------------------------------------------------------------------"
    process = subprocess.Popen(args)

    process.communicate()
    if process.returncode != 0:
        print "Error: %s did not terminate correctly. Return code: %i."%(name, process.returncode)
        sys.exit(1)



#convonlve the large image with the PSF
#this creates one image for each band of the image
run_process("jediconvolve", ['./jediconvolve', config['HST_image'], config['psf_file'], config['convolved_path']])

#combine each band into a single image
run_process("jedipaste (make convolved image)", ['./jedipaste', config['nx'], config['ny'], config['convolvedlist_file'], config['HST_convolved_image']])

#scale the image down from HST to LSST scale and trim the edges
run_process("jedirescale", ['./jedirescale', config['HST_convolved_image'], config['pix_scale'], config['final_pix_scale'], config['x_trim'], config['y_trim'], config['LSST_convolved_image']])

#simulate exposure time and add Poisson noise
run_process("jedinoise", ['./jedinoise', config['LSST_convolved_image'], config['exp_time'], config['noise_mean'], config['LSST_convolved_noise_image']])

print "jediconvolver successful! Exiting."