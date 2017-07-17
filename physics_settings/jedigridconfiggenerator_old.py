#! /usr/bin/env python
#last modified 8 July 2013
import os, sys, subprocess, math, re

#parse command line inputs
if(len(sys.argv) != 1):
    print("Usage: jedimaster config")
    sys.exit(1)

set = "production_1_grids/"
#masses = [2,5,10,20]
#concentrations = [4.0]

mass = 10
con = 4

#actual shears = [0.194147, 0.147361, 0.0850851, 0.0577599]
#pixel distance for desired shear
shears = [679.226, 1233.53, 3396.75, 5769.69]
directions=[0., 0.523599, 1.0472, 1.5708, 2.0944, 2.61799]
trials = 5

convolved_path = set
if not os.path.exists(set): 
    os.makedirs(set)

for s,shear in enumerate(shears):
		for d,direct in enumerate(directions):
				lens_file = "%slenses_%i_%i.txt"%(set,s,d)
				f = open("%slenses_%i_%i.txt"%(set,s,d),'w')
				f.write("6144 6144 2 %f %f"%(mass,con))
				f.close()
				for n in xrange(0,trials):
						r = s*len(directions)*trials + d*trials + n
						print "s: %f\td: %f\t n: %i\t r: %i"%(shear,direct,n,r)
						f = open("%sconfig_%i.txt"%(set,r),'w')
						config = """#---------------------physics settings-----------------------
pix_scale=0.03  		#pixel scale to use (arseconds per pixel)
final_pix_scale=0.2		#LSST pixscale (arcsecords per pixel)
exp_time=6000			#exposure time in seconds
noise_mean=10			#mean for poisson noise
nx=12288        		#number of pixels in the x direction
ny=12288        		#number of pixels in the y direction
x_border=301    		#must be large enough so that no image can overflow_
y_border=301			#must be large enough so that no image can overflow
x_trim=0				#larger than x_border to ensure no edge effects
y_trim=0				#larger than y_border to ensure no edge effects
num_galaxies=1000		#number of galaxies to simulate 138,000 default
min_mag=22      		#minimum magnitude galaxy to simulate (inclusive)
max_mag=26      		#maximum magnitude galaxy to simulate (inclusive)
single_redshift=1		#use a single source galaxy redshift? 0 = no, 1=yes
fixed_redshift=1.5		#the single source galaxy redshift to use
power=0.33      		#power for the power law galaxy distribution
lens_z=0.3				#the redshift of the lenses
lenses_file="physics_settings/%s"	#catalog of lenses to use
psf_file="physics_settings/psf_1024.fits"#the PSF to use	
#--------------------output settings--------------------------
output_folder="trial_%i/"
prefix="trial%i_"
HST_image="HST.fits"
HST_convolved_image="HST_convolved.fits"
LSST_convolved_image="LSST_convolved.fits"
LSST_convolved_noise_image="LSST_convolved_noise.fits"
grid_radius=%f
grid_angle=%f
#-----------database folders----------------------------------
#must contain files "n.txt" for n= min_mag to max_mag
radius_db_folder="simdatabase/radius_db/" 
red_db_folder="simdatabase/red_db/"
#-----------catalog file locations-----------------------------
catalog_file="catalog.txt"
dislist_file="dislist.txt"
convlist_file="toconvolvelist.txt"
dislist_grid_file="dislist_grid.txt"
distortedlist_file="distortedlist.txt"
convolvedlist_file="convolvedlist.txt"
#-----------source images-------------------------------------
num_source_images=126
#all postage stamp images should be on their own line, prefaced with image
#postage stamp images should be fits" file including the following header entries:
#MAG: magnitude of the postage stamp image
#MAG0: magnitude zeropoint of the postage stamp image
#PIXSCALE: pixel scale of the postage stamp image
#RADIUS: R50 radius of the image, in pixels
image="simdatabase/doneall2/scaled_finalnew_galaxy_10.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_11.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_12.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_13.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_14.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_15.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_16.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_17.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_18.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_19.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_1.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_20.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_21.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_22.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_23.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_24.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_25.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_26.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_27.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_28.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_29.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_2.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_30.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_31.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_32.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_33.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_34.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_35.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_36.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_37.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_38.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_39.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_3.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_40.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_41.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_42.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_43.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_44.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_45.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_46.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_47.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_48.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_49.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_4.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_50.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_51.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_52.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_53.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_54.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_55.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_56.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_57.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_58.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_59.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_5.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_60.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_61.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_62.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_63.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_64.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_65.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_66.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_67.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_68.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_69.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_6.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_70.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_71.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_72.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_73.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_74.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_75.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_76.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_77.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_78.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_79.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_7.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_80.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_81.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_82.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_83.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_84.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_85.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_86.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_87.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_88.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_89.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_8.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_90.fits"
image="simdatabase/doneall2/scaled_finalnew_galaxy_9.fits"
image="simdatabase/ivydoneall3/0.fits"
image="simdatabase/ivydoneall3/10.fits"
image="simdatabase/ivydoneall3/11.fits"
image="simdatabase/ivydoneall3/12.fits"
image="simdatabase/ivydoneall3/13.fits"
image="simdatabase/ivydoneall3/14.fits"
image="simdatabase/ivydoneall3/15.fits"
image="simdatabase/ivydoneall3/16.fits"
image="simdatabase/ivydoneall3/17.fits"
image="simdatabase/ivydoneall3/18.fits"
image="simdatabase/ivydoneall3/19.fits"
image="simdatabase/ivydoneall3/1.fits"
image="simdatabase/ivydoneall3/20.fits"
image="simdatabase/ivydoneall3/22.fits"
image="simdatabase/ivydoneall3/23.fits"
image="simdatabase/ivydoneall3/25.fits"
image="simdatabase/ivydoneall3/27.fits"
image="simdatabase/ivydoneall3/28.fits"
image="simdatabase/ivydoneall3/29.fits"
image="simdatabase/ivydoneall3/2.fits"
image="simdatabase/ivydoneall3/31.fits"
image="simdatabase/ivydoneall3/32.fits"
image="simdatabase/ivydoneall3/33.fits"
image="simdatabase/ivydoneall3/34.fits"
image="simdatabase/ivydoneall3/39.fits"
image="simdatabase/ivydoneall3/3.fits"
image="simdatabase/ivydoneall3/40.fits"
image="simdatabase/ivydoneall3/41.fits"
image="simdatabase/ivydoneall3/43.fits"
image="simdatabase/ivydoneall3/44.fits"
image="simdatabase/ivydoneall3/45.fits"
image="simdatabase/ivydoneall3/4.fits"
image="simdatabase/ivydoneall3/5.fits"
image="simdatabase/ivydoneall3/6.fits"
image="simdatabase/ivydoneall3/8.fits"
image="simdatabase/ivydoneall3/9.fits"
"""%(lens_file,r,r,shear,direct)
						
						f.write(config)
						f.close()

