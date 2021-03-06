![logo](https://github.com/rbliu/Memo-astro/blob/master/images/Brown_physics.jpg "logo of Brown Physics")

# jedisim

*jedisim* is a tool to generate galaxy cluster simulations from real galaxy observational images.

The goal of *jedisim* is to generate image simulations for galaxy clusters as a part of [ARCLETS (Analysis of Realistic Cluster Lensing through Extensive Training Simulations)](http://www.het.brown.edu/people/ian/ClustersChallenge/). It is also used to test and calibrate the shape measurements and deblanding in [LSST DM stack](https://pipelines.lsst.io/index.html), as well as mass measurements in [CFHT Clusters pipeline](https://github.com/nicolaschotard/Clusters).

> **Copyright**: jedisim is property of Ian Dell'Antonio's group. The original source code was developed by *Ian Dell'Antonio* and *Daniel Parker*, etc., maintained and updated by *Robert Binyang Liu* and *Shenming Fu*. For more details about the copyright, please refer to the [LICENSE](https://github.com/rbliu/jedisim/blob/master/LICENSE) file. <br>
> **Claim**: Due to the limitation of single file size on github, we cannot upload all the galaxy postage stamps here. We will upload them somewhere else.

**Updates:**

*2019-6-28*: Updated the Copyright section and the LICENSE file.

*2018-10-31*: *jedisim6.0* has created the new `jedidistortDC2` code to ingest DC2 alpha map as a single lens plane.

*2018-02-28*: *jedisim5.4* has updated `jedinoise` code to generate sky noise with a Gaussian distribution. The config parameter `noise_mean` is now replaced with `noise_std`. 

*2017-08-18*: *jedisim5.3* has scripts to add fake WCS information to the simulation images (at LSST pixel scale).

*2016-07-17*: *jedisim5.2* has added 201 images from *CANDELS GOODS-S sect23 F814W* corresponding to the *F606W* images. In the distortion code, type=3 generates the distortion table with a single shear, thus we suggest to use type=2 in grid simulations.

*2016-02-22*: *jedisim5.1* has included 201 images from *CANDELS GOODS-S sect23 F606W*.

-----------

## 1. Introduction and implementation

* Please refer to the documentations:

    * jedisim/documentation/jedisim4_jedisim.pdf

    * jedisim/documentation/jedisim4_lensing.pdf


## 2. Configuration

Modify your settings in *jedisim/physics_settings/* 

(The config files can be stored anywhere you like. By default, I put them under this folder. If you have them in different directories or with different file names, make sure you modify the simulation command in the next section accordingly).

### Lens parameter

* Modify the lens parameter file "lens*.txt".

* Modify the file name of "lens*.txt" in "config" or "config_grid" correspondingly.

* Format of lens parameter file "lens*.txt": x y type p1 p2
    * where x - x center of lens (in pixels)
        * y - y center of lens (in pixels)
        * type - type of mass profile

            1 = Singular isothermal sphere

            2 = Navarro-Frenk-White profile

            3 = NFW constant distortion profile for grid simulations

        * p1 - first profile parameter

            1 = sigma_v [km/s]

            2 = M200 parameter [10^14 solar masses]

            3 = Distance to center in pixels. M200 fixed at 20 default, which can be modified in case 3

        * p2 - second profile parameter

            1 = not applicable, can take any numerical

            2 = concentration parameter [dimentionless]

            3 = concentration parameter [dimentionless]

    * We recommend to use type=2 (NFW profile) for cluster simulations, and type=1 (SIS profile) for grid simulations.

### Configuration file

* Modify the configuration file "config" for cluster simulation or "config_grid" for grid simulation.

* Include the PSF file into the "physics_settings/" folder, and specify the PSF file name in the configuration file.

* The simulation output is rescaled from HST pixel scale to the final pixel scale set in the configuration file.

* For grid simulations:
    * Set the *grid_radius (in pixels)* and *grid_angle (in radian)* in the configuration file.
    * Set *sigma_v (in km/s)* in the lens parameter file.
    * Make sure the *grid_radius* is smaller than half of the image size.


## 3. Usage

Under "jedisim" directory:
```
make
python jedimaster.py physics_settings/config
```
Please note that you may need to modify the `makefile` according to your dependency library path.

## 4. Output

* Specify the directory of output images (eg. "trial_1/") in the output settings section of "config" file.

* Specify the output prefix (eg. "trial1_") in the output settings section of "config" file.

* There are corresponding simulation output in which the galaxies are 90 degree rotated in output directory (eg. "90_trial_1/").

## 5. WCS

* By default, cluster and grid simulations have no WCS information in fits headers. To add ra-dec coordinates, run
```
python jediWCS.py /path/to/LSST/simulation
```

## 6. Others

* The original images and radius/redshift data are in "jedisim/simdatabase/".

* The source codes are in "jedisim/sources/".

* Other parameters can be modified in "jedisim/sources/jedidistort.c" (eg. OMEGA_M, G, H0).

* If any modifications are made in "jedisim/source/" folder, remove the UNIX executable files in "jedisim/" (jedicatalog, jediconvolve, jedidistort, jedigrid_a, jedigrid_b, jedinoise, jedipaste, jedirescale, jeditransform) and regenerate them.
