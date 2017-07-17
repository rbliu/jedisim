#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

/*to compile, run:
 *gcc -O3 jeditransform.c -o jedi_imadjust -lm -lcfitsio
 */

#define PI  3.14159         //because C doesn't have it built-in for some silly reason
#define EPSILON 0.0000001   //a very small number
#define MAG0 30.0           //magnitude zero point

char *help[] = {
    "Takes in a catalog of images, and produces a FITS image for each entry, transformed to the correct specifications.",
    "Usage: jeditransform catalog distort_list",
    "Arguments: catalog - text file containing galaxy catalog",
    "           distort_list - file path to output instructions for jedidistort",
    "catalog file:  image x y angle redshift old_mag old_r50 new_mag new_r50 stamp1 stamp2 [tab separated]",
    "               image - file path for the base galaxy postage stamp image",
    "               x - x coordinate for the image center",
    "               y - y coordinate for the image center",
    "               angle - angle through which to rotate the input postage stamp",
    "               redshift - redshift for the galaxy",
    "               pixscale - pixel scale for the galaxy (arcseconds per pixel)",
    "               old_mag - magnitude of the base galaxy postage stamp image",
    "               old_r50 - r50-type radius of the base galaxy postage stamp image",
    "               new_mag - magnitude that the galaxy should have",
    "               new_r50 - r50-type radius that the galaxy should have",
    "               stamp1 - filepath for the galaxy's rotated, scaled, and photoscaled, postage stamp",
    "               stamp2 - filepath for the final distorted postage stamp",
    0};


float bilinear_interp(float x, float y, int xmax, int ymax, float *img);

//all the information needed to describe a single galaxy
typedef struct{
    char    *image;     //file path for the base galaxy postage stamp image",
    float   x;          //x coordinate for the image center",
    float   y;          //y coordinate for the image center",
    float   angle;      //angle through which to rotate the input postage stamp",
    float   redshift;   //redshift for the galaxy
    float   pixscale;   //pixelscale for the galaxy
    float   old_mag;    //magnitude of the base galaxy postage stamp image",
    float   old_r50;    //r50-type radius of the base galaxy postage stamp image",
    float   new_mag;    //magnitude that the galaxy should have",
    float   new_r50;    //r50-type radius that the galaxy should have",
    char    *stamp1;    //filepath for the galaxy's postage rotated, scaled, and photoscaled postage stamp",
    char    *stamp2;    //filepath for the final distorted postage stamp",
} galaxy;


int main(int argc, char *argv[]){
    //declare variables
    int     ngalaxies;          //the number of galaxies  
    char    gallist_path[1024]; //file path for the list of galaxies
    char    dislist_path[1024]; //file path for the distortions list
    FILE    *gallist_file;      //galaxy list file
    FILE    *dislist_file;      //distortion list file
    char    buffer1[1024], buffer2[1024], buffer3[1024]; //buffers for reading in strings
    galaxy  *galaxies;          //array of galaxies to read gallist into


    //print help
    if(argc != 3){
        int line;
        for(line = 0; help[line] !=0; line++)
            fprintf(stderr, "%s\n", help[line]);
        exit(1);
    }

    //parse command line input
    sscanf(argv[1], "%[^\t\n]", &gallist_path);
    sscanf(argv[2], "%[^\t\n]", &dislist_path);
    fprintf(stdout, "gallist_path: %s\n", gallist_path);

    //parse gallist file
    if((gallist_file = fopen(gallist_path, "r")) == NULL){
        fprintf(stderr,"Error: could not open galaxy list file \"%s\"\n", gallist_path);
        exit(1);
    }
    fprintf(stdout,"Opened galaxy list file \"%s\".\n", gallist_path);

    //open dislist file
    if((dislist_file = fopen(dislist_path, "w+")) == NULL){
        fprintf(stderr,"Error: could not open distortion instruction list file \"%s\"\n", dislist_path);
        exit(1);
    }
    fprintf(stdout,"Opened galaxy list file \"%s\".\n", gallist_path);

    //first count the number of galaxies
    ngalaxies = 0;
    galaxy temp;
    while(fscanf(gallist_file, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s", &buffer1, &temp.x, &temp.y, &temp.angle, &temp.redshift, &temp.pixscale, &temp.old_mag, &temp.old_r50, &temp.new_mag, &temp.new_r50, &buffer2, &buffer3) > 0)
        ngalaxies++;

    if(ngalaxies ==0){
        fprintf(stderr,"Error: found 0 galaxies in galaxy list \"%s\"\n", gallist_path);
        exit(1);
    }
    fprintf(stdout,"Found %i galaxies in galaxy list \"%s\"\n", ngalaxies, gallist_path);

    //allocate memory for the galaxies
    galaxies = (galaxy *) calloc(ngalaxies, sizeof(galaxy));
    if(galaxies == NULL){
        fprintf(stderr, "Error: could not allocate memory for galaxy list.");
        exit(1);
    }

    //read in the galaxies
    rewind(gallist_file);
    int     g;  //galaxies counter
    for(g = 0; g < ngalaxies; g++){
        fscanf(gallist_file, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s", &buffer1, &galaxies[g].x, &galaxies[g].y, &galaxies[g].angle, &galaxies[g].redshift, &galaxies[g].pixscale, &galaxies[g].old_mag, &galaxies[g].old_r50, &galaxies[g].new_mag, &galaxies[g].new_r50, &buffer2, &buffer3);
        galaxies[g].image = (char *) calloc(strlen(buffer1), sizeof(char));
        galaxies[g].stamp1 = (char *) calloc(strlen(buffer2)+1, sizeof(char));
        galaxies[g].stamp2 = (char *) calloc(strlen(buffer3)+1, sizeof(char));
        if(galaxies[g].image == NULL || galaxies[g].stamp1 == NULL || galaxies[g].stamp2 == NULL){
            fprintf(stderr,"Error: could not allocate memory for galaxy  %i.\n", g);
            exit(1);
        }
        strcpy(galaxies[g].image, buffer1);
        strcpy(galaxies[g].stamp1, buffer2);
        strcpy(galaxies[g].stamp2, buffer3);
        //fprintf(stdout,"Galaxy %i: %s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n",g, galaxies[g].image, galaxies[g].x, galaxies[g].y, galaxies[g].angle, galaxies[g].redshift, galaxies[g].pixscale, galaxies[g].old_mag, galaxies[g].old_r50, galaxies[g].new_mag, galaxies[g].new_r50, galaxies[g].stamp1, galaxies[g].stamp2);
    }
    fclose(gallist_file);
    fprintf(stdout,"Read in %i galaxies.\n", ngalaxies);



    //loop over each image, rotating, scaling, photoscaling, and "cutting out" each one
    for(g = 0; g < ngalaxies; g++){
        fitsfile    *galfptr, *outfptr;             //fits pointers for the input and output files
        int         status = 0, naxis = 0;          //cfitsio status counter & image dimensionality
        long        galnaxes[2], tgalnaxes[2];      //cfitsio image dimensions array
        long        fpixel[2] = {1,1}, lpixel[2];   //cfitsio first/last pixel to read in/out
        float       *galaxy, *tgal, *pgal;          //arrays for galaxy input image, transformed image, postage transformed image

        //open input image
        fits_open_file(&galfptr, galaxies[g].image, READONLY, &status);
        fits_get_img_dim(galfptr, &naxis, &status);
        fits_get_img_size(galfptr, 5, galnaxes, &status);

        //get image dimensions
        if(status){
            fits_report_error(stderr,status);
            exit(1);
        }

        //allocate enough memory for the galaxy and transformed galaxy images
        //fprintf(stdout, "Allocating memory for galaxy %i image.\n", g);
        galaxy = (float *) calloc (galnaxes[0] * galnaxes[1], sizeof(float));
        tgal = (float *) calloc (galnaxes[0] * galnaxes[1], sizeof(float));
        if (galaxy == NULL){
            fprintf(stderr, "Error: cannot allocate memory for galaxy %i image.\n", g);
            exit(1);
        }

        //read in the image
        //fprintf(stdout, "Reading in image.\n");
        if(fits_read_pix(galfptr, TFLOAT, fpixel, galnaxes[0]*galnaxes[1], NULL, galaxy, NULL, &status)){
            fprintf(stderr, "Error: cannot read in galaxy %i image.\n", g);
            fits_report_error(stderr, status);
            exit(1);
        }
        fits_close_file(galfptr, &status);

        //find magnifications
        float   scale = galaxies[g].old_r50/galaxies[g].new_r50;

        float c = cosf(galaxies[g].angle*PI/180);   //pre-compute cos
        float s = sinf(galaxies[g].angle*PI/180);   //and sin



        //variables for image transformation
        long int    xmin = galnaxes[0], xmax = -1, ymin = galnaxes[1], ymax = -1;   //information bounding box in transormed image
        long int    xc = 1+(galnaxes[0] >> 1), yc = 1+(galnaxes[0] >> 1);           //x and y center pixels
        float       newx, newy, oldx, oldy;                                         //x and y positions relative to the image center, pre- and post- transform
        long int    x, y;                                                           //counters

        float total = 0;        //the total value of all the pixels in the image
        //transform each pixel
        for(x = 0; x < galnaxes[0]; x++){
            for(y = 0; y < galnaxes[1]; y++){
                //find pixel relative to center
                oldx = (float) (x - xc);
                oldy = (float) (y - yc);

                //rotation and dilation matrix
                newx = scale * (oldx*c - oldy*s) + xc;
                newy = scale * (oldx*s + oldy*c) + yc;


                //do interpolation and find new magnitude
                float value = bilinear_interp(newx, newy, galnaxes[0], galnaxes[1], galaxy);
                if(value > EPSILON){
                    if(x < xmin)
                        xmin = x;
                    if(x > xmax)
                        xmax = x;
                    if(y < ymin)
                        ymin = y;
                    if(y > ymax)
                        ymax = y;
                }
                tgal[y*galnaxes[1]+x] = value;
                total += value;
            }
        }
        float mag = MAG0 - 2.5*log10f(total);
        float photoscale = pow(10.0,(mag - galaxies[g].new_mag)/2.5);

        //find the size of the smallest possible postage stamp needed
        tgalnaxes[0] = (xmax+1)-xmin;
        tgalnaxes[1] = (ymax+1)-ymin;

        //make smallest possible postage stamp galaxy, and adjust brightness to desired magnitude
        pgal = (float *) calloc(tgalnaxes[0]*tgalnaxes[1], sizeof(float));
        for(x = xmin; x < xmax + 1; x++)
            for(y = ymin; y < ymax + 1; y++)
                pgal[(y-ymin)*tgalnaxes[0] + (x-xmin)] = photoscale*tgal[y*galnaxes[1] + x];

        //specify the coordinates where the lower left-most pixel of this image should be embedded in the galaxy field
        long int    xembed = (long int) (0.5 + galaxies[g].x - (1+(float) tgalnaxes[0])/2);
        long int    yembed = (long int) (0.5 + galaxies[g].y - (1+(float) tgalnaxes[1])/2);


        //create and write output image
        fits_create_file(&outfptr, galaxies[g].stamp1, &status);
        fits_create_img(outfptr, FLOAT_IMG, naxis, tgalnaxes, &status);
        fits_update_key(outfptr, TLONG, "XEMBED", &xembed, "x pixel in targe image to embed lower left pixel.", &status);	
        fits_update_key(outfptr, TLONG, "YEMBED", &yembed, "y pixel in targe image to embed lower left pixel.", &status);
		fits_update_key(outfptr,TSTRING,"IMAGE NAME", galaxies[g].image,"",&status);	
        fits_write_pix(outfptr, TFLOAT, fpixel, tgalnaxes[0]*tgalnaxes[1], pgal, &status);

        //close the FITS file
        fits_close_file(outfptr, &status);
        fits_report_error(stderr, status);
        fprintf(stdout, "Transformed image %i.\n", g);

        //make a list of all the transformed galaxies so they can be distorted
        fprintf(dislist_file, "%li %li %li %li %f %s %s\n", xembed, yembed, tgalnaxes[0], tgalnaxes[1], galaxies[g].redshift, galaxies[g].stamp1, galaxies[g].stamp2);

        free(galaxy);
        free(tgal);
        free(pgal);
    }
    fclose(dislist_file);
    return 0;
}


//given a floating point location in an image, returns its pixel value by doing a bilinear interpolation
float bilinear_interp(float x, float y, int xmax, int ymax, float *img){
    int xi = (int) x, yi = (int) y;
    if(xi >= 0 && xi < xmax-1 && yi >= 0 && yi < ymax-1){
        float a = img[xi + ymax*yi];            //f(xi,yi)
        float b = img[xi+1+ymax*yi];      //f(xi+1,yi)
        float c = img[xi + ymax*(yi+1)];            //f(xi, yi+1)
        float d = img[xi+1+ymax*(yi+1)];      //f(xi+1, yi+1)
        float xf = x -xi, yf = y-yi;
        return a*(1-xf)*(1-yf) + b*xf*(1-yf) + c*(1-xf)*yf + d*xf*yf;
    } else
        return 0;        
}
