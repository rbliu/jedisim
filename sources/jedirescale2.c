#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"


#define NUMBANDS 2

char *help[] = {
    "Scales down a large image to a smaller pixelscale, and trims off a border. Determines values for new pixels by finding the box each pixel makes on the original image, integrating over the area of that box, and averaging. This avoids some of the possible pitfalls of bilinear interpolation: if you have a box over a mostly blank area with a few very bright pixels, the resulting pixel will be very bright. With bilinear interpolation, the resulting pixel may or may not be bright depending on which pixel inside the box was sampled.",
    "usage: jedirescale input_file pixscale new_pixscale trimx trimy output_file",
    "Arguments: input_file - 2D floating point FITS image to rescale",
    "           pixscale - the pixel scale of the input image in arcseconds per pixel",
    "           new_pixscale - the descired pixel scale of the output image in arcseconds per pixel",
    "           trimx - the number of pixels to trim on the left and right before rescaling",
    "           trimy - the number of pixels to trim on the top and bottom before rescaling",
    "           output_file - output FITS file that will be created",
    0};
int main(int argc, char *argv[]){
    //declare variables
    float       pixscale, new_pixscale, scale;  //input paramaters and the scale ratio
    long int    trimx, trimy;                   //input parameters
    char        infile[1024], outfile[1024];    //filenames

    //fits variables
    fitsfile    *infptr, *outfptr;                  //fits pointers for the input and output files
    int         status = 0, naxis = 0;              //cfitsio status counter & image dimensionality
    long int    inaxes[2], bnaxes[2], onaxes[2];    //cfitsio image dimensions array for input and output images
    long int    fpixel[2] = {1,1},bpixel[2] = {1,1};//cfitsio first/last pixel to read in/out
    float       *image, *oimage;                    //arrays for the input and output images


    //print help
    if(argc != 7){
        int line;
        for(line = 0; help[line] !=0; line++)
            fprintf(stderr, "%s\n", help[line]);
        exit(1);
    }

    //parse command line input
    sscanf(argv[1], "%s", infile);
    sscanf(argv[2], "%f", &pixscale);
    sscanf(argv[3], "%f", &new_pixscale);
    sscanf(argv[4], "%li", &trimx);
    sscanf(argv[5], "%li", &trimy);
    sscanf(argv[6], "%s", outfile);
    scale = new_pixscale/pixscale;
    if(scale < 1){
        fprintf(stderr, "Error: you are only allowed to scale down. This program cannot create information out of thin air.");
        exit(1);
    }

    //fprintf(stdout,"%s %f %f %li %li %s\n", infile, pixscale, new_pixscale, trimx, trimy, outfile);


    //open input image
    fits_open_file(&infptr, infile, READONLY, &status);
    fits_get_img_dim(infptr, &naxis, &status);
    fits_get_img_size(infptr, 5, inaxes, &status);
    if(status){
        fits_report_error(stderr,status);
        exit(1);
    }

    //we need to rescale the image in relatively small bands for memory reasons.
    //the dimensions of one band are given by the bnaxes array
    bnaxes[0] = inaxes[0];
    bnaxes[1] = inaxes[1]/NUMBANDS;

    //the output axes are given by the onaxes array
    onaxes[0] = ((float) (inaxes[0]-2*trimx))/scale;
    onaxes[1] = ((float) (inaxes[1]-2*trimy))/scale;

    fprintf(stdout,"input shape: (%li, %li)\tband shape: (%li, %li)\toutput shape: (%li, %li)\n", inaxes[0], inaxes[1], bnaxes[0], bnaxes[1], onaxes[0], onaxes[1]);

    //allocate space for the output image
    oimage = (float *) calloc(onaxes[0]*onaxes[1], sizeof(float));
    if(oimage == NULL){
        fprintf(stderr, "Error: cannot allocate memory for output image.\n");
        exit(1);
    }


    long int    band;   //counter
    for(band = 0; band < NUMBANDS; band++){

        //allocate enough memory for one band of the input image
        image = (float *) calloc (bnaxes[0] * bnaxes[1], sizeof(float));
        if (image == NULL){
            fprintf(stderr, "Error: cannot allocate memory for band %i.\n", band);
            exit(1);
        }
        //fprintf(stdout,"bpixel = {%li, %li}\n", bpixel[0], bpixel[1]);
        
        //read in the image
        fprintf(stdout, "Reading in %s band %i, size: %.0fMB. This may take a few minutes for large images.\n", infile, band, (float) inaxes[0]*inaxes[1]*sizeof(float)/(1000000*NUMBANDS));

        if(fits_read_pix(infptr, TFLOAT, bpixel, bnaxes[0]*bnaxes[1], NULL, image, NULL, &status)){
            fprintf(stderr, "Error: cannot read in input image.\n");
            fits_report_error(stderr, status);
            exit(1);
        }


        long int    rowmin = band*bnaxes[1];        //the first row of this band in the larger image
        long int    rowmax = (band+1)*bnaxes[1];    //the last row of this band in the larger image

        /*This next bit is very inefficient. it checks over all the rows for each band, even when
         * they're way outside the band we care about. This could be rectified by calculating what
         * the iteration bounds on "row" should be carefully. However, since the slowest part of the
         * program by far is the hard drive I/O for reading in large images, I haven't bothered to
         * optimize this part.*/


        fprintf(stdout, "Binning input image.\n");
        long int    row, col;   //counters
        for(col = 0; col < onaxes[0]; col++){
            fprintf(stdout, "Binning col %li/%li of band %li.\n", col, onaxes[0], band);
            for(row = 0; row < onaxes[1]; row++){
                float   total = 0;      //total value for this bin

                //borders of this bin
                float   xmin = trimx + scale*col;
                float   xmax = xmin + scale;
                float   ymin = trimy + scale*row;
                float   ymax = ymin + scale;

                //borders of this bin expanded to the nearest integer 
                long int    ixmin = (long int) floor(xmin);
                long int    ixmax = (long int) ceil(xmax);
                long int    iymin = (long int) floor(ymin);
                long int    iymax = (long int) ceil(ymax);

                //fraction of the edge row or column that is inside the bin
                float xmin_frac = (ixmin + 1) - xmin;
                float xmax_frac = xmax - (ixmax - 1);
                float ymin_frac = (iymin + 1) - ymin;
                float ymax_frac = ymax - (iymax - 1);  
                if(xmin_frac < 0 || xmax_frac < 0 || ymin_frac < 0 || ymax_frac < 0){
                    fprintf(stderr, "error: <0\n");
                    exit(1);
                }

                //loop over all the pixels that are partially inside the bin
                //adding their value time the fraction of the pixel which is inside the bin

                long int srow, scol;
                for(scol = ixmin; scol < ixmax; scol++){
                    for(srow = iymin; srow < iymax; srow++){
                        //only take the pixels inside the current band
                        if(srow >= rowmin && srow < rowmax){
                            //index of the pixel in question in the band
                            float pixel = image[scol+(srow-rowmin)*bnaxes[0]];
                            if(scol == ixmin)
                                pixel *= xmin_frac;
                            else if(scol == ixmax -1)
                                pixel *= xmax_frac;

                            if(srow == iymin)
                                pixel *= ymin_frac;
                            else if(srow == iymax -1)
                                pixel *= ymax_frac;
                            total += pixel;
                        }
                    }
                }
                oimage[col+onaxes[0]*row] += total;
            }
        }
        //release this band
        free(image);
        //and update the first pixel to read in
        bpixel[1] += bnaxes[1];
    }


    //write out the output image
    fprintf(stdout,"Writing output image.\n");
    fits_create_file(&outfptr, outfile, &status);
    fits_create_img(outfptr, FLOAT_IMG, naxis, onaxes, &status);
    fits_write_pix(outfptr, TFLOAT, fpixel, onaxes[0]*onaxes[1], oimage, &status);

    //clean up
    fits_close_file(infptr, &status);
    fits_close_file(outfptr, &status);
    fits_report_error(stderr, status);
    free(oimage);
    return 0;
}
