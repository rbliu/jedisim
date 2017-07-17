#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"


char *help[] = {
    "jedinoise simulates an exposure time by multiplying image intensities and adds poisson noise afterwards.",
    "usage: jedinoise input_file exp_time noise_mean output_file",
    "Arguments: input_file - 2D floating point image in FITS format to which noise will be added, interpretted as a 1 second exposure",
    "           exp_time - the exposure time for the image in seconds(multiplicative factor for the image before noise)",
    "           noise_mean - mean amount of Poisson noise, must be positive, not necessarily an integer, 0 -> no noise",
    "           output_file - output FITS file that will be created",
    0};


int main(int argc, char *argv[]){
    //declare variables
    float exp_time, noise_mean;             //the input parameters
    char infile[1024], outfile[1024];       //input filepaths

    //fits variables
    fitsfile    *infptr, *outfptr;          //fits pointers for the input and output files
    int         status = 0, naxis = 0;      //cfitsio status counter & image dimensionality
    long        inaxes[2];                  //cfitsio image dimensions array for iamges
    long        fpixel[2] = {1,1};          //cfitsio first pixel to read in/out
    float       *image;                     //array for the input and output images
    float       *CDF;                       //continuous distribution function 
    long int    CDF_bins;                   //number of bins for the CDF array

    //print help
    if(argc != 5){
        int line;   //counter
        for(line = 0; help[line] !=0; line++){
            fprintf(stderr, "%s\n", help[line]);
        }
        exit(1);
    }

    //parse command line input
    sscanf(argv[1], "%s", infile);
    sscanf(argv[2], "%f", &exp_time);
    sscanf(argv[3], "%f", &noise_mean);
    sscanf(argv[4], "%s", outfile);
    CDF_bins = (long int) 10*noise_mean;
    fprintf(stdout,"%s %f %f %s\n", infile, exp_time, noise_mean, outfile);
    if(noise_mean < 0){
        fprintf(stdout,"The mean noise must be positive. The mean noise is set as %f.\n", noise_mean);
        exit(1);
    }


    //open input image
    fits_open_file(&infptr, infile, READONLY, &status);
    fits_get_img_dim(infptr, &naxis, &status);
    fits_get_img_size(infptr, 2, inaxes, &status);
    if(status){
        fits_report_error(stderr,status);
        exit(1);
    }

    //allocate enough memory for the input image
    image = (float *) calloc (inaxes[0] * inaxes[1], sizeof(float));
    if (image == NULL){
        fprintf(stderr, "Error: cannot allocate memory for input image.\n");
        exit(1);
    }

    //read in the image
    //fprintf(stdout, "Reading in %s, size: %.0fMB. This may take a few minutes for large images.\n", infile, (float) inaxes[0]*inaxes[1]*sizeof(float)/1000000);
    if(fits_read_pix(infptr, TFLOAT, fpixel, inaxes[0]*inaxes[1], NULL, image, NULL, &status)){
        fprintf(stderr, "Error: cannot read in input image.\n");
        fits_report_error(stderr, status);
        exit(1);
    }
    fits_close_file(infptr, &status);


    //make CDF for the poisson noise
    CDF = (float *) calloc(CDF_bins, sizeof(float));
    if(CDF == NULL){
        fprintf(stderr, "Error: could not allocate memory for CDF.\n");
        exit(1);
    }
    {
        long int    i;                              //counter
        float       prefactor = expf(-noise_mean);  //exponentials are slow, so we only want to calculate it once
        float       partial_sum = 0;                //start the CDF at zero
        float       factorial;                      //the fastest way to calculate n! is to do it one multiplication at a time
        for(i = 0; i < CDF_bins; i++){
            factorial = (i == 0 ? 1 : factorial*i);     
            partial_sum += pow(noise_mean, i)/factorial;
            CDF[i] = prefactor*partial_sum;
            //fprintf(stdout,"CDF[%li]: %f\n", i, CDF[i]);
        }
    }


    fprintf(stdout, "Adding noise to input image.\n");
    long int    row, col;   //counters
    for(row = 0; row < inaxes[0]; row++){
        fprintf(stdout, "Noise-ing row %li/%li.\n", row, inaxes[0]);
        for(col = 0; col < inaxes[1]; col++){
            long int    index = row*inaxes[1]+col;      //index of this pixel in the image array  
            float       r = ((float) rand())/RAND_MAX;  //random integer between 0 and 1
            long int    j = 0;                          //counter
            //find the first bin above the random number, or get the last bin
            //while(CDF[j] <= r && r < CDF_bins)
            while(CDF[j] <= r && j < CDF_bins)
                j++;
            image[index] = exp_time*image[index]+exp_time*j;
        }
    }



    //write out the output image
    fprintf(stdout,"Writing output image.\n");
    fits_create_file(&outfptr, outfile, &status);
    fits_create_img(outfptr, FLOAT_IMG, naxis, inaxes, &status);
    fits_write_pix(outfptr, TFLOAT, fpixel, inaxes[0]*inaxes[1], image, &status);
    fits_close_file(outfptr, &status);
    fits_report_error(stderr, status);

    //release memory
    free(image);

    return 0;
}

