#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include <time.h>

#define NUMBANDS 2

char *help[] = {
    "Takes a list of 2D floating point FITS images and combines them into a single large image.",
    "usage: jedipaste x y imlist output_image",
    "Arguments: x - width of the final image",
    "           y - height of the final image",
    "           imlist - list of paths of images to insert in to the final image",
    "           output_image - path for the output image. Will overwrite other files.",
    "imlist description: text file of filepaths to valid FITS images to be inserted",
    "into the final image. There must be one filepath per line. The header of each",
    "FITS file must have the entries XEMBED and YEMBED, with integer values giving the",
    "x and y pixel where the lower left corner of the embedded image should be placed.",
    0};

int main(int argc, char *argv[]){

    //check command line input
    if(argc != 5){
        int line;
        for(line = 0; help[line] != 0; line++)
            fprintf(stderr, "%s\n", help[line]);
        exit(1);
    }

    FILE        *imlist_file;       //file object for the image list
    char        **imlist;           //array for the images to embed
    long int    nimages;            //number of images int embed
    char        buffer[1024];       //string buffer for reading imlist
    float       *fimage, *image;    //arrays for the final image and the current embedded image
    fitsfile    *ffptr, *efptr;     //final and embed CFITSIO file pointers
    int         status = 0,naxis=2; //CFITSIO status and number of axes parameters
    long        naxes[2], fnaxes[2], onaxes[2], fpixel[2] = {1,1}, opixel[2] = {1,1};    //CFITSIO axes lengths and first pixel to read in



    //parse command line arguments
    sscanf(argv[1], "%li", &fnaxes[0]);
    sscanf(argv[2], "%li", &fnaxes[1]);

    //dimensions of one band of the output image
    onaxes[0] = fnaxes[0];
    onaxes[1] = fnaxes[1]/NUMBANDS;
	
    fprintf(stdout,"(%li, %li)\n", fnaxes[0], fnaxes[1]);
    fprintf(stdout,"(%li, %li)\n", onaxes[0], onaxes[1]);

    //open imlist
    if((imlist_file = fopen(argv[3],"r"))==NULL){
        fprintf(stderr,"Error: could not open image list file \"%s\".\n",argv[3]);
        exit(1);
    }

    //count the number of files in imlist
    nimages = 0;
    while(fscanf(imlist_file, "%s\n", &buffer) >0)
        nimages++;
    fprintf(stdout,"%li images in \"%s\".\n", nimages, argv[3]);

    //initialize memory for the image names array
    imlist = (char **) calloc(nimages, sizeof(char*));
    if(imlist == NULL){
        fprintf(stdout, "Error: could not allocate memory for the list of image names.");
        exit(1);
    }

    //read in image names
    rewind(imlist_file);
    long int    im;     //image counter
    for(im = 0; im < nimages; im++){
        fscanf(imlist_file, "%s\n", &buffer);
        imlist[im] = (char *) calloc(strlen(buffer)+1, sizeof(char));
        if(imlist[im] == NULL){
            fprintf(stderr,"Error: could not allocate memory for image %li name.", im);
            exit(1);
        }
        strcpy(imlist[im], buffer);
    }

    fits_create_file(&ffptr, argv[4], &status);
    fits_create_img(ffptr, FLOAT_IMG, naxis, fnaxes, &status);


    int     band; //counter
    for(band = 0; band < NUMBANDS; band++){

        //we only want the parts of the images inside this band,
        //so we need the borders of this band
        long int    ymin = band*onaxes[1];
        long int    ymax = (band+1)*onaxes[1];
        long int    xmin = 0;
        long int    xmax = onaxes[0];


        fprintf(stdout,"(%li, %li)\n", opixel[0], opixel[1]);
        
        //make an array for this band
        fprintf(stdout, "Creating band %i array.\n", band);
        fimage = (float*) calloc(onaxes[0]*onaxes[1], sizeof(float));
        if(fimage == NULL){
            fprintf(stdout,"Error: could not allocate memory for the final image.");
            exit(1);
        }

        //loop over each image and embed them one at a time
        for(im = 0; im < nimages; im++){

            //read in galaxy image
            fits_open_file(&efptr, imlist[im], READONLY, &status);
            fits_get_img_dim(efptr, &naxis, &status);
            fits_get_img_size(efptr, 2, naxes, &status);
            if(status){
                fits_report_error(stderr, status);
                exit(1);
            }

			//fprintf(stdout,"naxes[0]: %li\nnaxes[1]: %li\n",naxes[0],naxes[1]);

            //read in where the embed this image
            char        xembedstr[32], yembedstr[32];   //strings of the embed coordinates
            long int    xembed, yembed;                 //integer embed coordinates
            fits_read_key_str(efptr, "XEMBED", xembedstr, NULL, &status);
            fits_read_key_str(efptr, "YEMBED", yembedstr, NULL, &status);
            sscanf(xembedstr, "%li", &xembed);
            sscanf(yembedstr, "%li", &yembed);



            //check if the image intersects this band
            if(yembed < ymax && yembed+naxes[1] >= ymin){
                //allocate enough memory for the image
                image = (float *) calloc (naxes[0]*naxes[1], sizeof(float));
                if (image == NULL){
                    fprintf(stderr, "Error allocating memory for galaxy image %li.\n", im);
                    exit(1);
                }

                //read in the image
                //fprintf(stdout, "Reading in image.\n");
                if(fits_read_pix(efptr, TFLOAT, fpixel, naxes[0]*naxes[1], NULL, image, NULL, &status)){
                    fprintf(stderr, "Can't read in embed image %li.\n", im);
                    fits_report_error(stderr, status);
                    exit(1);
                }
                fits_close_file(efptr, &status);

                //find the intersection
                //rowmin lowest row in im that's inside this band
                long int    rowmin = (yembed >= ymin ? 0 : ymin-yembed);
                //rowmax: highest row from im that's inside this band
                long int    rowmax = (yembed+naxes[1] < ymax ? naxes[1] : ymax-yembed);

                //same for columns
                long int    colmin = (xembed >= xmin ? 0 : xmin-xembed);
                long int    colmax = (xembed+naxes[0] < xmax ? naxes[0] : xmax-xembed);

                fprintf(stdout,"Embedding image %li in band %i at (%li, %li).\n", im, band, xembed, yembed);
                long int    row, col;   //counter to embed the image
                for(col = colmin; col < colmax; col++){
                    for(row = rowmin; row < rowmax; row++){
                        //fprintf(stderr, "index: %li\nmaxindex: %li\nrowmin: %li\nrowmax: %i\nymin: %li\nymax: %li\nxembed: %li\nymembed: %li\n(rowmax+yembed): %li\n", (row+yembed)*onaxes[0] + (col+xembed), onaxes[0]*onaxes[1], rowmin, rowmax, ymin, ymax, xembed, yembed, rowmax+yembed);
                        fimage[(col+xembed) + onaxes[0]*(row+yembed-ymin)] += image[col+naxes[0]*row];

                    } 

                }
                //release memory
                free(image);
            } else {
                //if it doesn't, move onto the next image
                fits_close_file(efptr, &status);
            } 
        }


		

        //write the band to disk
        fprintf(stdout,"Writing band %i in image %s. This may take several minutes.\n", band, argv[4]);
        fits_write_pix(ffptr, TFLOAT, opixel, onaxes[0]*onaxes[1], fimage, &status);

        //get ready for the next band
        opixel[1] += fnaxes[1]/NUMBANDS;
        free(fimage);
    }

    //close the output file
    fits_close_file(ffptr, &status);
    fits_report_error(stderr, status);

    return 0;
}
