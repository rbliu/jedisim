#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include <time.h>

#define NUMBANDS 2

char *help[] = {
	"Takes a list of 2D floating point FITS images and places them in a grid pattern in a single larger image.",
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
	fitsfile    *ffptr, *efptr;     //final and embed CFITSIO file pointers
	int         status = 0,naxis=2; //CFITSIO status and number of axes parameters
	long        naxes[2], fnaxes[2];    //CFITSIO axes lengths and first pixel to read in



	//parse command line arguments
	sscanf(argv[1], "%li", &fnaxes[0]);
	sscanf(argv[2], "%li", &fnaxes[1]);



	fprintf(stdout,"(%li, %li)\n", fnaxes[0], fnaxes[1]);
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


	long int max_width = 0, max_height=0;

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

		//fprintf(stdout,"naxes[0]: %li naxes[1]: %li\n",naxes[0],naxes[1]);

		//read in where the embed this image
		char        xembedstr[32], yembedstr[32];   //strings of the embed coordinates
		long int    xembed, yembed;                 //integer embed coordinates
		fits_read_key_str(efptr, "XEMBED", xembedstr, NULL, &status);
		fits_read_key_str(efptr, "YEMBED", yembedstr, NULL, &status);
		sscanf(xembedstr, "%li", &xembed);
		sscanf(yembedstr, "%li", &yembed);

		//calculate maximum rectangle
		//fprintf(stdout,"%li: size:(%li, %li); embed:(%li, %li)\n", im, naxes[0], naxes[1], xembed, yembed);
		if(naxes[0] > max_width)
			max_width = naxes[0];
		if(naxes[1] > max_height)
			max_height = naxes[1];


		fits_close_file(efptr,&status);
		fits_report_error(stderr, status);
	}

	fprintf(stdout,"Max size: (%li, %li).\n", max_width, max_height);

	//new embed coordinates
	long int    nxembed = 0, nyembed = 0; 
	long int	border = 150, border2 = 480;

	long int grid_x = 0,grid_y = 0;

	for(im = 0; im < nimages; im++){

		//read in galaxy image
		fits_open_file(&efptr, imlist[im], READWRITE, &status);
		fits_get_img_dim(efptr, &naxis, &status);
		fits_get_img_size(efptr, 2, naxes, &status);

		//read in where the embed this image
		char        xembedstr[32], yembedstr[32];   //strings of the embed coordinates
		long int    xembed, yembed;                 //integer embed coordinates
		fits_read_key_str(efptr, "XEMBED", xembedstr, NULL, &status);
		fits_read_key_str(efptr, "YEMBED", yembedstr, NULL, &status);
		sscanf(xembedstr, "%li", &xembed);
		sscanf(yembedstr, "%li", &yembed);

		if(status){
			fits_report_error(stderr, status);
			exit(1);
		}


		//check if we're at the end of a row
		if(grid_x*(border+max_width) + 2*border2 > fnaxes[0]){
			grid_y++;
			grid_x = 0;
		}

		//check if we're at the top
		if(grid_y*(border+max_height) + 2*border2 > fnaxes[1]){
			//if we are at the top, we don't want any image, so push offscreen
			nxembed = -1000;
			nyembed= -1000;

		} else {

			nxembed = border2+grid_x*(border+max_width) + (max_width/2)-(naxes[0]/2);
			nyembed = border2+grid_y*(border+max_height) + (max_height/2)-(naxes[1]/2);
		}
			fprintf(stdout, "im %li: grid(%li, %li) at (%li,%li).\n", im, grid_x, grid_y, nxembed, nyembed);
	
		


			fits_update_key(efptr, TLONG, "XEMBED", &nxembed , "x pixel in target image to embed lower left pixel.", &status);
			fits_update_key(efptr, TLONG, "YEMBED", &nyembed , "y pixel in target image to embed lower left pixel.", &status);

			grid_x++;

		fits_close_file(efptr,&status);
		fits_report_error(stderr, status);
	}

	return 0;
}
