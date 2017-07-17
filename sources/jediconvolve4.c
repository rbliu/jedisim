#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include <fftw3.h>

#define BANDHEIGHT 2048

char *help[] = {
    "Convolves an image with a specified PSF with a fourier-space convolution. The convolution is carried out in bands so that it can be performed on very large images.",
    "usage: jedirescale input_file PSF_file output_file",
    "Arguments: input_file - FITS file to rescale",
    "           PSF_file - the PSF with which to convolve",
    "           output_folder - output folder for images",
    0};

int main(int argc, char *argv[]){
    //declare variables
    char infile[1024], psffile[1024], outfolder[1024];

    //fits variables
    fitsfile    *infptr, *psffptr, *outfptr;    //fits pointers for the input, psf, and output files
    int         status = 0, naxis = 0;          //cfitsio status counter & image dimensionality
    long int    inaxes[2], bnaxes[2], pnaxes[2], onaxes[2];    //dimensions of input image, band, PSF, and output image
    long int    fpixel[2] = {1,1}, bpixel[2] = {1,1};//cfitsio first/last pixel to read in/out
    long int    lpixel[2];
    long int    inc[2] = {1,1};
    float       *bimage, *pimage, *pimage2, *oimage; //arrays for the input, psf, resized psf, and  output images
    int         num_bands, band;    //number of bands, and counter for the current band    
    long int    row, col;           //counters

    //print help
    if(argc != 4){
        int line;
        for(line = 0; help[line] !=0; line++)
            fprintf(stderr, "%s\n", help[line]);
        exit(1);
    }

    //parse command line input
    sscanf(argv[1], "%s", infile);
    sscanf(argv[2], "%s", psffile);
    sscanf(argv[3], "%s", outfolder);

    //open input image
    fprintf(stdout,"Getting size of the input image.\n");
    fits_open_file(&infptr, infile, READONLY, &status);
    fits_get_img_dim(infptr, &naxis, &status);
    fits_get_img_size(infptr, 2, inaxes, &status);
    if(status){
        fits_report_error(stderr,status);
        exit(1);
    }
    fprintf(stdout,"Input image has size: (%li, %li).\n", inaxes[0], inaxes[1]);

    //read in the PSF image
    fprintf(stdout,"Reading in PSF.\n");
    fits_open_file(&psffptr, psffile, READONLY, &status);
    fits_get_img_dim(psffptr, &naxis, &status);
    fits_get_img_size(psffptr, 2, pnaxes, &status);
    if(status){
        fits_report_error(stderr, status);
        exit(1);
    }

    //allocate enough memory for the PSF
    pimage = (float *) calloc (pnaxes[0]*pnaxes[1], sizeof(float));
    if (pimage == NULL){
        fprintf(stderr, "Error allocating memory for psf.\n"); 
        exit(1);
    }

    //read the PSF pixels 
    if(fits_read_pix(psffptr, TFLOAT, fpixel, pnaxes[0]*pnaxes[1], NULL, pimage, NULL, &status)){
        fprintf(stderr, "Can't read in image.\n");
        fits_report_error(stderr, status);
        exit(1);
    }
    fits_close_file(psffptr, &status);
    if(status){
        fprintf(stderr, "Error reading PSF file.\n");
        fits_report_error(stderr, status);
        exit(1);
    }


    //the bnaxes array specifies the dimension of one band of the input image
    bnaxes[0] = inaxes[0];
    bnaxes[1] = BANDHEIGHT;

    //the onaxes array specifies the dimension of one band of the output image
    //we need to put a border the size of the PSF around the band
    onaxes[0] = inaxes[0]+2*pnaxes[0];
    onaxes[1] = BANDHEIGHT+2*pnaxes[1];
    fprintf(stdout,"Band shape: (%li, %li).\n",bnaxes[0],bnaxes[1]);
    fprintf(stdout,"Output image shape: (%li, %li).\n", onaxes[0],onaxes[1]);
    num_bands = inaxes[1]/bnaxes[1];

    //the last pixel we want to read in for the first band
    lpixel[0] = inaxes[0];
    lpixel[1] = BANDHEIGHT;

    fprintf(stdout,"Resizing and wrapping PSF.\n");
    pimage2 = (float *) calloc(onaxes[0]*onaxes[1], sizeof(float));
    if(pimage2 == NULL){
        fprintf(stderr, "Error: could not allocate memory for PSF image 2.\n");
        exit(1);
    }

    //make PSF with the center pixel of the PSF at (0,0) and wrap around to the other sides
    //also change to row-major order
    for(row = 0; row < pnaxes[0]; row++){
        for(col = 0; col < pnaxes[1]; col++){
            //place in the new PSF array
            int new_row = (row - pnaxes[0]/2 + onaxes[0]) % onaxes[0]; 
            int new_col = (col - pnaxes[1]/2 + onaxes[1]) % onaxes[1];
            pimage2[new_row*onaxes[1]+new_col] = pimage[col*pnaxes[0]+row];
        }
    }
    //we don't need the old PSF array, so throw it away to save some memory
    free(pimage);


    //allocate memory to read in the one band
    bimage = (float *) calloc(bnaxes[0]*bnaxes[1], sizeof(float));
    if(bimage == NULL){
        fprintf(stderr, "Error allocating memory for band.\n");
        exit(1);
    }

    //allocate memory so there's a place to put one band
    oimage = (float *) calloc(onaxes[0]*onaxes[1], sizeof(float));
    if(oimage == NULL){
        fprintf(stderr, "Error allocating memory for output image.\n");
        exit(1);
    }
            

    //arrays for the input and output images in row-major order for FFTW3
    float *IMG, *OUT;
    IMG = (float *) calloc(onaxes[0]*onaxes[1], sizeof(float));
    if(IMG == NULL){
        fprintf(stderr,"Error allocating memory for FFTW3 input image.\n");
        exit(1);
    }
    OUT = (float *) calloc(onaxes[0]*onaxes[1], sizeof(float));
    if(OUT == NULL){
        fprintf(stderr,"Error allocating memory for FFTW3 input image.\n");
        exit(1);
    }




    fprintf(stdout,"Setting up FFT.\n");
    float   scale = 1.0/(onaxes[0]*onaxes[1]); //scaling factor for convolution

    //make arrays for the transformed input, PSF, and output images
    fftwf_complex *FIMG, *FPSF, *FOUT;
    
    //make FFTW3 plans for transforming them
    fftwf_plan pIMG, pPSF, pinv;

    //since this is a real-to-complex fourier transform, we don't need one complex number
    //for every pixel in the input array, but only about half that many
    long int    fcol = (onaxes[1]/2)+1;

    //allocate memory for the transformed images.
    //I think FFTW3 will give an error message if they can't be allocated
    //I'm also pretty sure that they shouldn't be de-allocated, because the FFTW3
    //plan depends on where exactly in memory they're stored
    FIMG = (fftwf_complex *) fftwf_malloc(onaxes[0]*fcol*sizeof(fftwf_complex));
    FPSF = (fftwf_complex *) fftwf_malloc(onaxes[0]*fcol*sizeof(fftwf_complex));
    FOUT = (fftwf_complex *) fftwf_malloc(onaxes[0]*fcol*sizeof(fftwf_complex));

    //make the FFTW3 plans
    pIMG = fftwf_plan_dft_r2c_2d(onaxes[0], onaxes[1], IMG, FIMG, FFTW_ESTIMATE);
    pPSF = fftwf_plan_dft_r2c_2d(onaxes[0], onaxes[1], pimage2, FPSF, FFTW_ESTIMATE);
    pinv = fftwf_plan_dft_c2r_2d(onaxes[0], onaxes[1], FOUT, OUT, FFTW_ESTIMATE);

    //transform the PSF because it's transform is the same for all bands
    fprintf(stdout,"Taking the FFT of the PSF.\n");
    fftwf_execute(pPSF);

    //we don't need the un-transformed PSF array anymore, so save some memory and get rid of it
    free(pimage2);

    //loop over all the bands
    for(band = 0; band < num_bands; band++){

        fprintf(stdout,"Band %i/%i: reading in image. Size: %.1fMB.\n", band+1, num_bands, (float) bnaxes[0]*bnaxes[1]*sizeof(float)/(1000*1000));

        //read in the band
        fprintf(stdout,"first: (%li, %li)\tlast:(%li, %li).\n", bpixel[0], bpixel[1], lpixel[0], lpixel[1]);
        fits_read_subset(infptr, TFLOAT, bpixel, lpixel, inc, NULL, bimage, NULL, &status);
        if(status){
            fprintf(stderr, "Error: can't read in band %i/%i of the input image.\n",band, num_bands);
            fits_report_error(stderr, status);
            exit(1);
        }
            
        //put image into row-major order
        fprintf(stdout,"Band %i/%i: transposing input to row-major order.\n",band+1,num_bands);
        for(row = 0; row < bnaxes[0]; row++){
            for(col = 0; col < bnaxes[1]; col++){
                long int col_m_index = row+col*bnaxes[0];
                long int row_m_index = (pnaxes[1] + col)+(row+pnaxes[0])*onaxes[1];
                IMG[row_m_index] = bimage[col_m_index];
            }
        }

        //transform the band
        fprintf(stdout,"Band %i/%i: taking FFT. \n", band+1, num_bands);
        fftwf_execute(pIMG);

        //perform convolution by multiplying in fourier space
        fprintf(stdout,"Band %i/%i: convolving.\n",band+1,num_bands);
        for(row = 0; row < onaxes[0]; row++){
            for(col = 0; col < fcol; col++){
                long int index = row*fcol + col; //the row-major index
                //do complex multiplication and get the scale right
                FOUT[index][0] = scale*(FIMG[index][0] * FPSF[index][0] - FIMG[index][1] * FPSF[index][1]);
                FOUT[index][1] = scale*(FIMG[index][0] * FPSF[index][1] + FIMG[index][1] * FPSF[index][0]);
            }
        }
    
        //transform back
        fprintf(stdout,"Band %i/%i: taking inverse FFT.\n", band+1,num_bands);
        fftwf_execute(pinv);

        fprintf(stdout,"Band %i/%i: transposing back to column-major order.\n",band+1,num_bands); 

        //put the output into column-major order
        for(row = 0; row < onaxes[0]; row++){
            for(col = 0; col < onaxes[1]; col++){
                long int col_m_index = row+col*onaxes[0];
                long int row_m_index = col+row*onaxes[1];
                oimage[col_m_index] = OUT[row_m_index];
            }
        }

        //write out the band to memory, each band as a separate file
        fprintf(stdout,"Band %i/%i: writing convolved band.\n",band+1,num_bands);
        char outfile[1024];
        sprintf(outfile,"%sconvolved_band_%i.fits",outfolder, band);
        long int xembed = bpixel[0]-pnaxes[0]-1;
        long int yembed = bpixel[1]-pnaxes[1]-1; 
        fits_create_file(&outfptr, outfile, &status);
        fits_create_img(outfptr, FLOAT_IMG, naxis, onaxes, &status);

        //keep track of where this band belongs so we can use jedipaste to recombine them
        fits_update_key(outfptr, TLONG, "XEMBED", &xembed , "x pixel in target image to embed lower left pixel.", &status);
        fits_update_key(outfptr, TLONG, "YEMBED", &yembed , "y pixel in target image to embed lower left pixel.", &status);
        fits_write_pix(outfptr, TFLOAT, fpixel, onaxes[0]*onaxes[1], oimage, &status);
        fits_close_file(outfptr, &status);

        //make sure nothing went wrong
        if(status){
            fprintf(stderr, "Error writing out band image.\n");
            fits_report_error(stderr, status);
            exit(1);
        }

        //step up to the next band for CFITSIO
        bpixel[1] += BANDHEIGHT;
        lpixel[1] += BANDHEIGHT;
    }

    //clean up
    fprintf(stdout,"Cleaning up.\n");

    fftwf_destroy_plan(pIMG);
    fftwf_destroy_plan(pPSF);
    fftwf_destroy_plan(pinv);

    fftwf_free(FIMG);
    fftwf_free(FPSF);
    fftwf_free(FOUT);

    return 0;
}   
