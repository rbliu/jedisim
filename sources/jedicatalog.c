#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fitsio.h"

/*jedicatalog
 *This program creates a realistic catalog of galaxies.
 */

#define BANDHEIGHT 2048
#define BINS_PER_MAG 100
#define DELTA_MAG 0.01
#define DATABASE_PIXSCALE 0.03

char *help[] = {
    "Creates a realistic catalog of galaxies.",
    "Usage: jedicatalog config_file",
    "Arguments: config_file - a text file with the configuration settings for this program.",
    0};

typedef struct {
    float   radius; //r50 radius of the stamp in arcseconds
    float   pixscale;   //pixel scale of the stamp in arcseconds per pixel
    float   magnitude;  //magnitude of the stamp
    char    *name;   //the file name for the postage stamp
} Stamp;


typedef struct {
    char    name[1024]; //the name of the galaxy
    float   x;  //the x position of the galaxy
    float   y;  //the y position of the galaxy
    float   angle;  //the angle to rotate the galaxy (in degrees)
    float   redshift;   //the redshift of the galaxy
    float   pixscale;   //the pixelscale of the galaxy in arcseconds per pixel
    float   old_mag;    //the original magnitude of the galaxy
    float   old_rad;    //the original radius (pixels) of the galaxy
    float   new_mag;    //the final magnitude of the galaxy
    float   new_rad;    //the final radius of the galaxy
    char    stamp_name[1024]; //the filepath for the postage stamp of this galaxy after its parameters are set
    char    dis_name[1024];   //the filepath for the postage stamp of this galaxy after it has been distorted
} Galaxy;


int compareInt(const void* p1, const void* p2);
int compareStamp(const void* p1, const void* p2);
void print_galaxy(Galaxy *gal, FILE *fptr);
float rand_float();
float get_mag(float *CDF, int min_mag, int num_bins);

int main(int argc, char *argv[]){
    FILE        *config_fp;      //filepointer for the configuration file
    char        buffer[1024], buffer3[1024];   //string buffers to read in the config file
    char        *buffer2;        //another string buffer

    //physics settings
    int         single_redshift = 0;    //0 -> get redshift from db; 1-> fixed redshift for all galaxies
    float       fixed_redshift;     //the fixed redshift to use if single_redshift=1
    float       pix_scale;      //pixel scale to use for the simulation in arcseconds per pixel
    long int    nx, ny;         //image size to simulate
    long int    xborder, yborder;   //borders on the image so we don't have galaxies going off the edges
    long int    ngalaxies;      //the total number of galaxies to simulate
    int         min_mag, max_mag;   //the minimum and maximum magnitude
    char        radius_db_folder[1024], red_db_folder[1024];    //the folders for the redshift and radius databases
    float       mag_power;      //the power for the power law distribution of galactic magnitudes
    
    //output settings
    char        output_folder[1024];    //the output folder for postage stamps
    char        prefix[1024];           //the prefix for all filenames associated with this trial

    //catalog file settings
    char        temp_catalog_file[1024], catalog_file[1024];//catalog of galaxy parameters
    FILE        *catalog_fp;
    char        temp_convlist_file[1024], convlist_file[1024];//list of postage stamps to be convolved
    FILE        *convlist_fp;
    char        temp_distortedlist_file[1024], distortedlist_file[1024];//list of distorted postage stamps
    FILE        *distortedlist_fp;
    char        temp_convolvedlist_file[1024], convolvedlist_file[1024];//list of convolved postage stamps
    FILE        *convolvedlist_fp;

    //source image settings
    long int    num_source_images;  //the number of source images
    Stamp       *source_images = NULL;  //array of source images
    long int    nimage=0;       //counts how many image filenames have been parsed so far

    //cosmic reality databases
    float       **radius_db;      //radius database: an array of radii for each integer bin of magnitudes, from min_mag to max_mag
    float       **red_db;         //redshift database: same format as radius_db
    long int    *radius_db_bin_sizes;   //array that lists the size of each magnitude bin for the radius database
    long int    *red_db_bin_sizes;   //array that lists the size of each magnitude bin for the redshift database
    float       *CDF;   //cumulative distribution function for the power law dist. of magnitudes  

    //parse command line input
    if(argc != 2){
        int line;
        for(line = 0; help[line] != 0; line++){
            fprintf(stderr, "%s\n", help[line]);
        }
        exit(1);
    }
    config_fp = fopen(argv[1],"r");
    if(config_fp == NULL){
        fprintf(stderr,"Error: cannot open config file.");
        exit(1);
    }

    //parse config file
    while(fgets(buffer, 1024, config_fp) != NULL){
        //read in whatever isn't a comment
        if(buffer[0] != '#'){
            buffer2 = strtok(buffer,"#\n");
            sscanf(buffer2,"%[a-z_]=",buffer3);
            //physics settings
            if(strcmp(buffer3,"pix_scale")==0) sscanf(buffer2,"pix_scale=%f",&pix_scale);
            else if(strcmp(buffer3,"nx")==0) sscanf(buffer2,"nx=%li", &nx);
            else if(strcmp(buffer3,"ny")==0) sscanf(buffer2,"ny=%li", &ny);
            else if(strcmp(buffer3,"x_border")==0) sscanf(buffer2,"x_border=%li", &xborder);
            else if(strcmp(buffer3,"y_border")==0) sscanf(buffer2,"y_border=%li", &yborder);
            else if(strcmp(buffer3,"num_galaxies")==0) sscanf(buffer2,"num_galaxies=%li", &ngalaxies);
            else if(strcmp(buffer3,"min_mag")==0) sscanf(buffer2,"min_mag=%i",&min_mag);
            else if(strcmp(buffer3,"max_mag")==0) sscanf(buffer2,"max_mag=%i",&max_mag);
            else if(strcmp(buffer3,"power")==0) sscanf(buffer2,"power=%f", &mag_power);
            else if(strcmp(buffer3,"single_redshift")==0) sscanf(buffer2,"single_redshift=%i", &single_redshift);
            else if(strcmp(buffer3,"fixed_redshift")==0) sscanf(buffer2,"fixed_redshift=%f", &fixed_redshift);

            //output settings
            else if(strcmp(buffer3,"output_folder")==0) sscanf(buffer2, "output_folder=\"%[^\"]", output_folder);
            else if(strcmp(buffer3,"prefix")==0) sscanf(buffer2, "prefix=\"%[^\"]", prefix);

            //catalog file settings
            else if(strcmp(buffer3,"catalog_file")==0) sscanf(buffer2, "catalog_file=\"%[^\"]", temp_catalog_file);

            else if(strcmp(buffer3,"convlist_file")==0) sscanf(buffer2, "convlist_file=\"%[^\"]", temp_convlist_file);
            else if(strcmp(buffer3,"distortedlist_file")==0) sscanf(buffer2, "distortedlist_file=\"%[^\"]", temp_distortedlist_file);
            else if(strcmp(buffer3,"convolvedlist_file")==0) sscanf(buffer2, "convolvedlist_file=\"%[^\"]", temp_convolvedlist_file);

            //cosmic reality databases
            else if(strcmp(buffer3,"radius_db_folder")==0) sscanf(buffer2, "radius_db_folder=\"%[^\"]", radius_db_folder);
            else if(strcmp(buffer3,"red_db_folder")==0) sscanf(buffer2, "red_db_folder=\"%[^\"]", red_db_folder);
            
            //image settings
            else if(strcmp(buffer3,"num_source_images")==0){
                sscanf(buffer2,"num_source_images=%li",&num_source_images);
                source_images = (Stamp*) calloc(num_source_images, sizeof(Stamp));
                if(source_images == NULL){
                    fprintf(stderr,"Error: could not allocate list of source images.\n");
                    exit(1);
                }
            } else if(strcmp(buffer3,"image")==0 && nimage < num_source_images){
                if(source_images == NULL){
                    fprintf(stderr, "Error: the parameter 'num_galaxies' must come before any 'image' parametersn");
                    exit(1);
                }
                char temp[1024];
                sscanf(buffer2, "image=\"%[^\"]", temp);
                source_images[nimage].name = (char*) calloc(strlen(temp)+1, sizeof(char));
                if(source_images[nimage].name == NULL){
                    fprintf(stderr,"Error: could not allocate source image struct %i.\n", nimage);
                    exit(1);
                }
                strcpy(source_images[nimage].name,temp);
                nimage++;
            }

        }
    }
    
    sprintf(catalog_file, "%s%s%s", output_folder, prefix, temp_catalog_file);
    sprintf(convlist_file, "%s%s%s", output_folder, prefix, temp_convlist_file);
    sprintf(distortedlist_file, "%s%s%s", output_folder, prefix, temp_distortedlist_file);
    sprintf(convolvedlist_file, "%s%s%s", output_folder, prefix, temp_convolvedlist_file);


    //print out what was just read in
    fprintf(stdout,"pixscale: %f\n", pix_scale);
    fprintf(stdout,"nx: %li\n", nx);
    fprintf(stdout,"ny: %li\n", ny);
    fprintf(stdout,"x_border: %li\n", xborder);
    fprintf(stdout,"y_border: %li\n", yborder);
    fprintf(stdout,"num_galaxies: %li\n",ngalaxies);
    fprintf(stdout,"min_mag: %li\n", min_mag);
    fprintf(stdout,"max_mag: %li\n", max_mag);
    fprintf(stdout,"radius_db_folder: %s\n", radius_db_folder);
    fprintf(stdout,"red_db_folder: %s\n", red_db_folder);
    fprintf(stdout,"output_folder: %s\n", output_folder);
    fprintf(stdout,"power: %f\n", mag_power);
    fprintf(stdout,"catalog_file: %s\n", catalog_file);
    fprintf(stdout,"convlist_file: %s\n", convlist_file);
    fprintf(stdout,"distortedlist_file: %s\n", distortedlist_file);
    fprintf(stdout,"convolvedlist_file: %s\n", convolvedlist_file);
    fprintf(stdout,"num_source_images: %li\n", num_source_images);

    //get magnitude and radius of each of the source galaxies with cfitsio
    {
        fitsfile    *sfptr;     //source image fits file pointer
        int         status;     //fits statusAN
        for(nimage = 0; nimage < num_source_images; nimage++){
            if(source_images[nimage].name == NULL){
                fprintf(stderr,"Error: image file path %i is null.", nimage);
                exit(1);
            }
            char temp[1024];
            fits_open_file(&sfptr, source_images[nimage].name, READONLY, &status);
            
            fits_read_key_str(sfptr, "MAG", temp, NULL, &status);
            sscanf(temp,"%f", &source_images[nimage].magnitude);

            fits_read_key_str(sfptr,"RADIUS", temp, NULL, &status);
            sscanf(temp,"%f", &source_images[nimage].radius);

            fits_read_key_str(sfptr,"PIXSCALE", temp, NULL, &status);
            sscanf(temp, "%f", &source_images[nimage].pixscale);
            source_images[nimage].radius *= source_images[nimage].pixscale;
            
            fits_close_file(sfptr,&status);
            if(status){
                fits_report_error(stderr,status);
                exit(1);
            }
        }
        qsort(source_images, num_source_images, sizeof(Stamp), compareStamp);
    }

    //print out the source images
    //for(nimage = 0; nimage < num_source_images; nimage++)
      //  fprintf(stdout,"image %i:\n\tName:%s\n\tMagnitude: %f\n\tRadius: %f (arcseconds)\n", nimage, source_images[nimage].name, source_images[nimage].magnitude, source_images[nimage].radius);

    //read in the radius and redshift databases
    {
        radius_db = (float **) calloc((max_mag-min_mag+1), sizeof(float*));
        radius_db_bin_sizes = (long int *) calloc((max_mag-min_mag+1), sizeof(long int));
        red_db_bin_sizes = (long int *) calloc((max_mag-min_mag+1), sizeof(long int));
        if(radius_db == NULL){
            fprintf(stderr,"Error: could not allocate radius database.\n");
            exit(1);
        } 
        fprintf(stdout,"Allocated radius db with %i bins.\n", (max_mag-min_mag+1));

        red_db = (float **) calloc((max_mag-min_mag+1), sizeof(float*));
        if(red_db == NULL){
            fprintf(stderr,"Error: could not allocate radius database.\n");
            exit(1);
        } 
        fprintf(stdout,"Allocated radius db with %i bins.\n", (max_mag-min_mag+1));

        int mag_bin;
        FILE    *radiusfp, *redfp;  //radius and redshift filepointers
        for(mag_bin=0; mag_bin <= (max_mag-min_mag); mag_bin++){
            char radius_filename[1024], redshift_filename[1024];    //names of the redshift and radius data files
            //using convention that the file names are radius_db_folder/n.dat where n from min_mag to max_mag
            sprintf(radius_filename,"%s%i.dat",radius_db_folder,mag_bin+min_mag);
            sprintf(redshift_filename,"%s%i.dat",red_db_folder,mag_bin+min_mag);

            //import radius bin
            radiusfp = fopen(radius_filename,"r");
            long int     nradius = 0;    //radius counter
            while(fgets(buffer, 1024, radiusfp) != NULL){
                nradius++;                
            }
            radius_db[mag_bin] = (float *) calloc(nradius, sizeof(float));
            radius_db_bin_sizes[mag_bin] = nradius;
            if(radius_db[mag_bin]==NULL){
                fprintf(stderr,"Error: could not allocate radius_db bin %i.\n", mag_bin);
                exit(1);
            }
            fprintf(stdout,"Allocated radius db bin %i of size %i.\n", mag_bin, nradius);
            rewind(radiusfp);
            nradius = 0;
            while(fgets(buffer, 1024, radiusfp) != NULL){
                if(sscanf(buffer,"%f",&radius_db[mag_bin][nradius]) != 1){
                    fprintf(stderr,"Error");
                    exit(1);
                }
                radius_db[mag_bin][nradius] *= DATABASE_PIXSCALE;
                //fprintf(stdout,"bin: %i\tnradius: %li\tsize: %li\tvalue: %.3f\tbuffer: %s\n", mag_bin+min_mag, nradius, size, radius_db[mag_bin][nradius], buffer);
                nradius++;
            }
            fclose(radiusfp);
            //sort the radii
            qsort(radius_db[mag_bin], nradius, sizeof(float), compareInt);

            /*int i;
              for(i = 0; i < nradius; i++)
              fprintf(stdout,"bin: %i\tnradius: %li\tsize: %li\tvalue: %.3f\n", mag_bin+min_mag, i, nradius, radius_db[mag_bin][i]);*/


            //import redshift bin

            redfp = fopen(redshift_filename,"r");
            long int     nred = 0;    //radius counter
            while(fgets(buffer, 1024, redfp) != NULL){
                nred++;                
            }
            red_db[mag_bin] = (float *) calloc(nred, sizeof(float));
            red_db_bin_sizes[mag_bin] = nred;
            if(radius_db[mag_bin]==NULL){
                fprintf(stderr,"Error: could not allocate radius_db bin %i.\n", mag_bin);
                exit(1);
            }
            fprintf(stdout,"Allocated red db bin %i of size %i.\n", mag_bin, nred);
            rewind(redfp);
            nred = 0;
            while(fgets(buffer, 1024, redfp) != NULL){
                if(sscanf(buffer,"%f",&red_db[mag_bin][nred]) != 1){
                    fprintf(stderr,"Error");
                    exit(1);
                }
                //fprintf(stdout,"bin: %i\tnradius: %li\tvalue: %.3f\tbuffer: %s\n", mag_bin+min_mag, nred, red_db[mag_bin][nred], buffer);
                nred++;
            }
            fclose(redfp);
        }
        for(mag_bin = 0; mag_bin <=(max_mag-min_mag); mag_bin++)
           fprintf(stdout,"mag bin: %i\t redshift size: %li\n", mag_bin, red_db_bin_sizes[mag_bin]); 
    }

    //set up the Cumulative distribution function for the magnitude distribution
    long int    nmag_dist_bins = BINS_PER_MAG*(max_mag-min_mag);    //number of bins in the CDF
    CDF = (float*) calloc(nmag_dist_bins, sizeof(float));
    if(CDF == NULL){
        fprintf(stdout,"Error: could not allocate CDF array.\n");
        exit(1);
    }
    {
        float total = 0;
        long int    i;  //counter
        for(i = 0; i < nmag_dist_bins; i++){
            total += expf(mag_power*logf(10)*((float) min_mag + (float) i*DELTA_MAG));
            CDF[i] = total;
        }
        for(i=0; i < nmag_dist_bins; i++)
            CDF[i]/=total;
        /*for(i =0; i < nmag_dist_bins; i++)
            fprintf(stdout, "%i: %f\n", i, CDF[i]);*/
    }

    //make and print the catalog
    long int    g;  //galaxy counter
    
    catalog_fp = fopen(catalog_file, "w");
    if(catalog_fp == NULL){
        fprintf(stderr, "Error: could not open catalog file.\n");
        exit(1);
    }
    distortedlist_fp = fopen(distortedlist_file, "w");

    for(g = 0; g < ngalaxies; g++){
        srand(time(NULL)+rand());   //the code runs too fast to use the time in miliseconds as the seed
        Galaxy gal;
        
        //set galaxy magnitude
        gal.new_mag = get_mag(CDF, min_mag, nmag_dist_bins);
        int     mag = (int) floor(gal.new_mag); //the magnitude bin we're in
        
        //fprintf(stdout,"mag bin: %i\tnew radius length: %li\told radius length: %i\n", mag, radius_db_bin_sizes[mag-min_mag], num_source_images);

        //make sure there is some pair (new_radius, old_radius) for this magnitude bin such that new_radius < old_radius
        //find the index of smallest old radius (i.e. from a source image) that is larger than the smallest new_radius and therefore doesn't cause problems
        int     i = 0;  //counter
        float   temp = radius_db[mag-min_mag][0];//smallest new radius for this magnitude
        while(temp > source_images[i].radius && i < num_source_images)
            i++;
        int     min_old_radius_index = i;   //the index of the smallest usable old radius

        //find the index of the largest new radius (i.e. from the database) that is smaller than the largest old radius and therefore doesn't cause problems
        i = radius_db_bin_sizes[mag-min_mag]-1;//last index of the new radius list for this mag bin
        temp = source_images[num_source_images-1].radius;//the largest old radius
        while(radius_db[mag-min_mag][i] > temp && i >= 0)
            i--;
        int max_new_radius_index = i;

        //fprintf(stdout,"min_old_radius_index: %i\t max_new_radius_index: %i\n", min_old_radius_index,  max_new_radius_index);

        //if there's no valid conbination, throw an error
        if(min_old_radius_index == num_source_images || max_new_radius_index == -1){
            fprintf(stderr,"Error: no possible combination of new and old radius can work. Magnitude bin: %i.\n", mag);
            exit(1);
        }
        //otherwise, choose some image
        int     image_num = min_old_radius_index + (rand() % (num_source_images - min_old_radius_index));
        if(image_num < min_old_radius_index || image_num >= num_source_images){
            fprintf(stderr,"Error: invalid old radius index.");
            exit(1);
        }
        
        //set all the parameters that depend on the old galaxy
        sprintf(gal.name,"%s", source_images[image_num].name);
        gal.pixscale = source_images[image_num].pixscale;
        gal.old_mag = source_images[image_num].magnitude;
        gal.old_rad = source_images[image_num].radius;

        //restrict the database radii to the compatible ones
        i = max_new_radius_index; 
        while(radius_db[mag-min_mag][i] > source_images[image_num].radius && i >= 0)
            i--;
        max_new_radius_index = i;

        //make sure that there are some compatible ones
        if(max_new_radius_index == -1){
            fprintf(stderr, "Error: no compatible new radii.");
            exit(1);
        }
        
        //choose a new radius
        int new_radius_num; 
        if(max_new_radius_index ==0)
            new_radius_num = 0;
        else
            new_radius_num = rand() % max_new_radius_index;
        
        gal.new_rad = radius_db[mag-min_mag][new_radius_num];

        //double-check, in case there's a bug
        if(gal.new_rad > gal.old_rad){
            fprintf(stderr,"Error: selected illegal new radius.");
            exit(1);
        }
        if(single_redshift == 0){
            //get a redshift from the database
            int redshift_num = rand() % red_db_bin_sizes[mag-min_mag];
            gal.redshift = red_db[mag-min_mag][redshift_num];
        } else if(single_redshift == 1){
            gal.redshift = fixed_redshift;
        } else {
            fprintf(stderr, "Error: single_redshift must be either 0 or 1.");
        }


        //set galaxy position and angle 
        gal.x = xborder + rand_float()*(nx-2*xborder);
        gal.y = yborder + rand_float()*(nx-2*yborder);
        gal.angle = rand_float()*360;

        //set galaxy filepaths
        sprintf(gal.stamp_name, "%sstamp_%i/stamp_%i.fits.gz", output_folder, g/1000, g);

		//can't write to .gz files with FITSIO, so this is disabled so jedigrid can work
        //sprintf(gal.dis_name, "%sdistorted_%i/distorted_%i.fits.gz", output_folder, g/1000, g);
        sprintf(gal.dis_name, "%sdistorted_%i/distorted_%i.fits", output_folder, g/1000, g);
        

        
        //print_galaxy(&gal, stdout);
        print_galaxy(&gal, catalog_fp);
        fprintf(distortedlist_fp, "%s\n", gal.dis_name);
    }


    convolvedlist_fp = fopen(convolvedlist_file, "w");
    int band, num_bands = ny/BANDHEIGHT;
    char conv_name[1024];
    for(band = 0; band < num_bands; band++){   
        sprintf(conv_name, "%sconvolved/convolved_band_%i.fits", output_folder, band);
        fprintf(convolvedlist_fp, "%s\n", conv_name);
    }

    fclose(catalog_fp);
    fclose(distortedlist_fp);
    fclose(convolvedlist_fp);
    fprintf(stdout,"Success! Catalog made with %i entries.\n", ngalaxies);
    return 0;
}

//comparator for sorting a list of integers
int compareInt(const void* p1, const void* p2){
    return (*(int*)p1 - *(int*)p2); 
}

//comparator for sorting a list of stamps
int compareStamp(const void* p1, const void* p2){
    float diff = ((Stamp*)p1)->radius - ((Stamp*)p2)->radius;
    if(diff > 0) return 1;
    else if(diff <0) return -1;
    else return 0; 
}

//prints out the given galaxy to the given file stream
void print_galaxy(Galaxy *gal, FILE *fptr){
    fprintf(fptr, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n",gal->name, gal->x, gal->y, gal->angle, gal->redshift, gal->pixscale, gal->old_mag, gal->old_rad, gal->new_mag, gal->new_rad, gal->stamp_name, gal->dis_name);
    //fprintf(fptr,"%f\t%f\n",gal->x, gal->y);
}

//returns a random float in [0,1)
float rand_float(){
    return ((float) rand())/RAND_MAX;
}

float get_mag(float *CDF, int min_mag, int num_bins){
    float r = rand_float();    
    int j=0;
    while(CDF[j]<= r && r < num_bins) 
        j++;
    //fprintf(stdout,"min_mag: %i\tnum_bins: %i\trand: %f\tj: %li\tCDF[j]: %f\n",min_mag, num_bins, r, j, CDF[j]); 
    return (float) min_mag + j*DELTA_MAG;
}
