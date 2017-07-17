#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include <stdbool.h>
/*to compile, run
 *gcc -O3 -g jedigrid_a.c -o jedigrid_a -lm -lcfitsio
 *to run,
 *./jedigrid_a 8192 8192 gallist.txt gallist_grid.txt
 */


char *help[] = {
	"Usage jedigrid_a gallist gallist_out",
	"Input galaxy parameter file: x y nx ny zs file",
	"           nx - width of the galaxy in pixels",
	"           ny - height of the galaxy in pixels",
	"			radius - radius of the circle on which to place the galaxies",
	"			angle - angle around the circle to place the galaxies at in radians",
	"           infile - filepath to the FITS file for this galaxy, 1024 chars max",
	"           outfile - filepath for the output FITS file for this galaxy, 1024 chars max",
	"gallist_out: output file for new galaxy list",
	0};


typedef struct {
	long int    x;         //x - x coord. of lower left pixel where galaxy should be embedded
	long int    y;         //y - y coord. of lower left pixel where galaxy should be embedded
	long int    nx;        //nx - width of the galaxy in pixels
	long int    ny;        //ny - height of the galaxy in pixels
	float       zs;        //zs - redshift of this galaxy
	char        *file;     //file - filepath to the FITS file for this galaxy
	char        *outfile;  //outfile - filepath for the output FITS file for this galaxy

} galaxy;


int main(int argc, char *argv[]){


	/*declare variables*/
	FILE        *galaxyfile;    //galaxy list file
	FILE        *galaxygridfile;    //galaxy grid list file
	FILE        *lensfile;      //lens list file
	long int	nx, ny;
	galaxy      *galaxies;      //array of galaxy structs
	char        buffer[1024], buffer2[1024];   //buffers for reading in strings
	galaxy      tempgal;        //temporary galaxy for reading in galaxies
	long int 	radius;
   	float 		angle;	

	/* Print help */
	if ((argc != 7) || (argv[1][0] == '^')) {
		int i;
		for (i = 0; help[i] != 0; i++)
			fprintf (stderr, "%s\n", help[i]);
		exit (1);
	}


	/* parse command line input*/
	sscanf(argv[1], "%li", &nx);
	sscanf(argv[2], "%li", &ny);
	sscanf(argv[3], "%li", &radius);
	sscanf(argv[4], "%f", &angle);

	/* parse galaxy list file*/
	if ((galaxyfile = fopen(argv[5], "r")) == NULL) {
		fprintf(stderr, "Error opening galaxy list file \"%s\".\n", argv[5]);
		exit(1);
	}
	fprintf(stdout, "Opened galaxy list file \"%s\".\n", argv[5]);

	/* parse galaxy list file*/
	if ((galaxygridfile = fopen(argv[6], "w")) == NULL) {
		fprintf(stderr, "Error opening galaxy grid list file \"%s\".\n", argv[6]);
		exit(1);
	}
	fprintf(stdout, "Opened galaxy grid list file \"%s\".\n", argv[6]);


	//count the number of galaxies
	long int ngalaxies = 0;
	while(fscanf(galaxyfile, "%li %li %li %li %f %s %s", &tempgal.x, &tempgal.y, &tempgal.nx, &tempgal.ny, &tempgal.zs, buffer, buffer2) > 0){
		ngalaxies++;
	}
	fprintf(stdout, "%li galaxies in \"%s\".\n", ngalaxies, argv[5]);

	//initialize memory for all the galaxies
	galaxies = (galaxy *) calloc(ngalaxies, sizeof(galaxy));
	if(galaxies == NULL){
		fprintf(stderr,"Error allocating memory for galaxy list.");
		exit(1);
	}

	//read in the galaxies
	rewind(galaxyfile);
	int ngal;
	for(ngal = 0; ngal < ngalaxies; ngal++){
		//fscanf(galaxyfile, "%li %li %li %li %f %s %s", &galaxies[ngal].x, &galaxies[ngal].y, &galaxies[ngal].nx, &galaxies[ngal].ny, &galaxies[ngal].zs, &galaxies[ngal].file, &galaxies[ngal].outfile);

		fscanf(galaxyfile, "%li %li %li %li %f %s %s", &galaxies[ngal].x, &galaxies[ngal].y, &galaxies[ngal].nx, &galaxies[ngal].ny, &galaxies[ngal].zs, &buffer, &buffer2);
		galaxies[ngal].file = calloc(strlen(buffer)+1, sizeof(char));
		if(galaxies[ngal].file == NULL){
			fprintf(stderr, "Error allocating memory for galaxy list.");
			exit(1);
		}
		strcpy(galaxies[ngal].file, buffer);
		galaxies[ngal].outfile = calloc(strlen(buffer2)+1, sizeof(char));
		if(galaxies[ngal].outfile == NULL){
			fprintf(stderr, "Error allocating memory for galaxy list.");
			exit(1);
		}
		strcpy(galaxies[ngal].outfile, buffer2);
		//fprintf(stdout,"Galaxy %li:\nx: %li\ny: %li\nnx: %li\nny: %li\nzs: %f\nfile: \"%s\"\noutfile: \"%s\"\n\n", ngal, galaxies[ngal].x,galaxies[ngal].y,galaxies[ngal].nx,galaxies[ngal].ny,galaxies[ngal].zs,galaxies[ngal].file,galaxies[ngal].outfile);
	
		
		galaxies[ngal].x = (nx/2) + radius*cos(angle) - (galaxies[ngal].nx/2);
		galaxies[ngal].y = (ny/2) + radius*sin(angle) - (galaxies[ngal].ny/2);

		//fprintf(stdout, "%li %li %li %li %f %s %s\n", galaxies[ngal].x, galaxies[ngal].y, galaxies[ngal].nx, galaxies[ngal].ny, galaxies[ngal].zs, galaxies[ngal].file, galaxies[ngal].outfile);
		fprintf(galaxygridfile, "%li %li %li %li %f %s %s\n", galaxies[ngal].x, galaxies[ngal].y, galaxies[ngal].nx, galaxies[ngal].ny, galaxies[ngal].zs, galaxies[ngal].file, galaxies[ngal].outfile);


	}
	fclose(galaxyfile);

    return 0;
}


