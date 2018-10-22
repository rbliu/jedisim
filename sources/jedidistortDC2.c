#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include <stdbool.h>
/*to compile, run
 *gcc -O3 -g jedidistortDC2.c -o jedidistortDC2 -lm -lcfitsio
 *to run,
 *./jedidistortDC2 8192 8192 gallist.txt lenses.txt 0.05 0.3
 */

#define DR      10          //lens table entries per pixel
#define PI      3.14159265  //surely this doesn't need an explanation
#define C       300000      //speed of light in km/s
#define NSUBPX  4           //number of subpixels in each dimension
#define BNSUBPX 2           //log base two of the number of subpixels for bitshift magic
#define OMEGA_M 0.315       //matter density of the universe today
#define OMEGA_D 0.685       //dark energy density of the universe
#define OPZEQ   3391        //1 + Z_eq
#define C_H_0   4424.778    //c [km/s] / Hubble constant today [km/(s*Mpc)] with H_0 = 67.80
#define G       4.302       //Gravitational constant [10^{-3} (pc/Solar mass) (km/s)^2]
#define H0      67.80       //Planck value for Hubble constant for present day [km/(Mps*s)]
#define EP      0           //epsilon

#define MAP_SIZE 4096				//the size and name of alpha map files are hard-coded!!!
#define AL_FILE1 "alpha1.txt"
#define AL_FILE2 "alpha2.txt"

char *help[] = {
	"Simulates gravitational lensing as realistically as possible. The gravitational lens is specified by a list of lenses and their parameters. The program takes in a list of galaxies and their parameters, to be distorted, and returned distorted versions of those images. The distortion was engineered to be as efficient as possible.",
	"Usage jedidistort x y gallist lenses scale zl",
	"Arguments: x - width of the image MUST BE AN INTEGER MULTIPLE OF 4096",
	"           y - height of the image MUST BE AN INTEGER MULTIPLE OF 4096",
	"           gallist - input galaxy parameter file",
	"           lenses - file containing lens parameters",
	"           scale - pixel scale in arcseconds per pixel",
	"           zl - lens redshift",
	"Input galaxy parameter file: x y nx ny zs file",
	"           x - x coord. of lower left pixel where galaxy should be embedded",
	"           y - y coord. of lower left pixel where galaxy should be embedded",
	"           nx - width of the galaxy in pixels",
	"           ny - height of the galaxy in pixels",
	"           zs - redshift of this galaxy",
	"           infile - filepath to the FITS file for this galaxy, 1024 chars max",
	"           outfile - filepath for the output FITS file for this galaxy, 1024 chars max",
	"Lens parameter file: x y rho0 Rs zl",
	"           x - x center of lens (in pixels)",
	"           y - y center of lens (in pixels)",
	"           type - type of mass profile",
	"                   1. Singular isothermal sphere",
	"                   2. Navarro-Frenk-White profile",
	"					3. NFW constant distortion profile for grid simulations",
	"           p1 - first profile parameter",
	"                   1. sigma_v [km/s]",
	"                   2. M200 parameter [10^14 solar masses]",
	"					3. Distance to center in px. M200 fixed at 10 default, which can be modified in case 3",
	"           p2 - second profile parameter",
	"                   1. not applicable, can take any numerical",
	"                   2. c parameter [unitless]",
	"                   3. c parameter [unitless]",
	0};

typedef struct {
	float       x;          // x coord of lens (pixels)
	float       y;          // y coord of lens (pixels)
	int         type;       // Identifies mass profile
	float       p1;         // mass profile parameter 1
	float       p2;         // mass profile parameter 2
	long int    nr;         // Number of entries in the distortion table (= max radius*DR)
	float       *table;     // Distortion table: M_enclosed(r) where r is measured in pixels
} lens;


typedef struct {
	long int    x;         //x - x coord. of lower left pixel where galaxy should be embedded
	long int    y;         //y - y coord. of lower left pixel where galaxy should be embedded
	long int    nx;        //nx - width of the galaxy in pixels
	long int    ny;        //ny - height of the galaxy in pixels
	float       zs;        //zs - redshift of this galaxy
	char        *file;     //file - filepath to the FITS file for this galaxy
	char        *outfile;  //outfile - filepath for the output FITS file for this galaxy

} galaxy;

typedef struct {
	float       xmin;        //smallest x value
	float       xmax;        //largest x value
	float       ymin;        //smallest y value
	float       ymax;        //largest y value
} rect;


// void get_alpha(long int x, long int y, int nlenses, lens* lenses, float* alphax, float* alphay);
void get_alpha(long int x, long int y, float scale, float* alpha_map1, float* alpha_map2, float* alphax, float* alphay);
float angular_di_dist(float z1, float z2);
void print_rect(rect *r);


int main(int argc, char *argv[]){
	/*declare variables*/
	int         ngr = 4;         //number of grids
	int         g = 8;          //exponential factor to grow by between grid sizes
	int         b[4] = {3, 6, 9, 12};       //array of log_2 of grid square sizes for bitshit ops
	int         ngx[4], ngy[4];         //width and height of grids.
	long int    nx, ny;         //width and height of the image
	long int    ngalaxies;      //number of galaxies
	long int    nlenses;        //number of lenses
	float       scale;          //pixel scale in arcseconds per pixel
	float       zl;             //lens redshift
	FILE        *galaxyfile;    //galaxy list file
	FILE        *lensfile;      //lens list file
	galaxy      *galaxies;      //array of galaxy structs
	lens        *lenses;        //array of lens structs
	char        buffer[1024], buffer2[1024];   //buffers for reading in strings
	galaxy      tempgal;        //temporary galaxy for reading in galaxies
	lens        templens;       //temporary lens for reading in lenses
	rect        *grids[ngr];    //array of pointers to the grids
	bool        *grids_t[ngr];  //array of pointers to a true table for each grid, true means that gridsquare contains information

	//fits variables
	fitsfile    *galfptr, *outfptr;         //FITS file pointers for the input and output images
	int         status = 0, naxis = 0;      //error reporting and # dimensions parameter for FITS
	long        naxes[2], fpixel[2] = {1,1};//variables for FITS to read in files
	float       *galimage, *outimage;       //input and output images as arrays

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
	sscanf(argv[5], "%f", &scale);
	sscanf(argv[6], "%f", &zl);

	//check that the width and height of the image are multiples of 4096
	if((nx & (4096-1)) != 0 || (ny & (4096-1)) != 0){
		fprintf(stderr,"Error: %i is not a multiple of 4096.\n",(nx & (4096-1))!=0 ? nx : ny);
		exit(1);
	}

	fprintf(stdout,"Running jedidistort on image of shape (%li, %li).\nPixel scale: %f.\nLens redshift: %f.\n\n", nx, ny, scale, zl);

	/* parse galaxy list file*/
	if ((galaxyfile = fopen(argv[3], "r")) == NULL) {
		fprintf(stderr, "Error opening galaxy list file \"%s\".\n", argv[3]);
		exit(1);
	}
	fprintf(stdout, "Opened galaxy list file \"%s\".\n", argv[3]);

	//count the number of galaxies
	ngalaxies = 0;
	while(fscanf(galaxyfile, "%li %li %li %li %f %s %s", &tempgal.x, &tempgal.y, &tempgal.nx, &tempgal.ny, &tempgal.zs, buffer, buffer2) > 0){
		ngalaxies++;
	}
	fprintf(stdout, "%li galaxies in \"%s\".\n", ngalaxies, argv[3]);

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

	}
	fclose(galaxyfile);

	ngal = 0;

	/*parse lens list file*/
	if ((lensfile = fopen(argv[4], "r")) == NULL) {
		fprintf(stderr, "Error opening lens list file \"%s\".\n", argv[4]);
		exit(1);
	}
	fprintf(stdout, "Opened lens list file \"%s\".\n", argv[4]);

	//count the number of lenses
	nlenses = 0;
	while(fscanf(lensfile, "%f %f %i %f %f", &templens.x, &templens.y, &templens.type, &templens.p1, &templens.p2) > 0){
		nlenses++;
	}
	fprintf(stdout, "%li lenses in \"%s\".\n", nlenses, argv[4]);

	//initialize memory for all the galaxies
	lenses = (lens *) calloc(nlenses, sizeof(lens));
	if(lenses == NULL){
		fprintf(stderr,"Error allocating memory for lens list.");
		exit(1);
	}

	//cluster mass profile variables
	float       px_per_rad = 180.0*3600.0/(PI*scale); //conversion factor from radians to pixels
	float       rad_per_px = 1/px_per_rad;          //conversion factor from pixels to radians

	//read in the lenses
	rewind(lensfile);
	int nlens;
	for(nlens = 0; nlens < nlenses; nlens++){
		fscanf(lensfile, "%f %f %i %f %f", &lenses[nlens].x, &lenses[nlens].y, &lenses[nlens].type, &lenses[nlens].p1, &lenses[nlens].p2);
		fprintf(stdout,"Lens %i:\nx: %f\ny: %f\ntype: %i\np1: %f\np2: %f\n\n",nlens,lenses[nlens].x, lenses[nlens].y, lenses[nlens].type, lenses[nlens].p1, lenses[nlens].p2);

		//find the maximum radius we need to put in the table
		//we need the further corner from the lens

		float dx = (lenses[nlens].x > nx/2) ? lenses[nlens].x : nx-lenses[nlens].x;
		float dy = (lenses[nlens].y > ny/2) ? lenses[nlens].y : nx-lenses[nlens].y;
		float dr = sqrt(dx*dx + dy*dy);
		lenses[nlens].nr = (long int) ceil(dr*DR);   //the number of table entries is the max radius * the increment
		fprintf(stdout,"dx: %f\ndy: %f\ndr: %f\nnr: %li\n", dx, dy, dr, lenses[nlens].nr);

		//allocate space for the table
		lenses[nlens].table = (float *) calloc(lenses[nlens].nr, sizeof(float));
		if(lenses[nlens].table == NULL){
			fprintf(stderr,"Error allocating memory for lens %i table.\n", nlens);
		}

		//special case at zero
		switch (lenses[nlens].type){
			case 1: //Singular Isothermal Sphere profile
				{
					long int     a; //counter in the table in units of [px/DR]
					lenses[nlens].table[0] = 0;     //special case at zero
					float alpha = (lenses[nlens].p1/C);
					alpha = px_per_rad*4*PI*alpha*alpha;
					for(a = 1; a < lenses[nlens].nr; a++){
						lenses[nlens].table[a] = alpha*((float) DR)/((float) a);
						//fprintf(stdout,"alpha: %f\n", alpha);
					}
				}
				break;
			case 2: //Navarro-Frenk-White profile
				{
					long int     a; //counter in the table in units of [px/DR]
					lenses[nlens].table[0] = 0;     //special case at zero
					//distance to the lens [Mpc]
					float   Dl = angular_di_dist(0, zl);

                    //calculate the concentration as a function of redshift zl
                    //lenses[nlens].p2 = (6.71/pow((1+zl),0.44))*pow((0.678*lenses[nlens].p1*1E14/(2E12)),-0.091);

					//rs parameter for converting from angle to unitless x variable
					//float   rs = cbrtf(G/(H0*H0))* (10.0/lenses[nlens].p2);
					//float   rs = 0.978146*cbrtf(lenses[nlens].p1)/lenses[nlens].p2;
                    //calculate the Hubble constant H as a function of zl, then substitute into the equation of rs(z); this is a correction to the above equation
                    float   rs = 0.978146*cbrtf(lenses[nlens].p1)* cbrtf(1/(OMEGA_M*(1+(1+zl)/OPZEQ)*(1+zl)*(1+zl)*(1+zl)+OMEGA_D))/ lenses[nlens].p2;

                    float   x_per_rad = Dl/rs;
					//slightly modified version of NFW delta_c parameter
					float rdelta_c = 1/(logf(1+lenses[nlens].p2) - (lenses[nlens].p2/(1+lenses[nlens].p2)));

					//constant prefactor for bending angle
					//float NFW_prefactor = px_per_rad*4*G*lenses[nlens].p1*rdelta_c/(9*Dl*10000);
					float NFW_prefactor = px_per_rad*4*G*lenses[nlens].p1*rdelta_c/(9*Dl*1E5);


					fprintf(stdout,"\nDl: %f\nrs: %f\nrdelta_c: %f\ndelta_c: %f\nNFW_prefactor: %e\n\n", Dl, rs, rdelta_c, (200.0/3.0)*lenses[nlens].p2*lenses[nlens].p2*lenses[nlens].p2*rdelta_c, NFW_prefactor);
					for(a = 1; a < lenses[nlens].nr; a++){
						//convert to radians
						double theta = (double) a *rad_per_px/((double) DR);
						//convert to unitless x variable
						double xvar = theta *( (double) x_per_rad);
						//calculate the enclosed mass for this distance, up to a multiplicative constant
						double Menc = log(xvar/2);
						if(xvar > 0 && xvar < 1)
							Menc += (2/sqrt(1-xvar*xvar))*atanh(sqrt((1-xvar)/(1+xvar)));
						else if(xvar == 1)
							Menc += 1;
						else if(xvar > 1)
							Menc += (2/sqrt(xvar*xvar-1))*atan(sqrt((xvar-1)/(1+xvar)));
						else{
							fprintf(stdout,"Error: could not calculate lens %i profile. xvar outside of range (0, infty). xvar: %f", nlens, xvar);
							exit(1);
						}

						//put it all together;
						double b = ((double) NFW_prefactor )*Menc*((double) DR)/(theta*a);
						float b_f = (float) b;
						lenses[nlens].table[a] = b_f;
						//if((a % 10)==0)
						//  fprintf(stdout,"table[%i]: alpha:%e\ttheta: %e\txvar: %e\tMenc: %e\n", a, b,theta, xvar, Menc);
					}
				}
				break;
			case 3: //Navarro-Frenk-White profile
				{

					//yes, this is hard-coded. since you can change the distance to center, overall mass doesn't really matter that much.
					float mass = 10;

					long int     a; //counter in the table in units of [px/DR]
					lenses[nlens].table[0] = 0;     //special case at zero
                    /**********************************************************************test with SIS profile for grid simulation*/
					//distance to the lens [Mpc]
					float   Dl = angular_di_dist(0, zl);

                    //calculate the concentration as a function of redshift zl
                    //lenses[nlens].p2 = (6.71/pow((1+zl),0.44))*pow((0.678*20*1E14/(2E12)),-0.091);

					//rs parameter for converting from angle to unitless x variable
					//float   rs = cbrtf(G/(H0*H0))* (10.0/lenses[nlens].p2);
					float   rs = 0.978146*cbrtf(mass)/lenses[nlens].p2;
					float   x_per_rad = Dl/rs;
					//slightly modified version of NFW delta_c parameter
					float rdelta_c = 1/(logf(1+lenses[nlens].p2) - (lenses[nlens].p2/(1+lenses[nlens].p2)));

					//constant prefactor for bending angle
					float NFW_prefactor = px_per_rad*4*G*mass*rdelta_c/(9*Dl*1E5);


					fprintf(stdout,"\nDl: %f\nrs: %f\nrdelta_c: %f\ndelta_c: %f\nNFW_prefactor: %e\n\n", Dl, rs, rdelta_c, (200.0/3.0)*lenses[nlens].p2*lenses[nlens].p2*lenses[nlens].p2*rdelta_c, NFW_prefactor);


					//convert to radians
					double theta = (double) lenses[nlens].p1 *rad_per_px;
					//convert to unitless x variable
					double xvar = theta *( (double) x_per_rad);
					//calculate the enclosed mass for this distance, up to a multiplicative constant
					double Menc = log(xvar/2);
					if(xvar > 0 && xvar < 1)
						Menc += (2/sqrt(1-xvar*xvar))*atanh(sqrt((1-xvar)/(1+xvar)));
					else if(xvar == 1)
						Menc += 1;
					else if(xvar > 1)
						Menc += (2/sqrt(xvar*xvar-1))*atan(sqrt((xvar-1)/(1+xvar)));
					else{
						fprintf(stdout,"Error: could not calculate lens %i profile. xvar outside of range (0, infty). xvar: %f", nlens, xvar);
						exit(1);
					}

					for(a = 1; a < lenses[nlens].nr; a++){
						//put it all together;
						double b = ((double) NFW_prefactor )*Menc*((double) DR)/(theta*a);
						float b_f = (float) b;
						lenses[nlens].table[a] = b_f;
						//if((a % 10)==0)
						//fprintf(stdout,"%i\t%e\t%e\t%e\t%e\n", a, b,theta, xvar, Menc);
					}
                     /***************************************************************************/
                     /***************************************************************************
                    float alpha = (1000/C); //Use sigma_v=200km/s as an example
					alpha = px_per_rad*4*PI*alpha*alpha;
					for(a = 1; a < lenses[nlens].nr; a++){
						lenses[nlens].table[a] = alpha*((float) DR)/((float) a);
						//fprintf(stdout,"alpha: %f\n", alpha);
					}
                      ***************************************************************************/
				}
				break;

		}
		fprintf(stdout,"Calculated profile for lens %i.\n\n", nlens);
	}
	fclose(lensfile);





	/**************************************************************************/
	//Variables to process DC2 alpha map files
	FILE        *alpha_file1, *alpha_file2;
	char        alpha_buffer01[1024], alpha_buffer02[1024];
	char        *buf_val1, *buf_val2;
	float       *alpha_arr1, *alpha_arr2;
	int         acounter_i, acounter_j;

	//Open DC2 alpha map files
	alpha_file1 = fopen(AL_FILE1, "rb");  //hard-coded file name!!!
	if (!alpha_file1){
		fprintf(stderr, "Unable to open alpha file1.");
		return 0;
	} else {
		fprintf(stdout,"Read in alpha file1 successfully! \n");
	}

	alpha_file2 = fopen(AL_FILE2, "rb");  //hard-coded file name!!!
	if (!alpha_file2){
		fprintf(stderr, "Unable to open alpha file2.");
		return 0;
	} else {
		fprintf(stdout,"Read in alpha file2 successfully! \n");
	}

	//Allocate memory
	alpha_arr1 = (float *) calloc(MAP_SIZE*MAP_SIZE, sizeof(float));
	alpha_arr2 = (float *) calloc(MAP_SIZE*MAP_SIZE, sizeof(float));
	if (!alpha_arr1){
		fprintf(stderr, "Memory error!");
		fclose(alpha_file1);
		fclose(alpha_file2);
		return 0;
	} else {
		fprintf(stdout,"alpha arr1 memory allocated successfully! \n");
	}

	if (!alpha_arr2){
		fprintf(stderr, "Memory error!");
		fclose(alpha_file1);
		fclose(alpha_file2);
		return 0;
	} else {
		fprintf(stdout,"alpha arr2 memory allocated successfully! \n");
	}

	//Read file contents into buffer
	acounter_i=0;
	acounter_j=0;

	while(fgets(alpha_buffer01, 1024, alpha_file1) != NULL){
		if(alpha_buffer01[0] != '#'){
			buf_val1 = strtok(alpha_buffer01, "#\n");
			sscanf(buf_val1, "%f", &(alpha_arr1[acounter_i]));
			acounter_i++;
		}
	}

	while(fgets(alpha_buffer02, 1024, alpha_file2) != NULL){
		if(alpha_buffer02[0] != '#'){
			buf_val2 = strtok(alpha_buffer02, "#\n");
			sscanf(buf_val2, "%f", &(alpha_arr2[acounter_j]));
			acounter_j++;
		}
	}

	fclose(alpha_file1);
	fclose(alpha_file2);










	//make the grids
	int         gr;         //the current grid
	long int    row, col;   //row and column in the current grid
	long int    srow, scol; //row and column in the sub grid
	for(gr = 0; gr < ngr; gr++){
		//allocate memory for the ith grid
		ngx[gr] = nx >> b[gr];
		ngy[gr] = ny >> b[gr];
		fprintf(stdout,"making grid %i with sieve %i: (%i,%i).\n",gr, 1 << b[gr],ngx[gr], ngy[gr]);
		grids[gr] = (rect *) calloc(ngx[gr] * ngy[gr], sizeof(rect));
		if(grids[gr] == NULL){
			fprintf(stderr,"Error: could not allocate memory for grid array %i.", gr);
			exit(1);
		}
		/*The first grid is different because we need to compute alpha directly instead of drawing on the other grids.*/
		//loop over grid
		for(row = 0; row < ngx[gr]; row++){
			for(col = 0; col < ngy[gr]; col++){
				//loop over subgrid
				rect bounding_box;
				for(srow = 0; srow < g; srow++){
					for(scol = 0; scol < g; scol++){
						if(gr==0){
							//get deflection angle
							float alphax = 0, alphay = 0;
							// get_alpha((row*g+srow) << BNSUBPX, (col*g+scol) << BNSUBPX, nlenses, lenses, &alphax, &alphay);
							get_alpha((row*g+srow) << BNSUBPX, (col*g+scol) << BNSUBPX, scale, alpha_arr1, alpha_arr2, &alphax, &alphay);
							// printf("x-y: %d %d; alpha: %f %f \n", (row*g+srow) << BNSUBPX, (col*g+scol) << BNSUBPX, alphax, alphay);
							//find maxima and minima
							if(srow==0 && scol==0){
								bounding_box.xmax = alphax;
								bounding_box.xmin = alphax;
								bounding_box.ymax = alphay;
								bounding_box.ymin = alphay;
							} else {
								if(alphay > bounding_box.ymax)
									bounding_box.ymax = alphay;
								if(alphay < bounding_box.ymin)
									bounding_box.ymin = alphay;

								if(alphax > bounding_box.xmax)
									bounding_box.xmax = alphax;
								if(alphax < bounding_box.xmin)
									bounding_box.xmin = alphax;
							}
						} else {
							//get index in finer array
							int index = (ngy[gr-1]*(col*g+scol))+row*g+srow;
							/*rect finer_bb = {grids[gr-1][index].xmin, grids[gr-1][index].xmax, grids[gr-1][index].ymin, grids[gr-1][index].ymax};*/
							rect finer_bb;
							finer_bb.xmin = grids[gr-1][index].xmin;
							finer_bb.xmax = grids[gr-1][index].xmax;
							finer_bb.ymin = grids[gr-1][index].ymin;
							finer_bb.ymax = grids[gr-1][index].ymax;
							//find maxima and minima
							if(srow ==0 && scol==0){
								bounding_box.xmax = finer_bb.xmin;
								bounding_box.xmin = finer_bb.xmax;
								bounding_box.ymax = finer_bb.ymax;
								bounding_box.ymin = finer_bb.ymin;
							} else {
								if(finer_bb.xmax > bounding_box.xmax)
									bounding_box.xmax = finer_bb.xmax;
								if(finer_bb.xmin < bounding_box.xmin)
									bounding_box.xmin = finer_bb.xmin;
								if(finer_bb.ymax > bounding_box.ymax)
									bounding_box.ymax = finer_bb.ymax;
								if(finer_bb.ymin < bounding_box.ymin)
									bounding_box.ymin = finer_bb.ymin;
							}
						}
					}
				}
				grids[gr][ngy[gr]*col + row].ymax = bounding_box.ymax;
				grids[gr][ngy[gr]*col + row].ymin = bounding_box.ymin;
				grids[gr][ngy[gr]*col + row].xmax = bounding_box.xmax;
				grids[gr][ngy[gr]*col + row].xmin = bounding_box.xmin;

				/*if(gr > 1){
				  fprintf(stdout,"gr: %i, (%i, %i)", gr, row, col);
				  print_rect(&grids[gr][ngy[gr]*col+row]);
				  }*/
			}
		}
	}



	//distort the galaxies one at a time
	long int gal;
	for(gal = 0; gal < ngalaxies; gal++){
		fprintf(stdout,"Distorting galaxy %i.\n", gal);
		//fprintf(stdout,"Dls: %f\t Ds: %f\n", angular_di_dist(zl,galaxies[gal].zs), angular_di_dist(0,galaxies[gal].zs));

		float     prefactor = angular_di_dist(zl,galaxies[gal].zs)/angular_di_dist(0,galaxies[gal].zs);          //the quantity D_ls/D_s
		//fprintf(stdout,"prefactor for galaxy %i: %f.\n", gal, prefactor);

		rect gal_box = {(float) galaxies[gal].x, (float) (galaxies[gal].x + galaxies[gal].nx), (float) galaxies[gal].y, (float) (galaxies[gal].y + galaxies[gal].ny)};
		//print_rect(&gal_box);

		//read in galaxy image
		fits_open_file(&galfptr, galaxies[gal].file, READONLY, &status);
		fits_get_img_dim(galfptr, &naxis, &status);
		fits_get_img_size(galfptr, 2, naxes, &status);
		if(status){
			fits_report_error(stderr, status);
			exit(1);
		}
		if(naxes[0] != galaxies[gal].nx || naxes[1] != galaxies[gal].ny){
			fprintf(stderr,"Error: galaxy %i is the wrong size. The galaxy list has shape (%i, %i) but the fits image has shape (%i, %i).\n", gal, galaxies[gal].nx, galaxies[gal].ny, naxes[0], naxes[1]);
			exit(1);
		}

		//fprintf(stdout,"Galaxy %i has size as follows. The galaxy list has shape (%i, %i) and the fits image has shape (%i, %i).\n", gal, galaxies[gal].nx, galaxies[gal].ny, naxes[0], naxes[1]);

		//allocate enough memory for the image
		//fprintf(stdout, "Allocating memory for image.\n");
		galimage = (float *) calloc (naxes[0]*naxes[1], sizeof(float));
		if (galimage == NULL){
			fprintf(stderr, "Error allocating memory for galaxy image %i.\n", gal);
			exit(1);
		}

		//read in the image
		//fprintf(stdout, "Reading in image.\n");
		if(fits_read_pix(galfptr, TFLOAT, fpixel, naxes[0]*naxes[1], NULL, galimage, NULL, &status)){
			fprintf(stderr, "Can't read in galaxy image %i.\n", gal);
			fits_report_error(stderr, status);
			exit(1);
		}
		fits_close_file(galfptr, &status);

		//allocate memory for the grid_t truth tables
		int gr;
		for(gr = 0; gr< ngr; gr++)
			grids_t[gr] = (bool *) calloc(ngx[gr] * ngy[gr], sizeof(bool));

		//we want to compute the smallest box that contains all the information for the finest grid
		int         ixmax = -1;          //convention: lowest row after all rows containing information
		int         ixmin = ngx[0];     //convetion: lowest row that contains information
		int         iymax = -1;          //convention: lowest column after all rows containing information
		int         iymin = ngy[0];     //convetion: lowest column that contains information

		for(gr = ngr; gr > 0; gr--){
			//loop over all grid squares
			long int     row, col, srow, scol;  //variables to loop over grids and sub-grids
			long int     gs_wd = (ngr == gr ? nx : g << (b[0]*gr));
			long int     gs_ht = (ngr == gr ? ny : g << (b[0]*gr));
			long int     s_gs_s= g << (b[0]*(gr-1));        //sub-gridsquare sieve length
			int subrows_per_row = (ngr == gr ? (nx >> b[gr-1]) : g);
			int subcols_per_col = (ngr == gr ? (nx >> b[gr-1]) : g);

			for(row = 0; row < (ngr == gr ? 1 : ngx[gr]); row++){
				for(col = 0; col < (ngr == gr ? 1 : ngy[gr]); col++){
					//check if we want ot look at this sub-gridsquare
					//always yes for the top level, otherwise look it up in the table
					if(ngr == gr || grids_t[gr][col*ngy[gr] + row]){
						//loop over all the sub-gridsquares
						//fprintf(stdout,"examining grid square (%i, %i)\n", row, col);
						for(srow = 0; srow < subrows_per_row; srow++){
							for(scol = 0; scol < subcols_per_col; scol++){
								//fprintf(stdout,"(%i, %i)\n", gs_wd*row+s_gs_s*srow, gs_ht*col+s_gs_s*scol);
								/*Now we want to see if the gal_box rectangle and the rectangle in the source plane from which it is possible for this grid square to draw information have an intersection. This is done using the efficient method from http://stackoverflow.com/questions/306316/determine-if-two-rectangles-overlap-each-other. An epsilon margin is added around gal_box, in case of rounding issues or other losses of precision.

								  The rect in the source plane is given by (note reversal of alpha min and max)
								  x'_min = x_min - prefactor * alpha.xmax
								  x'_max = x_max - prefactor * alpha.xmin */

								long int index = ((col << b[0]) + scol)*ngy[gr-1] + (row << b[0]) + srow;

								rect *alpha_box = &grids[gr-1][index];
								//print_rect(alpha_box);
								/*
								   fprintf(stdout,"CHECK: gr: %i\t row: %i\t col: %i\t x1: %i\t y1: %i\t x2: %i\t y2: %i\n", gr-1, (row << b[0]) + srow, (col << b[0]) + scol,((row << b[0]) + srow) << b[gr-1], ((col << b[0]) + scol) <<b[gr-1],((row << b[0]) + srow+1) << b[gr-1], ((col << b[0]) + scol+1) << b[gr-1]  );
								   print_rect(alpha_box);
								   fprintf(stdout,"xmin: %f\tymin: %f\txmax: %f\tymax: %f\n",gs_wd*row+s_gs_s*srow - prefactor*(*alpha_box).xmax, gs_ht*col+s_gs_s*scol - prefactor*(*alpha_box).ymax, gs_wd*row+s_gs_s*(srow+1) - prefactor*(*alpha_box).xmin,gs_ht*col+s_gs_s*(scol+1) - prefactor*(*alpha_box).ymin);
								   print_rect(&gal_box);
								   */
								//using 4 if's for readability, but one would suffice
								//let rect A be the source box and B be the gal box
								//check if A's Left Edge is to the left of B's right edge AND
								if(gs_wd*row+s_gs_s*srow - prefactor*(*alpha_box).xmax < gal_box.xmax){
									//check if A's right edge is to right of B's left edge AND

									if(gs_wd*row+s_gs_s*(srow+1) - prefactor*(*alpha_box).xmin > gal_box.xmin){
										//check if A's top is above B's bottom AND
										if(gs_ht*col+s_gs_s*(scol+1) - prefactor*(*alpha_box).ymin > gal_box.ymin){
											//check if A's bottom is below B's Top
											if(gs_ht*col+s_gs_s*scol - prefactor*(*alpha_box).ymax < gal_box.ymax){
												grids_t[gr-1][index] = TRUE;
												/*fprintf(stdout,"TRUE:gr: %i\t row: %i\t col: %i\t x1: %i\t y1: %i\t x2: %i\t y2: %i\n", gr-1, (row << b[0]) + srow, (col << b[0]) + scol,((row << b[0]) + srow) << b[gr-1], ((col << b[0]) + scol) <<b[gr-1],((row << b[0]) + srow+1) << b[gr-1], ((col << b[0]) + scol+1) << b[gr-1]  );
												  print_rect(alpha_box);
												  fprintf(stdout,"xmin: %f\tymin: %f\txmax: %f\tymax: %f\n",gs_wd*row+s_gs_s*srow - prefactor*(*alpha_box).xmax, gs_ht*col+s_gs_s*scol - prefactor*(*alpha_box).ymax, gs_wd*row+s_gs_s*(srow+1) - prefactor*(*alpha_box).xmin,gs_ht*col+s_gs_s*(scol+1) - prefactor*(*alpha_box).ymin);*/
												//find information bounding box for finest grid
												if(gr == 1){
													if((row << b[0]) + srow > ixmax)
														ixmax = (row << b[0]) + srow;
													if((row << b[0]) + srow < ixmin)
														ixmin = (row << b[0]) + srow;

													if((col << b[0]) + scol > iymax)
														iymax = (col << b[0]) + scol;
													if((col << b[0]) + scol < iymin)
														iymin = (col << b[0]) + scol;
												}
											}
										}
									}
								}
								//fprintf(stdout,"\n");
							}
						}
					}
				}
			}
		}

		ixmax++;    //we want one beyond the max
		iymax++;    //same
		//so many braces!
		/*j=0;
		  gr =0;
		  long int row1, col1;
		  for(row1 = 0; row1 < ngx[gr]; row1++)
		  for(col1 = 0; col1 < ngy[gr]; col1++)
		  if(grids_t[gr][col1*ngy[gr]+row1]){
		  j++;
		  fprintf(stdout,"%i: (%i, %i)\n", j, row1, col1);
		  }*/
		//fprintf(stdout,"Information box: xmin: %i; xmax: %i; ymin: %i; ymax: %i.\n",ixmin, ixmax, iymin, iymax);


		//if there's no grid squares that contain any information, set up the min and max so that the rest of the program does nothing
		if(ixmin > ixmax || iymin > iymax){
			//fprintf(stdout,"Warning: galaxy %i is completely of the screen.\n", gal);
			ixmin = 0;
			ixmax = 1;
			iymin = 0;
			iymax = 1;
		}


		/*Now we want to loop over each of the grid squares of the finest grid and see which ones contain information. We then determine the smallest box which contains all of them, and make a postage stamp output image which contains all information for the larger image. For each grid square that contains information, we want to loop over each individual pixel and break it into subpixels. For each subpixel, compute alpha and, via the lens equation beta = theta - alpha(theta), return the value of the pixel in the source plane which the lens equation points to. The value of the observed image at that pixel will be the average of all the subpixels.*/

		long int        onaxes[2] = {(ixmax - ixmin) << b[0],(iymax - iymin) << b[0]};      //(width,height) of output image in pixels.
		long int        ofx  = ixmin << b[0];                //x pixel of lower left corner of output image
		long int        ofy  = iymin << b[0];                //y pixel of lower left corner of output image


		//allocate memory for the output image
		outimage = (float *) calloc(onaxes[0]*onaxes[1], sizeof(float));
		if(outimage == NULL){
			fprintf(stderr,"Error allocating memory for output image %i.\n", gal);
			exit(1);
		}
		//fprintf(stdout,"onaxes: (%i, %i)\n", onaxes[0], onaxes[1]);

		//make output image WILL OVERWRITE EXISTING FILES
		char name[50];
		sprintf(name, "%s",galaxies[gal].outfile);
		//fprintf(stdout,"Making output image %i: %s.\n", gal, name);
		fits_create_file(&outfptr, name, &status);
		fits_create_img(outfptr, FLOAT_IMG, naxis, onaxes, &status);


		//set up some counters
		long int    row, col;       //row and column in the smallest grid
		long int    ox, oy;         //output image pixels
		long int    x, y;           //calculated source x and y pixels
		int         sx, sy;            //subpixels


		//Pretabulate some sub-pixel information to avoid some floating point calculations. dy would look the same as dx, so we don't need to bother.
		float dx[NSUBPX];
		for(sx = 0; sx < NSUBPX; sx++)
			dx[sx] = sx*(1.0/NSUBPX) + (1.0/(2*NSUBPX)) - 0.5;



		//loop over the finest grid
		for(row = ixmin; row < ixmax; row++){
			for(col = iymin; col < iymax; col++){
				if(grids_t[0][col*ngy[0]+row]){
					//calculate every pixel
					for(ox = (row << b[0]); ox < ((row+1) << b[0]); ox++){
						for(oy = (col << b[0]); oy < ((col+1) << b[0]); oy++){
							float     outvalue = 0;   //output value of this pixel
							//find the value for each subpixel

							for(sx = 0; sx < NSUBPX; sx++){
								for(sy = 0; sy < NSUBPX; sy++){
									//get alpha vector
									float alpha_x = 0, alpha_y = 0;
									// get_alpha((ox << BNSUBPX) + sy, (oy << BNSUBPX) +sy, nlenses, lenses, &alpha_x, &alpha_y);
									get_alpha((row*g+srow) << BNSUBPX, (col*g+scol) << BNSUBPX, scale, alpha_arr1, alpha_arr2, &alpha_x, &alpha_y);
									//calculate
									x = (long int) (ox + (((float) sx)/NSUBPX) - prefactor*alpha_x);
									y = (long int) (oy + (((float) sy)/NSUBPX) - prefactor*alpha_y);

									if(x >= gal_box.xmin && x < gal_box.xmax && y >= gal_box.ymin && y < gal_box.ymax){
										// the galimage array doesn't start at the origin, so we have to shift over
										x = x - (long int) gal_box.xmin;
										y = y - (long int) gal_box.ymin;
										// if(y > 295 && y < 305 && x > 295 && x < 305) fprintf(stdout,"(%i, %i): prefactor: %f\t x shift: %f\t y shift: %f\n", ox, oy, prefactor, prefactor*alpha_x, prefactor*alpha_y);
										//fprintf(stdout, "(%i, %i)\n", x, y);
										if(y*naxes[0]+x > naxes[0]*naxes[1]){
											fprintf(stderr,"Error: index too large. index: %li; max index: %li; ox: %li; oy: %li; x: %li; y: %li; naxes[0]: %li; naxes[1]: %li\n", y*naxes[0]+x, naxes[0]*naxes[1], ox, oy, x, y, naxes[0], naxes[1]);
											exit(1);
										} else
											outvalue += galimage[y*naxes[0] + x];
									}
								}
							}
							//outvalue = galimage[(ox-ofx)*naxes[0]+(oy-ofy)]; //replicate input image for debugging

							outimage[(oy -ofy)*onaxes[0] + (ox-ofx)] =  outvalue/(NSUBPX*NSUBPX);
						}
					}
				}
			}
		}

		//write header information for the output image
		fits_update_key(outfptr, TLONG, "XEMBED", &ofx , "x pixel in target image to embed lower left pixel.", &status);

		fits_update_key(outfptr, TLONG, "YEMBED", &ofy , "y pixel in target image to embed lower left pixel.", &status);

		//write output image
		fprintf(stdout, "Writing output image %i.\n", gal);
		fits_write_pix(outfptr, TFLOAT, fpixel, onaxes[0]*onaxes[1], outimage, &status);
		fits_close_file(outfptr, &status);
		fits_report_error(stderr, status);

		//release memory
		free(galimage);
		free(outimage);
		for(gr = 0; gr< ngr; gr++)
			free(grids_t[gr]);

	}

	free(alpha_arr1);
	free(alpha_arr2);
	return 0;
}



//given a  pixel (x,y) -- in subpixel, and the list of lenses,
//returns the vector (alpha_x,alpha_y) at that pixel
// void get_alpha(long int x, long int y, int nlenses, lens* lenses, float* alphax, float* alphay){
// 	long int         nlens;      //counter
// 	for(nlens = 0; nlens < nlenses; nlens++){
// 		float dx = ((float) x)/NSUBPX - lenses[nlens].x;
// 		float dy = ((float) y)/NSUBPX - lenses[nlens].y;
//
// 		//the table entry in the lens table is the dist from lens center to Point of Interest
// 		//times the table interval (DR per integer), rounded to an integer
// 		// (int) floating_point + 0.5 is a quick way to round a float to an int
// 		long int rad = (long int) (DR*sqrt(dx*dx +dy*dy)+0.5);
// 		//we want to add alpha(r) * (x/r) but the table stores alpha(r)/r to speed things up
// 		//where r is measured in px/NSUBPX
// 		*alphax += lenses[nlens].table[rad]*dx;
// 		*alphay += lenses[nlens].table[rad]*dy;
// 	}
// }
void get_alpha(long int x, long int y, float scale, float* alpha_map1, float* alpha_map2, float* alphax, float* alphay){
	int map_coord;
	map_coord = round(((float) x)/NSUBPX *MAP_SIZE + ((float) y)/NSUBPX);
	*alphax += alpha_map1[map_coord]/scale;
	*alphay += alpha_map2[map_coord]/scale;
}



//given two redshifts,
//returns the angular diameter distance between them in Mpc for a set cosmology
float angular_di_dist(float z1, float z2){
	float dist = 0, z;
	float dz = 0.001;
	for(z = 1+z1; z < 1+z2; z+= dz){
		dist += dz/sqrt(OMEGA_M*(1.0+(z/OPZEQ))*(z*z*z) + OMEGA_D);
	}
	dist = C_H_0*dist/z;
	return dist;
}

//for debugging
void print_rect(rect* r){
	fprintf(stdout,"xmin: %f\t ymin: %f\t xmax: %f\t  ymax: %f\n", (*r).xmin, (*r).ymin, (*r).xmax, (*r).ymax);
}
