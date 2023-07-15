#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BIN_NUMBER 300
#define TIME_GAP 100
#define	RADIUS_MAX 6

typedef struct{
	long double x;
	long double y;
	long double z;
	int atomId;
	int atomType;
} atom;

typedef struct{
	int step;
	atom particle[100];
} timestep;

typedef struct{
	long double x;
	long double y;
	long double z;
} vector;

vector a, b, c;
int numParticle;
int maxTimestep;
timestep trajectory[100000];
timestep *perTime;

/*
	Timestep starts from 1, atomId starts from 0.
	I know it's fucked up, but VASP did it.
	Don't blame me.
*/

long double getDistance(atom, atom);

int main(int argc, char* argv[])
{
	// Argument check
	if (argc != 2){
		printf("\tUSAGE : ./supercell_RDF XDATCAR_wrap\n");
		exit(1);
	}
	FILE *input = fopen(argv[1], "r");
	
	// File check
	if (input != NULL)
		printf("Successfully loaded %s.\n", argv[1]);
	else{
		printf("Something wrong with opening %s.\n", argv[1]);
		printf("Terminating.");
		exit(2);
	}

	// Reading Trajectory
	char line[200];
	int thisStep = 1;

	int numLi, numAl, numCl;

	while(fgets(line, sizeof(line), input)){
		// printf("line : %s", line);
		if (strncmp("Direct configuration=", line, strlen("Direct configuration=")) == 0){
			perTime = (timestep*)malloc(sizeof(timestep));
			perTime->step = thisStep;
			for (int i = 0; i < numParticle; i++){
				fscanf(input, "%Lf %Lf %Lf",  
					&perTime->particle[i].x,
					&perTime->particle[i].y,
					&perTime->particle[i].z);
				perTime->particle[i].atomId = i;
				if (i < numLi)
					perTime->particle[i].atomType = 0;
				else if (i >= numLi && i < numLi+numAl)
					perTime->particle[i].atomType = 1;
				else if (i >= numLi + numAl)
					perTime->particle[i].atomType = 2;
			}
			trajectory[thisStep++] = *perTime;
			free(perTime);
		}
		else if (strcmp("   Li   Al   Cl\n", line) == 0){
			fscanf(input, "%d %d %d", &numLi, &numAl, &numCl);
			numParticle = numLi + numAl + numCl;
			// printf("numParticle : %d\n", numParticle);
		}
		else if (strcmp("   Li   Cl   Al\n", line) == 0){
			fscanf(input, "%d %d %d", &numLi, &numCl, &numAl);
			numParticle = numLi + numAl + numCl;
			// printf("numParticle : %d\n", numParticle);
		}
		else if (strcmp("           1\n", line) == 0){
			fscanf(input, "%Lf %Lf %Lf",
					&a.x, &a.y, &a.z);
			fscanf(input, "%Lf %Lf %Lf",
					&b.x, &b.y, &b.z);
			fscanf(input, "%Lf %Lf %Lf",
					&c.x, &c.y, &c.z);
			/*
			printf("a : %Lf %Lf %Lf\n", a.x, a.y, a.z);
			printf("b : %Lf %Lf %Lf\n", b.x, b.y, b.z);
			printf("c : %Lf %Lf %Lf\n", c.x, c.y, c.z);
			*/
		}
	}

	maxTimestep = thisStep-1;
	
	printf("Reading trajectory from \"%s\" is complete.\n", argv[1]);
	printf("Number of Particles : %d\n", numParticle);
	printf("Maximum Timestep : %d\n", maxTimestep);

	// Calculating RDF
	FILE *output = fopen("./supercell_RDF.dat","w");
	if (output != NULL){
		fprintf(output, "#\ttimestep\tr\tg(r)\n");
	}
	else{
		printf("Error occured while opening \"supercell_RDF.cat\" file.\n");
		exit(4);
	}

	
	int hist[BIN_NUMBER] = { 0 };
	long double distance;
	long double binSize = (long double)RADIUS_MAX / (long double)BIN_NUMBER;
	int index;

	for (int i = 1; i <= maxTimestep; i++){
		// initializing histogram every 100 step
		if (i % TIME_GAP == 1)
			for(int j = 0; j < BIN_NUMBER; j++) hist[j] = 0;

		for (int j = 0; j < numParticle-1; j++){
			for (int k = j + 1; k < numParticle; k++){
				distance = getDistance(trajectory[i].particle[j],
						       trajectory[i].particle[k]);
				if (distance < (long double)RADIUS_MAX){
					index = (int)(distance / binSize);
					hist[index] = hist[index] + 2;
				}
			}
		}

		long double r, rdf;
		if (i % TIME_GAP == 0){
			for (int j = 0; j < BIN_NUMBER; j++){
				r = (j+0.5) * binSize;
				rdf = *(hist+j) / (long double)TIME_GAP;
				rdf /= (4 * M_PI * r * r * binSize);
				rdf /= (numParticle * (numParticle -1) / (a.x * b.y * c.z));
				// fprintf(stdout, "%d\t%Lf\t%10.8Lf\n" , i, r, rdf);
				fprintf(output, "%d\t%Lf\t%10.8Lf\n" , i, r, rdf);
			}
			printf("Now completed %d timestep...\n", i);
		}
	}

	fclose(input);
	fclose(output);

	return 0;
}

long double getDistance(atom a1, atom a2){
	long double da, db, dc;
	long double dx, dy, dz;
	long double distance;

	da = a1.x - a2.x;
	db = a1.y - a2.y;
	dc = a1.z - a2.z;

	da = da - round(da);
	db = db - round(db);
	dc = dc - round(dc);

	// printf("%Lf\t%Lf\t%Lf\n", da, db, dc);

	dx = a.x * da + b.x * db + c.x * dc;
	dy = a.y * da + b.y * db + c.y * dc;
	dz = a.z * da + b.z * db + c.z * dc;

	// printf("%Lf\t%Lf\t%Lf\n", dx, dy, dz);

	distance = dx*dx + dy*dy + dz*dz;
	distance = sqrtl(distance);

	// printf("%10.8Lf\n", distance);

	return distance;
}
