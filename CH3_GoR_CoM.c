#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_SIZE 100
#define MAX_NUM_MOLECULES 1000
#define INDEX_MAX 1000
// #define TOTAL_TRAJECTORY 10000000
#define TOTAL_TRAJECTORY 2000000
#define PLOTTING_GAP      500000

FILE *fp_traj;
FILE *fp_cluster;
FILE *fp_out;

int timestep;
int isCH3[MAX_NUM_MOLECULES] = {0};
int numCH3 = 0;
char line[MAX_LINE_SIZE];
int MAX_CLUSTER_ID;
int NUM_ATOMS;
double *x, *y, *z;
double pbc_x, pbc_y, pbc_z;
int *atom_id, *atom_type, *mol_id;
int histogram[INDEX_MAX];
double *x_com, *y_com, *z_com;

void checkSyntax(int, char **);
void readCH3();
void readTraj();
void perTraj();
void getHist();
double getDistance(char, double, double);
double getDistance_upgrade(char, double, double);
void getCoM();

int main(int argc, char *argv[])
{
	clock_t startT, endT, beginT;
	double perTime, totTime;

	checkSyntax(argc, argv);
	beginT = clock();
	for(; timestep < TOTAL_TRAJECTORY;)
	{
		startT = clock();
		perTraj();
		endT = clock();
		perTime = (double)(endT - startT) / CLOCKS_PER_SEC;
		totTime = (double)(endT - beginT) / CLOCKS_PER_SEC;
		if (timestep % 1000000 == 0)
			printf("Timestep : %10d | Per Trajectory : %.3lf | Total Time : %.3lf\n", timestep, perTime, totTime);
	}

	fclose(fp_traj); fclose(fp_cluster); fclose(fp_out);
	return 0;
}
void perTraj()
{
	readTraj();
	readCH3();
	getCoM();
	getHist();

	free(x); free(y); free(z);
	free(atom_id); free(atom_type); free(mol_id);
	free(x_com); free(y_com); free(z_com);

	x = NULL; y = NULL; z = NULL;
	atom_id = NULL; atom_type = NULL; mol_id = NULL;
	x_com = NULL; y_com = NULL; z_com = NULL;

	// if (timestep == TOTAL_TRAJECTORY) exit(0);
}
void checkSyntax(int argc, char *argv[])
{
	fp_traj = fopen("result_clustered.lammpstrj", "r");
	fp_cluster = fopen("result_molecule.dat", "r");
	fp_out = fopen("result_CH3_GoR.dat", "w");

	fprintf(fp_out, "# timestep r g(r)\n");
	fgets(line, sizeof(line), fp_cluster);
}
void readTraj()
{
	double min, max;
	while(fgets(line, sizeof(line), fp_traj))
	{
		if (strcmp(line, "ITEM: TIMESTEP\n") == 0)
		{
			fscanf(fp_traj, "%d", &timestep);
			// printf("TIMESTEP : %d\n", timestep);
		}
		else if (strcmp(line, "ITEM: NUMBER OF ATOMS\n") == 0)
		{
			fscanf(fp_traj, "%d", &NUM_ATOMS);
			// printf("NUM_ATOMS : %d\n", NUM_ATOMS);

		}
		else if (strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") == 0)
		{
			fscanf(fp_traj, "%lf %lf", &min, &max);
			pbc_x = max - min;
			fscanf(fp_traj, "%lf %lf", &min, &max);
			pbc_y = max - min;
			fscanf(fp_traj, "%lf %lf", &min, &max);
			pbc_z = max - min;
			// printf("PBCS : %lf %lf %lf\n", pbc_x, pbc_y, pbc_z);
		}
		else if (strcmp(line, "ITEM: ATOMS id type x y z mol_id\n") == 0)
		{
			MAX_CLUSTER_ID = 0;
			atom_id = malloc(sizeof(int) * NUM_ATOMS);
			atom_type = malloc(sizeof(int) * NUM_ATOMS);
			x = malloc(sizeof(double) * NUM_ATOMS);
			y = malloc(sizeof(double) * NUM_ATOMS);
			z = malloc(sizeof(double) * NUM_ATOMS);
			mol_id = malloc(sizeof(int) * NUM_ATOMS);
			for (int i = 0; i < NUM_ATOMS; i++)
			{
				fscanf(fp_traj, "%d %d %lf %lf %lf %d",
				       &atom_id[i], &atom_type[i], &x[i], &y[i], &z[i],
				       &mol_id[i]);
				if (mol_id[i] > MAX_CLUSTER_ID) 
					MAX_CLUSTER_ID = mol_id[i];
				// printf("%d %d %lf %lf %lf %d\n", 
				//        atom_id, atom_type, x[i], y[i], z[i],
				//        mol_id[i]);
			}
			return;
		}
	}
	return;
}
void readCH3()
{
	numCH3 = 0;
	int type_1, type_2, type_3;
	for (int i = 0; i < MAX_CLUSTER_ID; i++)
	{
		fgets(line, sizeof(line), fp_cluster);
		int lineLen = strlen(line);
		char *tempLine = (line+(lineLen-7));
		if (strcmp(tempLine, "1 3 0\n") == 22)
		{
			isCH3[i+1] = 1;
			numCH3++;
		}
	}
	
	/*
	for (int i = 0; i < MAX_CLUSTER_ID; i++)
		printf("%d %d %d\n", timestep, i, isCH3[i]);
	*/
	
	return;
}
void getCoM()
{
	x_com = malloc(sizeof(double) * MAX_CLUSTER_ID);
	y_com = malloc(sizeof(double) * MAX_CLUSTER_ID);
	z_com = malloc(sizeof(double) * MAX_CLUSTER_ID);

	int atomCounter = 0;
	for (int i = 1; i < MAX_CLUSTER_ID; i++)
	{
		x_com[i] = 0.0; y_com[i] = 0.0; z_com[i] = 0.0;
		atomCounter = 0;
		for (int j = 0; j < NUM_ATOMS; j++)
		{
			if (mol_id[j] != i) continue;
			x_com[i] += x[j]; y_com[i] += y[j]; z_com[i] += z[j];
			atomCounter++;
		}
		x_com[i] /= atomCounter;
		y_com[i] /= atomCounter;
		z_com[i] /= atomCounter;

		/*
		printf("timestep : %d | mol_id : %d\n", timestep, i);
		printf("x_com : %lf\ty_com : %lf\tz_com : %lf\n",
			x_com[i], y_com[i], z_com[i]);
		*/
	}
}
void getHist()
{
	double lenMax = sqrt((pbc_x * 0.5) * (pbc_x * 0.5)
			    +(pbc_y * 0.5) * (pbc_y * 0.5)
			    +(pbc_z * 0.5) * (pbc_z * 0.5));
	double binSize = lenMax / INDEX_MAX;
	double dx, dy, dz;

	// Initialize Histogram
	// 1000 step average
	if (timestep % PLOTTING_GAP == (PLOTTING_GAP - 100))
		for (int i = 0; i < INDEX_MAX; i++) histogram[i] = 0;

	for (int i = 1; i < MAX_CLUSTER_ID; i++)
	{
		if (!isCH3[mol_id[i]]) continue;
		double length = 0.0;
		for (int j = i + 1; j < MAX_CLUSTER_ID; j++)
		{
			if (!isCH3[mol_id[j]]) continue;
			dx = x_com[i] - x_com[j];
			dy = y_com[i] - y_com[j];
			dz = z_com[i] - z_com[j];

			dx = dx - round(dx / pbc_x) * pbc_x;
			dy = dy - round(dy / pbc_y) * pbc_y;
			dz = dz - round(dz / pbc_y) * pbc_z;

			length = sqrt(dx*dx + dy*dy + dz*dz);
			histogram[(int)(length / binSize)] += 2;
		}
	}

	double r = 0.0;
	double writeValue = 0.0;
	int normalize = 0;

	for (int i = 0; i < INDEX_MAX; i++)
		normalize += histogram[i];

	for (int i = 0; i < INDEX_MAX / 2; i++)
	{
		r = i * binSize;
		writeValue = (double)histogram[i];
		// dV = 4 * pi * r^2 * dr
		writeValue /= (4 * M_PI * r * r * binSize);
		// rho = N / V
		writeValue /= (MAX_CLUSTER_ID / (pbc_x * pbc_y * pbc_z)); 
		writeValue /= normalize;
		fprintf(fp_out, "%d %lf %lf\n",
			timestep, r, writeValue);
	}
}
