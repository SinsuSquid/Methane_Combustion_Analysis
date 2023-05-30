#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_SIZE 100
#define MAX_NUM_MOLECULES 1000
#define INDEX_MAX 2000
// #define TOTAL_TRAJECTORY 10000000
#define TOTAL_TRAJECTORY     2000000
#define PLOTTING_GAP          500000

int timestep;
int numCH3 = 0;
char line[MAX_LINE_SIZE];
int MAX_CLUSTER_ID;
int NUM_ATOMS;
long double pbc_x, pbc_y, pbc_z;
int histogram[INDEX_MAX];

int cumCH3;

typedef struct{
	int atom_id, atom_type, mol_id;
	long double x, y, z;
} Particle;

void checkSyntax(int argc, char *argv[], FILE **fp_traj, FILE **fp_out, FILE **fp_mol)
{
	*fp_traj = fopen("result_clustered.lammpstrj", "r");
	*fp_mol = fopen("result_molecule.dat", "r");
	fgets(line, sizeof(line), *fp_mol);
	*fp_out = fopen("result_GoR.dat", "w");


	fprintf(*fp_out, "# timestep r g(r)\n");
}

void getHist(FILE **fp_out, Particle **particles, int **isCH3)
{
	long double length, dx, dy, dz;
	long double lengthMax = sqrt((pbc_x * 0.5) * (pbc_x * 0.5)
		       	            +(pbc_y * 0.5) * (pbc_y * 0.5)
			            +(pbc_z * 0.5) * (pbc_z * 0.5));
	long double binSize = lengthMax / INDEX_MAX;

	// Initialize Histogram
	if (timestep % (PLOTTING_GAP) == (PLOTTING_GAP - 10000))
	{
		for (int i = 0; i < INDEX_MAX; i++) histogram[i] = 0;
		cumCH3 = 0;
	}

	for (int i = 0; i < NUM_ATOMS-1; i++)
	{
		if (*(*isCH3 + (*particles+i)->mol_id) != 1) continue;
		for (int j = i + 1; j < NUM_ATOMS; j++)
		{
			if (*(*isCH3 + (*particles+j)->mol_id) != 1) continue;
			dx = (*particles+i)->x -(*particles+j)->x;
			dy = (*particles+i)->y -(*particles+j)->y;
			dz = (*particles+i)->z -(*particles+j)->z;
			dx = dx - round(dx / pbc_x) * pbc_x;
			dy = dy - round(dy / pbc_y) * pbc_y;
			dz = dz - round(dz / pbc_z) * pbc_z;
			length = sqrt(dx*dx + dy*dy + dz*dz);
			histogram[(int)(length / binSize)] += 2;
		}
	}

	cumCH3 += numCH3;

	if (!(timestep % PLOTTING_GAP))
	{
		long double r = 0.0;
		long double writeValue = 0.0;

		for (int i = 0; i < INDEX_MAX / 2; i++)
		{
			r = i * binSize;
			writeValue = (long double)histogram[i];
			// dV = 4 * pi * r^2 * dr
			writeValue /= (4.0 / 3.0 * M_PI * pow(i+1,3) - pow(i,3) * pow(binSize,3));
			// rho = N / V
			writeValue /= (numCH3 * 4 / (pbc_x * pbc_y * pbc_z)); 
			writeValue *= 100;
			fprintf(*fp_out, "%d %Lf %Lf\n",
				timestep, r, writeValue);
		}
	}
}
void readTraj(FILE **fp_traj, Particle **particles)
{
	long double min, max;
	while(fgets(line, sizeof(line), *fp_traj))
	{
		if (strcmp(line, "ITEM: TIMESTEP\n") == 0)
		{
			fscanf(*fp_traj, "%d", &timestep);
			// printf("TIMESTEP : %d\n", timestep);
		}
		else if (strcmp(line, "ITEM: NUMBER OF ATOMS\n") == 0)
		{
			fscanf(*fp_traj, "%d", &NUM_ATOMS);
			// printf("NUM_ATOMS : %d\n", NUM_ATOMS);

		}
		else if (strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") == 0)
		{
			fscanf(*fp_traj, "%llf %llf", &min, &max);
			pbc_x = max - min;
			fscanf(*fp_traj, "%llf %llf", &min, &max);
			pbc_y = max - min;
			fscanf(*fp_traj, "%llf %llf", &min, &max);
			pbc_z = max - min;
			// printf("PBCS : %llf %llf %llf\n", pbc_x, pbc_y, pbc_z);
		}
		else if (strcmp(line, "ITEM: ATOMS id type x y z mol_id\n") == 0)
		{
			*particles = malloc(sizeof(Particle) * NUM_ATOMS);

			int tmp1, tmp2, tmp6;
			long double tmp3, tmp4, tmp5;
			MAX_CLUSTER_ID = 0;

			for (int i = 0; i < NUM_ATOMS; i++)
			{
				fscanf(*fp_traj, "%d %d %llf %llf %llf %d",
						&tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6);
				(*particles+i)->atom_id = tmp1;
				(*particles+i)->atom_type = tmp2;
				(*particles+i)->x = tmp3;
				(*particles+i)->y = tmp4;
				(*particles+i)->z = tmp5;
				(*particles+i)->mol_id = tmp6;

				if (tmp6 > MAX_CLUSTER_ID)
					MAX_CLUSTER_ID = tmp6;
			}
			return;
		}
	}
	return;
}

void checkCH3(FILE **fp_mol, int **isCH3)
{
	*isCH3 = malloc(sizeof(int) * (MAX_CLUSTER_ID + 1));
	numCH3 = 0;

	int t, mol, type_1, type_2, type_3;
	for (int i = 0; i < MAX_CLUSTER_ID; i++)
	{
		fscanf(*fp_mol, "%d %d %d %d %d",
			&t, &mol, &type_1, &type_2, &type_3);
		if (type_1 == 1 && type_2 == 3 && type_3 == 0)
		{
			*(*isCH3+mol) = 1;
			numCH3++;
		}
		else
			*(*isCH3+mol) = 0;

		// printf("%d %d %d\n", t, mol, *(*isCH3+mol)); 
	}
}

void perTraj(FILE **fp_traj, FILE **fp_out, FILE **fp_mol, Particle **particles, int **isCH3)
{
	readTraj(fp_traj, particles);
	checkCH3(fp_mol, isCH3);

	getHist(fp_out, particles, isCH3);
}
int main(int argc, char *argv[])
{
	FILE *fp_traj;
	FILE *fp_out;
	FILE *fp_mol;

	Particle *particles;
	int *isCH3;

	clock_t startT, endT, beginT;
	double perTime, totTime;

	checkSyntax(argc, argv, &fp_traj, &fp_out, &fp_mol);
	beginT = clock();
	for(; timestep < TOTAL_TRAJECTORY;)
	{
		startT = clock();
		perTraj(&fp_traj, &fp_out, &fp_mol, &particles, &isCH3);
		endT = clock();
		perTime = (double)(endT - startT) / CLOCKS_PER_SEC;
		totTime = (double)(endT - beginT) / CLOCKS_PER_SEC;
		if (timestep % (TOTAL_TRAJECTORY / 100) == 0)
		{
		printf("Timestep : %10d | Per Trajectory : %.3lf | Total Time : %.3lf ",
				timestep, perTime, totTime);
		printf("cumCH3 : %d\n", cumCH3);
		}
	}

	fclose(fp_traj); fclose(fp_out);
	return 0;
}
