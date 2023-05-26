#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_SIZE 100
#define INDEX_MAX 2000
#define TOTAL_TRAJECTORY 10000000

int timestep;
int numCH3 = 0;
char line[MAX_LINE_SIZE];
int MAX_CLUSTER_ID;
int NUM_ATOMS;
long double pbc_x, pbc_y, pbc_z;
long double xlo, xhi, ylo, yhi, zlo, zhi;

int cumCH3;

typedef struct{
	int atom_id, atom_type, mol_id;
	long double x, y, z;
} Particle;

void writeFile(FILE **fp_out, Particle **particles, int **isCH3)
{
	fprintf(*fp_out, "ITEM: TIMESTEP\n");
	fprintf(*fp_out, "%d\n", timestep);
	fprintf(*fp_out, "ITEM: NUMBER OF ATOMS\n");
	fprintf(*fp_out, "%d\n", NUM_ATOMS);
	fprintf(*fp_out, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(*fp_out, "%Lf %LF\n", xlo, xhi);
	fprintf(*fp_out, "%Lf %LF\n", ylo, yhi);
	fprintf(*fp_out, "%Lf %LF\n", zlo, zhi);
	fprintf(*fp_out, "ITEM: ATOMS id type x y z mol isCH3\n");

	for (int i = 0; i < NUM_ATOMS; i++)
	{
		fprintf(*fp_out, "%d %d %Lf %Lf %Lf %d %d\n",
				(*particles+i)->atom_id,
				(*particles+i)->atom_type,
				(*particles+i)->x,
				(*particles+i)->y,
				(*particles+i)->z,
				(*particles+i)->mol_id,
				*(*isCH3+(*particles+i)->mol_id));
	}
}

void checkSyntax(int argc, char *argv[], FILE **fp_traj, FILE **fp_out, FILE **fp_mol)
{
	*fp_traj = fopen("result_clustered.lammpstrj", "r");
	*fp_mol = fopen("result_molecule.dat", "r");
	fgets(line, sizeof(line), *fp_mol);
	*fp_out = fopen("result_new.lammpstrj", "w");
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
			fscanf(*fp_traj, "%llf %llf", &xlo, &xhi);
			pbc_x = xhi - xlo;
			fscanf(*fp_traj, "%llf %llf", &ylo, &yhi);
			pbc_y = yhi - ylo;
			fscanf(*fp_traj, "%llf %llf", &zlo, &zhi);
			pbc_z = zhi - zlo;
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
	writeFile(fp_out, particles, isCH3);
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
		}
	}

	fclose(fp_traj); fclose(fp_out);
	return 0;
}
