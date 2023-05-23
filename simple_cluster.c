#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define RCUT 2.0
#define SIZE_OF_LINE_MAX 100
#define MAX_TIMESTEP 10000000

int timestep = 0;
int *atomid, *atomtype, *mol_id;
double *x, *y, *z;
double pbc_x, pbc_y, pbc_z;
size_t NUM_ATOMS;

FILE *fp_in;
FILE *fp_out;

void checkSyntax(int argc, char *argv[])
{
	char *input_file_name = NULL;
	char *output_file_name = NULL;

	if (argc == 1){
		printf("USAGE : ./simple_cluster -i INPUT.lammpstrj");
		printf(" -o OUTPUT.lammpstrj\n");
		exit(1);
	}

	for (int i = 0; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			input_file_name = argv[i+1];
		else if (!strcmp(argv[i], "-o"))
			output_file_name = argv[i+1];
	}

	if (!input_file_name){
		printf("No input file found!\n");
		exit(1);
	}

	fp_in = fopen(input_file_name, "r");
	if (!fp_in)
	{
		printf("Error while opening file.\n");
		exit(1);
	}

	if(output_file_name)
		fp_out = fopen(output_file_name, "w");
	else
		fp_out = fopen("result.lammpstrj", "w");
}
void readInfos()
{
	char line[SIZE_OF_LINE_MAX];
	double max, min;

	while(fgets(line, sizeof(line), fp_in))
	{
		if (!strcmp(line, "ITEM: TIMESTEP\n"))
		{
			fputs(line, fp_out);
			fscanf(fp_in, "%d", &timestep);
			fprintf(fp_out, "%d\n", timestep);
		}
		else if (!strcmp(line, "ITEM: NUMBER OF ATOMS\n"))
		{
			fputs(line, fp_out);
			fscanf(fp_in, "%lu", &NUM_ATOMS);
			fprintf(fp_out, "%lu\n", NUM_ATOMS);

			atomtype = malloc(NUM_ATOMS * sizeof(int));
			atomid = malloc(NUM_ATOMS * sizeof(int));
			mol_id = malloc(NUM_ATOMS * sizeof(int));

			x = malloc(NUM_ATOMS * sizeof(double));
			y = malloc(NUM_ATOMS * sizeof(double));
			z = malloc(NUM_ATOMS * sizeof(double));
		}
		else if (!strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n"))
		{
			fputs(line, fp_out);
			fscanf(fp_in, "%lf %lf", &min, &max);
			fprintf(fp_out, "%lf %lf\n", min, max);
			pbc_x = max - min;
			fscanf(fp_in, "%lf %lf", &min, &max);
			fprintf(fp_out, "%lf %lf\n", min, max);
			pbc_y = max - min;
			fscanf(fp_in, "%lf %lf", &min, &max);
			fprintf(fp_out, "%lf %lf\n", min, max);
			pbc_z = max - min;
		}
		else if (!strcmp(line, "ITEM: ATOMS id type x y z\n"))
		{
			for (int i = 0; i < NUM_ATOMS; i++)
			{
				fscanf(fp_in,
				       "%d %d %lf %lf %lf",
				       &atomid[i], &atomtype[i],
				       &x[i], &y[i], &z[i]);
				mol_id[i] = 0;
			}
			return;
		}
	}
	exit(0);
}
void checkCluster()
{
	int id = 1;
	double length;
	double perAxis;

	for (int i = 0; i < NUM_ATOMS; i++)
	{
		if (mol_id[i] != 0) continue;
		mol_id[i] = id;
		for (int j = i + 1; j < NUM_ATOMS; j++)
		{
			// X axis
			if (mol_id[j] != 0) continue;
 			length = 0.0;

			perAxis = sqrt((x[i] - x[j]) * (x[i] - x[j]));
			if (perAxis > pbc_x * 0.5)
				length += (pbc_x - perAxis) * (pbc_x - perAxis);
			else
				length += perAxis * perAxis;

			// Y axis
			perAxis = sqrt((y[i] - y[j]) * (y[i] - y[j]));
			if (perAxis > pbc_y * 0.5)
				length += (pbc_y - perAxis) * (pbc_y - perAxis);
			else
				length += perAxis * perAxis;

			// Z axis
			perAxis = sqrt((z[i] - z[j]) * (z[i] - z[j]));
			if (perAxis > pbc_z * 0.5)
				length += (pbc_z - perAxis) * (pbc_z - perAxis);
			else
				length += perAxis * perAxis;

			// It's belong to a same cluster.
			if (length < RCUT * RCUT)
				mol_id[j] = id;
		}
		id++;
	}
}
void writeData()
{
	fprintf(fp_out, "ITEM: ATOMS id type x y z mol_id\n");
	for (int i = 0; i < NUM_ATOMS; i++)
		fprintf(fp_out,
		       "%d %d %lf %lf %lf %d\n",
		       atomid[i], atomtype[i],
		       x[i], y[i], z[i],
		       mol_id[i]);

	free(atomid); free(atomtype);
	free(x); free(y); free(z);
	free(mol_id);

	atomid = NULL; atomtype = NULL;
	x = NULL; y = NULL; z = NULL;
	mol_id = NULL;
}
void perTimestep()
{
	readInfos();
	checkCluster();
	writeData();
}
int main(int argc, char *argv[])
{
	checkSyntax(argc, argv);

	clock_t beginningT = clock();

	for (; timestep <= MAX_TIMESTEP ;)
	{
		clock_t startT = clock();
		perTimestep();
		clock_t endT = clock();
		float timePerTraj = (double)(endT - startT) / CLOCKS_PER_SEC;
		float totTime = (double)(endT - beginningT) / CLOCKS_PER_SEC;
		if (!(timestep % 1000000))
		{
			printf("TIMESTEP : %9d, TIME per Traj. : %.4f",
				timestep, timePerTraj);
			printf(", Total Time : %.4f\n", totTime);
		}
	}

	fclose(fp_in);
	fclose(fp_out);

	return 0;
}
