#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define RCUT 2.0
#define SIZE_OF_LINE_MAX 100
#define MAX_TIMESTEP 10000000

int timestep = 0;
int *atomid, *atomtype, *m;
double *x, *y, *z;
double pbc_x, pbc_y, pbc_z;
size_t NUM_ATOMS;
int MAX_NUM_TYPE;
int MAX_CLUSTER_ID;
int *mol_id;

FILE *fp_in;
FILE *fp_out_lmp;
FILE *fp_out_mol;

void checkSyntax(int argc, char *argv[])
{
	char *input_file_name = NULL;
	char *output_file_name = NULL;
	char *output_mol_file_name = NULL;

	if (argc == 1){
		printf("USAGE : ./simple_cluster -i INPUT.lammpstrj");
		printf(" -o OUTPUT_PREFIX\n");
		exit(1);
	}

	for (int i = 0; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			input_file_name = argv[i+1];
		else if (!strcmp(argv[i], "-o"))
		{
			char temp1[SIZE_OF_LINE_MAX], temp2[SIZE_OF_LINE_MAX];
			strcpy(temp1, argv[i+1]);
			strcpy(temp2, argv[i+1]);
			output_file_name = strcat(temp1, "_clustered.lammpstrj");
			output_mol_file_name = strcat(temp2, "_molecule.dat");
		}

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
	{
		fp_out_lmp = fopen(output_file_name, "w");
		fp_out_mol = fopen(output_mol_file_name, "w");
	}
	else
	{
		fp_out_lmp = fopen("result_clustered.lammpstrj", "w");
		fp_out_mol = fopen("result_molecule.dat", "w");
	}
	fprintf(fp_out_mol, "# timestep mol_id type_1 type_2 type_3\n");
}
void readInfos()
{
	char line[SIZE_OF_LINE_MAX];
	double max, min;

	while(fgets(line, sizeof(line), fp_in))
	{
		if (!strcmp(line, "ITEM: TIMESTEP\n"))
		{
			fputs(line, fp_out_lmp);
			fscanf(fp_in, "%d", &timestep);
			fprintf(fp_out_lmp, "%d\n", timestep);
		}
		else if (!strcmp(line, "ITEM: NUMBER OF ATOMS\n"))
		{
			fputs(line, fp_out_lmp);
			fscanf(fp_in, "%lu", &NUM_ATOMS);
			fprintf(fp_out_lmp, "%lu\n", NUM_ATOMS);

			atomtype = malloc(NUM_ATOMS * sizeof(int));
			atomid = malloc(NUM_ATOMS * sizeof(int));
			mol_id = malloc(NUM_ATOMS * sizeof(int));

			x = malloc(NUM_ATOMS * sizeof(double));
			y = malloc(NUM_ATOMS * sizeof(double));
			z = malloc(NUM_ATOMS * sizeof(double));
		}
		else if (!strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n"))
		{
			fputs(line, fp_out_lmp);
			fscanf(fp_in, "%lf %lf", &min, &max);
			fprintf(fp_out_lmp, "%lf %lf\n", min, max);
			pbc_x = max - min;
			fscanf(fp_in, "%lf %lf", &min, &max);
			fprintf(fp_out_lmp, "%lf %lf\n", min, max);
			pbc_y = max - min;
			fscanf(fp_in, "%lf %lf", &min, &max);
			fprintf(fp_out_lmp, "%lf %lf\n", min, max);
			pbc_z = max - min;
		}
		else if (!strcmp(line, "ITEM: ATOMS id type x y z\n"))
		{
			MAX_NUM_TYPE = 0;
			for (int i = 0; i < NUM_ATOMS; i++)
			{
				fscanf(fp_in,
				       "%d %d %lf %lf %lf",
				       &atomid[i], &atomtype[i],
				       &x[i], &y[i], &z[i]);
				mol_id[i] = 0;
				if (atomtype[i] > MAX_NUM_TYPE)
					MAX_NUM_TYPE = atomtype[i];
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
		if (mol_id[i] != 0 || !mol_id) continue;
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
	MAX_CLUSTER_ID = id;
}
void writeDataLammps()
{
	fprintf(fp_out_lmp, "ITEM: ATOMS id type x y z mol_id\n");
	for (int i = 0; i < NUM_ATOMS; i++)
		fprintf(fp_out_lmp,
		       "%d %d %lf %lf %lf %d\n",
		       atomid[i], atomtype[i],
		       x[i], y[i], z[i],
		       mol_id[i]);
}
void checkMolecule()
{
	int *type = malloc(MAX_NUM_TYPE * sizeof(int)); // array init
	for (int i = 1; i < MAX_CLUSTER_ID; i++)
	{
		for (int k = 0; k < MAX_NUM_TYPE; k++) *(type + k) = 0;
		// array init

		for (int j = 0; j < NUM_ATOMS; j++)
			if (mol_id[j] == i) *(type + atomtype[j] - 1) += 1; 

		fprintf(fp_out_mol, "%d %d ", timestep, i);
		for (int k = 0; k < MAX_NUM_TYPE; k++)
			fprintf(fp_out_mol, "%d ", *(type + k));
		fprintf(fp_out_mol, "\n");
	}
	free(type);
}
void perTimestep()
{
	readInfos();
	checkCluster();
	writeDataLammps();
	checkMolecule();

	free(atomid); free(atomtype);
	free(x); free(y); free(z);
	free(mol_id);

	atomid = NULL; atomtype = NULL;
	x = NULL; y = NULL; z = NULL;
	mol_id= NULL;
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
	fclose(fp_out_lmp);
	fclose(fp_out_mol);

	return 0;
}
