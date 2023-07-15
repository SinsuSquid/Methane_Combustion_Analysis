#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[])
{
	// Argument check
	if (argc != 2){
		printf("\tUSAGE : ./supercell_wraper XDATCAR\n");
		exit(1);
	}
	FILE *input = fopen(argv[1], "r");
	
	// File check
	if (input)
		printf("Successfully loaded %s.\n", argv[1]);
	else{
		printf("Something wrong with opening %s.\n", argv[1]);
		printf("Terminating.");
		exit(2);
	}

	FILE *output = fopen("./XDATCAR_wrap","w");
	char line[200];

	while(fgets(line, sizeof(line), input) != NULL){
		fprintf(output, "%s", line);
		if (strncmp(line, "Direct configuration=",
			    strlen("Direct configuration=")) == 0) break;
	}

	long double a, b, c;

	while(fgets(line, sizeof(line), input) != NULL){
		if (strncmp(line, "Direct configuration=",
			    strlen("Direct configuration=")) == 0){
			fprintf(output, "%s", line);
		}
		else{
			sscanf(line, "%LF %LF %LF", &a, &b, &c);
			if (a < 0) a = 1 + a;
			else if (a >= 1) a = a - 1;
			if (b < 0) b = 1 + b;
			else if (b >= 1) b = b - 1;
			if (c < 0) c = 1 + c;
			else if (c >= 1) c = c - 1;

			fprintf(output, "%10.8LF %10.8LF %10.8LF\n", a, b, c);
		}
	}

	printf("All task completed.\n");

	fclose(input);
	fclose(output);

	return 0;
}
