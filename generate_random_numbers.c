/*
Escrito por: Gabriel Borba
Data: 30/09/2022
Ultima atualiazacao: 03/10/2022

Powershell 7
clear && gcc ./generate_random_numbers.c -o ./generate_random_numbers.bin && ./generate_random_numbers.bin 40000 ./unsorted_numbers.txt
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(const int argc, const char **argv) {
	srand(time(0));
	if (argc!=3) {
		puts(
			"Expected to receive:\n"
			"\tan integer for the number of elements to generate\n"
			"\tthe path where to save the generated text file\n"
			"\n"
			"The output file will have:\n"
			"\ton the first line, the count of numbers on the file.\n"
			"\tcount lines, each with a integer.\n");
		exit(1);
	}
	FILE *fp = fopen(argv[2], "w");
	int n = atoi(argv[1]);
	fprintf(fp,"%d\n",n);
	for (size_t i=0;i<n;i++) {
		fprintf(fp,"%d\n", rand()%257);
	}
	fclose(fp);
	return 0;
}
