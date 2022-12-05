/*
Escrito por: Gabriel Borba
Data: 23/11/2022
Ultima atualiazacao: 30/11/2022

Powershell 7
clear && gcc -O2 bubble_sort_omp.c -o bubble_sort_omp.bin -fopenmp && [void]($Env:OMP_NUM_THREADS = 4) && .\bubble_sort_omp.bin unsorted_numbers.txt sorted_numbers.txt
*/


// =============================================================================
// SECTION - INCLUDES
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
// #include <mpi.h>
#include <omp.h>
// !SECTION


// =============================================================================
// SECTION - GLOBALS
typedef int ARRAY_DATA_TYPE;

// !SECTION


// =============================================================================
// SECTION - DECLARATIONS
/**
 * Return the start position for a slice (sId), among totalSlices, within an array
 * with nOriginal elements.
 * Considering load distribuition with slices at the start being bigger.
*/
size_t tArrayStart(
	const size_t sId,
	const size_t totalSlices,
	const size_t nOriginal);
/**
 * Return the size of a slice being sId, among totalSlices, within an array with
 *  nOriginal elements.
 * Considering load distribuition with slices at the start being bigger.
**/
size_t tArraySize(
	const size_t sId,
	const size_t totalSlices,
	const size_t nOriginal);
/**
 * Return the maximum number of elements for some slice, among totalSlices,
 * within nOriginal elements.
 */
size_t tArrayMaxSize(
	const size_t totalPieces,
	const size_t nOriginal);
// --------------------
/**
 * Sequential bubble sort algorithm with odd even transposition to aliviate data
 * dependency.
*/
void sequentialOddEvenBubbleSort(
	ARRAY_DATA_TYPE *anArr,
	const size_t n,
	const __compar_fn_t cmpFunction);
/**
 * Compare two data elements and exchange their position if the elements at i1
 * should be after i2.
*/
void compare_exchange(
	ARRAY_DATA_TYPE *anArr,
	const size_t i1,
	const size_t i2,
	const __compar_fn_t cmpFunction);
/**
 * Parallel bubble sort algorithm with odd even transposition.
*/
void parallelOddEvenBubbleSort(
	ARRAY_DATA_TYPE *anArr,
	const size_t nOriginal,
	const int nThreads,
	const __compar_fn_t cmpFunction);
/**
 * Compare function used to sort the data array in ascending order.
*/
int compareDataElements(const void *a, const void *b);
// --------------------
/**
 * Allocate an array on the heap with the desired size.
*/
void *newArray(int count, size_t elemSize);
/**
 * Output the data array to the stdout.
*/
void printDataArray(ARRAY_DATA_TYPE* arr, size_t n);
/**
 * Check whether or not the array is sorted according to the compare function.
*/
bool checkSorting(
	const void *arr,
	const size_t nmemb,
	size_t  sizememb,
	const __compar_fn_t cmpFunc);
/**
 * Return the value of x restricted between min and max.
*/
int clamp(int x, int min, int max);
// --------------------
/**
 * Program entry point.
*/
int main(int argc, char **argv);
/**
 * Read the numbers from the file.
*/
void loadData(
	char *filename,
	size_t *nOriginal,
	ARRAY_DATA_TYPE **theData);
/**
 * Write the given array to a text file on the given path.
*/
void writeArrayOnFile(
	const ARRAY_DATA_TYPE *arr,
	const size_t n,
	const char* newTextFilePath);
// !SECTION


// =============================================================================
// SECTION - IMPLEMENTATIONS
size_t tArrayStart(
	const size_t sId,
	const size_t totalSlices,
	const size_t nOriginal)
{
	const size_t leftOut = nOriginal%totalSlices;
	const size_t maxSize = tArrayMaxSize(totalSlices,nOriginal);
	if (!leftOut) {
		return sId*tArraySize(sId,totalSlices,nOriginal);
	}
	else if (sId<leftOut) {
		return sId*maxSize;
	}
	else {
		return leftOut*maxSize+(sId-leftOut)*(maxSize-1);
	}
}


size_t tArraySize(
	const size_t sId,
	const size_t totalSlices,
	const size_t nOriginal)
{
	const size_t leftOut  =nOriginal%totalSlices;
	const size_t nPerPiece=nOriginal/totalSlices;
	if (leftOut && sId<leftOut) {
		return nPerPiece+1;
	}
	// NOTE - !leftOut || pid>=leftOut
	else {
		return nPerPiece;
	}
}


size_t tArrayMaxSize(
	const size_t totalPieces,
	const size_t nOriginal)
{
	const size_t leftOut  =nOriginal%totalPieces;
	const size_t nPerPiece=nOriginal/totalPieces;
	if (leftOut) {
		return nPerPiece+1;
	}
	else {
		return nPerPiece;
	}
}


// --------------------
void sequentialOddEvenBubbleSort(
	int *anArr,
	const size_t n,
	const __compar_fn_t cmpFunction)
{
	size_t i, j;
	if (n<2) {
		return;
	}
	for (i=0;i<n;i++) {
		// NOTE - iteration x,_,..,.x,_
		if (i%2==0) {
			for (j=0;j<n/2;j++) {
				compare_exchange(anArr,j*2,j*2+1,cmpFunction);
			}
		}
		// NOTE - iteration _,x,...,_,x
		else {
			for (j=1;j<=n/2-1;j++) {
				compare_exchange(anArr,j*2-1,j*2,cmpFunction);
			}
			if(n%2==1) {
				compare_exchange(anArr,n-2,n-1,cmpFunction);
			}
		}
	}
}


void compare_exchange(
	ARRAY_DATA_TYPE *anArr,
	const size_t i1,
	const size_t i2,
	const __compar_fn_t cmpFunction)
{
	const ARRAY_DATA_TYPE e1=anArr[i1];
	const ARRAY_DATA_TYPE e2=anArr[i2];
	if (cmpFunction(&e1,&e2)>0) {
		anArr[i1]=e2;
		anArr[i2]=e1;
	}
	return;
}


void parallelOddEvenBubbleSort(
	ARRAY_DATA_TYPE *anArr,
	const size_t nOriginal,
	const int nThreads,
	const __compar_fn_t cmpFunction)
{
	if (nOriginal < 2*nThreads) {
		printf("Not enought elements to run the parallel bubble sort with %d threads\n"
				"number of threads should be at maximum (%ld), the number of data elements (%ld) integer divided by 2\n.", nThreads, nOriginal/2, nOriginal);
		exit(1);
	}

	/* NOTE - split the data on chunks  |  |  |...|  |  |
	 each chunk has at least 1 element */
	size_t nChunks = 2*nThreads;
	for (size_t i=0;i<nChunks;i++) {
		/* NOTE - iterations x,_,...,x,_
		 processing chunks ((0,1),(2,3),...,(nChunks-2,nChunks-1)) */
		if (i%2==0) {
			/* NOTE - explicitly set the num of threads [unnecessary because the
			 number of threads is rerad fom env] */
			omp_set_num_threads(nThreads);
			// NOTE - using parallel for to give a thread id to each thread
			#pragma omp parallel for schedule(static)
			for (size_t tId=0;tId<nThreads;tId++) {
				size_t chunkStart = tArrayStart(tId*2, nChunks, nOriginal);
				size_t chunkSize = tArraySize(tId*2, nChunks, nOriginal) + tArraySize(tId*2+1, nChunks, nOriginal);
				qsort(&anArr[chunkStart], chunkSize, sizeof(ARRAY_DATA_TYPE), cmpFunction);
			}
		}
		/* NOTE - iterations _,x,...,_,x
		 processing chunks ((1,2),(3,4),...,(nChunks-3,nChunks-2)) */
		else {
			/* NOTE - explicitly set the num of threads [unnecessary because the
			 number of threads is read fom env] */
			omp_set_num_threads(nThreads);
			// NOTE - using parallel for to give a thread id to each thread
			#pragma omp parallel for schedule(static)
			for (size_t tId=0;tId<nThreads+(nChunks%2-1);tId++) {
				size_t chunkStart = tArrayStart(tId*2+1, nChunks, nOriginal);
				size_t chunkSize = tArraySize(tId*2+1, nChunks, nOriginal) + tArraySize(tId*2+2, nChunks, nOriginal);
				qsort(&anArr[chunkStart], chunkSize, sizeof(ARRAY_DATA_TYPE), cmpFunction);
			}
		}
	}
	return;
}


int compareDataElements(const void *a, const void *b)
{
	return * (ARRAY_DATA_TYPE*) a - * (ARRAY_DATA_TYPE*) b;
}


// --------------------
void *newArray(int count, size_t elemSize)
{
	int *aux = malloc(elemSize * count);
	if (aux==NULL) {
		printf("Não foi possível alocar memória");
		exit(1);
	}
	return aux;
}


void printDataArray(ARRAY_DATA_TYPE* arr, size_t n)
{
	for (size_t i=0;i<n;i++) {
		printf("%d,", arr[i]);
	}
	printf("\n");
	return;
}


bool checkSorting(
	const void *arr,
	const size_t nmemb,
	size_t  sizememb,
	const __compar_fn_t cmpFunc)
{
	if (nmemb<2) {
		return true;
	}
	int cmpResult=0,cmpPrevious=0;
	const void *lastAddr=arr+sizememb*(nmemb-1);
	const void *nextAddr;
	for (const void *curAddr=arr;curAddr<lastAddr;curAddr+=sizememb) {
		nextAddr=curAddr+sizememb;
		cmpResult=clamp(cmpFunc(curAddr,nextAddr),-1,1);
		if ((cmpResult==-1 && cmpPrevious== 1)||
			(cmpResult== 1 && cmpPrevious==-1)){
			return false;
		}
		cmpPrevious = cmpResult;
	}
	return true;
}


int clamp(int x, int min, int max)
{
	if (x<min) return min;
	else if (x<max) return x;
	else return max;
}


// --------------------
int main(int argc, char **argv)
{
	int nThreads;
	size_t nOriginal;
	ARRAY_DATA_TYPE *theData;
	ARRAY_DATA_TYPE *sortedData;  // NOTE - only pMain use it
	char *inputFilePath=argv[1];
	char *outputFilePath=argv[2];
	// ----------
	// NOTE - getting from environment OMP_NUM_THREADS
	omp_set_dynamic(0);
	#pragma omp parallel
	{
	#pragma omp single
	nThreads = omp_get_num_threads();
	}
	if (argc<3) {
		puts(
			"Expected 2 arguments:\n"
			"\ta input file path\n"
			"\ta output file path\n"
			"\n"
			"The input file should have:\n"
			"\tfirst line with the count of numbers on the file.\n"
			"\tcount lines, each with a integer.\n"
			"\n"
			"The output file will have:\n"
			"\tfirst line with the count of numbers on the file.\n"
			"\tcount lines, now sorted, each with a integer.");
		return 1;
	}

	// NOTE - load the data
	loadData(inputFilePath, &nOriginal, &theData);

	// NOTE - sort the data
	if (nThreads==1) {
		sequentialOddEvenBubbleSort(theData,nOriginal,compareDataElements);
	}
	else if (nThreads>1) {
		parallelOddEvenBubbleSort(theData, nOriginal, nThreads, compareDataElements);
	}

	// NOTE - check the result
	bool sorted = checkSorting(theData,nOriginal,sizeof(ARRAY_DATA_TYPE),compareDataElements);
	puts(sorted? ">>> Array is now sorted.":">>> Array still not sorted.");

	// NOTE - output the result
	if (sorted) {
		writeArrayOnFile(theData,nOriginal,outputFilePath);
	}

	free(theData);
	// ----------
	return 0;
}


void loadData(
	char *filename,
	size_t *nOriginal,
	ARRAY_DATA_TYPE **theData)
{
	FILE *fp;
	fp = fopen(filename,"r");

	if (fp == NULL) {
		puts("file not found");
		exit(1);
	}

	// NOTE - number of elements in the original data
	// negating the returned value to bypass -Wunused-result
	!fscanf(fp,"%ld\n",nOriginal);
	// NOTE - allocating memory for the data
	*theData = (ARRAY_DATA_TYPE *) newArray(*nOriginal,sizeof(ARRAY_DATA_TYPE));
	// NOTE - load the data
	for (size_t i=0;i<*nOriginal;i++) {
		!fscanf(fp, "%d\n", &(*theData)[i]);
	}
	return;
}


void writeArrayOnFile(
	const ARRAY_DATA_TYPE *arr,
	const size_t n,
	const char* newTextFilePath)
{
	FILE *fp=fopen(newTextFilePath,"w");
	fprintf(fp,"%ld\n",n);
	for (size_t i=0;i<n;i++) {
		fprintf(fp,"%d\n",arr[i]);
	}
	fclose(fp);
	return;
}
// !SECTION
