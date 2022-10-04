/*
Escrito por: Gabriel Borba
Data: 28/09/2022
Ultima atualiazacao: 03/10/2022

Powershell 7
clear && mpicc main_bubble_sort.c -o main_bubble_sort.bin && mpirun -np 4 --oversubscribe main_bubble_sort.bin ./unsorted_numbers.txt ./sorted_numbers.txt
*/


// =============================================================================
// SECTION - INCLUDES
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
// !SECTION


// =============================================================================
// SECTION - GLOBALS
int pRank, pTotal, pMain=0;
MPI_Status status;
// NOTE - tag to send data
#define T_DATA 0

#define MPI_DATA_TYPE MPI_INT
typedef int ARRAY_DATA_TYPE;

// !SECTION


// =============================================================================
// SECTION - DECLARATIONS
/**
 * Return the start position for a slice, among totalSlices, within an array
 * with nOriginal elements.
 * Considering load distribuition with slices at the start being bigger.
*/
size_t pArrayStart(
	const int pid,
	const size_t totalPieces,
	const size_t nOriginal);
/**
 * Return the size of a slice, among totalSlices, within an array with nOriginal
 * elements.
 * Considering load distribuition with slices at the start being bigger.
**/
size_t pArraySize(
	const int pid,
	const size_t totalPieces,
	const size_t nOriginal);
/**
 * Return the maximum number of elements for some slice, among totalSlices,
 * within nOriginal elements.
 */
size_t pArrayMaxSize(
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
	const size_t nData,
	const size_t nOriginal,
	const int pidThis,
	const int nProcesses,
	const __compar_fn_t cmpFunction);
/**
 * Join, sort and split the array of this process with the one from the process
 * at the right.
 * Update the array on this process with the elements at the start of the
 * joined and sorted array.
 */
void compareSplitMin(
	ARRAY_DATA_TYPE *anArr,
	const size_t nArr,
	const int pidRight,
	const size_t nRight,
	const __compar_fn_t cmpFunction);
/**
 * Join, sort and split the array of this process with the one from the process
 * at the left.
 * Update the array on this process with the elements at the end of the
 * joined and sorted array.
 */
void compareSplitMax(
	ARRAY_DATA_TYPE *arr,
	const size_t nArr,
	const int pidLeft,
	const size_t nLeft,
	const __compar_fn_t cmpFunction);
/**
 * Join, sort and split the arrays from pidLeft process with the one from the
 * pidRight process.
 * Update the array on pidThis process with the elements at the start of the
 * joined and sorted array if pidThis is equal to pidLeft or with the elements
 * at the end.
 */
void compareSplit(
	ARRAY_DATA_TYPE *thisArr,
	const size_t nThisArr,
	const int pidThis,
	const int pidLeft,
	const size_t nLeft,
	const int pidRight,
	const size_t nRight,
	ARRAY_DATA_TYPE *commBuffer,
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
 * Preparations needed to start the parallel sorting.
*/
void sortStart(
	const int pidThis,
	const int pidMain,
	const int nProcesses,
	const char *filename,
	size_t *nOriginal,
	ARRAY_DATA_TYPE **theData,
	size_t *nData);
/**
 * Steps needed to finish the parallel sorting.
*/
void sortFinish(
	ARRAY_DATA_TYPE **sortedData,
	const int nProcesses,
	const size_t *nOriginal,
	ARRAY_DATA_TYPE *theData,
	const size_t *nData);
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
size_t pArrayStart(
	const int pid,
	const size_t totalPieces,
	const size_t nOriginal)
{
	const size_t leftOut = nOriginal%totalPieces;
	const size_t maxSize = pArrayMaxSize(totalPieces,nOriginal);
	if (!leftOut) {
		return pid*pArraySize(pid,totalPieces,nOriginal);
	}
	else if (pid<leftOut) {
		return pid*maxSize;
	}
	else {
		return leftOut*maxSize+(pid-leftOut)*(maxSize-1);
	}
}


size_t pArraySize(
	const int pid,
	const size_t totalPieces,
	const size_t nOriginal)
{
	const size_t leftOut  =nOriginal%totalPieces;
	const size_t nPerPiece=nOriginal/totalPieces;
	if (leftOut && pid<leftOut) {
		return nPerPiece+1;
	}
	// NOTE - !leftOut || pid>=leftOut
	else {
		return nPerPiece;
	}
}


size_t pArrayMaxSize(
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
	const size_t nArr,
	const size_t nOriginal,
	const int pidThis,
	const int nProcesses,
	const __compar_fn_t cmpFunction)
{
	int pidAux;
	ARRAY_DATA_TYPE *commBuffer=(ARRAY_DATA_TYPE *) newArray(pArrayMaxSize(nProcesses,nOriginal)*2,sizeof(ARRAY_DATA_TYPE));
	for (size_t i=0;i<nProcesses;i++) {
		if (i%2==0) { // NOTE - iterations x,_,...,x,_
		// processes ((0,1),(2,3),...,(nProcesses-2,nProcesses-1))
			if (pidThis%2==0) {  // NOTE - first process on the pair
			// (X,_)
				// compare-exchange with the other process on the pair
				if (pidThis<nProcesses-1) {
					pidAux=pidThis+1;
					// compareSplitMin(anArr,nArr,pidAux,pArraySize(pidAux,nProcesses,nOriginal),cmpFunction);
					compareSplit(anArr,nArr,pidThis,pidThis,nArr,pidAux,pArraySize(pidAux,nProcesses,nOriginal),commBuffer,cmpFunction);
				}
			}
			else {  // NOTE - second process on the pair
			// (_,Y)
				// compare-exchange with the other process on the paira
				if (pidThis>0) {
					pidAux=pidThis-1;
					// compareSplitMax(anArr,nArr,pidAux,pArraySize(pidAux,nProcesses,nOriginal),cmpFunction);
					compareSplit(anArr,nArr,pidThis,pidAux,pArraySize(pidAux,nProcesses,nOriginal),pidThis,nArr,commBuffer,cmpFunction);
				}
			}
		}
		else {  // NOTE - iterations _,x,...,_,x
		// processes ((1,2),(3,4),...,(nProcesses-3,nProcesses-2))
			if (pidThis%2==1) {  // NOTE - first process on the pair
			// (X,_)
				// compare-exchange with the other process on the pair
				if (pidThis<nProcesses-1) {
					pidAux=pidThis+1;
					// compareSplitMin(anArr,nArr,pidAux,pArraySize(pidAux,nProcesses,nOriginal),cmpFunction);
					compareSplit(anArr,nArr,pidThis,pidThis,nArr,pidAux,pArraySize(pidAux,nProcesses,nOriginal),commBuffer,cmpFunction);
				}
			}
			else {  // NOTE - second process on the pair
			// (_,Y)
				// compare-exchange with the other process on the paira
				if (pidThis>0) {
					pidAux=pidThis-1;
					// compareSplitMax(anArr,nArr,pidAux,pArraySize(pidAux,nProcesses,nOriginal),cmpFunction);
					compareSplit(anArr,nArr,pidThis,pidAux,pArraySize(pidAux,nProcesses,nOriginal),pidThis,nArr,commBuffer,cmpFunction);
				}
			}
		}
	}
	free(commBuffer);
}


void compareSplitMin(
	ARRAY_DATA_TYPE *anArr,
	const size_t nArr,
	const int pidRight,
	const size_t nRight,
	const __compar_fn_t cmpFunction)
{
	const int aUnionSize = nArr+nRight;
	ARRAY_DATA_TYPE *aUnionArr = (ARRAY_DATA_TYPE *) newArray(aUnionSize, sizeof(ARRAY_DATA_TYPE));

	// printf("(pid%d) data at the start: [%d,...]\n", pRank, anArr[0]);  // REVIEW

	// NOTE - send this process array to the other process on the pair
	MPI_Send(anArr, nArr, MPI_DATA_TYPE, pidRight, T_DATA, MPI_COMM_WORLD);
	// printf("pid%d sent to pid%d: [%d,...]\n", pRank, pidRight, anArr[0]);  // REVIEW
	// NOTE - receive the array from the other process on the pair
	MPI_Recv(aUnionArr, nRight, MPI_DATA_TYPE, pidRight, T_DATA, MPI_COMM_WORLD, &status);
	// printf("pid%d received from pid%d: [%d,...]\n", pRank, pidRight, aUnionArr[0]);  // REVIEW

	// NOTE - join the arrays and sort the generated array
	for ( int i=0;i<nArr;i++ ) {
		aUnionArr[nRight+i]=anArr[i];
	}

	qsort(aUnionArr, aUnionSize, sizeof(ARRAY_DATA_TYPE), cmpFunction);

	// NOTE - getting the initial part of the joined array
	for ( int i=0;i<nArr;i++ ) {
		anArr[i]=aUnionArr[i];
	}

	free(aUnionArr);
	// printf("(pid%d) data at the end: [%d,...]\n", pRank, anArr[0]);  // REVIEW
	return;
}


void compareSplitMax(
	ARRAY_DATA_TYPE *anArr,
	const size_t nArr,
	const int pidLeft,
	const size_t nLeft,
	const __compar_fn_t cmpFunction)
{
	const int aUnionSize = nArr+nLeft;
	ARRAY_DATA_TYPE *aUnionArr = (ARRAY_DATA_TYPE *) newArray(aUnionSize, sizeof(ARRAY_DATA_TYPE));

	// printf("(pid%d) data at the start: [%d,...]\n", pRank, anArr[0]);  // REVIEW

	// NOTE - send this process array to the other process on the pair
	MPI_Send(anArr, nArr, MPI_DATA_TYPE, pidLeft, T_DATA, MPI_COMM_WORLD);
	// printf("pid%d sent to pid%d: [%d,...]\n", pRank, pidLeft, anArr[0]);  // REVIEW
	// NOTE - receive the array from the other process on the pair
	MPI_Recv(aUnionArr, nLeft, MPI_DATA_TYPE, pidLeft, T_DATA, MPI_COMM_WORLD, &status);
	// printf("pid%d received from pid%d: [%d,...]\n", pRank, pidLeft, aUnionArr[0]);  // REVIEW

	// NOTE - join the arrays and sort the generated array
	for ( int i=0;i<nArr;i++ ) {
		aUnionArr[nLeft+i]=anArr[i];
	}

	qsort(aUnionArr, aUnionSize, sizeof(ARRAY_DATA_TYPE), cmpFunction);

	// NOTE - getting the last part of the joined array
	for ( int i=0;i<nArr;i++ ) {
		anArr[i]=aUnionArr[nLeft+i];
	}

	free(aUnionArr);
	// printf("(pid%d) data at the end: [%d,...]\n", pRank, anArr[0]);  // REVIEW
	return;
}


void compareSplit(
	ARRAY_DATA_TYPE *thisArr,
	const size_t nThisArr,
	const int pidThis,
	const int pidLeft,
	const size_t nLeft,
	const int pidRight,
	const size_t nRight,
	ARRAY_DATA_TYPE *commBuffer,
	const __compar_fn_t cmpFunction)
{
	int pidTarget;
	size_t nTarget;
	if (pidThis==pidLeft) {
		pidTarget = pidRight;
		nTarget = nRight;
	}
	else {
		pidTarget=pidLeft;
		nTarget=nLeft;
	}
	const int aUnionSize = nLeft+nRight;
	// ARRAY_DATA_TYPE *aUnionArr = (ARRAY_DATA_TYPE *) newArray(aUnionSize, sizeof(ARRAY_DATA_TYPE));
	ARRAY_DATA_TYPE *aUnionArr = commBuffer;

// --------------------
	// FIXME - large data are not sent
	// LINK - https://stackoverflow.com/a/37207677
	// NOTE - incorrect way of sending and receiving (accordingly to the link)
	// // NOTE - send this process array to the other process on the pair
	// MPI_Send(thisArr, (int) nThisArr, MPI_DATA_TYPE, pidTarget, T_DATA, MPI_COMM_WORLD);
	// printf("sending %ld elements from pid%d to pid%d\n", nThisArr, pidThis, pidTarget);  // REVIEW

	// // NOTE - receive the array from the other process on the pair
	// MPI_Recv(aUnionArr, (int) nTarget, MPI_DATA_TYPE, pidTarget, T_DATA, MPI_COMM_WORLD, &status);
	// printf("receiving %ld elements from pid%d to pid%d\n", nTarget, pidTarget, pidThis);  // REVIEW

	// NOTE - one of the correct way of sending and receiving (accordingly to the link)
	// NOTE - send and receive the processes data
	MPI_Sendrecv(
		thisArr,nThisArr,MPI_DATA_TYPE,pidTarget,T_DATA,
		aUnionArr,aUnionSize,MPI_DATA_TYPE,pidTarget,T_DATA,
		MPI_COMM_WORLD,&status);
// --------------------

	// NOTE - join the arrays and sort the generated array
	for ( int i=0;i<nThisArr;i++ ) {
		aUnionArr[nTarget+i]=thisArr[i];
	}

	qsort(aUnionArr, aUnionSize, sizeof(ARRAY_DATA_TYPE), cmpFunction);

	// NOTE - getting the initial part of the joined array if pidThis == pidLeft
	if (pidThis==pidLeft) {
		for ( int i=0;i<nThisArr;i++ ) {
			thisArr[i]=aUnionArr[i];
		}
	}
	// NOTE - getting the last part of the joined array if pidThis == pidRight
	else {
		for ( int i=0;i<nThisArr;i++ ) {
			thisArr[i]=aUnionArr[nLeft+i];
		}
	}

	// free(aUnionArr);
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
	size_t nOriginal;
	ARRAY_DATA_TYPE *theData;
	ARRAY_DATA_TYPE *sortedData;  // NOTE - only pMain use it
	size_t nData;
	// ----------
	if (argc<3) {
		if (pRank==pMain) {
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
		}
	}
	else {
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &pTotal);
		MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

		char *inputFilePath=argv[1];
		char *outputFilePath=argv[2];

		// NOTE - spliting the data among the processes
		sortStart(pRank,pMain,pTotal,inputFilePath,&nOriginal,&theData,&nData);

		if (pTotal==1) {
			if (pRank==pMain) {
				sequentialOddEvenBubbleSort(theData,nOriginal,compareDataElements);
			}
		}
		else if(pTotal>1) {
			parallelOddEvenBubbleSort(theData,nData,nOriginal,pRank,pTotal,compareDataElements);
		}

		// NOTE - gathering the data from the processes
		sortFinish(&sortedData,pTotal,&nOriginal,theData,&nData);

		MPI_Finalize();

		free(theData);

		// ----------
		// NOTE - output results
		if (!pRank) {
			bool sorted = checkSorting(sortedData,nOriginal,sizeof(ARRAY_DATA_TYPE),compareDataElements);
			puts(sorted? ">>> Array is now sorted.":">>> Array still not sorted.");
			writeArrayOnFile(sortedData,nOriginal,outputFilePath);
		}

		free(sortedData);
	}

	// ----------
	return 0;
}


void sortStart(
	const int pidThis,
	const int pidMain,
	const int nProcesses,
	const char *filename,
	size_t *nOriginal,
	ARRAY_DATA_TYPE **theData,
	size_t *nData)
{
	FILE *fp;  // NOTE - only pMain use it
	// NOTE - send the number of elements in the original data
	if (pRank==pMain) {
		fp = fopen(filename,"r");
		fscanf(fp,"%ld\n",nOriginal);
		for (size_t pid=0;pid<nProcesses;pid++) {
			if (pid!=pMain) {
				MPI_Send(nOriginal,1,MPI_UNSIGNED_LONG,pid,T_DATA,MPI_COMM_WORLD);
			}
		}
	}
	else {
		MPI_Recv(nOriginal,1,MPI_UNSIGNED_LONG,pMain,T_DATA,MPI_COMM_WORLD,&status);
	}
	*nData = pArraySize(pidThis,nProcesses,*nOriginal);
	*theData = (ARRAY_DATA_TYPE *) newArray(*nData,sizeof(ARRAY_DATA_TYPE));
	// NOTE - load and send the data array
	if (pRank==pMain) {
		// NOTE - load data for the main process
		for (size_t i=0;i<*nData;i++) {
			fscanf(fp, "%d\n", &(*theData)[i]);
		}
		// NOTE - send data to the other processes
		ARRAY_DATA_TYPE *aBuffer;
		aBuffer = (ARRAY_DATA_TYPE *) newArray(pArrayMaxSize(nProcesses,*nOriginal),sizeof(ARRAY_DATA_TYPE));
		size_t nBuffer;  // NOTE - number of elements currently in buffer
		for (size_t pid=0;pid<nProcesses;pid++) {
			if (pid!=pMain) {
				nBuffer = pArraySize(pid,nProcesses,*nOriginal);
				for (size_t i=0;i<*nData;i++) {
					fscanf(fp, "%d\n", &aBuffer[i]);
				}
				MPI_Send(aBuffer,nBuffer,MPI_DATA_TYPE,pid,T_DATA,MPI_COMM_WORLD);
			}
		}
		free(aBuffer);
		fclose(fp);
	}
	// NOTE - receive data from the main process
	else {
		MPI_Recv(*theData,*nData,MPI_DATA_TYPE,pMain,T_DATA,MPI_COMM_WORLD,&status);
	}
	return;
}


void sortFinish(
	ARRAY_DATA_TYPE **sortedData,
	const int nProcesses,
	const size_t *nOriginal,
	ARRAY_DATA_TYPE *theData,
	const size_t *nData)
{
	if (pRank!=pMain) {
	MPI_Send(
		theData,
		pArraySize(pRank,nProcesses,*nOriginal),
		MPI_INT, pMain, T_DATA,
		MPI_COMM_WORLD);
	}
	else {
		*sortedData = (ARRAY_DATA_TYPE *) newArray(*nOriginal,sizeof(ARRAY_DATA_TYPE));
		size_t iDataStartPid;
		ARRAY_DATA_TYPE *aBuffer;
		aBuffer = (ARRAY_DATA_TYPE *) newArray(pArrayMaxSize(nProcesses,*nOriginal),sizeof(ARRAY_DATA_TYPE));
		size_t nBuffer;  // NOTE - number of elements currently in the buffer

		for (size_t pid = 0;pid<nProcesses;pid++) {
			iDataStartPid=pArrayStart(pid,nProcesses,*nOriginal);
			if (pid!=pRank) {
				nBuffer = pArraySize(pid,nProcesses,*nOriginal);
				MPI_Recv(aBuffer,nBuffer,MPI_DATA_TYPE,pid,T_DATA,MPI_COMM_WORLD,&status);
				for (size_t i=0;i<nBuffer;i++) {
					(*sortedData)[iDataStartPid+i]=aBuffer[i];
				}
			}
			else {
				for (size_t i=0;i<(*nData);i++) {
					(*sortedData)[iDataStartPid+i]=theData[i];
				}
			}
		}
		free(aBuffer);
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
