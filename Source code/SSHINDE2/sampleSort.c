// SampleSort.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "stdlib.h"
#include "math.h"
#include "stdio.h"
#include "stdbool.h"
#include "string.h"
#include "time.h"
#include <sys/time.h>
#include "mpi.h"

// MPI Data structures
int    size;
int    my_rank;
MPI_Status  status;
MPI_Comm    comm;
//	

int * ReceiveInput(int *,int *);
void distribute(int *,int,int,int);
void sort(int *,int);
int * localSplitter(int *, int, int, bool);
void allToAllCalculation(int *, int *, int *,int *,int);
void countElements(int *, int , int *, int);
void merge(int *,int *, int, int, int *);
void prepareGatherv(int *, int *, int);

int main(int argc, char* argv[]) 
{
//
	//MPI_Comm    comm;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);

	int p,t,n; // processorts_count, thread_count, element_count
	int *input;
	int * inputList = (int *) malloc(n * sizeof(int));
	
	int * localList;
	// Time calculation - average of wtime() execution
	double startTime,endTime;
	double wtime_overhead = 0.0;
    for (int i = 0; i < 200; i++) {
        startTime = MPI_Wtime();
        endTime = MPI_Wtime();
        wtime_overhead = wtime_overhead + (startTime - endTime);
    }
    wtime_overhead = wtime_overhead/200.0;
	//
	//printf("i m started\n");
	
	if(my_rank == 0)
	{
		input = ReceiveInput(&t,&n);
		inputList = input;
		p = size;		//uncomment then
		printf("Value of p n t %d %d %d \n",p,n,t);
		
		//Start timer..
		startTime = MPI_Wtime();
		
		//Here we received all inputList array, processors and threadCount.
		/*for(int i=0;i<n;i++)
		{
			printf("%d ",*(inputList+i));
		}*/
		for(int i=1;i<p;i++)
		{
			MPI_Send(&n,1,MPI_INT,i,2,comm);
		}
		localList = (int *) malloc((n/p) * sizeof(int));
		for(int i=0;i<n/p;i++)
		{
			*(localList + i) = *(inputList + i);
		}
		
		distribute(inputList,p,t,n);
	}
	else
	{
		//receive n
		MPI_Recv(&n, 1, MPI_INT, 0, 2, comm,&status);
		p=size;
		
		localList = (int *) malloc((n/p) * sizeof(int));
		//recieve localList for all slave nodes.
		MPI_Recv(localList, (n/p), MPI_INT, 0, 1, comm,&status);
	}
	//Display received local copies..
	sort(localList,n/p);
	//MPI_Barrier(comm);
	//printf("\nlocalcopy sorted at %d \n",my_rank);
	//for(int i=0;i<(n/p);i++)
	//{
	//		printf("rank %d %d ",my_rank,*(localList+i));
	//}
	//printf("\n");
	//Find spliters at each node.
	bool isPsplitter = true;;
	int * splitterArray = (int *) malloc(p*sizeof(int));
	splitterArray = localSplitter(localList,(n/p),p,isPsplitter);
	//MPI_Barrier(comm);
	//printf("\nSplitters for %d are : \n",my_rank);
	//for(int i=0;i<p;i++)
	//{
	//		printf("rank %d %d ",my_rank,*(splitterArray+i));
	//}
	//All nodes should send splitterArray to root 0 for merging.
	int * sampleArray = (int *) malloc(p*p*sizeof(int));
MPI_Barrier(comm);	
	MPI_Allgather(splitterArray,p,MPI_INT,sampleArray,p,MPI_INT,comm);
MPI_Barrier(comm);
	//printf("\nSampleArray at %d :\n",my_rank);
	//for(int i=0;i<(p*p);i++)
	//{
	//	printf("%d ",*(sampleArray+i));
	//}
	//printf("\n");
	//Sort sample array.
	sort(sampleArray,(p*p));
	//Find final (p-1) splitter array for entire range of numbers.
	isPsplitter = false;
	int * FinalsplitterArray = (int *) malloc((p-1)*sizeof(int));
	FinalsplitterArray = localSplitter(sampleArray,(p*p),p,isPsplitter);
	//MPI_Barrier(comm);
	//printf("FinalsplitterArray at %d :\n",my_rank);
	//for(int i=0;i<(p-1);i++)
	//{
	//	printf("%d ",*(FinalsplitterArray+i));
	//}
	//printf("\n");
	//
	//**********************************************************
	//This will prepare scounts, sdispls, rcounts, rdispls for each process..
	int * scounts = (int *) malloc(p*sizeof(int));
	int * sdispls = (int *) malloc(p*sizeof(int));
	//Initialization
	for(int i=0;i<p;i++)
	{
		*(scounts + i) = 0;
		*(sdispls + i) = 0;
	}
	
	allToAllCalculation(localList,FinalsplitterArray,scounts,sdispls,n/p);
	//creating rcounts & rdispls with max n length
	int * rcounts = (int *)malloc(p*sizeof(int));
	int * rdispls = (int *)malloc(p*sizeof(int));
	//Initialization
	for(int i=0;i<p;i++)
	{
		*(rcounts + i) = 0;
		*(rdispls + i) = 0;
	}
	//This will move scount from all to all processed and each process will 
	// get rcounts updated.
MPI_Barrier(comm);	
	MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
MPI_Barrier(comm);
	//printf("\n");
	//for(int i=0;i<p;i++)
	//{
	//	printf("rcounts rank %d :  %d",my_rank,rcounts[i]);
	//}
	//printf("\n");
	//Generate rdispl list.
	rdispls[0] = 0;
	for(int i=1;i<p;i++)
	{
		//if(rcounts[i] != 0)
		{
			rdispls[i] = rdispls[i-1] + rcounts[i-1];
		}	
	}
	//for(int i=0;i<p;i++)
	//{
	//	printf("rdisp rank %d : %d ",my_rank,rdispls[i]);
	//}
	//printf("\n");
	//***********************************************************************
	//This will do MPI_Alltoallv communication with vectors
	int bucketSize = 0;
	for(int i=0;i<p;i++)
	{
		bucketSize += rcounts[i];
	}
	//printf("bucketSize rank %d : %d \n",my_rank,bucketSize);
	int * localBucket = (int *)malloc(bucketSize*sizeof(int));
MPI_Barrier(comm);	
	MPI_Alltoallv(localList,scounts,sdispls,MPI_INT,localBucket,rcounts,rdispls,MPI_INT,comm);
MPI_Barrier(comm);
	//printf("\n");
	//for(int i=0;i<bucketSize;i++)
	//{
	//	printf("bucket @ rank %d :  %d , ",my_rank,localBucket[i]);
	//}
	//printf("\n");
	//Sort bucket on every node. (Merge algorithm is used)
	int * sortedBucket = (int *) malloc(bucketSize * sizeof(int));
	merge(localBucket,rdispls,bucketSize,p,sortedBucket);
	//for(int i=0;i<bucketSize;i++)
	//{
	//	printf("sortedBucket %d: %d ",my_rank,sortedBucket[i]);
	//}
	//printf("\n");
	//***********************************************************************
	//Final result formation
	//send sortedBucket to root node 0.
	int * finalRcvCounts = (int *) malloc(p * sizeof(int));
	int * finalRcvDispls = (int *) malloc(p * sizeof(int));
	int * finalResults   = (int *) malloc(n * sizeof(int));
	//Gather counts of every process bucket.
MPI_Barrier(comm);	
	MPI_Gather(&bucketSize,1,MPI_INT,finalRcvCounts,1,MPI_INT,0,comm);
MPI_Barrier(comm);	
	
	if(my_rank == 0)
	{
		/*for(int i=0;i<p;i++)
		{
			printf("finalRcvCounts %d ",finalRcvCounts[i]);
		}
		printf("\n");
		finalRcvDispls[0] = 0;
		for(int i=1;i<p;i++)
		{
			if(finalRcvDispls[i] != 0)
			{
				finalRcvDispls[i] = finalRcvDispls[i-1] + finalRcvCounts[i-1];
			}	
		}
		for(int i=0;i<p;i++)
		{
			printf("finalRcvDispls %d ",finalRcvDispls[i]);
		}
		printf("\n");*/
		prepareGatherv(finalRcvCounts,finalRcvDispls,p);
	}
	//receive all to one gather vectors at root node 0.
MPI_Barrier(comm);	
	MPI_Gatherv(sortedBucket,bucketSize,MPI_INT,finalResults,finalRcvCounts,finalRcvDispls,MPI_INT,0,comm);
MPI_Barrier(comm);	
	
	if(my_rank == 0)
	{
		//Stop timer.
		endTime = MPI_Wtime();
		
		printf("\nTHE SORTED OUTPUT IS - \n");
		for(int i=0;i<n;i++)
		{
			printf(" %d ",finalResults[i]);
		}
		printf("\n");
		printf("Approx. Time for execution : %f ",(endTime - startTime - wtime_overhead));
	}
	
	MPI_Finalize();
	return 0;
}
int * ReceiveInput(int * threadCount,int *count)
{
	//Read input file..
	int i=0;int inputCount;int temp=0;
	FILE* f = fopen("numbers.txt","r");
	fscanf (f,"%d",count);
	fscanf (f,"%d",threadCount);
	int *receiveInput;
	receiveInput = (int *) malloc(*count * sizeof(int));

	while(i!= *count)
	{
		fscanf(f,"%d", &temp);
		*(receiveInput+i)=temp;
        i++;
	}
	//All values are set here,,,
	 
	return receiveInput;
}
void distribute(int *inputList,int p,int t,int n)
{
	int start=(n/p),finish=2*(n/p);
	for(int i=1;i<p;i++)
	{
		int *sendArray = (int *) malloc(n/p * sizeof(int));
		//printf("For processor id %d \n",i);
		for(int j=start,i=0;j<finish+1;j++,i++)
		{
			*(sendArray + i) = *(inputList + j);
		}
		start = finish ;
		finish = finish + n/p;

		//for(int k=0;k<n/p;k++)
		//{
		//	printf("%d ",*(sendArray+k));
		//}
		//printf("\n");
		// MPI_SEND 
		// send sendArray to processor i.
		MPI_Send(sendArray,(n/p),MPI_INT,i,1,comm);
		free(sendArray);
	}
	
	
	/*int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int *receiveBuff = (int *) malloc(n/p * sizeof(int));
	printf("scatter starts %d\n",my_rank);
	MPI_Scatter(inputList, 8, MPI_INT, receiveBuff, 8, MPI_INT, 0, comm);*/
}

//
void sort(int * localList,int n)
{
	//Insertion sort
	int key=0,index=0;
	for(int i=0;i<n;i++)
	{
		key = *(localList+i);
		index = i;
		for(int j=i+1;j<(n);j++)
		{
			if (key > *(localList + j))
			{
				index = j;
				key = *(localList + j);
			}
		}
		int temp = *(localList + i);
		*(localList + i) = *(localList + index);
		*(localList + index) = temp;
	}
}
//
int * localSplitter(int * localList, int n, int p,bool isPsplitter)
{
	int * splitArray = (int *) malloc(p*sizeof(int));
	int s,considerCount ;
	/*if(isPsplitter)
	{
		printf("\nisPsplitter true\n");
		float temp = (float) (n - p)/p;
		s = ceil(temp);
		considerCount = n - p;
	}
	else
	{
		printf("\nisPsplitter false %d %d %d %f\n",n,p,(n - p + 1),(n - p + 1)/p);
		float temp = (float)(n - p + 1)/p;
		printf("\n temp %f",temp);
		s = ceil(temp);
		considerCount = n - p + 1;
	}
	printf("\nsplitter will be multiples of %d \n",s);
	for(int i=s-1,j=0;i<n;i=i+s,j++)
	{
		//printf("splitter %d",*(localList+i));
		*(splitArray+j) = *(localList+i);
	}*/
	
	if(isPsplitter)
	{
		float temp = (float) n/p;
		s = floor(temp);
		for(int i=1,j=0;i<p+1;i++,j++)
		{
			*(splitArray+j) = *(localList+(i*s)-1);
		}
	}
	else
	{
		float temp = (float) n/(p-1);
		s = floor(temp);
		for(int i=1,j=0;i<p;i++,j++)
		{
			*(splitArray+j) = *(localList+(i*s)-1);
		}
	}
	return splitArray;
}

//
void allToAllCalculation(int * localList,int * FinalsplitterArray,int * scounts,int * sdispls,int n)
{
	int p,my_rank;
	int min;
	int max;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int maxCount,minCount;
	
	// Calculation of scounts
	//for(int i=0;i<p;i++)
	//{
	//	printf("before for rank %d : %d",my_rank,scounts[i]);
	//}
	for(int j=0;j<p;j++)	//Accessing FinalsplitterArray array
	{
		if(j==0)
		{
			min= -2147483648;
			max = *(FinalsplitterArray+0);
		}
		else if(j==p-1)
			{
				min = *(FinalsplitterArray+(p-2));
				max = 2147483647 ;
			}
			else
			{
				min = *(FinalsplitterArray+(j-1));
				max = *(FinalsplitterArray+(j));
			}
		//printf("\nRank %d : min max %d %d\n",my_rank,min,max);
		/*for(int i=0;i<n;i++)
		{
			if ((*(localList+i) > min) && (*(localList+i) <= max))
			{
				//*(scounts + j) = *(scounts + j)+1; 
				scounts[j]++;
				printf("%d , %d\n",*(localList+i),scounts[j]);
			}
		}*/
		countElements(localList,max,&maxCount,n);
		countElements(localList,min,&minCount,n) ;
		//printf("rank %d - Answer count is %d %d final count :%d\n",my_rank,minCount,maxCount,(maxCount-minCount));
		scounts[j] = maxCount - minCount;
	}
	//for(int i=0;i<p;i++)
	//{
	//	printf("for rank %d : %d ",my_rank,scounts[i]);
	//}
	//printf("\n");
	//Calculation of sdispls
	sdispls[0] = 0;
	for(int i=1;i<p;i++)
	{
		//if(scounts[i] != 0)
		{
			sdispls[i] = sdispls[i-1] + scounts[i-1];
		}	
	}
	//for(int i=0;i<p;i++)
	//{
	//	printf("disp rank %d : %d ",my_rank,sdispls[i]);
	//}
	//printf("\n");
	
}
//
void countElements(int *localList, int target, int *count, int n)
{
	int lo=0;
	//int n = 9;
	int hi= n -1;
	int mid;
	bool flag = false;

	while(lo<=hi)
	{
		float temp = (float)(lo+hi)/2;
		mid = floor(temp);
		if(localList[mid] > target)
		{
			if((mid-1) < 0)
			{
				*count = 0;
				flag = true;
				break;
			}
			else
			{
				if(localList[mid-1] <= target)
				{
					*count = mid ;
					flag = true;
					break;
				}
			}
		}
		if((mid+1) == n)
		{
			*count = n;
			flag = true;
			break;
		}
		//
		if(localList[mid] <= target)
		{
			lo = mid + 1;
		}
		else
			hi = mid - 1;
	}
	
	if(!flag)
	{
		*count = mid;
	}
}
//
void merge(int * localList, int * rdispls, int n,int p, int * sortedBucket)
{
		int min, minIndex=0;
		
		int sentinel = 2147483647; // Max int value.
			
		for (int k=0 ;k<n;k++) //localList
		{
			//min is initialized to max int value.
			min = sentinel;
			
			//printf("min to be compared %d\n",min);
			//Each processors top elemts is compared with min.
			//top elemts is identifed though rdispls array.
			for(int i=1;i<p+1;i++) //rdispls
			{
				if( rdispls[i-1] < n )
				{
					if(min > localList[rdispls[i-1]])
					{
						min = localList[rdispls[i-1]];
						minIndex = i-1 ;
					}
				}
			}
			
			
			//printf("Element : %d \n",min);
			*(sortedBucket + k) = min;
			//set min number found to sentilnel to avoid reconsidering it again..
			localList[rdispls[minIndex]] = sentinel; //sentinel
			rdispls[minIndex] += 1;
			/*for(int i=0;i<p;i++)
			{
				printf("%d ",rdispls[i]);
			}
			printf("\n");*/
		}

}
//
void prepareGatherv(int * finalRcvCounts,int * finalRcvDispls,int p)
{
		//for(int i=0;i<p;i++)
		//{
		//	printf("finalRcvCounts %d ",finalRcvCounts[i]);
		//}
		//printf("\n");
		finalRcvDispls[0] = 0;
		for(int i=1;i<p;i++)
		{
			//if(finalRcvCounts[i] != 0)
			{
				finalRcvDispls[i] = finalRcvDispls[i-1] + finalRcvCounts[i-1];
			}	
		}
		//for(int i=0;i<p;i++)
		//{
		//	printf("finalRcvDispls %d ",finalRcvDispls[i]);
		//}
		//printf("\n");
}