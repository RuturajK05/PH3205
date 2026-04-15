// C program for implementation of Heap Sort using array
#include <stdio.h>


void swap(int *a, int *b)
{
    int t;
    t = *a;
    *a = *b;
    *b = t;
}

void printArray(int arr[], int n)
{
    int i;
    printf("The array:\n\n");
	for (i = 0; i < n; ++i)
		printf("%d ", arr[i]);
	printf("\n");
}

void printArray_shape(int arr[], int n)
{
    int i, j;
    
    printf("\nThe array:\n\n");
	for (i = 0; i < n; ++i)
		{
			if(i == 0 || i == 1 || i == 3)
			{
				for(j = 1; j <= 3-i; j++)
				{
					printf("  ");
				}
			}
			
			printf("%d  ", arr[i]);
			if(i == 0 || i == 2 || i == 6)
			{
				printf("\n");
			}	
			
		}
		printf("\n\n\n\n");
	
}

void heapify(int arr[], int n, int i) /* max heap, ascending heap sort*/
{
	int largest = i; // Initialize largest as root (at arr[0])
	
	printf("\nlargest = %d\n", largest);
	
	int l = 2 * i + 1; // left = 2*i + 1
	int r = 2 * i + 2; // right = 2*i + 2

	// If left child is larger than root
	if (l < n && arr[l] > arr[largest])
		largest = l;

	// If right child is larger than largest so far
	if (r < n && arr[r] > arr[largest])
		largest = r;

	// If largest is not root
	printf("\nlargest = %d\n", largest);
	
	if (largest != i) {
		swap(&arr[i], &arr[largest]);

		printArray_shape(arr, n);
		// Recursively heapify the affected sub-tree
		heapify(arr, n, largest);
	}
}

void heapify_min(int arr[], int n, int i) //has descending heap sort
{
	int smallest = i; // Initialize smallest as root (at arr[0])
	int l = 2 * i + 1; // left = 2*i + 1
	int r = 2 * i + 2; // right = 2*i + 2

	// If left child is smaller than root
	if (l < n && arr[l] < arr[smallest])
		smallest = l;

	// If right child is smaller than largest so far
	if (r < n && arr[r] < arr[smallest])
		smallest = r;

	// If smallest is not root
	printf("smallest = %d\n", smallest);
	if (smallest != i) {
		swap(&arr[i], &arr[smallest]);

		// Recursively heapify the affected sub-tree
		heapify_min(arr, n, smallest);
	}
}
void buildHeap(int arr[], int n)
{
	int i;
	for (i = n / 2 - 1; i >= 0; i--)
		{
			printf("i = %d\n", i);
			heapify(arr, n, i);
		}
		

    printArray_shape(arr, n);
}
void heapSort(int arr[], int n) //ascending
{
    int i;
	
    printf("Build heap:\n\n");
	for (i = n / 2 - 1; i >= 0; i--)
		heapify(arr, n, i);

    printf("The max-heap:\n\n");
    	
    printArray(arr, n);
    
    

	// One by one extract an element from heap
	for (i = n - 1; i > 0; i--) {
		// Move current root to end
		swap(&arr[0], &arr[i]);

		// call heapify on the reduced heap
		heapify(arr, i, 0);
        printArray(arr, n);
	}
}

void heapSortDescending(int arr[], int n)
{
    int i;
	
    printf("Build heap:\n\n");
	for (i = n / 2 - 1; i >= 0; i--)
		heapify_min(arr, n, i);

    printArray(arr, n);

	// One by one extract an element from heap
	for (i = n - 1; i > 0; i--) {
		// Move current root to end
		swap(&arr[0], &arr[i]);

		// call heapify on the reduced heap
		heapify_min(arr, i, 0);
        printArray(arr, n);
	}
}




int main()
{
	int arr[] = { 4, 1, 3, 2, 16, 9, 10, 14, 8, 7 };
	int n = sizeof(arr) / sizeof(arr[0]);
	
	printArray_shape(arr, n);
	
	//heapify(arr, n, 2);
	
	//buildHeap(arr, n);
	
	//printArray_shape(arr, n);
	
	//buildHeap(arr, n);
	
	//heapify_min(arr, n, 4);
	heapSort(arr, n);
	//heapSortDescending(arr, n);

	printf("Sorted array is \n");
	printArray(arr, n);
}
