/*
Given the array:
A = {12, 3, 17, 8, 34, 25, 1}

(a) Construct a min heap from the array using the bottom-up heap construction method (show steps).
(b) Using the min heap, perform heap sort to obtain the array in descending order.
(c) Write a C function for min-heapify and use it to implement descending heap sort.

*/

#include <stdio.h>

// Swap function
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Min-heapify function
void heapify_min(int arr[], int n, int i) {
    int smallest = i;
    int left = 2*i + 1;
    int right = 2*i + 2;

    // Check left child
    if (left < n && arr[left] < arr[smallest])
        smallest = left;

    // Check right child
    if (right < n && arr[right] < arr[smallest])
        smallest = right;

    // If root is not smallest
    if (smallest != i) {
        swap(&arr[i], &arr[smallest]);
        heapify_min(arr, n, smallest);
    }
}

// Function to build min heap
void buildMinHeap(int arr[], int n) {
    for (int i = n/2 - 1; i >= 0; i--) {
        heapify_min(arr, n, i);
    }
}

// Descending Heap Sort using Min Heap
void heapSortDescending(int arr[], int n) {

    // Step 1: Build Min Heap
    buildMinHeap(arr, n);

    // Step 2: Extract elements one by one
    for (int i = n - 1; i > 0; i--) {
        // Move smallest element to end
        swap(&arr[0], &arr[i]);

        // Heapify reduced heap
        heapify_min(arr, i, 0);
    }
}

// Function to print array
void printArray(int arr[], int n) {
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

// Driver code
int main() {
    int n;

    printf("Enter number of elements: ");
    scanf("%d", &n);

    int arr[n];

    printf("Enter elements:\n");
    for (int i = 0; i < n; i++)
        scanf("%d", &arr[i]);

    printf("\nOriginal array:\n");
    printArray(arr, n);

    // Perform descending heap sort
    heapSortDescending(arr, n);

    printf("\nSorted array in DESCENDING order:\n");
    printArray(arr, n);

    return 0;
}