#include <stdio.h>

// Swap utility
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Max-heapify (downward)
void maxHeapify(int arr[], int n, int i) {
    int largest = i;
    int left = 2*i + 1;
    int right = 2*i + 2;

    // Compare with left child
    if (left < n && arr[left] > arr[largest])
        largest = left;

    // Compare with right child
    if (right < n && arr[right] > arr[largest])
        largest = right;

    // If root is not largest, swap and continue
    if (largest != i) {
        swap(&arr[i], &arr[largest]);
        maxHeapify(arr, n, largest);
    }
}

// Build max heap in O(n)
void buildMaxHeap(int arr[], int n) {
    for (int i = n/2 - 1; i >= 0; i--) {
        maxHeapify(arr, n, i);
    }
}

// Function to find k-th largest
int kthLargest(int arr[], int n, int k) {
    buildMaxHeap(arr, n);   // O(n)

    int heapSize = n;

    // Extract max (k-1) times
    for (int i = 0; i < k-1; i++) {
        swap(&arr[0], &arr[heapSize - 1]); // move max to end
        heapSize--;                        // reduce heap size
        maxHeapify(arr, heapSize, 0);      // O(log n)
    }

    return arr[0]; // k-th largest
}

// Driver code
int main() {
    int n, k;

    printf("Enter number of elements: ");
    scanf("%d", &n);

    int arr[n];
    printf("Enter elements:\n");
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }

    printf("Enter k: ");
    scanf("%d", &k);

    int result = kthLargest(arr, n, k);
    printf("%d-th largest element = %d\n", k, result);

    return 0;
}