#include <stdio.h>
#include <stdlib.h>

#define MAX 1000

// Patient structure
typedef struct {
    int age;
    int id;
} Patient;

Patient heap[MAX];
int size = 0;

// Swap two patients
void swap(Patient *a, Patient *b) {
    Patient temp = *a;
    *a = *b;
    *b = temp;
}

// Heapify UP (for insertion)
void heapifyUp(int i) {
    while (i > 0) {
        int parent = (i - 1) / 2;

        if (heap[parent].age < heap[i].age) {
            swap(&heap[parent], &heap[i]);
            i = parent;
        } else {
            break;
        }
    }
}

// Heapify DOWN (for deletion)
void heapifyDown(int i) {
    while (1) {
        int largest = i;
        int left = 2*i + 1;
        int right = 2*i + 2;

        if (left < size && heap[left].age > heap[largest].age)
            largest = left;

        if (right < size && heap[right].age > heap[largest].age)
            largest = right;

        if (largest != i) {
            swap(&heap[i], &heap[largest]);
            i = largest;
        } else {
            break;
        }
    }
}

// Insert new patient
void insertPatient(int id, int age) {
    if (size == MAX) {
        printf("Heap is full!\n");
        return;
    }

    heap[size].id = id;
    heap[size].age = age;
    size++;

    heapifyUp(size - 1); // O(log n)
}

// Extract oldest patient
Patient extractMax() {
    if (size == 0) {
        printf("No patients!\n");
        Patient p = {-1, -1};
        return p;
    }

    Patient root = heap[0];
    heap[0] = heap[size - 1];
    size--;

    heapifyDown(0); // O(log n)

    return root;
}

// Display heap (for debugging)
void displayHeap() {
    printf("Current Patients (Heap):\n");
    for (int i = 0; i < size; i++) {
        printf("(ID:%d Age:%d) ", heap[i].id, heap[i].age);
    }
    printf("\n");
}

// Driver
int main() {
    int choice, id, age;

    while (1) {
        printf("\n1. Add Patient\n2. Call Next Patient\n3. Display Heap\n4. Exit\n");
        printf("Enter choice: ");
        scanf("%d", &choice);

        switch (choice) {
            case 1:
                printf("Enter Patient ID and Age: ");
                scanf("%d %d", &id, &age);
                insertPatient(id, age);
                break;

            case 2: {
                Patient p = extractMax();
                if (p.id != -1)
                    printf("Next patient: ID=%d Age=%d\n", p.id, p.age);
                break;
            }

            case 3:
                displayHeap();
                break;

            case 4:
                return 0;

            default:
                printf("Invalid choice!\n");
        }
    }
}