/*
Given a hash table of size 7 and hash function:
h(k)=kmod7

Insert the keys:
50, 700, 76, 85, 92, 73, 101
Resolve collisions using:

1. Separate Chaining
2. Linear Probing
3. Quadratic Probing
4. Double Hashing

Write C code for each method and display the final table.
*/

/*
//SEPARATE CHAINING
#include <stdio.h>
#include <stdlib.h>

#define SIZE 7

typedef struct Node {
    int key;
    struct Node* next;
} Node;

Node* table[SIZE] = {NULL};

void insert(int key) {
    int index = key % SIZE;

    Node* newNode = malloc(sizeof(Node));
    newNode->key = key;
    newNode->next = table[index];
    table[index] = newNode;
}

void display() {
    printf("Chaining:\n");
    for(int i=0;i<SIZE;i++) {
        printf("%d: ", i);
        Node* temp = table[i];
        while(temp) {
            printf("%d -> ", temp->key);
            temp = temp->next;
        }
        printf("NULL\n");
    }
}


// Time Complexity:
// Insertion: O(1) avg, O(n) worst


int main() {
    int keys[] = {50,700,76,85,92,73,101};

    for(int i=0;i<7;i++) insert(keys[i]);

    display();
}
*/
 /*
 //LINEAR PROBING

#include <stdio.h>

#define SIZE 7

int table[SIZE];

void insert(int key) {
    int index = key % SIZE;

    while(table[index] != -1) {
        index = (index + 1) % SIZE;
    }
    table[index] = key;
}

void display() {
    printf("Linear Probing:\n");
    for(int i=0;i<SIZE;i++)
        printf("%d ", table[i]);
}


//Time Complexity:
//Insertion: O(1) avg, O(n) worst
//(Primary clustering)


int main() {
    int keys[] = {50,700,76,85,92,73,101};

    for(int i=0;i<SIZE;i++) table[i] = -1;

    for(int i=0;i<7;i++) insert(keys[i]);

    display();
}
 */

/*
//QUADRATIC PROBING
#include <stdio.h>

#define SIZE 7

int table[SIZE];

void insert(int key) {
    int index = key % SIZE;
    int i = 0;

    while(table[(index + i*i) % SIZE] != -1) {
        i++;
    }
    table[(index + i*i) % SIZE] = key;
}

void display() {
    printf("Quadratic Probing:\n");
    for(int i=0;i<SIZE;i++)
        printf("%d ", table[i]);
}


//Time Complexity:
//Insertion: O(1) avg, O(n) worst
//Reduces primary clustering


int main() {
    int keys[] = {50,700,76,85,92,73,101};

    for(int i=0;i<SIZE;i++) table[i] = -1;

    for(int i=0;i<7;i++) insert(keys[i]);

    display();
}

*/

/*
//DOUBLE HASHING
#include <stdio.h>

#define SIZE 7

int table[SIZE];

int hash2(int key) {
    return 1 + (key % (SIZE - 1));
}

void insert(int key) {
    int index = key % SIZE;
    int step = hash2(key);
    int i = 0;

    while(table[(index + i*step) % SIZE] != -1) {
        i++;
    }
    table[(index + i*step) % SIZE] = key;
}

void display() {
    printf("Double Hashing:\n");
    for(int i=0;i<SIZE;i++)
        printf("%d ", table[i]);
}


//Time Complexity:
//Insertion: O(1) avg, O(n) worst
//Best among open addressing methods


int main() {
    int keys[] = {50,700,76,85,92,73,101};

    for(int i=0;i<SIZE;i++) table[i] = -1;

    for(int i=0;i<7;i++) insert(keys[i]);

    display();
}
*/


