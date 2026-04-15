/* 📄 CS3201 – CLASS TEST (25 MARKS)
🔴 SECTION A — Red-Black Trees (12.5 marks)
Q1 (7.5 marks)

Insert the following keys into an empty red-black tree in order:

20, 15, 25, 10, 5
Show tree after each insertion
Clearly indicate rotations and recoloring
State which case (1/2/3) is applied
Q2 (5 marks)

Delete the node 15 from the red-black tree obtained above.

Show all steps
Explain how violations are fixed
🪵 SECTION B — Heaps (7.5 marks)
Q3 (5 marks)

Given:

A = {5, 77, 2, 3, 44, 22, 1, 65}

(a) Build a max heap using bottom-up method (show steps)
(b) Draw the final heap as a tree

Q4 (2.5 marks)

Write an algorithm to find the k-th largest element using heaps.

State time complexity
Explain why it is efficient */
/*
#include <stdio.h>
#include <stdlib.h>

typedef struct Node {
    int key;
    char color; // 'R' or 'B'
    struct Node *left, *right, *parent;
} Node;

Node *root = NULL;

// Create node
Node* newNode(int key) {
    Node* n = malloc(sizeof(Node));
    n->key = key;
    n->color = 'R';
    n->left = n->right = n->parent = NULL;
    return n;
}

// Left rotate
void leftRotate(Node *x) {
    Node *y = x->right;
    x->right = y->left;

    if (y->left)
        y->left->parent = x;

    y->parent = x->parent;

    if (!x->parent)
        root = y;
    else if (x == x->parent->left)
        x->parent->left = y;
    else
        x->parent->right = y;

    y->left = x;
    x->parent = y;
}

// Right rotate
void rightRotate(Node *y) {
    Node *x = y->left;
    y->left = x->right;

    if (x->right)
        x->right->parent = y;

    x->parent = y->parent;

    if (!y->parent)
        root = x;
    else if (y == y->parent->left)
        y->parent->left = x;
    else
        y->parent->right = x;

    x->right = y;
    y->parent = x;
}

// Fix insert
void fixInsert(Node *z) {
    while (z->parent && z->parent->color == 'R') {

        if (z->parent == z->parent->parent->left) {
            Node *y = z->parent->parent->right;

            if (y && y->color == 'R') {
                z->parent->color = 'B';
                y->color = 'B';
                z->parent->parent->color = 'R';
                z = z->parent->parent;
            } else {
                if (z == z->parent->right) {
                    z = z->parent;
                    leftRotate(z);
                }
                z->parent->color = 'B';
                z->parent->parent->color = 'R';
                rightRotate(z->parent->parent);
            }
        } else {
            Node *y = z->parent->parent->left;

            if (y && y->color == 'R') {
                z->parent->color = 'B';
                y->color = 'B';
                z->parent->parent->color = 'R';
                z = z->parent->parent;
            } else {
                if (z == z->parent->left) {
                    z = z->parent;
                    rightRotate(z);
                }
                z->parent->color = 'B';
                z->parent->parent->color = 'R';
                leftRotate(z->parent->parent);
            }
        }
    }
    root->color = 'B';
}

// Insert
void insert(int key) {
    Node *z = newNode(key);
    Node *y = NULL;
    Node *x = root;

    while (x) {
        y = x;
        if (key < x->key)
            x = x->left;
        else
            x = x->right;
    }

    z->parent = y;

    if (!y)
        root = z;
    else if (key < y->key)
        y->left = z;
    else
        y->right = z;

    fixInsert(z);
}

// Find minimum
Node* minimum(Node* node) {
    while (node->left)
        node = node->left;
    return node;
}

// Transplant
void transplant(Node *u, Node *v) {
    if (!u->parent)
        root = v;
    else if (u == u->parent->left)
        u->parent->left = v;
    else
        u->parent->right = v;

    if (v)
        v->parent = u->parent;
}

// Fix delete (core logic)
void fixDelete(Node *x) {
    while (x != root && (!x || x->color == 'B')) {

        if (x == x->parent->left) {
            Node *w = x->parent->right;

            if (w && w->color == 'R') {
                w->color = 'B';
                x->parent->color = 'R';
                leftRotate(x->parent);
                w = x->parent->right;
            }

            if ((!w->left || w->left->color == 'B') &&
                (!w->right || w->right->color == 'B')) {
                w->color = 'R';
                x = x->parent;
            } else {
                if (!w->right || w->right->color == 'B') {
                    if (w->left) w->left->color = 'B';
                    w->color = 'R';
                    rightRotate(w);
                    w = x->parent->right;
                }
                w->color = x->parent->color;
                x->parent->color = 'B';
                if (w->right) w->right->color = 'B';
                leftRotate(x->parent);
                x = root;
            }
        } else {
            Node *w = x->parent->left;

            if (w && w->color == 'R') {
                w->color = 'B';
                x->parent->color = 'R';
                rightRotate(x->parent);
                w = x->parent->left;
            }

            if ((!w->right || w->right->color == 'B') &&
                (!w->left || w->left->color == 'B')) {
                w->color = 'R';
                x = x->parent;
            } else {
                if (!w->left || w->left->color == 'B') {
                    if (w->right) w->right->color = 'B';
                    w->color = 'R';
                    leftRotate(w);
                    w = x->parent->left;
                }
                w->color = x->parent->color;
                x->parent->color = 'B';
                if (w->left) w->left->color = 'B';
                rightRotate(x->parent);
                x = root;
            }
        }
    }
    if (x) x->color = 'B';
}

// Delete node
void deleteNode(int key) {
    Node *z = root;

    while (z && z->key != key)
        z = (key < z->key) ? z->left : z->right;

    if (!z) return;

    Node *y = z;
    Node *x;
    char y_original_color = y->color;

    if (!z->left) {
        x = z->right;
        transplant(z, z->right);
    } else if (!z->right) {
        x = z->left;
        transplant(z, z->left);
    } else {
        y = minimum(z->right);
        y_original_color = y->color;
        x = y->right;

        if (y->parent == z) {
            if (x) x->parent = y;
        } else {
            transplant(y, y->right);
            y->right = z->right;
            y->right->parent = y;
        }

        transplant(z, y);
        y->left = z->left;
        y->left->parent = y;
        y->color = z->color;
    }

    if (y_original_color == 'B')
        fixDelete(x);
}

// Inorder print
void inorder(Node *root) {
    if (root) {
        inorder(root->left);
        printf("%d(%c) ", root->key, root->color);
        inorder(root->right);
    }
}

// Driver
int main() {
    insert(20);
    insert(15);
    insert(25);
    insert(10);
    insert(5);

    printf("After insertions:\n");
    inorder(root);

    deleteNode(15);

    printf("\nAfter deletion:\n");
    inorder(root);

    return 0;
}

*/

#include <stdio.h>

void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
}

void heapify(int arr[], int n, int i) {
    int largest = i;
    int l = 2*i + 1, r = 2*i + 2;

    if (l < n && arr[l] > arr[largest])
        largest = l;

    if (r < n && arr[r] > arr[largest])
        largest = r;

    if (largest != i) {
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

void buildMaxHeap(int arr[], int n) {
    for (int i = n/2 - 1; i >= 0; i--)
        heapify(arr, n, i);
}

int kthLargest(int arr[], int n, int k) {
    buildMaxHeap(arr, n);
    int size = n;

    for (int i = 0; i < k-1; i++) {
        swap(&arr[0], &arr[size-1]);
        size--;
        heapify(arr, size, 0);
    }
    return arr[0];
}

int main() {
    int arr[] = {5,77,2,3,44,22,1,65};
    int n = 8;

    buildMaxHeap(arr, n);

    printf("Max Heap:\n");
    for(int i=0;i<n;i++) printf("%d ", arr[i]);

    printf("\n3rd largest: %d\n", kthLargest(arr,n,3));

    return 0;
}