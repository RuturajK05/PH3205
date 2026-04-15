/*
 * ============================================================
 *  CS3201 – Class Test 3 Code Bank: HEAPS + HASHING
 *  Covers every PPT pseudocode + predicted question answers
 * ============================================================
 *
 *  INDEXING CONVENTION:
 *    All heap arrays are 1-indexed (A[1] = root, as in CLRS).
 *    So declare:  int A[MAX+1];  and set A[0] unused.
 *
 *  COMPILE:   gcc -o bank heap_hash_codebank.c
 * ============================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>   /* INT_MIN, INT_MAX */

/* ============================================================
 *  SECTION 1 – HEAP PRIMITIVES  (1-indexed)
 * ============================================================ */

#define MAX_HEAP 200
int heap[MAX_HEAP];     /* 1-indexed global heap */
int heap_size = 0;

/* --- index helpers --- */
int left(int i)   { return 2 * i; }
int right(int i)  { return 2 * i + 1; }
int parent(int i) { return i / 2; }

void swap(int *a, int *b) { int t = *a; *a = *b; *b = t; }

/* ============================================================
 *  1a. MAX-HEAPIFY  (sift-down)
 *      Restores heap property at node i.
 *      Assumption: left/right subtrees of i are already max-heaps.
 *      Time: O(log n)
 * ============================================================ */
void max_heapify(int A[], int n, int i)
{
    int l = left(i), r = right(i), largest = i;

    if (l <= n && A[l] > A[largest]) largest = l;
    if (r <= n && A[r] > A[largest]) largest = r;

    if (largest != i) {
        swap(&A[i], &A[largest]);
        max_heapify(A, n, largest);   /* recurse on affected subtree */
    }
}

/* ============================================================
 *  1b. BUILD-MAX-HEAP
 *      Converts arbitrary array A[1..n] into a max-heap.
 *      Time: O(n)  — tight bound proven via geometric series
 * ============================================================ */
void build_max_heap(int A[], int n)
{
    for (int i = n / 2; i >= 1; i--)
        max_heapify(A, n, i);
}

/* ============================================================
 *  1c. HEAP-EXTRACT-MAX
 *      Removes and returns the max element.
 *      Time: O(log n)
 * ============================================================ */
int heap_extract_max(int A[], int *n)
{
    if (*n < 1) { fprintf(stderr, "heap underflow\n"); return INT_MIN; }
    int max = A[1];
    A[1] = A[*n];
    (*n)--;
    max_heapify(A, *n, 1);
    return max;
}

/* ============================================================
 *  1d. HEAP-INCREASE-KEY
 *      Increases A[i] to key, then bubbles up.
 *      Time: O(log n)
 * ============================================================ */
void heap_increase_key(int A[], int n, int i, int key)
{
    if (key < A[i]) { fprintf(stderr, "new key is smaller\n"); return; }
    A[i] = key;
    while (i > 1 && A[parent(i)] < A[i]) {
        swap(&A[i], &A[parent(i)]);
        i = parent(i);
    }
}

/* ============================================================
 *  1e. MAX-HEAP-INSERT
 *      Inserts a new key into the heap.
 *      Time: O(log n)
 * ============================================================ */
void max_heap_insert(int A[], int *n, int key)
{
    (*n)++;
    A[*n] = INT_MIN;                    /* sentinel: -inf */
    heap_increase_key(A, *n, *n, key);
}

/* ============================================================
 *  1f. HEAP-SORT
 *      Sorts A[1..n] in ascending order using heap.
 *      Time: O(n log n)
 * ============================================================ */
void heap_sort(int A[], int n)
{
    build_max_heap(A, n);
    for (int i = n; i >= 2; i--) {
        swap(&A[1], &A[i]);
        max_heapify(A, i - 1, 1);
    }
    /* A[1..n] is now sorted ascending */
}

/* ============================================================
 *  1g. MIN-HEAP variant  (swap comparisons)
 *      Just for quick reference / predicted min-heap question.
 * ============================================================ */
void min_heapify(int A[], int n, int i)
{
    int l = left(i), r = right(i), smallest = i;
    if (l <= n && A[l] < A[smallest]) smallest = l;
    if (r <= n && A[r] < A[smallest]) smallest = r;
    if (smallest != i) {
        swap(&A[i], &A[smallest]);
        min_heapify(A, n, smallest);
    }
}

void build_min_heap(int A[], int n)
{
    for (int i = n / 2; i >= 1; i--)
        min_heapify(A, n, i);
}

int heap_extract_min(int A[], int *n)
{
    if (*n < 1) { fprintf(stderr, "heap underflow\n"); return INT_MAX; }
    int min = A[1];
    A[1] = A[*n];
    (*n)--;
    min_heapify(A, *n, 1);
    return min;
}

/* ============================================================
 *  SECTION 2 – PREDICTED QUESTION SOLUTIONS
 * ============================================================ */

/* --- Q2: K-th largest in O(n) + (k-1)*O(log n) ----------- */
/*
 *  Algorithm:
 *    1. build_max_heap → O(n)
 *    2. Call extract_max (k-1) times → each O(log n)
 *    3. A[1] is now the k-th largest
 */
int kth_largest(int A[], int n, int k)
{
    /* A is 1-indexed, n elements */
    build_max_heap(A, n);
    for (int i = 0; i < k - 1; i++)
        heap_extract_max(A, &n);
    return A[1];
}

/* --- Q3: Priority-Queue Clinic (SWASTH) ------------------- */
/*
 *  Each patient: age (priority) + patient_id.
 *  Max-heap on age → eldest patient first.
 *  insert  → max_heap_insert   O(log n)
 *  serve   → heap_extract_max  O(log n)
 */
typedef struct {
    int age;
    int patient_id;
} Patient;

#define MAX_PATIENTS 100
Patient clinic[MAX_PATIENTS + 1]; /* 1-indexed */
int clinic_size = 0;

void clinic_swap(int i, int j)
{
    Patient tmp = clinic[i];
    clinic[i] = clinic[j];
    clinic[j] = tmp;
}

void clinic_heapify_up(int i)
{
    while (i > 1 && clinic[parent(i)].age < clinic[i].age) {
        clinic_swap(i, parent(i));
        i = parent(i);
    }
}

void clinic_heapify_down(int i)
{
    int l = left(i), r = right(i), largest = i;
    if (l <= clinic_size && clinic[l].age > clinic[largest].age) largest = l;
    if (r <= clinic_size && clinic[r].age > clinic[largest].age) largest = r;
    if (largest != i) {
        clinic_swap(i, largest);
        clinic_heapify_down(largest);
    }
}

/* Patient arrives → insert */
void patient_arrive(int age, int pid)
{
    clinic_size++;
    clinic[clinic_size].age = age;
    clinic[clinic_size].patient_id = pid;
    clinic_heapify_up(clinic_size);   /* O(log n) */
}

/* Patient is served → extract eldest */
Patient patient_serve()
{
    if (clinic_size == 0) { fprintf(stderr, "No patients\n"); Patient p = {-1,-1}; return p; }
    Patient eldest = clinic[1];
    clinic[1] = clinic[clinic_size];
    clinic_size--;
    clinic_heapify_down(1);           /* O(log n) */
    return eldest;
}

/* ============================================================
 *  SECTION 3 – HASHING
 * ============================================================ */

/* --- Q1 (Hashing): h_sum ---------------------------------- */
/*
 *  h_sum(x) = sum of ASCII values of characters of x.
 *  Obvious flaw: anagrams ("stop","pots") → same hash → collision.
 */
int h_sum(const char *s)
{
    int sum = 0;
    while (*s) sum += (unsigned char)(*s++);
    return sum;
}

/* --- Q2 (Hashing): h_prime (DJB2 variant) ----------------- */
/*
 *  Recurrence:  h[0] = 7
 *               h[i] = 31*h[i-1] + ascii(s[i-1])   for i=1..n
 *  Returns h[n].
 *  Better than h_sum: positional sensitivity → "stop" ≠ "pots".
 */
int h_prime(const char *s)
{
    int h = 7;
    int n = strlen(s);
    for (int i = 1; i <= n; i++)
        h = 31 * h + (unsigned char)s[i - 1];
    return h;
}

/* ============================================================
 *  3a. HASH TABLE WITH CHAINING
 *      Collision resolution: separate chaining with linked lists.
 *      Operations: insert, search, delete
 *      Avg case O(1+α) where α = n/m (load factor)
 * ============================================================ */
#define CHAIN_M 11   /* table size (prime recommended) */

typedef struct ChainNode {
    int key;
    struct ChainNode *next;
} ChainNode;

ChainNode *chain_table[CHAIN_M]; /* initialise all to NULL */

int chain_hash(int k) { return k % CHAIN_M; }

/* Insert at head of chain */
void chain_insert(int key)
{
    int h = chain_hash(key);
    ChainNode *node = malloc(sizeof(ChainNode));
    node->key = key;
    node->next = chain_table[h];
    chain_table[h] = node;          /* O(1) */
}

/* Search; returns 1 if found */
int chain_search(int key)
{
    int h = chain_hash(key);
    ChainNode *cur = chain_table[h];
    while (cur) {
        if (cur->key == key) return 1;
        cur = cur->next;
    }
    return 0;
}

/* Delete first occurrence */
void chain_delete(int key)
{
    int h = chain_hash(key);
    ChainNode **pp = &chain_table[h];
    while (*pp) {
        if ((*pp)->key == key) {
            ChainNode *tmp = *pp;
            *pp = tmp->next;
            free(tmp);
            return;
        }
        pp = &(*pp)->next;
    }
}

void chain_print()
{
    for (int i = 0; i < CHAIN_M; i++) {
        printf("[%2d]: ", i);
        for (ChainNode *c = chain_table[i]; c; c = c->next)
            printf("%d -> ", c->key);
        printf("NULL\n");
    }
}

/* ============================================================
 *  3b. OPEN ADDRESSING: LINEAR PROBING
 *      H(k,i) = (h(k) + i) mod m
 *      Suffers from PRIMARY CLUSTERING.
 *      Deletion: mark deleted slots with a DELETED sentinel.
 * ============================================================ */
#define OA_M 11
#define EMPTY   -1
#define DELETED -2  /* tombstone for deletion */

int lp_table[OA_M];

void lp_init() { for (int i=0;i<OA_M;i++) lp_table[i]=EMPTY; }

int lp_h(int k) { return k % OA_M; }

void lp_insert(int key)
{
    for (int i = 0; i < OA_M; i++) {
        int idx = (lp_h(key) + i) % OA_M;
        if (lp_table[idx] == EMPTY || lp_table[idx] == DELETED) {
            lp_table[idx] = key; return;
        }
    }
    printf("Table full!\n");
}

int lp_search(int key)
{
    for (int i = 0; i < OA_M; i++) {
        int idx = (lp_h(key) + i) % OA_M;
        if (lp_table[idx] == EMPTY) return -1;         /* not found */
        if (lp_table[idx] == key)   return idx;        /* found */
        /* DELETED → skip and continue probing */
    }
    return -1;
}

void lp_delete(int key)
{
    int idx = lp_search(key);
    if (idx != -1) lp_table[idx] = DELETED;
}

void lp_print()
{
    for (int i = 0; i < OA_M; i++) {
        if (lp_table[i] == EMPTY)   printf("[%2d]: EMPTY\n", i);
        else if (lp_table[i]==DELETED) printf("[%2d]: DELETED\n",i);
        else                         printf("[%2d]: %d\n", i, lp_table[i]);
    }
}

/* ============================================================
 *  3c. OPEN ADDRESSING: QUADRATIC PROBING
 *      H(k,i) = (h(k) + i²) mod m
 *      Suffers from SECONDARY CLUSTERING.
 *      Use prime m to ensure all slots reachable.
 * ============================================================ */
int qp_table[OA_M];

void qp_init() { for (int i=0;i<OA_M;i++) qp_table[i]=EMPTY; }

void qp_insert(int key)
{
    for (int i = 0; i < OA_M; i++) {
        int idx = (lp_h(key) + i*i) % OA_M;
        if (qp_table[idx] == EMPTY || qp_table[idx] == DELETED) {
            qp_table[idx] = key; return;
        }
    }
    printf("Table full or insert failed!\n");
}

int qp_search(int key)
{
    for (int i = 0; i < OA_M; i++) {
        int idx = (lp_h(key) + i*i) % OA_M;
        if (qp_table[idx] == EMPTY) return -1;
        if (qp_table[idx] == key)   return idx;
    }
    return -1;
}

/* ============================================================
 *  3d. OPEN ADDRESSING: DOUBLE HASHING
 *      H(k,i) = (h(k) + i * h'(k)) mod m
 *      h(k)  = k mod m
 *      h'(k) = p - (k mod p)   where p < m, p prime
 *              (h'(k) must never be 0)
 *      No clustering issues — best of open addressing.
 * ============================================================ */
#define DH_M  11
#define DH_P   7   /* prime < DH_M */

int dh_table[DH_M];
void dh_init() { for(int i=0;i<DH_M;i++) dh_table[i]=EMPTY; }

int dh_h1(int k)  { return k % DH_M; }
int dh_h2(int k)  { return DH_P - (k % DH_P); }   /* always 1..DH_P */

void dh_insert(int key)
{
    int h1 = dh_h1(key), h2 = dh_h2(key);
    for (int i = 0; i < DH_M; i++) {
        int idx = (h1 + i * h2) % DH_M;
        if (dh_table[idx] == EMPTY || dh_table[idx] == DELETED) {
            dh_table[idx] = key; return;
        }
    }
    printf("Table full!\n");
}

int dh_search(int key)
{
    int h1 = dh_h1(key), h2 = dh_h2(key);
    for (int i = 0; i < DH_M; i++) {
        int idx = (h1 + i * h2) % DH_M;
        if (dh_table[idx] == EMPTY) return -1;
        if (dh_table[idx] == key)   return idx;
    }
    return -1;
}

/* ============================================================
 *  SECTION 4 – UTILITY / PRINT HELPERS
 * ============================================================ */

void print_array(int A[], int n)   /* 1-indexed */
{
    for (int i = 1; i <= n; i++) printf("%d ", A[i]);
    printf("\n");
}

void print_array0(int A[], int n)  /* 0-indexed */
{
    for (int i = 0; i < n; i++) printf("%d ", A[i]);
    printf("\n");
}

/* ============================================================
 *  SECTION 5 – MAIN: DEMOS FOR EVERY FUNCTION
 *  Read this to verify your understanding, then adapt.
 * ============================================================ */

int main()
{
    /* ------ HEAP DEMOS ------ */
    printf("=== BUILD-MAX-HEAP demo ===\n");
    /* Q1: array {5,77,2,3,44,22,1,65} step by step */
    int A[20] = {0, 5, 77, 2, 3, 44, 22, 1, 65};  /* 1-indexed, n=8 */
    int n = 8;
    printf("Before: "); print_array(A, n);
    build_max_heap(A, n);
    printf("After BUILD-MAX-HEAP: "); print_array(A, n);
    /* Should be: 77 65 22 3 44 2 1 5 (or similar valid max-heap) */
    /* Verify at https://heapsortvisualizer.web.app/ */

    printf("\n=== HEAP-SORT demo ===\n");
    int B[10] = {0, 5, 77, 2, 3, 44, 22, 1, 65};
    int nb = 8;
    heap_sort(B, nb);
    printf("Sorted: "); print_array(B, nb);

    printf("\n=== K-TH LARGEST (Q2) ===\n");
    int C[10] = {0, 10, 4, 3, 50, 8, 25};  /* n=6 */
    int nc = 6;
    printf("3rd largest: %d\n", kth_largest(C, nc, 3));  /* expected 10 */

    printf("\n=== INSERT + EXTRACT-MAX ===\n");
    int D[20]; int nd = 0;
    max_heap_insert(D, &nd, 5);
    max_heap_insert(D, &nd, 77);
    max_heap_insert(D, &nd, 44);
    max_heap_insert(D, &nd, 3);
    printf("Heap after inserts: "); print_array(D, nd);
    printf("Extracted: %d\n", heap_extract_max(D, &nd));
    printf("Heap after extract: "); print_array(D, nd);

    printf("\n=== CLINIC / PRIORITY QUEUE (Q3) ===\n");
    patient_arrive(35, 101);
    patient_arrive(72, 102);
    patient_arrive(60, 103);
    patient_arrive(45, 104);
    printf("Serving (eldest first):\n");
    while (clinic_size > 0) {
        Patient p = patient_serve();
        printf("  Patient %d, age %d\n", p.patient_id, p.age);
    }

    /* ------ HASHING DEMOS ------ */
    printf("\n=== h_sum and h_prime (Q1,Q2 Hashing) ===\n");
    printf("h_sum(\"stop\") = %d\n", h_sum("stop"));
    printf("h_sum(\"pots\") = %d\n", h_sum("pots"));   /* same → collision! */
    printf("h_prime(\"stop\") = %d\n", h_prime("stop"));
    printf("h_prime(\"pots\") = %d\n", h_prime("pots")); /* different → better */

    printf("\n=== CHAINING ===\n");
    memset(chain_table, 0, sizeof(chain_table));
    int keys[] = {45, 34, 22, 60, 52, 13};
    for (int i = 0; i < 6; i++) chain_insert(keys[i]);
    chain_print();
    printf("Search 52: %s\n", chain_search(52) ? "FOUND" : "NOT FOUND");
    chain_delete(52);
    printf("After delete 52:\n"); chain_print();

    printf("\n=== LINEAR PROBING ===\n");
    lp_init();
    int lpkeys[] = {29, 46, 18, 36, 43, 21, 24, 23};
    for (int i = 0; i < 8; i++) lp_insert(lpkeys[i]);
    lp_print();

    printf("\n=== QUADRATIC PROBING ===\n");
    qp_init();
    int qpkeys[] = {46, 28, 21, 35, 57, 39, 19, 50};
    for (int i = 0; i < 8; i++) qp_insert(qpkeys[i]);
    for (int i = 0; i < OA_M; i++) {
        if (qp_table[i]==EMPTY)   printf("[%2d]: EMPTY\n",i);
        else if(qp_table[i]==DELETED) printf("[%2d]: DELETED\n",i);
        else printf("[%2d]: %d\n", i, qp_table[i]);
    }

    printf("\n=== DOUBLE HASHING ===\n");
    dh_init();
    int dhkeys[] = {46, 18, 36, 43, 21, 24, 54};
    for (int i = 0; i < 7; i++) dh_insert(dhkeys[i]);
    for (int i = 0; i < DH_M; i++) {
        if (dh_table[i]==EMPTY)    printf("[%2d]: EMPTY\n",i);
        else if(dh_table[i]==DELETED) printf("[%2d]: DELETED\n",i);
        else printf("[%2d]: %d\n", i, dh_table[i]);
    }

    return 0;
}
