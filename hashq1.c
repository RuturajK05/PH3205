#include <stdio.h>

// Function to compute hash value
int h_sum(char *x) {
    int sum = 0;

    for (int i = 0; x[i] != '\0'; i++) {
        sum += (int)x[i];   // add ASCII value of character
    }

    return sum;
}

// Driver code
int main() {
    char str[100];

    printf("Enter a string: ");
    scanf("%s", str);

    int hash = h_sum(str);

    printf("Hash value (h_sum) = %d\n", hash);

    return 0;
}