#include <stdio.h>

// h_prime function
int h_prime(char *s) {
    int h = 7;   // h[0]

    for (int i = 0; s[i] != '\0'; i++) {
        h = 31 * h + (int)s[i];  // ascii(s[i])
    }

    return h;  // h[n]
}

// Driver code
int main() {
    char str[100];

    printf("Enter string: ");
    scanf("%s", str);

    printf("Hash value = %d\n", h_prime(str));

    return 0;
}