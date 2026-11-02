/* glibc strtod/atoi over every line of a file (C locale): consumed length and bits. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "rb"); char* line = NULL; size_t cap = 0; ssize_t n;
    while ((n = getline(&line, &cap, f)) >= 0) {
        if (n > 0 && line[n-1] == '\n') line[--n] = 0;
        char* end; double v = strtod(line, &end); uint64_t u; memcpy(&u, &v, 8);
        printf("%ld %016llx %d\n", (long)(end - line), (unsigned long long)u, atoi(line));
    }
}
