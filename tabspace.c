#include <stdio.h>

int main(int argc, char**argv) {
    int c;
    do {
        c = getc(stdin);
        if (c > 0) {
            if (c == '\t') {
                printf("    ");
            } else {
                printf("%c",c);
            }
        }
    } while (c > 0);
}
