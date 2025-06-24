#include <stdio.h>
#include <math.h>


int main(int argc, char* argv[]) {
    // ------------------------------------ stdout(1) and stderr(2) -------------------------------------
    /* printf("This goes to stdout.\n"); */
    /* fprintf(stdout, "This is stdout using fprintf.\n"); */
    /* fprintf(stderr, "Error: something went wrong!\n"); */
    /* // ./testc 1> output.txt */
    // ./testc 2> error.txt
    // ./testc > all.txt 2>&1
    // ./testc &> all.txt
    // ------------------------------------ M_PI --------------------------------------------
    printf("%lf\n", M_PI);

    return 0;
}
