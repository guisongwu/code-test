#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>




bool get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (true) {
        memset(token, 0, 100);
        if (fscanf(fp, "%s", token) != 1)  // everytime 'fscanf' is called, the pointer will move to next string.
            return false;
        if (token[0] != '#')  
            break;
        /* skip to newline */
        do {
            if ((c = fgetc(fp)) == EOF)
                return false;
        } while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
        *p = '\0';
    return true;
}


int main(int argc, char* argv[]) {
    
    // ------------------------------- isalpha(int c) ---------------------------------
    // instead of isalpha(char c)
    /* printf("%d\n", isalpha(9)); */
    /* printf("%s\n", isalpha(9) ? "TRUE" : "FALSE"); */
    /* printf("%d\n", isalpha('a')); */
    /* printf("%s\n", isalpha('a') ? "TRUE" : "FALSE"); */
    /* printf("%d\n", isalpha('A')); */
    /* printf("%s\n", isalpha('A') ? "TRUE" : "FALSE"); */
    /* printf("%d\n", isalpha((int)'A')); */
    /* printf("%s\n", isalpha((int)'A') ? "TRUE" : "FALSE"); */

    // ------------------------------- atoi and atof -----------------------------------
    /* char *str = "98.8"; */
    /* printf("%d\n", atoi(str)); // Skips leading whitespace (spaces, tabs, etc.). */
    /*                            // Reads optional sign (+ or -). */
    /*                            // Converts subsequent digits until a non-digit character is encountered. */
    /*                            // Ignores the rest of the string. */
    /* printf("%d\n", atoi("abc123")); */
    /* printf("%d\n", atoi("-123abc")); */
    /* printf("%d\n", atoi("   -\n12abc123")); */
    /* printf("%d\n", atoi("   -  12abc123")); */
    /* printf("%d\n", atoi("     -12abc123")); */
    /* printf("%f\n", atof("0.234")); */
    /* printf("%f\n", atof("0.234abc")); */
    /* printf("%f\n", atof("abc0.234")); */

    // ------------------------------------ goto ----------------------------------------
    /* goto error; */
    /* printf("Test if this sentence can be executed or not\n"); // this will be ignored. */
/* error: */
    /* printf("This is a test for goto sentence\n"); */

    /* goto errorinif; */
    /* if (0) { */
/* errorinif: */
    /*     printf("This is a test for goto sentence in a if condition.\n"); */
    /* } */

    // ---------------------------------- __LINE__ --------------------------------------
    /* printf("%s\n", __FILE__); */
    /* printf("%d\n", __LINE__); */
    /* printf("%d\n", __LINE__); */
    /* printf("%d\n", __LINE__); */

    // ----------------------------- strcmp and strcasecmp -------------------------------
    /* char *str1 = "helloworld"; */
    /* char *str2 = "HelloWorld"; */
    /* printf("%d\n", strcmp(str1, str2));  // case-sensitive */
    /*                                      // >0: str1[0] > str2[0] */
    /*                                      // <0: str1[0] < str2[0] */
    /*                                      // =0: str1 == str2 */
    /* printf("%d\n", strcasecmp(str1, str2));  // case-insensitive */

    // ------------------------------------ get_token -------------------------------------
    /* FILE *file = fopen("quad.msh", "r"); */
    /* if (file == NULL) { */
    /*     perror("Failed to open file"); */
    /*     return 1; */
    /* } */
    /* char buffer[256]; */ 
    /* bool flag = get_token(file, buffer); */
    /* printf("%s\n", buffer); */
    /* get_token(file, buffer); */
    /* printf("%s\n", buffer); */
    
    /* fclose(file);  // Always close the file! */

    // ------------------------------------ strchr ----------------------------------------
    /* char *str1 = "hello world";  // This is a read-only str. */
    /* char str2[] = "hello world";  // This is changable. */
    /* char *num = strchr(str2, 'o');  // char *strchr(const char *str, int c); */
    /*                                // Return pointer to the first occurrence of c in str, */
    /*                                // or NULL if the character is not found. */
    /* printf("%c\n", *num); */
    /* *num = 'e'; */
    /* printf("%s\n", str2); */

    // ----------------------------------- Write File --------------------------------------
    /* FILE *file = fopen("example.txt", "w");  // Open in write mode ("w") */
    /* if (file == NULL) { */
    /*     perror("Failed to open file");  // Error handling */
    /*     return 1; */
    /* } */
    /* // Write data to the file */
    /* fprintf(file, "Hello, World!\n");  // Formatted write (like printf) */
    /* fputs("This is a line of text.\n", file);  // Writes a string */
    /* fputc('A', file);  // Writes a single character */
    /* fclose(file);  // Always close the file! */

    // ------------------------------------- Read File ---------------------------------------
    /* FILE *file = fopen("quad.msh", "r"); */
    /* if (file == NULL) { */
    /*     perror("Failed to open file"); */
    /*     return 1; */
    /* } */
    /* char buffer[256];  // Buffer to store each line */
    /* // Read line by line until EOF (End of File) */
    /* while (fgets(buffer, sizeof(buffer), file) != NULL) {  // char *fgets(char *str, int size, FILE *file); */
    /*     printf("%s", buffer);  // Print each line */
    /* } */
    /* // Read structured data */
    /* while (fscanf(file, "%s", buffer) == 1) {  // Read strings */
    /*     printf("string: %s\n", buffer); */
    /* } */
    /* fclose(file);  // Always close the file! */

    // -------------------------------------- memset ------------------------------------------
    /* char token[10]; */
    /* memset(token, 0, 10); */
    /* printf("%s\n", token); */
    /* memcpy(token, "hello world", 10); // void *memcpy(void *dest, const void *src, size_t n); */
    /* printf("%s\n", token); */
    /* memset(token, 0, 10); */
    /* printf("%s\n", token); */
    /* strcpy(token, "hello shi"); // char *strcpy(char *dest, const char *src); */
    /* printf("%s\n", token); */

    // --------------------------------------- char --------------------------------------------
    /* char a = 65; */
    /* printf("%c\n", a); */

    // ------------------------- COORD *verts and COORD *verts[2] -----------------------------
    double points[4][2] = {
        {1,2},
        {2,1},
        {3,4},
        {4,5}
    };

    typedef double Coord[2];
    Coord *ref_coord = &points[2];
    printf("%f\t%f\n", (*ref_coord)[0], (*ref_coord)[1]);

    /* Coord *vert; */
    /* vert = (Coord *) calloc (5, sizeof(Coord)); */
    /* for (int i = 0; i < 5; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         vert[i][j] = (double)i*j; */
    /*         printf("%lf\t", vert[i][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */

    /* Coord vert1 = {1.0, 2.0}; */
    /* Coord vert2 = {4.0, 5.0}; */
    /* Coord *verts[2] = {&vert1, &vert2}; */
    /* for (int i = 0; i < 3; i++) { */
    /*     for (int j = 0; j < 2; j++) { */
    /*         printf("%lf\t", (*verts[i])[j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */



    // --------------------------- stdout(1) and stderr(2) ---------------------------------
    /* printf("This goes to stdout.\n"); */
    /* fprintf(stdout, "This is stdout using fprintf.\n"); */
    /* fprintf(stderr, "Error: something went wrong!\n"); */
    /* // ./testc 1> output.txt */
    // ./testc 2> error.txt
    // ./testc > all.txt 2>&1
    // ./testc &> all.txt

    // ------------------------------------ M_PI --------------------------------------------
    /* printf("%lf\n", M_PI); */

    return 0;
}
