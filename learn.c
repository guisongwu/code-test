#include <stdio.h>
#include "learn.h"
//extern int Add(int a, int b);
//extern int a;

void print(Student* ps){
    printf("%s %d %s %s\n", (*ps).name, (*ps).age, (*ps).sex, (*ps).tele);
    printf("%s %d %s %s\n", ps->name, ps->age, ps->sex, ps->tele);
}



int main(int argc, char* argv[]){
    int i = 0;
    int j = 0;
    for (i=0,j=0; i<2&&j<4; i++,j++)
    {
        printf("hehe");
    }
    return 0;
}


