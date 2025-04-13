#include <stdio.h>
#include "learn.h"
#include <stdbool.h>
//extern int Add(int a, int b);
//extern int a;

void print(Student* ps){
    printf("%s %d %s %s\n", (*ps).name, (*ps).age, (*ps).sex, (*ps).tele);
    printf("%s %d %s %s\n", ps->name, ps->age, ps->sex, ps->tele);
}

void swap(int* p1, int* p2) {
    int c = 0;
    c = *p1;
    *p1 = *p2;
    *p2 = c;
}

int binary_search(int arr[], int target, int sz) {
    int left = 0;
    int right = sz - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] < target) {
            left = mid + 1;
        } else if (arr[mid] > target) {
            right = mid - 1;
        } else {
            return mid;
        }
    }
    return -1;
}

int main(int argc, char* argv[]) {
    //int x = 10;
    //int y = 20;
    //printf("%d %d\n", x, y);
    //swap(&x, &y);
    //printf("%d %d\n", x, y);
    int arr[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    int target = 8;
    int sz = sizeof(arr) / sizeof(arr[0]);
    int position = binary_search(arr, target, sz);
    if (position == -1) {
        printf("target doesn't exsit\n");
    } else {
        printf("done, index is %d\n", position);
    }
    printf("%zu\n",sizeof(bool));
    printf("%zu\n",sizeof(char));
    printf("%zu\n",sizeof(size_t));


    return 0;
}


