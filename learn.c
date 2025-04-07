#include <stdio.h>
//extern int a;

enum trafficLight{
    RED,
    GREEN,
    YELLOW
};


int main()
{
    // enum trafficLight light = RED;
    // printf("%d\n",light);
    // char arr[3] = "abc";
    // printf("%s\n",arr);
    //printf("%zu\n",sizeof(char));
    //printf("%zu\n",sizeof(short));

    //printf("%zu\n",sizeof(int));

    //printf("%zu\n",sizeof(long));
    //printf("%zu\n",sizeof(long long));
    //
    //printf("%zu\n",sizeof(float));
    //printf("%zu\n",sizeof(double));
    //
    //
    // int arr[10] = {1,2,3,4,5,6,7,8,9,10};
    // int i = 0;
    // while (i < 10)
    // {
    //     printf("%d\n",arr[i]);
    //     i++;
    // }
    // printf("%d\n",sizeof(arr));



    // int a = 10;
    // int b = a++;
    // int b = ++a;
    // printf("%d\n",b);

    // int arr[10] = {1,2,3,4,5,6,7,8,9,10};
    enum Fruit{
        APPLE,
        ORANGE,
        BANANA
    };
    enum Fruit fruit = APPLE;
    switch (fruit)
    {
        case APPLE:
            printf("Apple\n");
            break;
        case ORANGE:
            printf("Orange\n");
            break;
        case BANANA:
            printf("Banana\n");
            break;
        default:
            printf("Unknown fruit\n");
            break;
    }


    return 0;
}
