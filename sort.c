#include "sort.h"
#include <stdio.h>

void swap(double *x, double *y){
   double temp;
   temp = *x;
   *x = *y;
   *y = temp;
}

void bublesort(int n, double *a){
   int i,j;
   for(i=0;i<(n-1);i++)
      for(j=0;j<(n-(i+1));j++)
             if(a[j] > a[j+1])
                    swap(&a[j],&a[j+1]);
}

