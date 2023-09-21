#include "InitArrPointer.h"
#include <malloc.h>
// 初始化二维数组指针，并赋初值为0
double **InitArrPointer(int row, int col)
{
    double **p;
    p = (double **)malloc(sizeof(double *) * row); // 开辟行
    for (int i = 0; i < row; i++)
    {
        *(p + i) = (double *)malloc(sizeof(double) * col); // 开辟列
    }
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            p[i][j] = 0.0;
        }
    }
    return p;
}

// 初始化二维数组指针，并赋初值为init
double **InitArrPointer(int row, int col, int init)
{
    double **p;
    p = (double **)malloc(sizeof(double *) * row); // 开辟行
    for (int i = 0; i < row; i++)
    {
        *(p + i) = (double *)malloc(sizeof(double) * col); // 开辟列
    }
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            p[i][j] = init;
        }
    }
    return p;
}

int RowsOfArray(double **array)
{
    int a = sizeof(array) / sizeof(array[0]);
    return a;
}

int ColsOfArray(double **array)
{
    return sizeof(array[0]) / sizeof(double);
}