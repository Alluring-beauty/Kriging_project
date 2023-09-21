#include <math.h>
#include <float.h>
#include <iostream>

// int compare_doubles(const void *a, const void *b)
// {
//     double x = *(double *)a;
//     double y = *(double *)b;
//     if (x > y)
//     {
//         return 1;
//     }
//     else if (x < y)
//     {
//         return -1;
//     }
//     else
//     {
//         return 0;
//     }
// }

// double normal_cdf(double x)
// {
//     return 0.5 * (1.0 + erf(x / sqrt(2.0)));
// }

// double ks_test(double *data, int data_size)
// {
//     double max_diff = 0;
//     double temp[data_size];
//     double sum = 0;
//     double ave;
//     for (int i = 0; i < data_size; i++)
//     {
//         sum += data[i];
//     }
//     ave = sum / data_size;
//     double var = 0;
//     for (int i = 0; i < data_size; i++)
//     {
//         var += pow(data[i] - ave, 2) / data_size;
//     }
//     var = sqrt(var);
//     for (int i = 0; i < data_size; i++)
//     {
//         temp[i] = data[i];
//     }
//     qsort(temp, data_size, sizeof(double), compare_doubles);
//     for (int i = 0; i < data_size; i++)
//     {
//         double diff = fabs(normal_cdf((temp[i] - ave) / var) - (double)(i + 1) / data_size);
//         if (diff > max_diff)
//         {
//             max_diff = diff;
//         }
//     }
//     return max_diff;
// }
// 对数据进行BoxCox变换，使数据满足正态
// double BoxCox_before(double *data, int data_size)
// {
//     double lamda = 0;
//     //
//     double temp[data_size];
//     for (int i = 0; i < data_size; i++)
//     {
//         temp[i] = data[i];
//         std::cout << temp[i] << std::endl;
//     }

//     double max_diff = 0;
//     double min = 1;
//     double result_lamda = 0;
//     double pValue = 0;
//     for (lamda; lamda < 40; lamda += 0.01)
//     {
//         double change[data_size];
//         if (lamda == 0)
//         {
//             for (int i = 0; i < data_size; i++)
//             {
//                 change[i] = log(data[i]);
//             }
//         }
//         else if (lamda > 0)
//         {
//             for (int i = 0; i < data_size; i++)
//             {
//                 change[i] = (pow(data[i], lamda) - 1) / lamda;
//             }
//         }
//         // for (int i = 0; i < data_size; i++)
//         // {
//         //     std::cout << "change:" << change[i] << std::endl;
//         // }
//         max_diff = ks_test(change, data_size);
//         std::cout << "max_diff:" << max_diff << std::endl;
//         if (max_diff < min)
//         {
//             min = max_diff;
//             result_lamda = lamda;
//         }
//         for (int i = 0; i < data_size; i++)
//         {
//             data[i] = temp[i];
//         }
//     }
//     return result_lamda;
// }

// 根据设计的3sigma原则，返回正态评价标准
double sigma_3(double *data, int data_size)
{
    // 保存原始数据
    double temp[data_size];
    for (int i = 0; i < data_size; i++)
    {
        temp[i] = data[i];
    }
    double ave = 0;
    double sigma = 0;
    for (int i = 0; i < data_size; i++)
    {
        ave += temp[i];
    }
    ave /= data_size;
    for (int i = 0; i < data_size; i++)
    {
        sigma += pow(data[i] - ave, 2) / data_size;
    }
    sigma = sqrt(sigma);
    double p_1, p_2, p_3;
    int count_1 = 0, count_2 = 0, count_3 = 0;
    for (int i = 0; i < data_size; i++)
    {
        if (temp[i] > ave - sigma && temp[i] < ave + sigma)
        {
            count_1++;
        }
        if (temp[i] > ave - 2 * sigma && temp[i] < ave + 2 * sigma)
        {
            count_2++;
        }
        if (temp[i] > ave - 3 * sigma && temp[i] < ave + 3 * sigma)
        {
            count_3++;
        }
    }
    p_1 = (double)count_1 / data_size;
    p_2 = (double)count_2 / data_size;
    p_3 = (double)count_3 / data_size;
    std::cout << "count_1:" << count_1 << "\tcount_1:" << count_2 << "\tcount_3:" << count_3 << std::endl;
    std::cout << "p_1:" << p_1 << "\tp_2:" << p_2 << "\tp_3:" << p_3 << std::endl;
    return 3 * (p_1 - 0.683) + 2 * (p_2 - 0.954) + 1 * (p_3 - 0.997);
}

// 返回对数据的卡方检验的值
double kafang(double *data, int data_size)
{
    // 保存原始数据,确保每次调用的原始数据不会被改变
    double temp[data_size];
    for (int i = 0; i < data_size; i++)
    {
        temp[i] = data[i];
    }
    double ave = 0;
    double sigma = 0;
    for (int i = 0; i < data_size; i++)
    {
        ave += temp[i];
    }
    ave /= data_size;
    for (int i = 0; i < data_size; i++)
    {
        sigma += pow(data[i] - ave, 2) / data_size;
    }
    sigma = sqrt(sigma);
    double e_1, e_2, e_3, e_4, e_5, e_6;
    e_1 = data_size * (0.49865 - 0.4772);
    e_2 = data_size * (0.4772 - 0.3413);
    e_3 = data_size * 0.3413;
    e_4 = e_3;
    e_5 = e_2;
    e_6 = e_1;
    double o_1 = 0, o_2 = 0, o_3 = 0, o_4 = 0, o_5 = 0, o_6 = 0, o_7 = 0, o_8 = 0, o_9 = 0, o_10 = 0;
    for (int i = 0; i < data_size; i++)
    {
        if (temp[i] > ave - 3 * sigma && temp[i] < ave - 2 * sigma)
            o_1++;
        if (temp[i] > ave - 2 * sigma && temp[i] < ave - sigma)
            o_2++;
        if (temp[i] > ave - sigma && temp[i] < ave)
            o_3++;
        if (temp[i] > ave && temp[i] < ave + sigma)
            o_4++;
        if (temp[i] > ave + sigma && temp[i] < ave + 2 * sigma)
            o_5++;
        if (temp[i] > ave + 2 * sigma && temp[i] < ave + 3 * sigma)
            o_6++;
    }
    printf("sum of o:%lf\n", o_1 + o_2 + o_3 + o_4 + o_5 + o_6 + o_7 + o_8 + o_9 + o_10);
    return pow(o_1 - e_1, 2) / e_1 + pow(o_2 - e_2, 2) / e_2 + pow(o_3 - e_3, 2) / e_3 +
           pow(o_4 - e_4, 2) / e_4 + pow(o_5 - e_5, 2) / e_5 + pow(o_6 - e_6, 2) / e_6;
}

// 根据步长，遍历所有的lamda，找到最优的lamda
double BoxCox(double *data, int data_size)
{
    double lamda = 0;
    double temp[data_size];
    for (int i = 0; i < data_size; i++)
    {
        temp[i] = data[i];
    }
    double min = 10000;
    double result = 0;
    for (lamda; lamda < 10; lamda += 0.01)
    {
        if (lamda == 0)
        {
            for (int i = 0; i < data_size; i++)
            {
                temp[i] = log(temp[i]);
            }
        }
        else if (lamda > 0)
        {
            for (int i = 0; i < data_size; i++)
            {
                temp[i] = (pow(temp[i], lamda) - 1) / lamda;
            }
        }
        double sigma_weight = kafang(temp, data_size);
        std::cout << "kafang:" << sigma_weight << std::endl;
        if (sigma_weight < min)
        {
            min = sigma_weight;
            result = lamda;
        }
        for (int i = 0; i < data_size; i++)
        {
            temp[i] = data[i];
        }
    }
    return result;
}