#include <iostream>
#include <iomanip>
#include <math.h>
#include "InitArrPointer.h"
#include "LevenbergMarquardt.h"
#include "DataPreprocessing.h"
#include <eigen3/Eigen/Dense>
#include <fstream>
#define MAX_DISTANCE 17300
#define MAX_LENGTH 10000
using namespace Eigen;
using namespace std;

// 原始数据点结构体
struct sourceData
{
    int length;
    double x[MAX_LENGTH];
    double y[MAX_LENGTH];
    double z[MAX_LENGTH];
} source_data;

// 从文本中读入原始数据点
void LoadSourceData(struct sourceData *source_data, const char *doc)
{
    // 从文本中读如浮点数数据，分别存放入x,y,z数组中
    FILE *fp;
    fp = fopen(doc, "r");
    if (!fp)
        cout << "No such file!" << endl;
    char str[100];
    fscanf(fp, "%s", str);
    fscanf(fp, "%d", &source_data->length);
    for (int i = 0; i < source_data->length; i++)
    {
        fscanf(fp, "%lf", &source_data->x[i]);
        fscanf(fp, "%lf", &source_data->y[i]);
        fscanf(fp, "%lf", &source_data->z[i]);
    }
    fclose(fp);
}

// 预测数据点结构体
struct predictData
{
    int length;
    double x[MAX_LENGTH][3];
} predict_data;

// 从文本中读入预测数据点
void LoadPredictData(struct predictData *predict_data, const char *doc)
{
    // 从文本中读如浮点数数据，分别存放入x,y,z数组中
    FILE *fp;
    fp = fopen(doc, "r");
    if (!fp)
        cout << "No such file!" << endl;
    char str[100];
    fscanf(fp, "%s", str);
    fscanf(fp, "%d", &predict_data->length);
    for (int i = 0; i < predict_data->length; i++)
    {
        fscanf(fp, "%lf", &predict_data->x[i][0]);
        fscanf(fp, "%lf", &predict_data->x[i][1]);
    }
    fclose(fp);
}

// 可选参数结构体
struct optionalParameter
{
    /* data */
    int mode_choice = 2;
    double mode_parameter = 5;
    int section_choice = 2;
    int section_parameter = 10;
    int data_scale = 100;
    int norm = 1;
} optional_parameter;

// 设置可选参数
void SetOptionalParameter(struct optionalParameter *optional_parameter, int mode_choice, double mode_parameter, int section_choice, int section_parameter, int data_scale, int norm)
{
    optional_parameter->mode_choice = mode_choice;
    optional_parameter->mode_parameter = mode_parameter;
    optional_parameter->section_choice = section_choice;
    optional_parameter->section_parameter = section_parameter;
    optional_parameter->data_scale = data_scale;
    optional_parameter->norm = norm;
}

// 计算各点坐标之间的距离
double **Distance(double x[], double y[], int length)
{
    // int length = sizeof(x) / sizeof(x[0]);
    // int length = 4;
    double **p = InitArrPointer(length, length);
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            double temp_1 = x[i] - x[j];
            double temp_2 = y[i] - y[j];
            double temp = temp_1 * temp_1 + temp_2 * temp_2;
            p[i][j] = sqrt(temp);
        }
    }
    // FILE *fp;
    // fp = fopen("distance.txt", "w");
    // for (int i = 0; i < length; i++)
    // {
    // for (int j = 0; j < length; j++)
    // {
    // fprintf(fp, "%lf\t", p[i][j]);
    // }
    // fprintf(fp, "\n");
    // }
    return p;
}

// 计算半差方
double **Omega(double z[], int length)
{
    double **p = InitArrPointer(length, length);
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            double temp_1 = z[i] - z[j];
            double temp = temp_1 * temp_1 / 2;
            p[i][j] = temp;
        }
    }
    return p;
}

// 将距离数组和半差方数组进行归一化处理，处理成(n*n）*2的二维数组
double **Integration(double **array_1, double **array_2, int length)
{
    double **integration = InitArrPointer(length * length, 2);
    int k = 0;
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            integration[k][0] = array_1[i][j];
            integration[k][1] = array_2[i][j];
            k++;
        }
    }
    return integration;
}

// 将归一化得到的二维n*n数组进行插入排序
double **InsertSort(double **array, int length)
{
    for (int i = 1; i < length; i++)
    {
        double temp_1 = array[i][0];
        double temp_2 = array[i][1];
        int insert = i - 1;
        while (insert >= 0 && temp_1 < array[insert][0])
        {
            array[insert + 1][0] = array[insert][0];
            array[insert + 1][1] = array[insert][1];
            insert--;
        }
        array[insert + 1][0] = temp_1;
        array[insert + 1][1] = temp_2;
    }
    return array;
}

// 分组求平均
double **Average(double **array, int length, optionalParameter optional_parameter)
{
    //
    if (optional_parameter.section_choice == 1)
    {
        int section_num = ceil(length * length / optional_parameter.section_parameter);
        double **result = InitArrPointer(section_num, 2);
        // // 第一组全0,舍去
        // for (int i = 0; i < length; i++)
        // {
        //     double sum_distance = 0;
        //     double sum_omega = 0;
        //     for (int j = i * optional_parameter.section_parameter; j < (i + 1) * optional_parameter.section_parameter; j++)
        //     {
        //         sum_distance += array[j][0];
        //         sum_omega += array[j][1];
        //     }

        //     result[i][0] = sum_distance / length;
        //     result[i][1] = sum_omega / length;
        // }
        int index = 0;
        double sum_distance;
        double sum_omega;
        for (int i = 0; i < section_num; i++)
        {
            sum_distance = 0;
            sum_omega = 0;
            int count = 0;
            for (int j = index; j < length * length; j++)
            {
                sum_distance += array[j][0];
                sum_omega += array[j][1];
                count++;
                if (count >= optional_parameter.section_parameter)
                {
                    count = 0;
                    index = j + 1;
                    result[i][0] = sum_distance / optional_parameter.section_parameter;
                    result[i][1] = sum_omega / optional_parameter.section_parameter;
                    break;
                }
            }
        }
        result[section_num - 1][0] = sum_distance / optional_parameter.section_parameter;
        result[section_num - 1][1] = sum_omega / optional_parameter.section_parameter;
        FILE *fp;
        fp = fopen("average_1.txt", "w");
        if (!fp)
            cout << "No such file!" << endl;
        for (int i = 0; i < section_num; i++)
        {
            fprintf(fp, "%lf\t%lf\n", result[i][0], result[i][1]);
        }
        fclose(fp);
        return result;
    }
    else if (optional_parameter.section_choice == 2)
    {
        int every_distance = MAX_DISTANCE / optional_parameter.section_parameter;
        double **result = InitArrPointer(optional_parameter.section_parameter, 2);
        int index = 0;
        int num = 0;
        double sum_distance;
        double sum_omega;
        // 等间隔将array分成section_parameter组，间隔距离为every_distance，结果存入result

        for (int i = 0; i < optional_parameter.section_parameter; i++)
        {
            num = 0;
            sum_distance = 0;
            sum_omega = 0;
            for (int j = index; j < length * length; j++)
            {
                if (array[j][0] <= (i + 1) * every_distance)
                {
                    sum_distance += array[j][0];
                    sum_omega += array[j][1];
                    num++;
                }
                if (array[j][0] > (i + 1) * every_distance)
                {
                    index = j;
                    result[i][0] = sum_distance / num;
                    result[i][1] = sum_omega / num;
                    break;
                }
            }
        }
        result[optional_parameter.section_parameter - 1][0] = sum_distance / num;
        result[optional_parameter.section_parameter - 1][1] = sum_omega / num;
        FILE *fp;
        fp = fopen("average_2.txt", "w");
        if (!fp)
            cout << "No such file!" << endl;
        for (int i = 0; i < optional_parameter.section_parameter; i++)
        {
            fprintf(fp, "%lf\t%lf\n", result[i][0], result[i][1]);
        }
        fclose(fp);
        return result;
        return result;
    }
    // // 将数据输出到文本，便于比对验证
    // FILE *fp;
    // fp = fopen("average.txt", "w");
    // if (!fp)
    //     cout << "No such file!" << endl;
    // for (int i = 0; i < length; i++)
    // {
    //     fprintf(fp, "%lf\t%lf\n", result[i][0], result[i][1]);
    // }
    // fclose(fp);
    // return result;
    return NULL;
}

// 计算代求点的临近点，返回二维数组，一个存放距离，一个存放索引
double **NextToPosition(double *x, double *y, double *position, int length, optionalParameter optional_parameter)
{
    // choice=1,表示通过距离求临近点信息
    if (optional_parameter.mode_choice == 1)
    {
        double position_x = position[0];
        double position_y = position[1];
        int count = 0;
        int k = 0;
        double **nextTo = InitArrPointer(length + 1, 2, -1);
        for (int i = 0; i < length; i++)
        {
            double distance = sqrt((position_x - x[i]) * (position_x - x[i]) + (position_y - y[i]) * (position_y - y[i]));
            if (distance < optional_parameter.mode_parameter)
            {
                nextTo[k + 1][0] = distance;
                nextTo[k + 1][1] = i;
                k++;
                count++;
            }
        }
        // 临近点个数
        nextTo[0][0] = count;
        return nextTo;
    }
    // choice=2,表示通过点的个数求临近点信息
    else if (optional_parameter.mode_choice == 2)
    {
        int count = (int)optional_parameter.mode_parameter;
        double position_x = position[0];
        double position_y = position[1];
        // 代求点和各个已知点之间的距离，和索引
        double **toBeSorted = InitArrPointer(length, 2);
        for (int i = 0; i < length; i++)
        {
            double distance = sqrt((position_x - x[i]) * (position_x - x[i]) + (position_y - y[i]) * (position_y - y[i]));
            toBeSorted[i][0] = distance;
            toBeSorted[i][1] = i;
        }
        // 通过距离作为标准，进行排序
        double **sorted = InsertSort(toBeSorted, length);
        double **nextTo = InitArrPointer(count + 1, 2);
        // 选取距离最近的count个点
        for (int i = 1; i <= count; i++)
        {
            nextTo[i][0] = sorted[i - 1][0];
            nextTo[i][1] = sorted[i - 1][1];
        }
        nextTo[0][0] = count;
        return nextTo;
    }
    return NULL;
}

// 计算临近点的预测半方差，通过拟合的高斯模型参数，代入求解
double **CalOmega(double **parameter, double **nextTo)
{
    int count = nextTo[0][0];
    int k = 1;
    double **omega = InitArrPointer(count, 1);
    for (int i = 1; i <= count; i++)
    {
        omega[i - 1][0] = parameter[0][0] + parameter[1][0] * (1 - exp(-1 * nextTo[i][0] * nextTo[i][0] / parameter[2][0] / parameter[2][0]));
        if (omega[i - 1][0] < 0)
            omega[i - 1][0] = 0;
    }
    return omega;
}

// 计算权重
double **CalWeight(double *x, double *y, double **nextTo, double **omega_0, int count, double **parameter)
{
    // double position[count];
    // for (int i = 0; i < count; i++)
    // {
    //     position[i] = z[(int)nextTo[i + 1][1]];
    // }
    double position_x[count];
    for (int i = 0; i < count; i++)
    {
        position_x[i] = x[(int)nextTo[i + 1][1]];
    }
    double position_y[count];
    for (int i = 0; i < count; i++)
    {
        position_y[i] = y[(int)nextTo[i + 1][1]];
    }
    double **omega_1 = InitArrPointer(count, count);
    // omega_1 = Omega(position, count);
    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
            double dis_x = (position_x[i] - position_x[j]) * (position_x[i] - position_x[j]);
            double dis_y = (position_y[i] - position_y[j]) * (position_y[i] - position_y[j]);
            double dis = sqrt(dis_x + dis_y);
            omega_1[i][j] = parameter[0][0] + parameter[1][0] * (1 - exp(-1 * dis * dis / parameter[2][0] / parameter[2][0]));
        }
    }
    MatrixXd matrix_0;
    matrix_0.resize(count + 1, count + 1);
    for (int i = 0; i < count + 1; i++)
    {
        for (int j = 0; j < count + 1; j++)
        {
            matrix_0(i, j) = 1;
        }
    }
    matrix_0(count, count) = 0;
    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
            matrix_0(i, j) = omega_1[i][j];
        }
    }
    for (int i = 0; i < count + 1; i++)
    {
        matrix_0(i, i) += 1e-10;
    }
    cout << "Matrix_0" << endl
         << matrix_0 << endl;
    matrix_0 = matrix_0.inverse();
    cout << "Matrix_0.inverse()" << endl
         << matrix_0 << endl;
    MatrixXd matrix_1;
    matrix_1.resize(count + 1, 1);
    for (int i = 0; i < count; i++)
    {
        matrix_1(i, 0) = omega_0[i][0];
    }
    matrix_1(count, 0) = 1;
    cout << "Matrix_1" << endl
         << matrix_1 << endl;
    MatrixXd weight;
    weight.resize(count + 1, 1);
    weight = matrix_0 * matrix_1;
    cout << "Weight:" << endl
         << weight << endl;
    double **result = InitArrPointer(count + 1, 1);
    double add = 0;
    for (int i = 0; i < count; i++)
    {
        if (weight(i, 0) < 0)
            weight(i, 0) = 0;
    }
    for (int i = 0; i < count; i++)
    {
        add += weight(i, 0);
    }
    for (int i = 0; i < count; i++)
    {
        weight(i, 0) /= add;
    }
    for (int i = 0; i < count + 1; i++)
    {
        result[i][0] = weight(i, 0);
        // add += weight(i, 0);
    }
    cout << "sub of weight: " << add << endl;
    return result;
}

// 代入临近点和权重求得最终插值结果
double CalResult(double *z, double **nextTo, double **weight, int count)
{
    double result = 0.0;
    for (int i = 0; i < count; i++)
    {
        cout << "z[i]:" << z[(int)nextTo[i + 1][1]] << "\tweight[i]:" << weight[i][0] << endl;
        result += z[(int)nextTo[i + 1][1]] * weight[i][0];
    }
    return result;
}

// 求解克里金插值结果，根据参数选择进行不同的计算过程
/*
    1、optional_parameter.mode_choice：选择临近点的方式，1选择距离，2选择点数
    2、optional_parameter.mode_parameter：确定准确的距离和点数
    3、optional_parameter.norm：是否进行正态化
    4、optional_parameter.data_scale：数据的尺度
    5、optional_parameter.section_choice：选择分段的方式，1选择距离，2选择点数
    6、optional_parameter.section_parameter：确定准确的距离和点数
*/
void Kriging(sourceData *source_data, predictData *predict_data, optionalParameter optional_parameter)
{
    int length = source_data->length;
    double save_z[length];
    // 是否正态化
    if (optional_parameter.norm == 1)
    {
        for (int i = 0; i < length; i++)
        {
            save_z[i] = source_data->z[i];
            source_data->z[i] = (pow(source_data->z[i], 3.5) - 1) / 3.5 / pow(10, 10);
        }
    }
    else if (optional_parameter.norm == 0)
    {
        for (int i = 0; i < length; i++)
        {
            save_z[i] = source_data->z[i];
        }
    }
    // 求解各点之间的距离
    double **distance = InitArrPointer(length, length);
    distance = Distance(source_data->x, source_data->y, length);
    // 求解各点之间的半方差
    double **omega_1 = InitArrPointer(length, length);
    omega_1 = Omega(source_data->z, length);
    for (int i = 0; i < 89; i++)
    {
        for (int j = 0; j < 89; j++)
        {
            cout << omega_1[i][j] << " ";
        }
        cout << endl;
    }
    // 排序行化
    double **dis_omega = InitArrPointer(length * length, 2);
    dis_omega = Integration(distance, omega_1, length);
    dis_omega = InsertSort(dis_omega, length * length);
    // 分组求平均
    double **average = InitArrPointer(length - 1, 2);
    average = Average(dis_omega, length, optional_parameter);
    // 拟合高斯曲线，求得拟合曲线参数
    double a = 50.0, b = 100.0, c = 10000.0; // 初值
    LevenbergMarquardt lm(&a, &b, &c);
    lm.setParameters(1e-10, 1e-10, 1000, true);
    // two choices(length/para)
    for (int i = 0; i < (int)optional_parameter.section_parameter * optional_parameter.data_scale * 0.01; i++)
    {
        lm.addObservation(average[i][0], average[i][1]);
    }
    double **parameter = InitArrPointer(3, 1);
    parameter = lm.solve();
    parameter[0][0] = 0;
    cout << "Result:" << parameter[0][0] << " + " << parameter[1][0]
         << " * (1 - exp(-x^2 / " << parameter[2][0] << "^2))" << endl;

    // 循环计算待求点的插值结果，并存入文本中
    FILE *fp = fopen("out.txt", "w");
    fprintf(fp, "predict_result\n");
    for (int times = 0; times < predict_data->length; times++)
    {
        // 计算待求点的邻接点的索引，和距离
        double **nextTo = InitArrPointer(source_data->length + 1, 2);
        nextTo = NextToPosition(source_data->x, source_data->y, predict_data->x[times], length, optional_parameter);
        int count = nextTo[0][0];
        cout << "count" << count << endl;
        for (int i = 1; i <= count; i++)
        {
            cout << nextTo[i][0] << "/" << nextTo[i][1] << endl;
        }
        // 带入邻近距离求得邻近半方差
        double **omega_0 = InitArrPointer(count, 1);
        omega_0 = CalOmega(parameter, nextTo);
        for (int i = 0; i < count; i++)
        {
            cout << omega_0[i][0] << "~~~~" << endl;
        }
        // 计算权重
        double **weight = InitArrPointer(count + 1, 1);
        weight = CalWeight(source_data->x, source_data->y, nextTo, omega_0, count, parameter);
        // for (int i = 0; i < count; i++)
        // {
        //     weight[i][0] /= weight[count][0];
        //     cout
        //         << "weight/lamuda:\t" << weight[i][0] << endl;
        // }
        // 计算插值结果
        double result = CalResult(source_data->z, nextTo, weight, count);
        if (optional_parameter.norm == 1)
            result = pow(result * pow(10, 10) * 3.5 + 1, (double)1 / 3.5);
        fprintf(fp, "%lf\n", result);
    }
    fclose(fp);
    // 还原z数组
    for (int i = 0; i < length; i++)
    {
        source_data->z[i] = save_z[i];
    }
    // return result;
}