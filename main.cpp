#include "Kriging.h"
int main(int argc, char **argv)
{
    char *source_data_file = "source_data.txt";
    char *predict_data_file = "predict_data.txt";
    // 命令汉参数传递文件名称
    // source_data_file = argv[1];
    // predict_data_file = argv[2];
    // 加载原始数据点
    LoadSourceData(&source_data, source_data_file);
    // 加载待预测数据点
    LoadPredictData(&predict_data, predict_data_file);
    // 设置可选参数
    SetOptionalParameter(&optional_parameter, 2, 5, 2, 10, 100, 1);
    // 调用Kriging函数，将插值结果保存在out.txt中
    Kriging(&source_data, &predict_data, optional_parameter);
    return 0;
}
