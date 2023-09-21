# Kriging_project 说明文档

# source_data.txt格式说明
line 1：字符串=source_data  数据点个数=n
line 2:数据点x坐标 数据点y坐标 数据点属性值
line 3:数据点x坐标 数据点y坐标 数据点属性值
......
line n+1:数据点x坐标 数据点y坐标 数据点属性值


# predict_data.txt格式说明

line 1：字符串=predict_data  数据点个数=n
line 2:待预测数据点x坐标 待预测数据点y坐标 
line 3:待预测数据点x坐标 待预测数据点y坐标
......
line n+1:待预测数据点x坐标 待预测数据点y坐标

项目启动方式：运行main.cpp文件即可，保证上述两个txt文件和main.cpp在同一目录下。
