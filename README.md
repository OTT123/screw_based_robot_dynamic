# 基于旋量理论的固定基座开链机械臂动力学解算
基于《现代机器人学》编写。实现了从urdf文件解析旋量理论参数的解析器以及迭代动力学算法。
该工程仅能实现urdf的解析，以及动力学参数的数值解算。不具备其他任何功能。

# 编译方法
## Windows
将该工程下载到工作空间文件夹，在工作空间文件夹打开shell，（假设你的eigen路径为DIR）

`g++ symrobotic.cpp pugixml.cpp my_dynamic.cpp main.cpp -o main -IDIR`

运行

`./main`

# screw_based_robot_dynamic
## 1. main.cpp
提供了一些测试方法。给出了如何计算样例机器人动力学参数的例子。

## 2. my_dynamic.h
动力学算法头文件，注释都写在这个文件

## 3. my_dynamic.cpp
动力学算法

## 4.pugixml.cpp\pugixml.hpp\pugiconfig.hpp
解析xml语言的轻量级库

## 5.robot_basic_params.h
留一个接口来设置机器人的一些基本参数，目前主要是设置机器人自由度DOF，默认是6

## 6.symrobotic.cpp\symrobotic.h
计划的符号运算功能，还未实现，请忽略。

# 使用方法
在main.cpp中可以看到，可以切换一些机器人型号，设置关节角度，角速度，角加速度等。

如果想要新增机器人，必须获得该机器人的urdf文件，并按照模板的格式设置urdf解析器的参数。
