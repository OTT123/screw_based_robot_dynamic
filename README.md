# 编译方法
## Windows
将该工程下载到工作空间文件夹，假设你的eigen路径为DIR，打开shell

`g++ symrobotic.cpp pugixml.cpp my_dynamic.cpp main.cpp -o main -IDIR`

运行

`./main`

# screw_based_robot_dynamic
## 1. main.cpp
提供了一些测试方法

## 2. my_dynamic.h
动力学算法头文件，注释都写在这个文件

## 3. my_dynamic.cpp
动力学算法

## 4.pugixml.cpp\pugixml.hpp\pugiconfig.hpp
解析xml语言的轻量级库

## 5.robot_basic_params.h
留一个接口来设置机器人的一些基本参数，目前主要是设置机器人自由度DOF，默认是6

## 6.symrobotic.cpp\symrobotic.h
测试文件，忽略

# 使用方法
在main.cpp中可以看到
