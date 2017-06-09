#include <opencv2/core/core.hpp>//这两个头文件时支持opencv版本的头文件，包含数据结构，矩阵运算，数据变换内存管理数学文本等
#include <opencv2/highgui/highgui.hpp>//包含图形界面和视频图像处理的头文件
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>//文件写操作，内存写入存储设备
#include <assert.h>//表达式结果正确性测试并可使程序终止
#include <time.h>
#include <cmath>//数学函数

using namespace cv;
using namespace std;

//default filter and iteration number
int ItNum = 10;
int Type = 2;

#include "DM.h"

int main(int argc, char** argv)
{
    
    DM DualMesh;
    if (argc!=4)
    {
       cout<<"usage: main filename filterType Iterations.\n For example: ./cf lena.bmp m 30, where m means Mean Curvature Filter. \n";
       return -1;
    }
    DualMesh.read(argv[1]);//读取图像帧
    ItNum = atoi(argv[3]);//把字符串转换成整型

    char * filterType = argv[2];//区分滤波器种类
    if (*filterType == 't') Type = 0;
    if (*filterType == 'm') Type = 1;
    if (*filterType == 'd') Type = 3;


    DualMesh.split();//分割函数，本函数将原图像分割成四类像素点，白园黑圆白三角黑三角，并将原图像缩小一半
    double mytime;
    DualMesh.Filter(Type, mytime, ItNum);//根据输入参数选择滤波器种类，并将图像分割为四类像素点
    cout<<"runtime is "<<mytime<<" milliseconds."<<endl;//显示算法执行时间

    DualMesh.merge();//将四类像素值组合为一幅完整的图像
    DualMesh.write();//将Mat类型的矩阵保存到图像

    DualMesh.read(argv[1]);
    DualMesh.FilterNoSplit(Type, mytime, ItNum);//不对图像进行分割，整体采用曲率滤波
    cout<<"runtime (noSplit) is "<<mytime<<" milliseconds."<<endl;
    DualMesh.write("CF_NoSplit_result.png");
    
    return 0;
}


