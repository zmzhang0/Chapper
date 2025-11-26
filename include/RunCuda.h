#pragma once
#include <cuda_runtime.h>
#include <cufft.h>
#include <device_launch_parameters.h>
#include <cstdio>
#include <iostream>
#include <vector>


#if 1   // 3D

// 定义体素网格的数量
#define NX 128
#define NY 128
#define NZ 128
#define TX NX * NY * NZ

#define ContainerSize 128// 容器大小

#define CutterLength 90 // 刀具长度 
#define CutterWidth 4 //刀具宽度

#define VOXELSIZE 1 // 体素大小

//定义blockSize
#define BSX 8 
#define BSY 8
#define BSZ 8

#define OriNum 1 //物体方向数量

#define CutterNum 6 //刀具方向数量

#define SweepNum 4 // 扫除方向个数

#else     // 2D


// 定义体素网格的数量
#define NX 128
#define NY 128
#define NZ 1
#define TS NX * NY * NZ

#define ContainerSize 256  // 容器大小

//定义blockSize
#define BSX 8
#define BSY 8
#define BSZ 1

#define OriNum 4 //方向采样数量
#define CutterNum 8 //物体数量

#define OutPut2DResult false //是否输出中间结果

#endif



