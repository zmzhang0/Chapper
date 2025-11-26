#pragma once
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define WINDOWS_LEAN_AND_MEAN // Please, not too much windows shenanigans
#endif

// Standard libs
#include <string>
#include <cstdio>
#include "RunCuda.h"

// Trimesh for model importing
#include "TriMesh.h"
// Util
#include "util.h"
#include "util_io.h"
#include "util_cuda.h"
#include "timer.h"
// CPU voxelizer fallback
#include "cpu_voxelizer.h"


// CUDA voxelizer
// 体素化函数，输入mesh、x/y/z网格尺寸，输出体素bool三维数组
void cuda_voxelizeMesh(
	trimesh::TriMesh* mesh,
	const size_t gridX,
	const size_t gridY,
	const size_t gridZ,
	const double voxel_size,
	std::vector<bool>& voxelGrid);
