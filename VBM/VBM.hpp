#ifndef VBM_HPP
#define VBM_HPP

#include <iostream>
#include "VBM/solver.hpp"

namespace vbm
{

    // 计算距离场
    void DistanceFunction(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, std::vector<float>& unsignedDis);

    void DistanceFunction(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, std::vector<int>& unsignedDis);

    // 计算梯度方向敏感的距离场
    void DistanceFunction(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, 
        const std::vector<float>& gradientX, const std::vector<float>& gradientY, const std::vector<float>& gradientZ, std::vector<float>& unsignedDis);


    // 计算sdf同时得到source点信息
    void DistanceFunctionWithSources(const size_t nx, const size_t ny, const size_t nz,
        const std::vector<bool>& visibilityField, std::vector<int>& unsignedDis, std::vector<int>& sources, std::vector<std::vector<int>>& sourcePoints);

    // 计算internal visaul hull
    void computeIVH(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, std::vector<int>& internalCH);


    // 将2二维矩阵或者三维矩阵的第一层转换为图像
    void Vector2Image(const size_t nx, const size_t ny, const size_t nz, const size_t ind, const std::vector<int>& gScore, const std::string& outputPath);


} // namespace vbm

#endif // VBM_HPP