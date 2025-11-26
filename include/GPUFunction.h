#pragma once
#include "Runcuda.h"


__global__ void MultiplyInFrequencyDomain(cufftDoubleComplex* A, cufftDoubleComplex* B, cufftDoubleComplex* C);

__global__ void MultiplyInFrequencyDomain(cufftDoubleComplex* A, cufftDoubleComplex* B, cufftDoubleComplex* C, const int nx, const int ny, const int nz);



void MultiplyInFrequencyDomainWithLoop(cufftDoubleComplex* A, cufftDoubleComplex* B, cufftDoubleComplex* C);

inline void Bool2Complex(const std::vector<bool>& h_, cufftDoubleComplex* h_complex);
inline void Complex2Bool(const cufftDoubleComplex* h_complex, std::vector<bool>& h_);

inline void Bool2DoubleRealWithFlip(const std::vector<bool>& h_, cufftDoubleReal* h_complex);
inline void Bool2DoubleRealWithFlip(const std::vector<bool>& h_, cufftDoubleReal* h_complex, const int nx, const int ny, const int nz);


inline void Bool2DoubleReal(const std::vector<bool>& h_, cufftDoubleReal* h_complex);
inline void Bool2DoubleReal(const std::vector<bool>& h_, cufftDoubleReal* h_complex, const int nx, const int ny, const int nz);

inline void DoubleReal2Bool(const cufftDoubleReal* h_complex, std::vector<bool>& h_);
inline void DoubleReal2Bool(const cufftDoubleReal* h_complex, std::vector<bool>& h_, const int nx, const int ny, const int nz);


inline void DoubleReal2BoolDirectly(const cufftDoubleReal* h_complex, std::vector<bool>& h_);


inline void Int2DoubleReal(const std::vector<int>& h_, cufftDoubleReal* h_complex);
inline void DoubleReal2Int(const cufftDoubleReal* h_complex, std::vector<int>& h_);

inline void DoubleReal2IntDirectly(const cufftDoubleReal* h_complex, std::vector<int>& h_);

void convolution(const std::vector<bool>& h_A, const std::vector<bool>& h_B, std::vector<bool>& h_C);


__global__ void ConvolutionPointwise(const cufftDoubleReal* A, cufftDoubleReal* B, cufftDoubleReal* C, int* box);

void ConvolutionWithoutFFT(const std::vector<int>& h_A, const std::vector<bool>& h_B, std::vector<int>& h_C);


void ConvolutionOnReal(const std::vector<bool>& h_A, std::vector<bool>& h_B, std::vector<bool>& h_C);

void ConvolutionOnReal(const std::vector<int>& h_A, const std::vector<bool>& h_B, std::vector<int>& h_C);

void ConvolutionOnReal(const std::vector<bool>& h_A, const std::vector<bool>& h_B, std::vector<bool>& h_C, const int nx, const int ny, const int nz);

void ConvolutionOnRealMulti(std::vector<cufftDoubleReal>& h_A, 
    std::vector<std::vector<cufftDoubleReal>>& h_B, 
    std::vector<std::vector<cufftDoubleReal>>& h_C, 
    const int nx, const int ny, const int nz);




__global__ void ComputeVisibility(const int* Obs, int* Vis);

void VisibilityAnalysisByGPU(const std::vector<bool>& h_C, std::vector<int>& visibility);



__global__ void computeKernel(
    int* Box, int* kernel,
    int loopNum, int kx, int ky, int kz,
    int offsetx, int offsety, int offsetz,
    int kykz, int newBoxYZ, int newBoxZ);

void ExtractKernel(std::vector<std::vector<bool>>& OffsetBox,
    std::vector<int>& kernel,
    std::vector<int>& newBox,
    std::vector<int>& offset,
    int kx, int ky, int kz,
    int LoopNum);


