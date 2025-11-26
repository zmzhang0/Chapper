#include "GPUFunction.h"


__global__ void MultiplyInFrequencyDomain(cufftDoubleComplex* A, cufftDoubleComplex* B, cufftDoubleComplex* C) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;


    if (x < NX && y < NY && z < NZ) {
        int idx = x * NY * NZ + y * NZ + z;
        cufftDoubleComplex a = A[idx];
        cufftDoubleComplex b = B[idx];

        C[idx].x = a.x * b.x - a.y * b.y;
        C[idx].y = a.x * b.y + a.y * b.x;
    }
}

__global__ void MultiplyInFrequencyDomain(cufftDoubleComplex* A, cufftDoubleComplex* B, cufftDoubleComplex* C,const int nx, const int ny, const int nz) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;


    if (x < nx && y < ny && z < nz) {
        int idx = x * ny * nz + y * nz + z;
        cufftDoubleComplex a = A[idx];
        cufftDoubleComplex b = B[idx];

        C[idx].x = a.x * b.x - a.y * b.y;
        C[idx].y = a.x * b.y + a.y * b.x;
    }
}


void MultiplyInFrequencyDomainWithLoop(cufftDoubleComplex* A, cufftDoubleComplex* B, cufftDoubleComplex* C) {

    for (int i = 0; i < TX; i++) {
        cufftDoubleComplex a = A[i];
        cufftDoubleComplex b = B[i];

        C[i].x = a.x * b.x - a.y * b.y;
        C[i].y = a.x * b.y + a.y * b.x;
    }
}


inline void Bool2ComplexWithFlip(const std::vector<bool>& h_, cufftDoubleComplex* h_complex) {

    for (int i = 0; i < TX; ++i) {
        h_complex[TX - i - 1].y = 0.0f;
        if (h_[i])
            h_complex[TX - i - 1].x = 1.0f;
        else
            h_complex[TX - i - 1].x = 0.0f;
    }

}

inline void Bool2Complex(const std::vector<bool>& h_, cufftDoubleComplex* h_complex) {

    for (int i = 0; i < TX; ++i) {
        h_complex[i].y = 0.0f;
        if (h_[i])
            h_complex[i].x = 1.0f;
        else
            h_complex[i].x = 0.0f;
    }

}




inline void Complex2Bool(const cufftDoubleComplex* h_complex, std::vector<bool>& h_) {

	int index = 0;
    for (int i = 0; i < TX; ++i)
    {
        if (h_complex[i].x + 0.5 > 1.0f)
        {
            index++;
            h_[i] = true;
        }
    }
	printf("Voxel Num: %d\n", index);
}



void convolution(const std::vector<bool>& h_A, const std::vector<bool>& h_B, std::vector<bool>& h_C) {

    cufftDoubleComplex* d_A, * d_B, * d_C;
    cudaMalloc((void**)&d_A, sizeof(cufftDoubleComplex) * NX * NY * NZ);
    cudaMalloc((void**)&d_B, sizeof(cufftDoubleComplex) * NX * NY * NZ);
    cudaMalloc((void**)&d_C, sizeof(cufftDoubleComplex) * NX * NY * NZ);

    cufftDoubleComplex* h_A_complex = new cufftDoubleComplex[NX * NY * NZ]();
    cufftDoubleComplex* h_B_complex = new cufftDoubleComplex[NX * NY * NZ]();
    cufftDoubleComplex* h_C_complex = new cufftDoubleComplex[NX * NY * NZ]();

    Bool2Complex(h_A, h_A_complex);
    Bool2ComplexWithFlip(h_B, h_B_complex);

    cudaMemcpy(d_A, h_A_complex, sizeof(cufftDoubleComplex) * NX * NY * NZ, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B_complex, sizeof(cufftDoubleComplex) * NX * NY * NZ, cudaMemcpyHostToDevice);

    cufftHandle planA, planB, planC;
    cufftPlan3d(&planA, NX, NY, NZ, CUFFT_Z2Z);
    cufftPlan3d(&planB, NX, NY, NZ, CUFFT_Z2Z);
    cufftPlan3d(&planC, NX, NY, NZ, CUFFT_Z2Z);

    cufftExecZ2Z(planA, d_A, d_A, CUFFT_FORWARD);
    cufftExecZ2Z(planB, d_B, d_B, CUFFT_FORWARD);

    dim3 blockSize(BSX, BSY, BSZ);   
    dim3 gridSize((NX + blockSize.x - 1) / blockSize.x,
        (NY + blockSize.y - 1) / blockSize.y,
        (NZ + blockSize.z - 1) / blockSize.z); 

    MultiplyInFrequencyDomain <<< gridSize, blockSize >>> (d_A, d_B, d_C); 
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }

    cufftExecZ2Z(planC, d_C, d_C, CUFFT_INVERSE);

    cudaMemcpy(h_C_complex, d_C, sizeof(cufftDoubleComplex) * NX * NY * NZ, cudaMemcpyDeviceToHost);
    
   Complex2Bool(h_C_complex, h_C);

   cufftDestroy(planA);
   cufftDestroy(planB);
   cufftDestroy(planC);
   cudaFree(d_A);
   cudaFree(d_B);
   cudaFree(d_C);
   delete[] h_A_complex;
   delete[] h_B_complex;
   delete[] h_C_complex;
}


inline void Bool2DoubleRealWithFlip(const std::vector<bool>& h_, cufftDoubleReal* h_complex) {

    for (int i = 0; i < TX; ++i) {
        if (h_[i])
            h_complex[TX - i - 1] = 1.0f;
        else
            h_complex[TX - i - 1] = 0.0f;
    }

}



inline void Bool2DoubleReal(const std::vector<bool>& h_, cufftDoubleReal* h_complex) {

    for (int i = 0; i < TX; ++i) {
        if (h_[i])
            h_complex[i] = 1.0f;
        else
            h_complex[i] = 0.0f;
    }

}



inline void DoubleReal2Bool(const cufftDoubleReal* h_complex, std::vector<bool>& h_) {

    for(int i = 0; i < NX; ++i)
        for(int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k)
            {
                int idx = ((i + NX - 1) % NX) * NY * NZ + ((j + NY - 1) % NY) * NZ + ((k + NZ - 1) % NZ);
                int idx_ = i * NY * NZ + j * NZ + k;
                if (h_complex[idx] + 0.5 > 1.0f)
                    h_[idx_] = true;
                else
                    h_[idx_] = false;
            }

}



inline void DoubleReal2BoolDirectly(const cufftDoubleReal* h_complex, std::vector<bool>& h_) {

    for (int i = 0; i < TX; ++i)
    {
        if (h_complex[i] > 0.5f)
            h_[i] = true;
        else
            h_[i] = false;
    }
}


__global__ void ConvolutionPointwise(const cufftDoubleReal* A, cufftDoubleReal* B, cufftDoubleReal* C, int* box) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;


    if (x < NX && y < NY && z < NZ) {
        int idx_c = x * NY * NZ + y * NZ + z;
        C[idx_c] = 0.0f;

        for (int i = 0; i < box[0] + 1; i++)
            for (int j = 0; j < box[1] + 1; j++)
                for (int k = 0; k < box[2] + 1; k++)
                {
                    int idx_b = i * NY * NZ + j * NZ + k;
                    int idx_a = (i + x) % NX * NY * NZ + (j + y) % NY * NZ + (k + z) % NZ;
                    C[idx_c] += A[idx_a] * B[idx_b];
                }
    }
}

void ConvolutionWithoutFFT(const std::vector<bool>& h_A, const std::vector<bool>& h_B, std::vector<bool>& h_C)
{
    int c_box[3] = { 0, 0, 0 };
    
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k)
            {
                int idx = i * NY * NZ + j * NZ + k;
                if (h_B[idx])
                {
                    c_box[0] = std::max(c_box[0], i);
                    c_box[1] = std::max(c_box[1], j);
                    c_box[2] = std::max(c_box[2], k);
                }
            }
          
    int* d_box;
    cudaMalloc((void**)&d_box, sizeof(int) * 3);
    cudaMemcpy(d_box, c_box, sizeof(int) * 3, cudaMemcpyHostToDevice);

    cufftDoubleReal* d_A, * d_B, * d_C;
    cudaMalloc((void**)&d_A, sizeof(cufftDoubleReal) * NX * NY * NZ);
    cudaMalloc((void**)&d_B, sizeof(cufftDoubleReal) * NX * NY * NZ);
    cudaMalloc((void**)&d_C, sizeof(cufftDoubleReal) * NX * NY * NZ);


    cufftDoubleReal* h_A_double = new cufftDoubleReal[NX * NY * NZ]();
    cufftDoubleReal* h_B_double = new cufftDoubleReal[NX * NY * NZ]();
    cufftDoubleReal* h_C_double = new cufftDoubleReal[NX * NY * NZ]();

    Bool2DoubleReal(h_A, h_A_double);
    Bool2DoubleReal(h_B, h_B_double);

    cudaMemcpy(d_A, h_A_double, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B_double, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyHostToDevice);


    dim3 blockSize(BSX, BSY, BSZ);   
    dim3 gridSize((NX + blockSize.x - 1) / blockSize.x,
        (NY + blockSize.y - 1) / blockSize.y,
        (NZ + blockSize.z - 1) / blockSize.z); 

    ConvolutionPointwise << < gridSize, blockSize >> > (d_A, d_B, d_C, d_box);
     
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }


    cudaMemcpy(h_C_double, d_C, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyDeviceToHost);

    DoubleReal2BoolDirectly(h_C_double, h_C);
    

    cudaFree(d_box);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    delete[] h_A_double;
    delete[] h_B_double;
    delete[] h_C_double;
}


void ConvolutionOnReal(const std::vector<bool>& h_A, std::vector<bool>& h_B, std::vector<bool>& h_C) {
    cufftDoubleReal* d_A, * d_B, * d_C;
    cudaMalloc((void**)&d_A, sizeof(cufftDoubleReal) * NX * NY * NZ);
    cudaMalloc((void**)&d_B, sizeof(cufftDoubleReal) * NX * NY * NZ);
    cudaMalloc((void**)&d_C, sizeof(cufftDoubleReal) * NX * NY * NZ);

    cufftDoubleComplex* c_A, * c_B, * c_C;
    cudaMalloc((void**)&c_A, sizeof(cufftDoubleComplex) * NX * NY * NZ);
    cudaMalloc((void**)&c_B, sizeof(cufftDoubleComplex) * NX * NY * NZ);
    cudaMalloc((void**)&c_C, sizeof(cufftDoubleComplex) * NX * NY * NZ);


    cufftDoubleReal* h_A_double = new cufftDoubleReal[NX * NY * NZ]();
    cufftDoubleReal* h_B_double = new cufftDoubleReal[NX * NY * NZ]();
    cufftDoubleReal* h_C_double = new cufftDoubleReal[NX * NY * NZ]();

    Bool2DoubleReal(h_A, h_A_double);
    Bool2DoubleRealWithFlip(h_B, h_B_double);


    cudaMemcpy(d_A, h_A_double, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B_double, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyHostToDevice);

    cufftHandle planA, planB, planC;
    cufftPlan3d(&planA, NX, NY, NZ, CUFFT_D2Z);
    cufftPlan3d(&planB, NX, NY, NZ, CUFFT_D2Z);
    cufftPlan3d(&planC, NX, NY, NZ, CUFFT_Z2D);

    cufftExecD2Z(planA, d_A, c_A);
    cufftExecD2Z(planB, d_B, c_B);

    dim3 blockSize(BSX, BSY, BSZ);   
    dim3 gridSize((NX + blockSize.x - 1) / blockSize.x,
        (NY + blockSize.y - 1) / blockSize.y,
        (NZ + blockSize.z - 1) / blockSize.z); 

    MultiplyInFrequencyDomain << < gridSize, blockSize >> > (c_A, c_B, c_C);

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }

    cufftExecZ2D(planC, c_C, d_C);

    cudaMemcpy(h_C_double, d_C, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyDeviceToHost);


    DoubleReal2Bool(h_C_double, h_C);
  

    cufftDestroy(planA);
    cufftDestroy(planB);
    cufftDestroy(planC);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(c_A);
    cudaFree(c_B);
    cudaFree(c_C);

    delete[] h_A_double;
    delete[] h_B_double;
    delete[] h_C_double;
}


inline void Int2DoubleReal(const std::vector<int>& h_, cufftDoubleReal* h_complex) {

    for (int i = 0; i < TX; ++i) {
		double tmp = h_[i];
		h_complex[i] = tmp;
    }

}


inline void DoubleReal2Int(const cufftDoubleReal* h_complex, std::vector<int>& h_) {

    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k)
            {
                int idx = ((i + NX - 1) % NX) * NY * NZ + ((j + NY - 1) % NY) * NZ + ((k + NZ - 1) % NZ);
                int idx_ = i * NY * NZ + j * NZ + k;
                
                double tmp = h_complex[idx];
                double total = TX;
                tmp = tmp / total;
                h_[idx_] = int(tmp);
            }
}

inline void DoubleReal2IntDirectly(const cufftDoubleReal* h_complex, std::vector<int>& h_) {
    
    for (int i = 0; i < TX; ++i)
    {
        double tmp = h_complex[i];
        h_[i] = int(tmp);
    }
}




void ConvolutionOnReal(const std::vector<int>& h_A, const std::vector<bool>& h_B, std::vector<int>& h_C) {
    cufftDoubleReal* d_A, * d_B, * d_C;
    cudaMalloc((void**)&d_A, sizeof(cufftDoubleReal) * NX * NY * NZ);
    cudaMalloc((void**)&d_B, sizeof(cufftDoubleReal) * NX * NY * NZ);
    cudaMalloc((void**)&d_C, sizeof(cufftDoubleReal) * NX * NY * NZ);

    cufftDoubleComplex* c_A, * c_B, * c_C;
    cudaMalloc((void**)&c_A, sizeof(cufftDoubleComplex) * NX * NY * NZ);
    cudaMalloc((void**)&c_B, sizeof(cufftDoubleComplex) * NX * NY * NZ);
    cudaMalloc((void**)&c_C, sizeof(cufftDoubleComplex) * NX * NY * NZ);


    cufftDoubleReal* h_A_double = new cufftDoubleReal[NX * NY * NZ]();
    cufftDoubleReal* h_B_double = new cufftDoubleReal[NX * NY * NZ]();
    cufftDoubleReal* h_C_double = new cufftDoubleReal[NX * NY * NZ]();

    Int2DoubleReal(h_A, h_A_double);
    Bool2DoubleRealWithFlip(h_B, h_B_double);


    cudaMemcpy(d_A, h_A_double, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B_double, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyHostToDevice);

    cufftHandle planA, planB, planC;
    cufftPlan3d(&planA, NX, NY, NZ, CUFFT_D2Z);
    cufftPlan3d(&planB, NX, NY, NZ, CUFFT_D2Z);
    cufftPlan3d(&planC, NX, NY, NZ, CUFFT_Z2D);

    cufftExecD2Z(planA, d_A, c_A);
    cufftExecD2Z(planB, d_B, c_B);

    dim3 blockSize(BSX, BSY, BSZ);   
    dim3 gridSize((NX + blockSize.x - 1) / blockSize.x,
        (NY + blockSize.y - 1) / blockSize.y,
        (NZ + blockSize.z - 1) / blockSize.z); 

    MultiplyInFrequencyDomain << < gridSize, blockSize >> > (c_A, c_B, c_C);

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }

    cufftExecZ2D(planC, c_C, d_C);

    cudaMemcpy(h_C_double, d_C, sizeof(cufftDoubleReal) * NX * NY * NZ, cudaMemcpyDeviceToHost);

    DoubleReal2Int(h_C_double, h_C);

    cufftDestroy(planA);
    cufftDestroy(planB);
    cufftDestroy(planC);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(c_A);
    cudaFree(c_B);
    cudaFree(c_C);

    delete[] h_A_double;
    delete[] h_B_double;
    delete[] h_C_double;
}


inline void Bool2DoubleRealWithFlip(const std::vector<bool>& h_, cufftDoubleReal* h_complex, const int nx, const int ny, const int nz) {

    for (int i = 0; i < nx * ny * nz; ++i) {
        if (h_[i])
            h_complex[nx * ny * nz - i - 1] = 1.0f;
        else
            h_complex[nx * ny * nz - i - 1] = 0.0f;
    }
}

inline void Bool2DoubleReal(const std::vector<bool>& h_, cufftDoubleReal* h_complex, const int nx, const int ny, const int nz) {

    for (int i = 0; i < nx * ny * nz; ++i) {
        if (h_[i])
            h_complex[i] = 1.0f;
        else
            h_complex[i] = 0.0f;
    }

}

inline void DoubleReal2Bool(const cufftDoubleReal* h_complex, std::vector<bool>& h_, const int nx, const int ny, const int nz) {

    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
            {
                int idx = ((i + nx - 1) % nx) * ny * nz + ((j + ny - 1) % ny) * nz + ((k + nz - 1) % nz);
                int idx_ = i * ny * nz + j * nz + k;

                if (h_complex[idx] + 0.5 > 1.0f)
                    h_[idx_] = true;
                else
                    h_[idx_] = false;
            }

}


void ConvolutionOnReal(const std::vector<bool>& h_A, const std::vector<bool>& h_B, std::vector<bool>& h_C, const int nx, const int ny, const int nz)
{
   
    cufftDoubleReal* d_A, * d_B, * d_C;
    cudaMalloc((void**)&d_A, sizeof(cufftDoubleReal) * nx * ny * nz);
    cudaMalloc((void**)&d_B, sizeof(cufftDoubleReal) * nx * ny * nz);
    cudaMalloc((void**)&d_C, sizeof(cufftDoubleReal) * nx * ny * nz);

    cufftDoubleComplex* c_A, * c_B, * c_C;
    cudaMalloc((void**)&c_A, sizeof(cufftDoubleComplex) * nx * ny * nz);
    cudaMalloc((void**)&c_B, sizeof(cufftDoubleComplex) * nx * ny * nz);
    cudaMalloc((void**)&c_C, sizeof(cufftDoubleComplex) * nx * ny * nz);


    cufftDoubleReal* h_A_double = new cufftDoubleReal[nx * ny * nz]();
    cufftDoubleReal* h_B_double = new cufftDoubleReal[nx * ny * nz]();
    cufftDoubleReal* h_C_double = new cufftDoubleReal[nx * ny * nz]();


    Bool2DoubleReal(h_A, h_A_double, nx, ny, nz);
    Bool2DoubleRealWithFlip(h_B, h_B_double, nx, ny, nz);


    cudaMemcpy(d_A, h_A_double, sizeof(cufftDoubleReal) * nx * ny * nz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B_double, sizeof(cufftDoubleReal) * nx * ny * nz, cudaMemcpyHostToDevice);

    
    cufftHandle planA, planB, planC;
    cufftPlan3d(&planA, nx, ny, nz, CUFFT_D2Z);
    cufftPlan3d(&planB, nx, ny, nz, CUFFT_D2Z);
    cufftPlan3d(&planC, nx, ny, nz, CUFFT_Z2D);


    cufftExecD2Z(planA, d_A, c_A);
    cufftExecD2Z(planB, d_B, c_B);


    dim3 blockSize(BSX, BSY, BSZ);   
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x,
        (ny + blockSize.y - 1) / blockSize.y,
        (nz + blockSize.z - 1) / blockSize.z); 

    MultiplyInFrequencyDomain << < gridSize, blockSize >> > (c_A, c_B, c_C, nx, ny, nz);

   
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }
   

    cufftExecZ2D(planC, c_C, d_C);


    cudaMemcpy(h_C_double, d_C, sizeof(cufftDoubleReal) * nx * ny * nz, cudaMemcpyDeviceToHost);
    
    clock_t start, end;
    start = clock();

    
    DoubleReal2Bool(h_C_double, h_C, nx, ny, nz);
   
    cufftDestroy(planA);
    cufftDestroy(planB);
    cufftDestroy(planC);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(c_A);
    cudaFree(c_B);
    cudaFree(c_C);


    end = clock();
    delete[] h_A_double;
    delete[] h_B_double;
    delete[] h_C_double;
}



void ConvolutionOnRealMulti(std::vector<cufftDoubleReal>& h_A, 
    std::vector<std::vector<cufftDoubleReal>>& h_B, 
    std::vector<std::vector<cufftDoubleReal>>& h_C, 
    const int nx, const int ny, const int nz)
{
    cufftDoubleReal* d_A, * d_B, * d_C;
    cudaMalloc((void**)&d_A, sizeof(cufftDoubleReal) * nx * ny * nz);
    cudaMalloc((void**)&d_B, sizeof(cufftDoubleReal) * nx * ny * nz);
    cudaMalloc((void**)&d_C, sizeof(cufftDoubleReal) * nx * ny * nz);

    cufftDoubleComplex* c_A, * c_B, * c_C;
    cudaMalloc((void**)&c_A, sizeof(cufftDoubleComplex) * nx * ny * nz);
    cudaMalloc((void**)&c_B, sizeof(cufftDoubleComplex) * nx * ny * nz);
    cudaMalloc((void**)&c_C, sizeof(cufftDoubleComplex) * nx * ny * nz);


    cufftHandle planA, planB, planC;
    cufftPlan3d(&planA, nx, ny, nz, CUFFT_D2Z);
    cufftPlan3d(&planB, nx, ny, nz, CUFFT_D2Z);
    cufftPlan3d(&planC, nx, ny, nz, CUFFT_Z2D);


    for (int i = 0; i < CutterNum; ++i) {
        cufftDoubleReal* h_A_double = h_A.data();
        cufftDoubleReal* h_B_double = h_B[i].data();
        cufftDoubleReal* h_C_double = h_C[i].data();
        clock_t start, end;
        start = clock();
        cudaMemcpy(d_A, h_A_double, sizeof(cufftDoubleReal) * nx * ny * nz, cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, h_B_double, sizeof(cufftDoubleReal) * nx * ny * nz, cudaMemcpyHostToDevice);

        
        cufftExecD2Z(planA, d_A, c_A);
        cufftExecD2Z(planB, d_B, c_B);

        end = clock();
        dim3 blockSize(BSX, BSY, BSZ);   
        dim3 gridSize((nx + blockSize.x - 1) / blockSize.x,
            (ny + blockSize.y - 1) / blockSize.y,
            (nz + blockSize.z - 1) / blockSize.z); 

        MultiplyInFrequencyDomain << < gridSize, blockSize >> > (c_A, c_B, c_C, nx, ny, nz);

        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            printf("CUDA error: %s\n", cudaGetErrorString(error));
        }


        cufftExecZ2D(planC, c_C, d_C);

        cudaMemcpy(h_C_double, d_C, sizeof(cufftDoubleReal) * nx * ny * nz, cudaMemcpyDeviceToHost);       
    }

    cufftDestroy(planA);
    cufftDestroy(planB);
    cufftDestroy(planC);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(c_A);
    cudaFree(c_B);
    cudaFree(c_C);
}


__global__ void ComputeVisibility(const int* Obs, int * Vis) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;


    if (x < NX && y < NY && z < NZ) {
        int Index = x * NY * NZ + y * NZ + z;
        if (Vis[Index] != 2) 
            return;

        bool visible = true;
        for (int iz = z + 1; iz < NZ; ++iz)
        {
            int idx = x * NY * NZ + y * NZ + iz;
            if (Obs[idx] == 1)
            {
                visible = false;
                break;
            }
        }
        if (visible)
        {
            Vis[Index] = 0;
            return;
        }

        visible = true;
        for (int iz = z - 1; iz >= 0; --iz)
        {
            int idx = x * NY * NZ + y * NZ + iz;
            if (Obs[idx] == 1)
            {
                visible = false;
                break;
            }
        }
        if (visible)
        {
            Vis[Index] = 0;
            return;
        }


        visible = true;
        for (int iy = y + 1; iy < NY; ++iy)
        {
            int idx = x * NY * NZ + iy * NZ + z;
            if (Obs[idx] == 1)
            {
                visible = false;
                break;
            }
        }
        if (visible)
        {
            Vis[Index] = 0;
            return;
        }

        visible = true;
        for (int iy = y - 1; iy >= 0; --iy)
        {
            int idx = x * NY * NZ + iy * NZ + z;
            if (Obs[idx] == 1)
            {
                visible = false;
                break;
            }
        }
        if (visible)
        {
            Vis[Index] = 0;
            return;
        }

        visible = true;
        for (int ix = x + 1; ix < NX; ++ix)
        {
            int idx = ix * NY * NZ + y * NZ + z;
            if (Obs[idx] == 1)
            {
                visible = false;
                break;
            }
        }
        if (visible)
        {
            Vis[Index] = 0;
            return;
        }

        visible = true;
        for (int ix = x - 1; ix >= 0; --ix)
        {
            int idx = ix * NY * NZ + y * NZ + z;
            if (Obs[idx] == 1)
            {
                visible = false;
                break;
            }
        }
        if (visible)
        {
            Vis[Index] = 0;
            return;
        }

    }
}


void VisibilityAnalysisByGPU(const std::vector<bool>& h_C, std::vector<int>& visibility)
{
    for (int i = 0; i < TX; i++) {
        if (h_C[i]) {
            visibility[i] = 1; 
        }
        else {
            visibility[i] = 2;
        }
    }

    
    int* Obs, * Vis;
    cudaMalloc((void**)&Obs, sizeof(int) * NX * NY * NZ);
    cudaMalloc((void**)&Vis, sizeof(int) * NX * NY * NZ);

    int* cpu_Obs = new int[NX * NY * NZ]();
    int* cpu_Vis = new int[NX * NY * NZ]();

    for (int i = 0; i < NX * NY * NZ; ++i) {
        cpu_Obs[i] = visibility[i];
        cpu_Vis[i] = visibility[i];
    }

    cudaMemcpy(Obs, cpu_Obs, sizeof(int) * NX * NY * NZ, cudaMemcpyHostToDevice);
    cudaMemcpy(Vis, cpu_Vis, sizeof(int) * NX * NY * NZ, cudaMemcpyHostToDevice);

    dim3 blockSize(BSX, BSY, BSZ);   
    dim3 gridSize((NX + blockSize.x - 1) / blockSize.x,
        (NY + blockSize.y - 1) / blockSize.y,
        (NZ + blockSize.z - 1) / blockSize.z); 

    ComputeVisibility << < gridSize, blockSize >> > (Obs, Vis);


    cudaMemcpy(cpu_Vis, Vis, sizeof(int) * NX * NY * NZ, cudaMemcpyDeviceToHost);

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }

    for (int i = 0; i < NX * NY * NZ; ++i) {
        visibility[i] = cpu_Vis[i];
    }

    cudaFree(Obs);
    cudaFree(Vis);

    delete[] cpu_Obs;
    delete[] cpu_Vis;
}




__global__ void computeKernel(
    int* Box, int* kernel,
    int loopNum, int kx, int ky, int kz,
    int offsetx, int offsety, int offsetz,
    int kykz, int newBoxYZ, int newBoxZ)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= offsetx && ix - offsetx < kx &&
        iy >= offsety && iy - offsety < ky &&
        iz >= offsetz && iz - offsetz < kz)
    {
        int idx = ix * newBoxYZ + iy * newBoxZ + iz;

        int idx_ = (ix - offsetx) * kykz +
            (iy - offsety) * kz +
            (iz - offsetz);

        if (!Box[idx]) {
            kernel[idx_] += 1;  
        }
    }
}


void ExtractKernel(std::vector<std::vector<bool>>& OffsetBox, 
    std::vector<int>& kernel, 
    std::vector<int>& newBox,
    std::vector<int>& offset, 
    int kx, int ky, int kz,
    int LoopNum) {
    int kykz = ky * kz;
    int newBoxYZ = newBox[1] * newBox[2];

    int* d_kernel;
    int* d_Box;

    size_t BoxSize = newBox[0] * newBox[1] * newBox[2] * sizeof(int);
    size_t kernelSize = kx * ky * kz * sizeof(int);

    cudaMalloc((void**)&d_Box, BoxSize);
    
    cudaMalloc((void**)&d_kernel, kernelSize);

    cudaMemcpy(d_kernel, &kernel[0], kernelSize, cudaMemcpyHostToDevice);



    for (int i = 0; i < LoopNum; i++) {
        std::vector<int> tempBox(newBox[0] * newBox[1] * newBox[2], 0);
        for (int j = 0; j < newBox[0] * newBox[1] * newBox[2]; j++) {
            tempBox[j] = OffsetBox[i][j];
        }

        cudaMemcpy(d_Box, &tempBox[0], BoxSize, cudaMemcpyHostToDevice);

        dim3 blockSize(BSX, BSY, BSZ);   
        dim3 gridSize((newBox[0] + blockSize.x - 1) / blockSize.x,
            (newBox[1] + blockSize.y - 1) / blockSize.y,
            (newBox[2] + blockSize.z - 1) / blockSize.z); 

        computeKernel << <gridSize, blockSize >> > (
            d_Box, d_kernel,
            LoopNum, kx, ky, kz,
            offset[0], offset[1], offset[2],
            kykz, newBoxYZ, newBox[2]);


    }

    cudaDeviceSynchronize();

    cudaMemcpy(&kernel[0], d_kernel, kernelSize, cudaMemcpyDeviceToHost);

    cudaFree(d_Box);
    cudaFree(d_kernel);
    return;
}
