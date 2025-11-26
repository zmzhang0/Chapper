#ifndef CUTTER_H
#define CUTTER_H
#include "CGALib.h"
#include "CUDA_VOXEL/cuda_voxel.h"

using namespace std;
class Cutter {

public:
    Cutter() {
        Rnum = CutterNum;

        if (NZ == 1)
            CutterFile = "data/cutter/cutter-2d-" + to_string(ContainerSize) + ".off";
        else
            CutterFile = "data/cutter/cutter-3d-" + to_string(ContainerSize) + "-" + to_string(CutterLength) + "-" + to_string(CutterWidth) + ".off";


        VoxelCutter.resize(Rnum);
        VoxelCutter_.resize(Rnum);
        Feasibility.resize(Rnum);
        Feasibility_.resize(Rnum);
        TipPoints.resize(Rnum);
        maxBox.resize(3, 0);
        NewSize.resize(3, 0);

    }

    const vector<int>& getMaxBox() {
        return maxBox;
    }

    const vector<int>& getNewSize() {
        return NewSize;
    }

    void InitVoxelCutter() {
        std::cout << "Init Cutter..." << std::endl;
        vector<Polyhedron> mesh_;

        VoxelSize vs_(VOXELSIZE);
        VoxelColor vc_{ 0.6, 0.9, 0.7 };

        mesh_.clear();
        mesh_.resize(Rnum);


        std::ifstream input(CutterFile);
        if (!input) {
            std::cout << "Cannot open file!" << std::endl;
            return;
        }
        mesh_[0].clear();
        input >> mesh_[0];

        for (int j = 1; j < Rnum; j++) {
            mesh_[j].clear();
            mesh_[j] = mesh_[0];
        }

        input.clear();
        input.close();

        SAMPLE_ON_BALL sampleOri;
        sampleOri.num_ori_sample = Rnum;
        sampleOri.EnumerationOritation();
        for (int j = 0; j < Rnum; j++) {
            std::cout << "Sample point " << j << ": " << sampleOri.sample_points[j][0] << " " << sampleOri.sample_points[j][1] << " " << sampleOri.sample_points[j][2] << std::endl;
        }

        trimesh::TriMesh* tri_mesh = new trimesh::TriMesh();
        for (int j = 0; j < Rnum; j++) {
            rotateMesh(mesh_[j], sampleOri.sample_points[j]);
            VoxelCutter[j].resize(NX * NY * NZ);

            tri_mesh->clear();

            Polyhedron2Trimesh(mesh_[j], tri_mesh);

            cuda_voxelizeMesh(tri_mesh, NX, NY, NZ, VOXELSIZE, VoxelCutter[j]);

        }
        delete tri_mesh;

        maxBox.resize(3, 0);
        fill(maxBox.begin(), maxBox.end(), -1);

        for (int j = 0; j < Rnum; j++) {
            for (int k = 0; k < NX * NY * NZ; k++) {
                if (VoxelCutter[j][k]) {
                    int ix = k / (NY * NZ);
                    int iy = (k / NZ) % NY;
                    int iz = k % NZ;
                    maxBox[0] = max(maxBox[0], ix);
                    maxBox[1] = max(maxBox[1], iy);
                    maxBox[2] = max(maxBox[2], iz);
                }
            }
        }

        maxBox[0] += 1;
        maxBox[1] += 1;
        maxBox[2] += 1;


        NewSize.resize(3, 0);
        NewSize[0] = NX + 2 * maxBox[0];
        NewSize[1] = NY + 2 * maxBox[1];

        if (NZ == 1)   
            NewSize[2] = NZ;
        else   
            NewSize[2] = NZ + 2 * maxBox[2];


        for (int j = 0; j < Rnum; j++)
            VoxelCutter_[j].resize(NewSize[0] * NewSize[1] * NewSize[2], 0);


        for (int j = 0; j < Rnum; j++) {
            vector<bool> tmp(maxBox[0] * maxBox[1] * maxBox[2], false);

            for (int ix = 0; ix < maxBox[0]; ix++)
                for (int iy = 0; iy < maxBox[1]; iy++)
                    for (int iz = 0; iz < maxBox[2]; iz++) {
                        int idx = ix * maxBox[1] * maxBox[2] + iy * maxBox[2] + iz;
                        int idx_ = ix * NY * NZ + iy * NZ + iz;
                        tmp[idx] = VoxelCutter[j][idx_];
                    }

            VoxelCutter[j].resize(NewSize[0] * NewSize[1] * NewSize[2], false);
            std::fill(VoxelCutter[j].begin(), VoxelCutter[j].end(), false);

            for (int ix = 0; ix < maxBox[0]; ix++)
                for (int iy = 0; iy < maxBox[1]; iy++)
                    for (int iz = 0; iz < maxBox[2]; iz++) {
                        int idx = ix * maxBox[1] * maxBox[2] + iy * maxBox[2] + iz;
                        int idx_ = ix * NewSize[1] * NewSize[2] + iy * NewSize[2] + iz;
                        VoxelCutter[j][idx_] = tmp[idx];
                    }

            for (int k = 0; k < NewSize[0] * NewSize[1] * NewSize[2]; ++k) {
                if (VoxelCutter[j][k])
                    VoxelCutter_[j][NewSize[0] * NewSize[1] * NewSize[2] - k - 1] = 1.0f;
                else
                    VoxelCutter_[j][NewSize[0] * NewSize[1] * NewSize[2] - k - 1] = 0.0f;
            }

        }

        TipPoints.resize(Rnum);
        for (int j = 0; j < Rnum; j++) {
            TipPoints[j] = findTip(VoxelCutter[j], sampleOri.sample_points[j]);

        }

        if (NZ == 1)   
            maxBox[2] = 0;


        VoxelCutter.clear();
        Feasibility.clear();
    }

    void ConvoluteWith(const vector<bool>& Container, vector<int>& feasibility_) {

        cout << "Convolute with container..." << endl;

        if (NZ == 1)   
            maxBox[2] = 0;

        VoxelSize vs_(VOXELSIZE);
        VoxelColor vc_{ 0.6, 0.9, 0.7 };


        vector<bool> Container_(NewSize[0] * NewSize[1] * NewSize[2], false);

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                for (int k = 0; k < NZ; k++) {
                    int idx = i * NY * NZ + j * NZ + k;
                    int idx_ = (i + maxBox[0]) * NewSize[1] * NewSize[2] + (j + maxBox[1]) * NewSize[2] + k + maxBox[2];

                    Container_[idx_] = Container[idx];
                }


        for (int j = 0; j < Rnum; j++)
            fill(Feasibility[j].begin(), Feasibility[j].end(), false);


        for (int j = 0; j < Rnum; j++) {
            ConvolutionOnReal(Container_, VoxelCutter[j], Feasibility[j], NewSize[0], NewSize[1], NewSize[2]);
        }

        for (int j = 0; j < Rnum; j++) {
            for (int ix = 0; ix < NewSize[0]; ix++)
                for (int iy = 0; iy < NewSize[1]; iy++)
                    for (int iz = 0; iz < NewSize[2]; iz++) {
                        int idx = ix * NewSize[1] * NewSize[2] + iy * NewSize[2] + iz;
                        int idx_ = (ix + TipPoints[j].x() - maxBox[0]) * NY * NZ +
                            (iy + TipPoints[j].y() - maxBox[1]) * NZ + (iz + TipPoints[j].z() - maxBox[2]);

                        if (ix + TipPoints[j].x() - maxBox[0]<0 || ix + TipPoints[j].x() - maxBox[0]>NX - 1 ||
                            iy + TipPoints[j].y() - maxBox[1]<0 || iy + TipPoints[j].y() - maxBox[1]>NY - 1 ||
                            iz + TipPoints[j].z() - maxBox[2]<0 || iz + TipPoints[j].z() - maxBox[2]>NZ - 1)
                            continue;

                        if (idx_ < 0 || idx_ >= NX * NY * NZ)
                            std::cout << "idx_ out of range!" << std::endl;

                        if (!Feasibility[j][idx])
                            feasibility_[idx_]++;
                    }
        }

    }


    void ConvoluteWith(const vector<bool>& object, vector<vector<bool>>& NFV) {

        if (NZ == 1)   
            maxBox[2] = 0;

        VoxelSize vs_(VOXELSIZE);
        VoxelColor vc_{ 0.6, 0.9, 0.7 };


        vector<cufftDoubleReal> object_(NewSize[0] * NewSize[1] * NewSize[2], false);

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                for (int k = 0; k < NZ; k++) {
                    int idx = i * NY * NZ + j * NZ + k;
                    int idx_ = (i + maxBox[0]) * NewSize[1] * NewSize[2] + (j + maxBox[1]) * NewSize[2] + k + maxBox[2];

                    if (object[idx])
                        object_[idx_] = 1.0f;
                    else
                        object_[idx_] = 0.0f;
                }


#pragma omp parallel for
        for (int j = 0; j < Rnum; j++) {
            Feasibility[j].resize(NewSize[0] * NewSize[1] * NewSize[2], false);
            Feasibility_[j].resize(NewSize[0] * NewSize[1] * NewSize[2], 0.0f);
        }

        clock_t start_time, end_time;
        start_time = clock();

        ConvolutionOnRealMulti(object_, VoxelCutter_, Feasibility_, NewSize[0], NewSize[1], NewSize[2]);

        
#pragma omp parallel for
        for (int cn = 0; cn < Rnum; cn++) {
            for (int i = 0; i < NewSize[0]; ++i)
                for (int j = 0; j < NewSize[1]; ++j)
                    for (int k = 0; k < NewSize[2]; ++k)
                    {
                        int idx = ((i + NewSize[0] - 1) % NewSize[0]) * NewSize[1] * NewSize[2] + ((j + NewSize[1] - 1) % NewSize[1]) * NewSize[2] + ((k + NewSize[2] - 1) % NewSize[2]);
                        int idx_ = i * NewSize[1] * NewSize[2] + j * NewSize[2] + k;

                        if (Feasibility_[cn][idx] + 0.5 > 1.0f)
                            Feasibility[cn][idx_] = true;
                        else
                            Feasibility[cn][idx_] = false;
                    }
        }
    
#pragma omp parallel for
        for (int j = 0; j < Rnum; j++) {
            for (int ix = 0; ix < NewSize[0]; ix++)
                for (int iy = 0; iy < NewSize[1]; iy++)
                    for (int iz = 0; iz < NewSize[2]; iz++) {
                        int idx = ix * NewSize[1] * NewSize[2] + iy * NewSize[2] + iz;

                        if (ix + TipPoints[j].x() >= 0 && ix + TipPoints[j].x() < NewSize[0] &&
                            iy + TipPoints[j].y() >= 0 && iy + TipPoints[j].y() < NewSize[1] &&
                            iz + TipPoints[j].z() >= 0 && iz + TipPoints[j].z() < NewSize[2])
                        {
                            int idx_ = (ix + TipPoints[j].x()) * NewSize[1] * NewSize[2] +
                                (iy + TipPoints[j].y()) * NewSize[2] + (iz + TipPoints[j].z());
                            if (idx_ < 0 || idx_ >= NewSize[0] * NewSize[1] * NewSize[2])
                                std::cout << "idx_ out of range!" << std::endl;
                            NFV[j][idx_] = Feasibility[j][idx];
                        }
                    }
        }

        end_time = clock();
        std::cout << "Compute NFV-object time: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

    }
    
    ~Cutter() = default;

private:
    string CutterFile;
    vector<vector<bool>> VoxelCutter;    
    vector<vector<cufftDoubleReal>> VoxelCutter_;    

    vector<vector<bool>> Feasibility;   

    vector<vector<cufftDoubleReal>> Feasibility_;   

    vector<Point> TipPoints;   
    vector<int> maxBox;   
    vector<int> NewSize;   
    int Rnum = 1;   

    Point findTip(vector<bool>&VoxelCutter, const Vector & sampleOri) {
        int x, y, z;
        double max_val = -1e8;
        for (int i = 0; i < maxBox[0]; i++)
            for (int j = 0; j < maxBox[1]; j++)
                for (int k = 0; k < maxBox[2]; k++)
                    if (VoxelCutter[i * NewSize[1] * NewSize[2] + j * NewSize[2] + k])
                    {
                        double val = i * sampleOri[0] + j * sampleOri[1] + k * sampleOri[2];
                        if (val > max_val)
                        {
                            max_val = val;
                            x = i;
                            y = j;
                            z = k;
                        }
                    }
        return Point(x, y, z);
    }

    };

#endif