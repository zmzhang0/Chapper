#include "CUDA_VOXEL/cuda_voxel.h"
#include "Evaluation.h"
#include "VBM/VBM.hpp"
#include <filesystem>

using namespace std;

int ToPacking() {
    clock_t total_time;
    total_time = clock();
    PackingInfo packInfo;
    clock_t start, end;

    int PackedNum = 0; 
    int ObjNum = 0;  
    int layerNum = 0; 


    VoxelSize voxelS;
    voxelS.x = voxelS.y = voxelS.z = VOXELSIZE;
    VoxelColor voxelcolor(0.6, 0.9, 0.7);

    namespace fs = std::filesystem;
    std::string folderPath = "./data";
    int offFileCount = 0;
    vector<string> offFiles;

    try {
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.is_regular_file() && entry.path().extension() == ".off") {
                ++offFileCount;
                offFiles.push_back(entry.path().string());
            }
        }
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << std::endl << std::endl << std::endl;
    std::cout << "Found " << offFileCount << " OFF files in " << folderPath << std::endl;
    ObjNum = offFileCount;  


    vector<bool> container(NX * NY * NZ, false), containerNoBorder(NX * NY * NZ, false); 

    vector<vector<bool>> ConV(OriNum, vector<bool>(NX * NY * NZ, false)), ConVwithBorder(OriNum, vector<bool>(NX * NY * NZ, false));
    vector<vector<vector<bool>>> toPack(ObjNum, vector<vector<bool>>(OriNum, vector<bool>(NX * NY * NZ, false)));

    vector<int> unsigndis(NX * NY * NZ, 0), point_sdf(NX * NY * NZ, 0);
    vector<vector<int>>  evalu(OriNum, vector<int>(NX * NY * NZ, 0)), evaluBBX(OriNum, vector<int>(NX * NY * NZ, 0));


    AddBorder(container);
    
    
    vector<vector<trimesh::TriMesh>> mesh(ObjNum, vector<trimesh::TriMesh>(OriNum));


    std::cout << std::endl << std::endl << std::endl;
    std::cout << "*****************Reading and sorting Mesh!******************" << std::endl;
    start = clock();

    int meshes_num = ObjNum;
    vector<trimesh::TriMesh> meshes(meshes_num);
    vector<Polyhedron> Cgalmeshes(meshes_num);

    for (int i = 0; i < meshes_num; i++) {
        std::cout << std::endl << i << "th object, " <<
            "Reading " << offFiles[i] << std::endl;

        trimesh::TriMesh* the_mesh = trimesh::TriMesh::read(offFiles[i]);
        meshes[i] = *the_mesh;

        {
            std::ifstream input(offFiles[i]);
            if (!input) {
                std::cout << "Cannot open file!" << std::endl;
                return EXIT_FAILURE;
            }
            input >> Cgalmeshes[i];

            if (Cgalmeshes[i].is_empty()) {
                std::cerr << "Error: Loaded mesh is empty." << std::endl;
                return EXIT_FAILURE;
            }

            input.clear();
            input.close();
        }

        std::cout << "Mesh has " << meshes[i].vertices.size() << " vertices and " <<
            meshes[i].faces.size() << " facets." << std::endl << std::endl << std::endl;

    }

    sort(meshes.begin(), meshes.end(), compare_bbx_volume_trimesh);
    sort(Cgalmeshes.begin(), Cgalmeshes.end(), compare_bbx_volume);


    for (int i = 0; i < ObjNum; i++) {
        for (int j = 0; j < OriNum; j++) {
            mesh[i][j].clear();
            mesh[i][j] = meshes[i];
        }
    }

    SAMPLE_ON_BALL sampleOri;
    sampleOri.num_ori_sample = OriNum;
    sampleOri.EnumerationOritation();

    meshes.clear();

    end = clock();
    packInfo.timeRecord.read_and_sort = double(end - start) / CLOCKS_PER_SEC; 
    std::cout << std::endl << std::endl << std::endl <<
        "*****************Starting Voxelization!******************" << std::endl;
    start = clock();


    for (int i = 0; i < ObjNum; i++)
    {
        for (int j = 0; j < OriNum; j++) {
            std::cout << std::endl
                << "######################################## " << i << "th object, " << offFiles[i] << std::endl;

            trimesh::TriMesh* tri_mesh;
            if (NZ != 1)
            {
                rotateMesh(mesh[i][j], sampleOri.sample_points[j]);

                tri_mesh = new trimesh::TriMesh();
                tri_mesh->vertices = mesh[i][j].vertices;
                tri_mesh->faces = mesh[i][j].faces;

                cuda_voxelizeMesh(tri_mesh, NX, NY, NZ, VOXELSIZE, toPack[i][j]);
            } else {
                tri_mesh = new trimesh::TriMesh();
                tri_mesh->vertices = mesh[i][j].vertices;
                tri_mesh->faces = mesh[i][j].faces;

                cuda_voxelizeMesh(tri_mesh, NX, NY, NZ, VOXELSIZE, toPack[i][j]);
            }
            delete tri_mesh;
        }
    }


    end = clock();
    packInfo.timeRecord.voxelization = double(end - start) / CLOCKS_PER_SEC; 
    std::cout << "Voxelization time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  





    std::cout << "*****************Voxelization completed!*****************" << std::endl;
    std::cout << std::endl << std::endl << std::endl;


    for (int i = 0; i < ObjNum; i++) {
        meshes[i].clear();
        for (int j = 0; j < OriNum; j++) {
            mesh[i][j].clear();
        }
    }

    std::cout << "****************************Start To Pack!***************************" << std::endl;
    clock_t start_part, end_part;

    packInfo.RecordObjecs(toPack);

    vector<bool> containter_temp(NX * NY * NZ, false);
    containter_temp[0] = true;

    for (int LoopInd = 0; LoopInd < ObjNum; LoopInd++) {
        int cnt = LoopInd;

        std::cout << "*********************Computing Collision Metric!*********************" << std::endl;

        start_part = clock();

        for (int j = 0; j < OriNum; j++) {
            ConvolutionOnReal(container, toPack[cnt][j], ConVwithBorder[j]);
            ConvolutionOnReal(containerNoBorder, toPack[cnt][j], ConV[j]);
        }


        int offsetWidth = CutterWidth / VOXELSIZE + 1;
        cout << "offsetWidth = " << offsetWidth << endl;


        if (cnt != 0) { 
#pragma omp parallel for
            for (int j = 0; j < OriNum; j++)
                SpaceOffset(ConV[j], offsetWidth);
#pragma omp barrier
        }

#pragma omp parallel for
        for (int j = 0; j < OriNum; j++)
            for (int k = 0; k < NX * NY * NZ; k++)
                ConV[j][k] = ConV[j][k] || ConVwithBorder[j][k];
#pragma omp barrier



        end_part = clock();
        packInfo.timeRecord.collision_detection += double(end_part - start_part) / CLOCKS_PER_SEC; 
        std::cout << "Convolution time = " << double(end_part - start_part) / CLOCKS_PER_SEC << "s" << std::endl;  



        std::cout << "***************Collision Metric Convolution completed!***************" << std::endl << std::endl;

        std::cout << "**********************Computing Proximity Metric!********************" << std::endl;
        start_part = clock();


        if (cnt == 0) {
            vbm::DistanceFunction(NX, NY, NZ, containter_temp, unsigndis);
        }
        
        

        for (int j = 0; j < OriNum; j++) {
            ConvolutionOnReal(unsigndis, toPack[cnt][j], evalu[j]);
        }


        for (int j = 0; j < OriNum; j++)
            HightPenalTerm(evalu[j]);

        end_part = clock();
        packInfo.timeRecord.compute_evaluation += double(end_part - start_part) / CLOCKS_PER_SEC; 
        std::cout << "Evaluation time = " << double(end_part - start_part) / CLOCKS_PER_SEC << "s" << std::endl;  
        std::cout << "****************Proximity Metric Computation completed!**************" << std::endl << std::endl;

        pair<int, ValuePoint> bestOriPlace;
        bestOriPlace.first = -1;
        bestOriPlace.second.value = 1000000000;

        start = clock();

        start_part = clock();

        vector<bool> need_conv(OriNum, false);

        for (int j = 0; j < OriNum; j++) {
            for (int k = 0; k < NX * NY * NZ; k++) {
                if (!ConV[j][k]) {
                    need_conv[j] = true;
                    break;
                }
            }
        }
        packInfo.ComputeNFV_OBJ(toPack[cnt], need_conv);
        
        end_part = clock();
        packInfo.timeRecord.compute_NFV += double(end_part - start_part) / CLOCKS_PER_SEC; 
        std::cout << "Compute NFV-object time = " << double(end_part - start_part) / CLOCKS_PER_SEC << "s" << std::endl << std::endl;  

        std::cout << "***********************Trying to Pack object " << cnt << "!*********************" << std::endl;
        int check_cnt = 0; 

        for (int j = 0; j < OriNum; j++) {

            std::vector<ValuePoint> valuePoints;
            int Pointsnum = 0;

            start_part = clock();
            sortEvaluatedPoints(ConV[j], evalu[j], evaluBBX[j], valuePoints, Pointsnum);
            end_part = clock();
            packInfo.timeRecord.sort_evaluated_points += double(end_part - start_part) / CLOCKS_PER_SEC; 

            start_part = clock(); 
            
            atomic<bool> stop(false);   
            int first_success_index = -1;    

#pragma omp parallel for schedule(dynamic, 1) shared(stop)
            for (int k = 0; k < Pointsnum; k++) {
                if (stop.load()) break;

                bool PackSuccess = PackNewObject(packInfo, containerNoBorder, toPack[cnt][j], j, valuePoints[k].x, valuePoints[k].y, valuePoints[k].z);

                if (PackSuccess) {
#pragma omp critical
                    {
                        if (first_success_index == -1 || k < first_success_index) {
                            cout << "-------------------------------------" << k << endl;
                            first_success_index = k;
                            stop.store(true);
                        }
                    }
                }
            }
#pragma omp barrier


            if (first_success_index != -1) {
                if (bestOriPlace.first == -1 || compare_two_value(valuePoints[first_success_index], bestOriPlace.second)) {
                    bestOriPlace.first = j;
                    bestOriPlace.second = valuePoints[first_success_index];
                }

                std::cout << "Orientation: " << j << ". The Trying times is " << first_success_index + 1 << " / Total Times is " << Pointsnum << "." << std::endl;

                check_cnt += first_success_index + 1;

                std::cout << "Orientation: " << j
                    << ". The best position is (" << valuePoints[first_success_index].x << "," << valuePoints[first_success_index].y << "," << valuePoints[first_success_index].z
                    << "), the value is " << valuePoints[first_success_index].value << " the bbox value is " << valuePoints[first_success_index].value_bbx << std::endl;
            } else {
                std::cout << "Orientation: " << j << ". The Trying times is " << Pointsnum << " / Total Times is " << Pointsnum << "." << std::endl;
                check_cnt += Pointsnum;
            }

            end_part = clock();
            packInfo.timeRecord.check_disassemble += double(end_part - start_part) / CLOCKS_PER_SEC; 
            std::cout << "The " << j << "-th orientation trying time = " << double(end_part - start_part) / CLOCKS_PER_SEC << "s" << std::endl << std::endl;  
        }
        packInfo.timeRecord.check_disassemble_num.push_back(check_cnt);


        if (bestOriPlace.first != -1) {

            Pack(containerNoBorder, toPack[cnt][bestOriPlace.first], bestOriPlace.second.x, bestOriPlace.second.y, bestOriPlace.second.z);



            packInfo.addObject(toPack[cnt][bestOriPlace.first], cnt, bestOriPlace.first, Vector(bestOriPlace.second.x, bestOriPlace.second.y, bestOriPlace.second.z));
            PackedNum++;


            std::cout << "******************Object " << cnt << " is successfully placed!******************" << std::endl;

            std::cout << "The best Orientation: " << bestOriPlace.first
                << ". The best position is (" << bestOriPlace.second.x << "," << bestOriPlace.second.y << "," << bestOriPlace.second.z
                << "), the value is " << bestOriPlace.second.value << " the bbox value is " << bestOriPlace.second.value_bbx << std::endl;
        } else {
            packInfo.addObject(cnt, -1, Vector(-1, -1, -1), Vector(-1, -1, -1));

            PackedNum++;

            std::cout << "****************New object " << cnt << " is failed to be placed!****************" << std::endl;
        }


        end = clock();
        packInfo.timeRecord.placement += double(end - start) / CLOCKS_PER_SEC; 
        std::cout << "Placement time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  


    }


    std::cout << "**************************Packing completed!*************************" << std::endl;
    std::cout << std::endl << std::endl << std::endl;



    ofstream out("output/Density.txt");

    out << "************************Density of Container!************************" << std::endl;

    out << "Density Of container is " << double((packInfo.GetDensityOf(Cgalmeshes, ContainerSize - VOXELSIZE, ContainerSize - VOXELSIZE, ContainerSize - VOXELSIZE) / 10.0)) << "%" << std::endl;
    
    out << "Density (based voxel) Of container is " << double((DensityOf(containerNoBorder) / 10.0)) << "%" << std::endl;

    Bool2Voxel("output/Pack-Result.obj", containerNoBorder, voxelcolor, voxelS);

    packInfo.ExportIVHT(out);

    packInfo.ExportMeshes(out, Cgalmeshes);


    start_part = clock();
    packInfo.ComputeDisassembleOrder();
    end_part = clock();
    packInfo.timeRecord.compute_disassemble_order += double(end_part - start_part) / CLOCKS_PER_SEC; 

    packInfo.ExportDOT();

    packInfo.ExportDisassembleMeshes(out, Cgalmeshes);


    out.close();

    total_time = clock() - total_time;
    packInfo.timeRecord.total = double(total_time) / CLOCKS_PER_SEC; 

    packInfo.printTimeRecord();

    return 0;
}

int main() {
    clock_t start, end;
    start = clock();

    ToPacking();

    end = clock();
    std::cout << "Total time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  
    return 0;
}