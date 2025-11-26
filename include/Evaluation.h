#pragma once
#include "RunCuda.h"
#include <stack>
#include <algorithm>
#include "CGALib.h"   

#include "GPUFunction.h" 
#include "VBM/VBM.hpp"
#include <execution>
#include "LimitedSizeMap.h"

#include <random>

#include "cutter.h"

#include <omp.h>

#include "Meshcut.h"

struct Plank {
    Point P0;        
    Point P1;        
    Vector Dir;       
    Vector Sweep;     
    Vector normal;    

    float Width;       
    float Height;      
    float score;      
};


struct ValuePoint {
    int value; 
    int value_bbx;  
    int x, y, z;   
};

struct ObjectInfo {
    Vector center;  
    Vector pos;  
    int ObjId;  
    int OriId;  
    bool packed = false;  
    ObjectInfo(const int objId,const int OriId,const Vector& center,const Vector& pos) : center(center), pos(pos), ObjId(objId), OriId(OriId) {}
    ObjectInfo() : center(Vector(0, 0, 0)), pos(Vector(-1, -1, -1)), ObjId(-1), OriId(-1) {}
};


struct TreeNode {
    int id;  
    vector<int> data;
    string imagePath;
    vector<shared_ptr<TreeNode>> children;
};

struct LegalPlane {
    Point point;              
    Vector n1, n2;            
    Vector normal;            
    int side;                 
};

struct LegalCone {
    Point point;              
    Vector n1, n2, n3;        
    double angle;             
    bool isInCone;            
    bool isCone;                 

    LegalCone() : isCone(false), isInCone(false) {}

    LegalCone(const Point& p, const Vector& n1, const Vector& n2, const Vector& n3, const bool& isInCone)
        : point(p), n1(n1), n2(n2), n3(n3), isInCone(isInCone) {
        isCone = true; angle = 0.0;
    }

    LegalCone(const Point& p, const Vector& n1, const Vector& n2, const bool& isInCone)
        : point(p), n1(n1), n2(n2), isInCone(isInCone) {
        isCone = false; angle = 0.0;
    }
};

struct TimeRecord {
    double total = 0.01;
    double read_and_sort = 0.01;
    double voxelization = 0.01;
    double collision_detection = 0.01;
    double compute_evaluation = 0.01;

    double sort_evaluated_points = 0.01;

    double placement = 0.01;

    double compute_NFV = 0.01; 

    double check_disassemble = 0.01; 

    double Recursion_time = 0.01;
    int Recursion_num = 0;
    double compute_NFV_Recursion = 0.01;
    double compute_visibility = 0.01;
    double compute_connectedness = 0.01;


    double compute_disassemble_order = 0.01;  

    double export_result = 0.01;
    double tmp = 0.01;

    vector<int> check_disassemble_num;

    double alpha_shape = 0.01;   
    double export_alpha_shape = 0.01;   
};

class PackingInfo {
public:
    PackingInfo();

    ~PackingInfo() = default;

    int GetDensityOf(const vector<Polyhedron>& meshes, const double x, const double y, const double z)
    {

        double tot_volume = 0.0;
        for (int i = 0; i < visitedNum; i++) {
            if (positions[i].packed) {
                tot_volume += CGAL::Polygon_mesh_processing::volume(meshes[i]);
            }
        }
        double ConSize = x * y * z;

        int density = tot_volume * 1000 / ConSize;

        return density;
    }

    void ResetUniqueId() {
        unique_id = 0;
    }

    void RecordObjecs(vector<vector<vector<bool>>>& toPack) {
        objects = make_shared<vector<vector<vector<bool>>>>(toPack);
    }
    

    void ComputeNFV_OBJ(const vector<vector<bool>>& object, const vector<bool>& need_conv);

    void RecordCombinations(const int& objId, const int& oriId, const int& x, const int& y, const int& z);


    void ComputeVis(const int& oriId, const int& comb_num, const int& x, const int& y, const int& z, vector<int>& visibility);

    void ComputeVis(const vector<vector<bool>>& NFV, vector<int>& visibility);

    void addObject(const vector<bool>& object, const int& objId, const int& oriId, const Vector& pos);

    void addObject(const int& objId, const int& oriId, const Vector& center, const Vector& pos) {
        positions.push_back(ObjectInfo(objId, oriId, center, pos));
        visitedNum++;
        if (visitedNum - 1 != objId)
            std::cout << "OrderError: visitedNum is not equal to objId" << std::endl;
    }

    void ExportMeshes(ofstream& txtout, vector<Polyhedron>& meshes, string label = "");
    void ExportDisassembleMeshes(ofstream& txtout, vector<Polyhedron>& meshes);


    void ExportMeshes(ofstream& txtout, vector<vector<trimesh::TriMesh>>& meshes);



    void ExportIVHT(ofstream& txtout);

    void ExportIVHTInfo(ofstream& txtout, const vector<int> connected_objects, const int layer, shared_ptr<TreeNode>& current);

    void ComputeDisassembleOrder();

    void ExportDOT();

    void ExportIVHTNode(const shared_ptr<TreeNode>& node, std::ofstream& file);


    bool CheckDisassembleWithIVHT(
        const int oriId,
        const int x, const int y, const int z,
        vector<int> connected_objects, int depth);
    bool CheckDisassembleWithIVHTLoop(const int oriId, const int x, const int y, const int z, vector<int> connected_objects, int depth);

    int computeConnectedComponent(const vector<int>& visibility, const vector<int>& connected_objects, vector<int>& visited);


    void GetObjIds(vector<int>& objIds) {
        for (int i = 0; i < visitedNum; i++)
            if (positions[i].packed)
                objIds.push_back(i);
    }

    bool isSlientPrint() const {
        return isSlient;
    }

    void printTimeRecord();

    TimeRecord timeRecord;  
private:
    vector<ObjectInfo> positions;
    Cutter cutter;  

    vector<int> newBox = { 0, 0, 0 };

    vector<vector<bool>> NFV_container; 
    vector<vector<vector<bool>>> NFV_object;
    vector<vector<vector<Vector>>> NFV_pos;
    int visitedNum = 0;
    int objectNum = 0;


    map<vector<int>, pair<int, int>> CombinationMap;

    vector<vector<vector<bool>>> NFV_objects_comb;
    int ConbinationNum = 0;  

    shared_ptr<vector<vector<vector<bool>>>> objects;  

    bool isSlient = true; 

    int unique_id = 0;

    shared_ptr<TreeNode> root;

    vector<int> disassembleOrder;  
    vector<vector<LegalCone>> disassembleCones;  

    int currentDepth = 0;
};



inline bool isSubsequence(const vector<int>& toMatch, const vector<int>& sub);
inline void DiffSequence(const vector<int>& toMatch, const vector<int>& sub, vector<int>& difference);

void SpaceOffset(std::vector<bool>& space, const int offset);

static bool compare(const ValuePoint& a, const ValuePoint& b);

bool compare_two_value(const ValuePoint& a, const ValuePoint& b);


void sortEvaluatedPoints(const std::vector<bool>& ConV, const std::vector<int>& evalu, const std::vector<int>& evalu_bbx, std::vector<ValuePoint>& points, int& num);


int DensityOf(const std::vector<bool>& s_omiga);


void HightPenalTerm(std::vector<int>& evalu);


bool PackNewObject(PackingInfo& packInfo, const std::vector<bool>& h_A, std::vector<bool>& h_B, const int oriId, const int x, const int y, const int z);

void Pack(std::vector<bool>& container, std::vector<bool>& object, const int x, const int y, const int z);

void InternalVisualHull(const std::vector<bool>& h_C, std::vector<int>& visibility);






void AddBorder(std::vector<bool>& container);

void RemoveBorder(const std::vector<bool>& container, std::vector<bool>& containerNoBorder);

void RemoveBorder(std::vector<bool>& container);


