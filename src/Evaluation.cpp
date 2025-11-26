#include "Evaluation.h"
#include "gco-v3.0/GCoptimization.h"
#include <igl/marching_cubes.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/remove_unreferenced.h>


void SpaceOffset(std::vector<bool>& space, const int offset) {
	std::vector<int> unsigndis(NX * NY * NZ, 0);
	vbm::DistanceFunction(NX, NY, NZ, space, unsigndis);  
	for (int x = 0; x < NX; x++)
		for (int y = 0; y < NY; y++)
			for (int z = 0; z < NZ; z++)
			{
				if (!space[x * NY * NZ + y * NZ + z]) { 
					if (unsigndis[x * NY * NZ + y * NZ + z] <= offset)  
						space[x * NY * NZ + y * NZ + z] = true;
				}
			}
}

static bool compare(const ValuePoint& a, const ValuePoint& b) {
	return a.value < b.value;
}
bool compare_two_value(const ValuePoint& a, const ValuePoint& b) {
	if (a.value_bbx < b.value_bbx)
		return true;
	else if (a.value_bbx > b.value_bbx)
		return false;
	else
		return a.value < b.value;
}

void sortEvaluatedPoints(const std::vector<bool>& ConV, const std::vector<int>& evalu, const std::vector<int>& evalu_bbx, std::vector<ValuePoint>& points, int& num) {

	points.clear();
	num = 0;

	for (int iz = 0; iz < NZ; iz++)
		for (int iy = 0; iy < NY; iy++)
			for (int ix = 0; ix < NX; ix++) {
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (!ConV[idx]) {   
					points.push_back(ValuePoint{ evalu[idx],evalu_bbx[idx], ix, iy, iz });
					num++;
				}
			}

	std::sort(points.begin(), points.end(), compare);
}





int DensityOf(const std::vector<bool>& s_omiga)
{
	double ConSize = TX;
	double cnt = 0;
	for (int i = 0; i < TX; i++)
	{
		if (s_omiga[i])
			cnt++;
	}
	int density = cnt * 1000 / ConSize;

	return density;
}





void HightPenalTerm(std::vector<int>& evalu) {
	for (int ix = 0; ix < NX; ix++) {
		for (int iy = 0; iy < NY; iy++) {
			for (int iz = 0; iz < NZ; iz++) {
				int idx = ix * NY * NZ + iy * NZ + iz;
				double NumOfN = NZ;
				double penalty = (iz + 1) / NumOfN;
				penalty = penalty * penalty * penalty;
				penalty = 1e6 * penalty;
				evalu[idx] += penalty;

			}
		}
	}
}




bool PackNewObject(PackingInfo& packInfo, const std::vector<bool>& h_A, std::vector<bool>& h_B, const int oriId, const int x, const int y, const int z) {
	std::vector<bool> h_C(TX, false);
	for (int i = 0; i < TX; i++) {
		h_C[i] = h_A[i];
	}

	for (int ix = 0; ix < NX; ++ix)
		for (int iy = 0; iy < NY; ++iy)
			for (int iz = 0; iz < NZ; ++iz)
			{
				int idx_b = ix * NY * NZ + iy * NZ + iz;
				if (h_B[idx_b]) {
					int idx_c = ((x + ix) % NX) * NY * NZ + ((y + iy) % NY) * NZ + ((z + iz) % NZ);
					if (h_C[idx_c]) {
						std::cout << "PackError: The new object will cause collision with the existing objects." << std::endl;
						return false;
					}
					h_C[idx_c] = true;
				}
			}

	clock_t start, end;
	start = std::clock();

	vector<int> connectedObjIds;
	packInfo.GetObjIds(connectedObjIds);
	if (connectedObjIds.size() == 0)
		return true;

	if (!packInfo.isSlientPrint())
		cout << "############### Internally Visible Hull Tree ##################" << endl;

	packInfo.ResetUniqueId();  
	bool isDisassemble = packInfo.CheckDisassembleWithIVHTLoop(oriId, x, y, z, connectedObjIds, 0);

	if (!packInfo.isSlientPrint())
		cout << "###############################################################" << endl << endl;

	end = std::clock();
	packInfo.timeRecord.tmp += (double)(end - start) / CLOCKS_PER_SEC;
	if (isDisassemble)
		return true;
	else
		return false;
}


void Pack(std::vector<bool>& container, std::vector<bool>& object, const int x, const int y, const int z) {

	for (int ix = 0; ix < NX; ++ix)
		for (int iy = 0; iy < NY; ++iy)
			for (int iz = 0; iz < NZ; ++iz)
			{
				int idx_b = ix * NY * NZ + iy * NZ + iz;
				if (object[idx_b]) {
					int idx_c = ((x + ix) % NX) * NY * NZ + ((y + iy) % NY) * NZ + ((z + iz) % NZ);
					if (container[idx_c]) {
						std::cout << "PackError: The new object will cause collision with the existing objects." << std::endl;
						return;
					}
					container[idx_c] = true;
				}
			}
	return;
}

void InternalVisualHull(const std::vector<bool>& h_C, std::vector<int>& visibility) {
	for (int i = 0; i < TX; i++) {
		if (h_C[i]) {
			visibility[i] = 1; 
		}
		else {
			visibility[i] = 2;
		}
	}



	for (int ix = 0; ix < NX; ++ix)
		for (int iy = 0; iy < NY; ++iy) {
			for (int iz = 0; iz < NZ - 1; ++iz)
			{
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (visibility[idx] != 1)
					visibility[idx] = 0;
				else
					break;
			}
			for (int iz = NZ - 1; iz >= 1; --iz) {
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (visibility[idx] != 1)
					visibility[idx] = 0;
				else
					break;
			}
		}

	for (int ix = 0; ix < NX; ++ix)
		for (int iz = 0; iz < NZ; ++iz) {
			for (int iy = 0; iy < NY - 1; ++iy)
			{
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (visibility[idx] != 1)
					visibility[idx] = 0;
				else
					break;
			}
			for (int iy = NY - 1; iy >= 1; --iy) {
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (visibility[idx] != 1)
					visibility[idx] = 0;
				else
					break;
			}
		}

	for (int iy = 0; iy < NY; ++iy)
		for (int iz = 0; iz < NZ; ++iz) {
			for (int ix = 0; ix < NX - 1; ++ix)
			{
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (visibility[idx] != 1)
					visibility[idx] = 0;
				else
					break;
			}
			for (int ix = NX - 1; ix >= 1; --ix) {
				int idx = ix * NY * NZ + iy * NZ + iz;
				if (visibility[idx] != 1)
					visibility[idx] = 0;
				else
					break;
			}
		}

	if (OutPut2DResult)
		Int2VoxelUseColorBar("output/midOutput/2d-Visiblity.obj", visibility, VoxelSize(VOXELSIZE));
}

void AddBorder(std::vector<bool>& container) {
	for (int ix = 0; ix < NX; ix++)
		for (int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++)
			{
				if (NZ != 1) { 
					if((ix== 0 && iy==0)||(ix==NX-1 && iy==0)||(ix==0 && iy==NY-1)||(ix==NX-1 && iy==NY-1)||
						(iz==0 && iy==0)||(iz==NZ-1 && iy==0)||(iz==0 && iy==NY-1)||(iz==NZ-1 && iy==NY-1)||
						(ix==0 && iz==0)||(ix==NX-1 && iz==0)||(ix==0 && iz==NZ-1)||(ix==NX-1 && iz==NZ-1))
						container[ix * NY * NZ + iy * NZ + iz] = 1;
				}
				else 
				{
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1)
						container[ix * NY * NZ + iy * NZ + 0] = 1;
				}
			}

	Bool2Voxel("output/midOutput/x-border.obj", container, VoxelColor(0.6, 0.6, 0.7), VoxelSize(VOXELSIZE));

	for (int ix = 0; ix < NX; ix++)
		for (int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++)
			{
				if (NZ != 1) { 
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1 || iz == 0 || iz == NZ - 1)
						container[ix * NY * NZ + iy * NZ + iz] = 1;
				}
				else 
				{
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1)
						container[ix * NY * NZ + iy * NZ + 0] = 1;
					
				}
			}

	if (ForFabrication) {
		for (int ix = 0; ix < NX; ix++)
			for (int iy = 0; iy < NY; iy++)
				for (int iz = 0; iz < NZ; iz++) {
					if (iz >= 47)
						container[ix * NY * NZ + iy * NZ + iz] = 1;
				}
	}
}

void RemoveBorder(const std::vector<bool>& container, std::vector<bool>& containerNoBorder) {
	for (int ix = 0; ix < NX; ix++)
		for (int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++)
			{
				if (NZ != 1) { 
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1 || iz == 0 || iz == NZ - 1)
						containerNoBorder[ix * NY * NZ + iy * NZ + iz] = 0;
					else
						containerNoBorder[ix * NY * NZ + iy * NZ + iz] = container[ix * NY * NZ + iy * NZ + iz];
				}
				else 
				{
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1)
						containerNoBorder[ix * NY * NZ + iy * NZ + 0] = 0;
					else
						containerNoBorder[ix * NY * NZ + iy * NZ + 0] = container[ix * NY * NZ + iy * NZ + 0];
				}
			}
}

void RemoveBorder(std::vector<bool>& container) {
	for (int ix = 0; ix < NX; ix++)
		for (int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++)
			{
				if (NZ != 1) { 
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1 || iz == 0 || iz == NZ - 1)
						container[ix * NY * NZ + iy * NZ + iz] = 0;
				}
				else 
				{
					if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1)
						container[ix * NY * NZ + iy * NZ + 0] = 0;
				}

			}
}


PackingInfo::PackingInfo()
{
	visitedNum = 0;
	cutter.InitVoxelCutter();

	std::cout << "Cutter Initialized!" << std::endl;

	newBox[0] = cutter.getNewSize()[0];
	newBox[1] = cutter.getNewSize()[1];
	newBox[2] = cutter.getNewSize()[2];
	NFV_container.resize(CutterNum);

	for (int j = 0; j < CutterNum; j++) {
		NFV_container[j].resize(newBox[0] * newBox[1] * newBox[2], false);
	}

	NFV_object.resize(OriNum);
	NFV_pos.resize(OriNum);
	for (int i = 0; i < OriNum; i++) {
		NFV_object[i].resize(CutterNum, vector<bool>(newBox[0] * newBox[1] * newBox[2], false));
		NFV_pos[i].resize(CutterNum, vector<Vector>(0));
	}

	positions.clear();
	NFV_objects_comb.clear();
	CombinationMap.clear();

}

void PackingInfo::ComputeNFV_OBJ(const vector<vector<bool>>& object, const vector<bool>& need_conv) {
	for (int i = 0; i < OriNum; i++)
#pragma omp parallel for
		for (int j = 0; j < CutterNum; j++) {
			fill(NFV_object[i][j].begin(), NFV_object[i][j].end(), false);
			NFV_pos[i][j].clear();
		}

	clock_t start_time, end_time;
	start_time = clock();

	for (int i = 0; i < OriNum; i++)
	{
		if (need_conv[i])
			cutter.ConvoluteWith(object[i], NFV_object[i]);
	}

	end_time = clock();

	for (int i = 0; i < OriNum; i++)
#pragma omp parallel for
		for (int j = 0; j < CutterNum; j++) {
			for (int k = 0; k < newBox[0] * newBox[1] * newBox[2]; k++)
				if (NFV_object[i][j][k]) {
					int x = k / newBox[1] / newBox[2];
					int y = (k / newBox[2]) % newBox[1];
					int z = k % newBox[2];
					NFV_pos[i][j].push_back(Vector(x, y, z));
				}
		}
}

void PackingInfo::RecordCombinations(const int& objId, const int& oriId, const int& x, const int& y, const int& z)
{
	vector<vector<bool>> NFV_objects_comb_tmp(CutterNum, vector<bool>(newBox[0] * newBox[1] * newBox[2], false));
	vector<int> combination_tmp;

	int dx, dy, dz; 
	for (int i = 0; i < TX; i++)
		if ((*objects)[visitedNum - 1][oriId][i]) {

			int ix = i / NY / NZ;
			int iy = (i / NZ) % NY;
			int iz = i % NZ;

			dx = x;
			dy = y;
			dz = z;
			if (ix + x >= NX)
				dx = x - NX;
			if (iy + y >= NY)
				dy = y - NY;
			if (iz + z >= NZ)
				dz = z - NZ;
			break;
		}

	for (int i = 0; i < CutterNum; i++)
		for (int j = 0; j < NFV_pos[oriId][i].size(); j++)
		{
			int cx = NFV_pos[oriId][i][j].x() + dx;

			int cy = NFV_pos[oriId][i][j].y() + dy;

			int cz = NFV_pos[oriId][i][j].z() + dz;

			if (cx >= 0 && cx < newBox[0] && cy >= 0 && cy < newBox[1] && cz >= 0 && cz < newBox[2])
			{
				NFV_container[i][cx * newBox[1] * newBox[2] + cy * newBox[2] + cz] = true;
				NFV_objects_comb_tmp[i][cx * newBox[1] * newBox[2] + cy * newBox[2] + cz] = true;
			}
		}

	combination_tmp.clear();
	combination_tmp.push_back(objId);

	vector<vector<bool>> kernel(CutterNum, vector<bool>(NX * NY * NZ, false));

	for (int j = 0; j < CutterNum; j++)
		for (int ix = 0; ix < NX; ix++)
			for (int iy = 0; iy < NY; iy++)
				for (int iz = 0; iz < NZ; iz++)
				{
					int idx = ix * NY * NZ + iy * NZ + iz;

					int idx_ = (ix + cutter.getMaxBox()[0]) * newBox[1] * newBox[2] + (iy + cutter.getMaxBox()[1]) * newBox[2] + (iz + cutter.getMaxBox()[2]);

					if (idx_ < 0 || idx_ >= newBox[0] * newBox[1] * newBox[2])
						std::cout << "idx_ out of range!" << std::endl;

					if (NFV_objects_comb_tmp[j][idx_])
						kernel[j][idx] = true;
					else
						kernel[j][idx] = false;
				}


	if (CombinationMap.find(combination_tmp) == CombinationMap.end())
	{
		NFV_objects_comb.push_back(kernel);
		CombinationMap[combination_tmp] = make_pair(ConbinationNum, 1);   
		ConbinationNum++;
	}

	for (int i = 0; i < CutterNum; i++)
		fill(kernel[i].begin(), kernel[i].end(), false);

	for (int j = 0; j < CutterNum; j++)
		for (int ix = 0; ix < NX; ix++)
			for (int iy = 0; iy < NY; iy++)
				for (int iz = 0; iz < NZ; iz++)
				{
					int idx = ix * NY * NZ + iy * NZ + iz;

					int idx_ = (ix + cutter.getMaxBox()[0]) * newBox[1] * newBox[2] + (iy + cutter.getMaxBox()[1]) * newBox[2] + (iz + cutter.getMaxBox()[2]);

					if (idx_ < 0 || idx_ >= newBox[0] * newBox[1] * newBox[2])
						std::cout << "idx_ out of range!" << std::endl;

					if (NFV_container[j][idx_])
						kernel[j][idx] = true;
					else
						kernel[j][idx] = false;
				}

	vector<int> visibility_tmp(TX, 0);
	ComputeVis(kernel, visibility_tmp);


	for (int i = 0; i < CutterNum; i++) {
		vector<int> kernel_tmp(NX * NY * NZ, 1);
		for (int j = 0; j < NX * NY * NZ; j++)
		{
			if (kernel[i][j] == 1)
				kernel_tmp[j] = 0;
		}
	}

	for (int i = 0; i < CutterNum; i++) {
		vector<int> kernel_tmp(NX * NY * NZ, 1);
		for (int j = 0; j < NX * NY * NZ; j++)
		{
			if (visibility_tmp[j] >= i)
				kernel_tmp[j] = 0;
		}
	}


	combination_tmp.clear();
	GetObjIds(combination_tmp);
	vector<int> visited_tmp(TX, 0);

	int comb_num_tmp = computeConnectedComponent(visibility_tmp, combination_tmp, visited_tmp);


	if (CombinationMap.find(combination_tmp) == CombinationMap.end())
	{
		NFV_objects_comb.push_back(kernel);
		CombinationMap[combination_tmp] = make_pair(ConbinationNum, comb_num_tmp);   
		ConbinationNum++;
	}
}

void PackingInfo::ComputeVis(const int& oriId, const int& comb_num, const int& x, const int& y, const int& z, vector<int>& visibility)
{
	vector<vector<bool>> Vis(CutterNum, vector<bool>(NX * NY * NZ, false));

#pragma omp critical
	{
		for (int i = 0; i < CutterNum; i++)
			copy(NFV_objects_comb[comb_num][i].begin(), NFV_objects_comb[comb_num][i].end(), Vis[i].begin());
	}

	int dx, dy, dz; 
	for (int i = 0; i < TX; i++)
		if ((*objects)[visitedNum][oriId][i]) {
			int ix = i / NY / NZ;
			int iy = (i / NZ) % NY;
			int iz = i % NZ;

			dx = x;
			dy = y;
			dz = z;
			if (ix + x >= NX)
				dx = x - NX;
			if (iy + y >= NY)
				dy = y - NY;
			if (iz + z >= NZ)
				dz = z - NZ;
			break;
		}

	for (int i = 0; i < CutterNum; i++)
		for (int j = 0; j < NFV_pos[oriId][i].size(); j++)
		{
			int cx = NFV_pos[oriId][i][j].x() + dx;

			int cy = NFV_pos[oriId][i][j].y() + dy;

			int cz = NFV_pos[oriId][i][j].z() + dz;

			if ((cx - cutter.getMaxBox()[0]) >= 0 && (cx - cutter.getMaxBox()[0]) < NX &&
				(cy - cutter.getMaxBox()[1]) >= 0 && (cy - cutter.getMaxBox()[1]) < NY &&
				(cz - cutter.getMaxBox()[2]) >= 0 && (cz - cutter.getMaxBox()[2]) < NZ) {
				if (NZ != 1)
					Vis[i][(cx - cutter.getMaxBox()[0]) * NY * NZ + (cy - cutter.getMaxBox()[1]) * NZ + (cz - cutter.getMaxBox()[2])] = true;
				else
					Vis[i][(cx - cutter.getMaxBox()[0]) * NY * NZ + (cy - cutter.getMaxBox()[1]) * NZ + 0] = true;
			}
		}


	clock_t start, end;
	start = clock();

	for (int j = 0; j < CutterNum; j++) {
		for (int ix = 0; ix < NX; ix++)
			for (int iy = 0; iy < NY; iy++)
				for (int iz = 0; iz < NZ; iz++)
				{
					int idx = ix * NY * NZ + iy * NZ + iz;

					if (!Vis[j][idx])
						visibility[idx]++;
				}
	}

	end = clock();
}

void PackingInfo::ComputeVis(const vector<vector<bool>>& NFV, vector<int>& visibility)
{
	for (int j = 0; j < CutterNum; j++) {
		for (int ix = 0; ix < NX; ix++)
			for (int iy = 0; iy < NY; iy++)
				for (int iz = 0; iz < NZ; iz++)
				{
					int idx = ix * NY * NZ + iy * NZ + iz;

					if (!NFV[j][idx])
						visibility[idx]++;
				}
	}
}

void PackingInfo::addObject(const vector<bool>& object, const int& objId, const int& oriId, const Vector& pos)
{
	int x, y, z;
	for (int i = 0; i < TX; i++)
		if (object[i]) {
			x = i / NY / NZ;
			y = (i / NZ) % NY;
			z = i % NZ;
			break;
		}

	x += pos.x();
	x %= NX;
	y += pos.y();
	y %= NY;
	z += pos.z();
	z %= NZ;

	positions.push_back(ObjectInfo(objId, oriId, Vector(x, y, z), pos));

	positions[visitedNum].packed = true;

	visitedNum++;
	objectNum++;

	clock_t start, end;
	start = clock();
	RecordCombinations(objId, oriId, pos.x(), pos.y(), pos.z());
	end = clock();
	timeRecord.compute_NFV += (double)(end - start) / CLOCKS_PER_SEC;

	if (visitedNum - 1 != objId)
		std::cout << "OrderError: visitedNum is not equal to objId" << std::endl;
}

void PackingInfo::ExportMeshes(ofstream& txtout, vector<Polyhedron>& meshes, string label)
{
	SAMPLE_ON_BALL sampleOri;
	sampleOri.num_ori_sample = OriNum;
	sampleOri.EnumerationOritation();

	for (int i = 0; i < visitedNum; i++) {
		txtout << "The position of object " << i << " is (" << positions[i].pos[0] << "," << positions[i].pos[1] << "," << positions[i].pos[2] << ")" << std::endl;

		if (positions[i].packed) {
			rotateMesh(meshes[positions[i].ObjId], sampleOri.sample_points[positions[i].OriId]);

			translateMesh(meshes[positions[i].ObjId], positions[i].pos);
			std::string filename = "output/midOutput/Object/Object_" + label + std::to_string(i) + ".obj";
			std::ofstream output(filename);

			CGAL::IO::write_OBJ(output, meshes[positions[i].ObjId]);

			output.close();

		}
		else {
			std::string filename = "output/midOutput/Object/fail_Object_" + label + std::to_string(i) + ".obj";
			std::ofstream output(filename);

			CGAL::IO::write_OBJ(output, meshes[positions[i].ObjId]);

			output.close();
		}
	}

}

void PackingInfo::ExportDisassembleMeshes(ofstream& txtout, vector<Polyhedron>& meshes)
{

	for (int i = 0; i < disassembleOrder.size(); i++) {
		int idx = disassembleOrder[i];
		txtout << "The position of object " << idx << " is (" << positions[idx].pos[0] << "," << positions[idx].pos[1] << "," << positions[idx].pos[2] << ")" << std::endl;

		if (positions[idx].packed) {
			std::string filename = "output/midOutput/Object/Disassemble_Object_" + std::to_string(i) + ".obj";
			std::ofstream output(filename);

			CGAL::IO::write_OBJ(output, meshes[positions[idx].ObjId]);

			output.close();
		}
	}

}

void PackingInfo::ExportMeshes(ofstream& txtout, vector<vector<trimesh::TriMesh>>& meshes)
{
	for (int i = 0; i < visitedNum; i++) {
		txtout << "The position of object " << i << " is (" << positions[i].pos[0] << "," << positions[i].pos[1] << "," << positions[i].pos[2] << ")" << std::endl;

		if (positions[i].packed) {
			translateMesh(meshes[positions[i].ObjId][positions[i].OriId], positions[i].pos);
			std::string filename = "output/midOutput/Object/Object_" + std::to_string(i) + ".off";
			meshes[positions[i].ObjId][positions[i].OriId].write(filename);
		}
	}

	for (int i = 0; i < disassembleOrder.size(); i++) {
		int idx = disassembleOrder[i];
		txtout << "The position of object " << idx << " is (" << positions[idx].pos[0] << "," << positions[idx].pos[1] << "," << positions[idx].pos[2] << ")" << std::endl;

		if (positions[idx].packed) {
			std::string filename = "output/midOutput/Object/Disassemble_Object_" + std::to_string(i) + ".off";

			meshes[positions[idx].ObjId][positions[idx].OriId].write(filename);
		}
	}

}

int PackingInfo::computeConnectedComponent(const vector<int>& visibility, const vector<int>& connected_objects, vector<int>& visited)
{

	const int sDir[26][3] = {
		{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
		{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
		{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
		{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
		{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
		{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };

	fill(visited.begin(), visited.end(), 0);
	std::stack<int> dfs_stack;
	std::stack<int> init_stack;

	int component_num = 0;

	for (int i = 0; i < connected_objects.size(); i++) {
		int x = positions[connected_objects[i]].center[0];
		int y = positions[connected_objects[i]].center[1];
		int z = positions[connected_objects[i]].center[2];
		init_stack.push(x * NY * NZ + y * NZ + z);
		if (x * NY * NZ + y * NZ + z < 0 || x * NY * NZ + y * NZ + z >= NX * NY * NZ)
			std::cout << "3Init stack Error: The start point is out of range!" << std::endl;
		if (visibility[x * NY * NZ + y * NZ + z] != 0)
		{
			cout << x << "," << y << "," << z << endl;
			std::cout << "3Init stack Error: The start point is Unvisible!" << std::endl;
		}
	}


	while (!init_stack.empty()) {
		int idx_ = init_stack.top();
		init_stack.pop();

		if (visited[idx_] != 0)
			continue;

		dfs_stack.push(idx_);
		visited[idx_] = ++component_num;

		while (!dfs_stack.empty()) {
			int idx = dfs_stack.top();
			dfs_stack.pop();

			int ix, iy, iz;
			ix = idx / (NY * NZ);
			iy = (idx / NZ) % NY;
			iz = idx % NZ;
			for (int j = 0; j < 26; j++)
			{
				int jx = ix + sDir[j][0];
				int jy = iy + sDir[j][1];
				int jz = iz + sDir[j][2];
				if (jx < 0 || jx >= NX || jy < 0 || jy >= NY || jz < 0 || jz >= NZ)
					continue;

				int jidx = jx * NY * NZ + jy * NZ + jz;
				if (visibility[jidx] == 0 && visited[jidx] == 0) {
					visited[jidx] = component_num;
					dfs_stack.push(jidx);
				}
			}
		}
	}

	return component_num;
}

void PackingInfo::printTimeRecord()
{
	ofstream out("output/midOutput/TimeRecord.txt");
	out << "Time Record: " << std::endl;
	out << "Total Time: " << timeRecord.total << "s" << std::endl;
	out << " >> Read and Sort Objects time: " << timeRecord.read_and_sort << "s("
		<< timeRecord.read_and_sort / timeRecord.total * 100 << "%)" << std::endl;
	out << " >> Voxelization time: " << timeRecord.voxelization << "s("
		<< timeRecord.voxelization / timeRecord.total * 100 << "%)" << std::endl;
	out << " >> Collision Detection time: " << timeRecord.collision_detection << "s("
		<< timeRecord.collision_detection / timeRecord.total * 100 << "%)" << std::endl;
	out << " >> Compute Evaluation time: " << timeRecord.compute_evaluation << "s("
		<< timeRecord.compute_evaluation / timeRecord.total * 100 << "%)" << std::endl;
	out << " >> Placement time: " << timeRecord.placement << "s("
		<< timeRecord.placement / timeRecord.total * 100 << "%)" << std::endl;

	out << " >> Compute disassemble time: "<< timeRecord.compute_disassemble_order << "s("
		<< timeRecord.compute_disassemble_order / timeRecord.total * 100 << "%)" << std::endl;
	out << " >>>> Compute alpha shape: " << timeRecord.alpha_shape << "s("
		<< timeRecord.alpha_shape / timeRecord.compute_disassemble_order * 100 << "%)" << std::endl;
	out << " >>>> Export alpha shape: " << timeRecord.export_alpha_shape << "s("
		<< timeRecord.export_alpha_shape / timeRecord.compute_disassemble_order * 100 << "%)" << std::endl;

	out << " >>>> Sort Evaluated points time: " << timeRecord.sort_evaluated_points << "s("
		<< timeRecord.sort_evaluated_points / timeRecord.placement * 100 << "%)" << std::endl;

	out << " >>>> Compute NFV time: " << timeRecord.compute_NFV << "s("
		<< timeRecord.compute_NFV / timeRecord.placement * 100 << "%)" << std::endl;
	out << " >>>> Check Disassemble time: " << timeRecord.check_disassemble << "s("
		<< timeRecord.check_disassemble / timeRecord.placement * 100 << "%)" << std::endl;

	out << " >>>>>> Compute visibility time: " << timeRecord.compute_visibility << "s("
		<< timeRecord.compute_visibility / timeRecord.check_disassemble * 100 << "%)" << std::endl;
	out << " >>>>>> Compute connectedness time: " << timeRecord.compute_connectedness << "s("
		<< timeRecord.compute_connectedness / timeRecord.check_disassemble * 100 << "%)" << std::endl;

	out << " >>>>>> Compute recursion time: " << timeRecord.Recursion_time << "s("
		<< timeRecord.Recursion_time / timeRecord.check_disassemble * 100 << "%)" << std::endl;
	out << " >>>>>>>> Compute recursion NFV time: " << timeRecord.compute_NFV_Recursion << "s("
		<< timeRecord.compute_NFV_Recursion / timeRecord.Recursion_time * 100 << "%)" << std::endl;


	out << " >> Export Visible Result time: " << timeRecord.export_result << "s" << std::endl;
	out << " >> Tmp Time (equal to CheckDisassembleWithIVHT): " << timeRecord.tmp << "s" << std::endl << std::endl << std::endl;

	int Check_disassemble_tot_num = 0;
	out << "Check Points Time: ";
	for (int i = 0; i < timeRecord.check_disassemble_num.size(); i++) {
		out << timeRecord.check_disassemble_num[i] << " ";
		Check_disassemble_tot_num += timeRecord.check_disassemble_num[i];
	}
	out << std::endl;

	out << "Total Check Points Times: " << Check_disassemble_tot_num << std::endl;

	out << "Total Recursion Times: " << timeRecord.Recursion_num << std::endl;

	out.close();
}


void PackingInfo::ExportIVHT(ofstream& txtout)
{
	vector<int> connected_objects;
	for (int i = 0; i < visitedNum; i++) {
		if (positions[i].packed)
			connected_objects.push_back(i);
	}

	root = make_shared<TreeNode>();

	ResetUniqueId();

	ExportIVHTInfo(txtout, connected_objects, 0, root);

}

inline bool isSubsequence(const vector<int>& toMatch, const vector<int>& sub) {
	int n = toMatch.size();
	int m = sub.size();

	int i = 0, j = 0;

	while (i < n && j < m) {
		if (toMatch[i] == sub[j]) {
			j++;
		}
		i++;
	}

	return j == m;
}
inline void DiffSequence(const vector<int>& toMatch, const vector<int>& sub, vector<int>& difference) {

	int n = toMatch.size();
	int m = sub.size();

	int i = 0, j = 0;

	while (i < n && j < m) {
		if (toMatch[i] == sub[j]) {
			j++;
		}
		else {
			difference.push_back(toMatch[i]);
		}
		i++;
	}
	while (i < n) {
		difference.push_back(toMatch[i]);
		i++;
	}

	return;
}

void PackingInfo::ExportIVHTInfo(ofstream& txtout, const vector<int> connected_objects, const int layer, shared_ptr<TreeNode>& current)
{
	vector<int> visibility(TX, 0);
	vector<int> visited(TX, 0);

	txtout << "Layer " << layer << ", Unique ID " << unique_id << " : ";
	for (int i = 0; i < connected_objects.size(); i++) {
		txtout << connected_objects[i] << ", ";
	}
	txtout << std::endl;

	if (CombinationMap.find(connected_objects) == CombinationMap.end())   
	{
		txtout << endl << "The combination has NOT been calculated!" << std::endl;

		vector<vector<bool>> NFV_c(CutterNum, vector<bool>(NX * NY * NZ, false));

		vector<int> bestMatchKey;
		int maxMatchLength = 0;

		for (const auto& entry : CombinationMap) {
			const vector<int>& currentKey = entry.first;
			if (currentKey.size() >= connected_objects.size())
				continue;

			if (isSubsequence(connected_objects, currentKey))
			{
				int matchLength = currentKey.size();
				if (matchLength > maxMatchLength) {
					maxMatchLength = matchLength;
					bestMatchKey = currentKey;
				}
			}
		}
		auto key = CombinationMap[bestMatchKey];

		for (int i = 0; i < CutterNum; i++)
			copy(NFV_objects_comb[key.first][i].begin(), NFV_objects_comb[key.first][i].end(), NFV_c[i].begin());

		vector<int> difference;
		DiffSequence(connected_objects, bestMatchKey, difference);

		for (int i = 0; i < difference.size(); i++) {
			int objid = difference[i];

			if (!positions[objid].packed)
				txtout << "Combination Error: The object is not packed!" << std::endl;
			if (positions[objid].OriId == -1)
				txtout << "Combination Error: The object has not been placed!" << std::endl;

			vector<int> tmp;
			tmp.push_back(objid);

			key = CombinationMap[tmp];
			for (int j = 0; j < CutterNum; j++) {
				for (int k = 0; k < NFV_objects_comb[key.first][j].size(); k++) {
					if (NFV_objects_comb[key.first][j][k])
						NFV_c[j][k] = true;
				}
			}
		}

		NFV_objects_comb.push_back(NFV_c);
		CombinationMap[connected_objects] = make_pair(ConbinationNum, -1);   
		ConbinationNum++;
	}

	auto key = CombinationMap[connected_objects];
	int comb_num = key.first;

	if (!isSlient)
	{
		txtout << "objects: [";
		for (int i = 0; i < connected_objects.size(); i++)
			txtout << connected_objects[i] << ", ";
		txtout << "], compoent_num: " << key.second << " ;  ";
	}

	fill(visibility.begin(), visibility.end(), 0);
	fill(visited.begin(), visited.end(), 0);

	ComputeVis(NFV_objects_comb[comb_num], visibility);

	if (OutPut2DResult)
	{
		vector<bool> material(TX, false);  
		for (int m = 0; m < TX; m++)
			if (visibility[m] == 0)
				material[m] = true;
		Bool2Voxel("output/midOutput/IVHT/Material_" + std::to_string(layer) + "_" + std::to_string(unique_id) + "_visibility.obj",
			material, VoxelColor(0.6, 0.9, 0.7), VoxelSize(VOXELSIZE));
	}

	int comp_tmp = computeConnectedComponent(visibility, connected_objects, visited);


	current->id = unique_id;
	for (int i = 0; i < connected_objects.size(); i++) {
		current->data.push_back(connected_objects[i]);
	}
	current->imagePath = "F:/Project/Cuda_Pack/output/midOutput/IVHT/IVHT_" + std::to_string(layer) + "_" + std::to_string(unique_id) + "_visited.png";

	unique_id++;

	if (comp_tmp == 1)  
		return;

	vector<vector<int>> components(comp_tmp);
	for (int i = 0; i < visitedNum; i++) {
		if (!positions[i].packed) continue;

		int idx = positions[i].center[0] * NY * NZ + positions[i].center[1] * NZ + positions[i].center[2];
		if (visited[idx] != 0)
			components[visited[idx] - 1].push_back(i);
	}

	vector<int> UnremovedObjects;
	for (int i = 0; i < comp_tmp; i++) {
		if (components[i].empty())
			continue;

		if (components[i].size() == 1)  
		{
			shared_ptr<TreeNode> child = make_shared<TreeNode>();
			current->children.push_back(child);

			ExportIVHTInfo(txtout, components[i], layer + 1, child);
		}
		else
		{
			for (int j = 0; j < components[i].size(); j++)
				UnremovedObjects.push_back(components[i][j]);
		}
	}

	if (!UnremovedObjects.empty()) {
		if (UnremovedObjects.size() == connected_objects.size())
			return;  

		shared_ptr<TreeNode> child = make_shared<TreeNode>();
		current->children.push_back(child);
		sort(UnremovedObjects.begin(), UnremovedObjects.end());

		ExportIVHTInfo(txtout, UnremovedObjects, layer + 1, child);
	}
}

double mixedProduct(const Vector& v1, const Vector& v2, const Vector& v3) {
	return CGAL::cross_product(v2, v3) * v1;
}


float ExpHalfDecay(float x, float k = 0.1) {
	return 0.5 + 0.5 * std::exp(-k * x);
}

inline int index(int x, int y, int z, int ny, int nz) {
	return x * ny * nz + y * nz + z;
}

std::vector<std::vector<std::vector<float>>> create_gaussian_kernel_3d(int radius, float sigma, float center_scale = 1.0f) {
	int size = 2 * radius + 1;
	std::vector<std::vector<std::vector<float>>> kernel(size,
		std::vector<std::vector<float>>(size, std::vector<float>(size)));

	float sum = 0.0f;
	float sigma2 = sigma * sigma;
	int center = radius;

	for (int dx = -radius; dx <= radius; ++dx) {
		for (int dy = -radius; dy <= radius; ++dy) {
			for (int dz = -radius; dz <= radius; ++dz) {
				float value = std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0f * sigma2));

				if (dx == 0 && dy == 0 && dz == 0) {
					value *= center_scale;
				}

				kernel[dx + radius][dy + radius][dz + radius] = value;
				sum += value;
			}
		}
	}

	for (int dx = 0; dx < size; ++dx)
		for (int dy = 0; dy < size; ++dy)
			for (int dz = 0; dz < size; ++dz)
				kernel[dx][dy][dz] /= sum;

	return kernel;
}

std::vector<float> gaussian_filter_3d(
	const std::vector<float>& voxel_grid,
	int nx, int ny, int nz,
	float sigma = 1.0f,
	float center_scale = 1.0f     
) {
	int radius = static_cast<int>(std::ceil(2 * sigma));
	auto kernel = create_gaussian_kernel_3d(radius, sigma, center_scale);
	std::vector<float> output(voxel_grid.size(), 0.0f);

	for (int x = 0; x < nx; ++x) {
		for (int y = 0; y < ny; ++y) {
			for (int z = 0; z < nz; ++z) {
				float sum = 0.0f;

				for (int dx = -radius; dx <= radius; ++dx) {
					for (int dy = -radius; dy <= radius; ++dy) {
						for (int dz = -radius; dz <= radius; ++dz) {
							int xx = std::clamp(x + dx, 0, nx - 1);
							int yy = std::clamp(y + dy, 0, ny - 1);
							int zz = std::clamp(z + dz, 0, nz - 1);

							float weight = kernel[dx + radius][dy + radius][dz + radius];
							float value = voxel_grid[index(xx, yy, zz, ny, nz)];
							sum += weight * value;
						}
					}
				}

				output[index(x, y, z, ny, nz)] = sum;
			}
		}
	}

	return output;
}


void extractIosSurface(const std::vector<float>& voxel_grid, const string &filename, 
	double iso_value = 0.0, float sigma = 1.0f) {
	Eigen::MatrixXd V;   
	Eigen::MatrixXi F;   

	Eigen::MatrixXd S(NX * NY * NZ, 3);     
	Eigen::VectorXd D(NX * NY * NZ);      

	std::vector<float> blurred = gaussian_filter_3d(voxel_grid, NX, NY, NZ, sigma);

	for (int x = 0; x < NX; x++) {
		for (int y = 0; y < NY; y++) {
			for (int z = 0; z < NZ; z++) {
				int idx = x * NY * NZ + y * NZ + z;
				S(idx, 0) = x;
				S(idx, 1) = y;
				S(idx, 2) = z;
				D(idx) = blurred[idx];   
			}
		}
	}
	
	igl::marching_cubes(D, S, NX, NY, NZ, iso_value, V, F);

	Eigen::MatrixXd V_scaled = VOXELSIZE * V;

	if (!igl::writeOBJ(filename, V_scaled, F)) {
		std::cerr << "Failed to write mesh to " << filename << std::endl;
		return;
	}

	std::cout << "Extracted and saved mesh to " << filename << std::endl;

}


void ToolPathSimulation(const vector<int>& label_input, vector<int>& Result, 
	const vector<Vector>& sample_ori,
	const vector<vector<bool>>& non_reachability_map
) {
	for (int i = 0; i < TX; i++) {
		if(!label_input[i])
			continue;

		int maching_label = label_input[i] - 1;
		if (!non_reachability_map[maching_label][i])  
		{
			int source_idx = i;

			int sx, sy, sz;
			sx = source_idx / (NY * NZ);
			sy = (source_idx / NZ) % NY;
			sz = source_idx % NZ;

			double dir_x = sample_ori[maching_label][0];
			double dir_y = sample_ori[maching_label][1];
			double dir_z = sample_ori[maching_label][2];

			int step_x = (dir_x > 0) ? 1 : (dir_x < 0) ? -1 : 0;
			int step_y = (dir_y > 0) ? 1 : (dir_y < 0) ? -1 : 0;
			int step_z = (dir_z > 0) ? 1 : (dir_z < 0) ? -1 : 0;

			double t_max_x = (step_x == 0)
				? std::numeric_limits<double>::infinity()
				: (step_x > 0
					? 1.0 / dir_x
					: 0.0);

			double t_max_y = (step_y == 0)
				? std::numeric_limits<double>::infinity()
				: (step_y > 0
					? 1.0 / dir_y
					: 0.0);

			double t_max_z = (step_z == 0)
				? std::numeric_limits<double>::infinity()
				: (step_z > 0
					? 1.0 / dir_z
					: 0.0);

			double t_delta_x = (dir_x == 0)
				? std::numeric_limits<double>::infinity()
				: 1.0 / std::abs(dir_x);

			double t_delta_y = (dir_y == 0)
				? std::numeric_limits<double>::infinity()
				: 1.0 / std::abs(dir_y);

			double t_delta_z = (dir_z == 0)
				? std::numeric_limits<double>::infinity()
				: 1.0 / std::abs(dir_z);


			int ox = sx;
			int oy = sy;
			int oz = sz;

			while (ox >= 0 && ox < NX && oy >= 0 && oy < NY && oz >= 0 && oz < NZ) {
				int idx = ox * NY * NZ + oy * NZ + oz;

				if (non_reachability_map[maching_label][idx])
					break;

				Result[idx] = maching_label + 1;

				if (t_max_x < t_max_y && t_max_x < t_max_z) {
					ox += step_x;
					t_max_x += t_delta_x;
				}
				else if (t_max_y < t_max_z) {
					oy += step_y;
					t_max_y += t_delta_y;
				}
				else {
					oz += step_z;
					t_max_z += t_delta_z;
				}
			}
		}
	}
}



void FloodFill(
	const std::vector<bool>& UnAccessible,  
	const std::vector<int>& obstacle,  
	std::vector<bool>& kernal, 
	vector<int> &inner_surface) {

	const int sDir[26][3] = {
			{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
			{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
			{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
			{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
			{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
			{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };
	
	stack<int> flood_stack;
	vector<bool> flood_visited(TX, false);

	for (int i = 0; i < TX; i++) {
		if (kernal[i])
		{
			flood_stack.push(i);
			flood_visited[i] = 1;  	
		}
	}

	while (!flood_stack.empty()) {
		int idx = flood_stack.top();
		flood_stack.pop(); 

		int ix, iy, iz;
		ix = idx / (NY * NZ);
		iy = (idx / NZ) % NY;
		iz = idx % NZ;
		for (int j = 0; j < 6; j++)
		{
			int jx = ix + sDir[j][0];
			int jy = iy + sDir[j][1];
			int jz = iz + sDir[j][2];
			if (jx < 0 || jx >= NX || jy < 0 || jy >= NY || jz < 0 || jz >= NZ)
				continue;

			int jidx = jx * NY * NZ + jy * NZ + jz;

			if (obstacle[jidx])
				inner_surface[jidx] = obstacle[jidx];  

			if (!obstacle[jidx] && !flood_visited[jidx] && !UnAccessible[jidx]) {
				flood_visited[jidx] = 1;  
				flood_stack.push(jidx);
				kernal[jidx] = true;   
			}
		}
	}
}

void FloodFill(
	const std::vector<bool>& UnAccessible,  
	const std::vector<int>& obstacle,  
	const std::vector<bool>& remain_material,  
	std::vector<bool>& kernal,
	vector<int>& inner_surface) {

	const int sDir[26][3] = {
			{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
			{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
			{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
			{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
			{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
			{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };

	stack<int> flood_stack;
	vector<bool> flood_visited(TX, false);

	for (int i = 0; i < TX; i++) {
		if (kernal[i])
		{
			flood_stack.push(i);
			flood_visited[i] = 1;  	
		}
	}


	while (!flood_stack.empty()) {
		int idx = flood_stack.top();
		flood_stack.pop(); 

		int ix, iy, iz;
		ix = idx / (NY * NZ);
		iy = (idx / NZ) % NY;
		iz = idx % NZ;
		for (int j = 0; j < 6; j++)
		{
			int jx = ix + sDir[j][0];
			int jy = iy + sDir[j][1];
			int jz = iz + sDir[j][2];
			if (jx < 0 || jx >= NX || jy < 0 || jy >= NY || jz < 0 || jz >= NZ)
				continue;

			int jidx = jx * NY * NZ + jy * NZ + jz;

			if (obstacle[jidx] && obstacle[jidx] != CutterNum + 2)
				inner_surface[jidx] = obstacle[jidx];  

			if (remain_material[jidx] && !obstacle[jidx] && !flood_visited[jidx] && !UnAccessible[jidx]) {
				flood_visited[jidx] = 1;  
				flood_stack.push(jidx);
				kernal[jidx] = true;   
			}
		}
	}
}

	void SetMachiningDirection(const vector<bool>& Input, vector<int>& ResultWithLabel, 
		const vector<vector<bool>>& non_reachability_map, const vector<Vector>& sample_ori,
		const vector<bool>& remain_material, double weight_d = 1.0, double weight_s = 1.0) {

		const int sDir[26][3] = {
				{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
				{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
				{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
				{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
				{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
				{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };

		int num_pixels = 0;
		int num_labels = CutterNum;  

		vector<int> pixels;  
		map<int, int> pixel2index;  

		for (int i = 0; i < TX; i++)
			if (Input[i])   
			{
				num_pixels++;
				pixels.push_back(i);
			}
		for (int i = 0; i < num_pixels; i++)
			pixel2index[pixels[i]] = i;

		int* result = new int[num_pixels];       

		int* data = new int[num_pixels * num_labels];
		for (int i = 0; i < num_pixels; i++) {
			int pixel_idx = pixels[i];
			for (int l = 0; l < num_labels; l++) {
				if (non_reachability_map[l][pixel_idx])  
				{
					data[i * num_labels + l] = 1000000;
					continue;
				}

				float penal_count = 0, tot_length = 0;
				int source_idx = pixel_idx;

				int sx, sy, sz;
				sx = source_idx / (NY * NZ);
				sy = (source_idx / NZ) % NY;
				sz = source_idx % NZ;

				double dir_x = sample_ori[l][0];
				double dir_y = sample_ori[l][1];
				double dir_z = sample_ori[l][2];

				int step_x = (dir_x > 0) ? 1 : (dir_x < 0) ? -1 : 0;
				int step_y = (dir_y > 0) ? 1 : (dir_y < 0) ? -1 : 0;
				int step_z = (dir_z > 0) ? 1 : (dir_z < 0) ? -1 : 0;

				double t_max_x = (step_x == 0)
					? std::numeric_limits<double>::infinity()
					: (step_x > 0
						? 1.0 / dir_x
						: 0.0);

				double t_max_y = (step_y == 0)
					? std::numeric_limits<double>::infinity()
					: (step_y > 0
						? 1.0 / dir_y
						: 0.0);

				double t_max_z = (step_z == 0)
					? std::numeric_limits<double>::infinity()
					: (step_z > 0
						? 1.0 / dir_z
						: 0.0);

				double t_delta_x = (dir_x == 0)
					? std::numeric_limits<double>::infinity()
					: 1.0 / std::abs(dir_x);

				double t_delta_y = (dir_y == 0)
					? std::numeric_limits<double>::infinity()
					: 1.0 / std::abs(dir_y);

				double t_delta_z = (dir_z == 0)
					? std::numeric_limits<double>::infinity()
					: 1.0 / std::abs(dir_z);


				int ox = sx;
				int oy = sy;
				int oz = sz;

				while (ox >= 0 && ox < NX && oy >= 0 && oy < NY && oz >= 0 && oz < NZ) {
					int idx = ox * NY * NZ + oy * NZ + oz;

					tot_length += 1.0;
					if (!Input[idx])
						penal_count += 1.0;

					if (t_max_x < t_max_y && t_max_x < t_max_z) {
						ox += step_x;
						t_max_x += t_delta_x;
					}
					else if (t_max_y < t_max_z) {
						oy += step_y;
						t_max_y += t_delta_y;
					}
					else {
						oz += step_z;
						t_max_z += t_delta_z;
					}
				}
				data[i * num_labels + l] = (penal_count + 1.0) / (tot_length + 1.0) * 10;   
			}
		}

		int* smooth = new int[num_labels * num_labels];
		for (int l1 = 0; l1 < num_labels; l1++)
			for (int l2 = 0; l2 < num_labels; l2++)
				if (abs(l2 - l1) == 0)
					smooth[l1 + l2 * num_labels] = 0;
				else
					smooth[l1 + l2 * num_labels] = 1;

		for (int i = 0; i < num_pixels * num_labels; i++)
			data[i] = data[i] * (int)weight_d;

		for(int i = 0;i< num_labels * num_labels;i++)
			smooth[i] = smooth[i] * (int)weight_s;

		try {
			GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(num_pixels, num_labels);
			gc->setDataCost(data);
			gc->setSmoothCost(smooth);
			for (int p = 0; p < num_pixels; p++) {
				int cur_idx = pixels[p];
				for (int ni = 0; ni < 26; ni++) {
					float dirx = sDir[ni][0];
					float diry = sDir[ni][1];
					float dirz = sDir[ni][2];

					int cur_x = cur_idx / (NY * NZ);
					int cur_y = (cur_idx / NZ) % NY;
					int cur_z = cur_idx % NZ;

					int nei_x = cur_x + dirx;
					int nei_y = cur_y + diry;
					int nei_z = cur_z + dirz;

					if (nei_x < 0 || nei_x >= NX || nei_y < 0 || nei_y >= NY || nei_z < 0 || nei_z >= NZ)
						continue;

					int nei_idx = nei_x * NY * NZ + nei_y * NZ + nei_z;

					if (!Input[nei_idx])
						continue;

					if (cur_idx < nei_idx)  
					{
						int weight = 100;  

						Vector cur_dir(dirx, diry, dirz);
						for (int c = 0; c < CutterNum; c++) {
							Vector sample_dir = sample_ori[c];
							Vector cross = CGAL::cross_product(sample_dir, cur_dir);
							if (sqrt(cross.squared_length()) < 1e-6 && !non_reachability_map[c][cur_idx] && !non_reachability_map[c][nei_idx]) {  
								weight = 100;  
								break;
							}
						}
						int nei = pixel2index[nei_idx];
						gc->setNeighbors(p, nei, weight);
					}
				}
			}

			printf("\nBefore optimization energy is %d", gc->compute_energy());
			gc->expansion(2);         
			printf("\nAfter optimization energy is %d\n", gc->compute_energy());

			for (int i = 0; i < num_pixels; i++)
				result[i] = gc->whatLabel(i);

			delete gc;
		}
		catch (GCException e) {
			e.Report();
		}

		for (int i = 0; i < num_pixels; i++) {
			int label = result[i];
			ResultWithLabel[pixels[i]] = label + 1;
		}

		delete[] data;
		delete[] smooth;
		delete[] result;
	}

void ToolPathExtraction(const vector<int>& outer_surface, const vector<Vector>& sample_ori,const int& file_index) {

	vector<vector<int>> direction_count(CutterNum);
	for (int i = 0; i < TX; i++) {
		if (outer_surface[i] > 0) {
			direction_count[outer_surface[i] - 1].push_back(i);
		}
	}
	for (int c = 0; c < CutterNum; c++) {
		if (direction_count[c].size() > 0) {
			vector<bool> is_machining(TX, 0);
			vector<int> machining_path(TX, 0);
			for (int i = 0; i < direction_count[c].size(); i++) {
				is_machining[direction_count[c][i]] = 1;
				machining_path[direction_count[c][i]] = c + 1;
			}

			Int2VoxelUseColorBar("output/midOutput/IVHT/0-label-" + to_string(file_index) + "-" + to_string(c) + ".obj", machining_path, VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));

			vector<float> iosField(TX, 0);

			vbm::DistanceFunction(NX, NY, NZ, is_machining, iosField);

			for(int i = 0; i < TX; i++)
				iosField[i] = -iosField[i];

			extractIosSurface(iosField, "output/midOutput/IVHT/11-surface-" + to_string(file_index) + "-" + to_string(c) + ".obj", -1.0, 0.5);
		}
	}
}

void ToolPathExtraction(const vector<bool> & remaing_material, const vector<int>& outer_surface, const vector<Vector>& sample_ori, const int& file_index) {

	vector<vector<int>> direction_count(CutterNum);
	for (int i = 0; i < TX; i++) {
		if (outer_surface[i] > 0) {
			direction_count[outer_surface[i] - 1].push_back(i);
		}
	}
	for (int c = 0; c < CutterNum; c++) {
		if (direction_count[c].size() > 0) {
			vector<bool> is_machining(TX, 0);
			vector<int> machining_path(TX, 0);
			for (int i = 0; i < direction_count[c].size(); i++) {
				is_machining[direction_count[c][i]] = 1;
				if (remaing_material[direction_count[c][i]])
					machining_path[direction_count[c][i]] = c + 1;
			}

			Int2VoxelUseColorBar("output/midOutput/IVHT/0-label-" + to_string(file_index) + "-" + to_string(c) + ".obj", machining_path, VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));

			vector<float> iosField(TX, 0);

			vbm::DistanceFunction(NX, NY, NZ, is_machining, iosField);

			for (int i = 0; i < TX; i++)
				iosField[i] = -iosField[i];

			extractIosSurface(iosField, "output/midOutput/IVHT/11-surface-" + to_string(file_index) + "-" + to_string(c) + ".obj", -1.0, 0.5);
		}
	}
}





void VoronoiBasedInitSolverDetailed(
	const std::vector<bool>& object_A,
	const std::vector<bool>& object_B,
	const std::vector<std::vector<bool>>& non_reachability_map,
	const vector<Vector>& sample_ori,
	vector<bool>& remain_material,  
	int file_index 
) {

	vector<bool> skeleton(TX, false);

	for(int ix = 0; ix < NX; ix++)
		for(int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++)
			{
				if (((ix == 0 && iy == 0) || (ix == NX - 1 && iy == 0) || (ix == 0 && iy == NY - 1) || (ix == NX - 1 && iy == NY - 1) ||
					(iz == 0 && iy == 0) || (iz == NZ - 1 && iy == 0) || (iz == 0 && iy == NY - 1) || (iz == NZ - 1 && iy == NY - 1) ||
					(ix == 0 && iz == 0) || (ix == NX - 1 && iz == 0) || (ix == 0 && iz == NZ - 1) || (ix == NX - 1 && iz == NZ - 1))
					&& remain_material[ix * NY * NZ + iy * NZ + iz])
					skeleton[ix * NY * NZ + iy * NZ + iz] = true;
					
			}

	Bool2Voxel("output/midOutput/IVHT/0-BESTskeleton-" + to_string(file_index) + ".obj", skeleton, VoxelColor(0.9, 0.6, 0.8), VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));
	vector<float> a_field(TX, 1);
	for (int i = 0; i < TX; i++) {
		if (!remain_material[i])
			a_field[i] = -1;
	}
	extractIosSurface(a_field, "output/midOutput/IVHT/0-remain material-" + to_string(file_index) + ".obj");

	const int sDir[26][3] = {
			{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
			{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
			{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
			{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
			{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
			{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };
	vector<float> gradientX_map(TX, 0), gradientY_map(TX, 0), gradientZ_map(TX, 0);
	for (int i = 0; i < TX; i++) {
		float grad_x = 0, grad_y = 0, grad_z = 0;
		for (int c = 0; c < CutterNum; c++) {
			if (!non_reachability_map[c][i]) {
				grad_x += sample_ori[c][0];
				grad_y += sample_ori[c][1];
				grad_z += sample_ori[c][2];
			}
		}
		gradientX_map[i] = grad_x;
		gradientY_map[i] = grad_y;
		gradientZ_map[i] = grad_z;
	}

	vector<float> dist_to_A(TX, 0), dist_to_B(TX, 0);
	vbm::DistanceFunction(NX, NY, NZ, object_A, dist_to_A);
	vbm::DistanceFunction(NX, NY, NZ, object_B, dist_to_B);

	
	vector<bool> is_machining(TX, false);
	vector<float> DistanceField(TX, 0);
	for (int i = 0; i < TX; i++) {
		DistanceField[i] = dist_to_A[i] - dist_to_B[i];
		int x = i / NY / NZ;
		int y = (i / NZ) % NY;
		int z = i % NZ;
	}

	extractIosSurface(DistanceField, "output/midOutput/IVHT/1-VoronoiBoundarySurface-" + to_string(file_index) + ".obj");

	cout<< "extracted iso surface" << endl;

	vector<bool> BelongToVor(TX, false);   
	for (int i = 0; i < TX; i++) {
		float cur_x = min(dist_to_A[i], dist_to_B[i]);
		float err = 1.0;
		if (abs(dist_to_A[i] - dist_to_B[i]) <= err)
			BelongToVor[i] = true;
	}

	vector<bool> tmp_voronoi_boundary(TX, false);
	for (int i = 0; i < TX; i++) {
		if (BelongToVor[i] && remain_material[i])
			tmp_voronoi_boundary[i] = true;
	}
	Bool2Voxel("output/midOutput/IVHT/1-VoronoiBoundary-" + to_string(file_index) + ".obj", tmp_voronoi_boundary, VoxelColor(0.9, 0.6, 0.8), VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));

	vector<int> label_distribution(TX, 0);  
	SetMachiningDirection(BelongToVor, label_distribution, non_reachability_map, sample_ori, remain_material, 10, 1);

	Int2VoxelUseColorBar("output/midOutput/IVHT/2-label_distribution-" + std::to_string(file_index) + ".obj", label_distribution, VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));

	int IDX = 101;
	ToolPathExtraction(remain_material, label_distribution, sample_ori, IDX);
	IDX++;	cout << "extracted tool path" << endl;


	vector<int> ToolPathResult(TX, 0);  

	ToolPathSimulation(label_distribution, ToolPathResult, sample_ori, non_reachability_map);

	ToolPathExtraction(remain_material, ToolPathResult, sample_ori, IDX);
	IDX++;
	cout << "extracted tool path" << endl;
	Int2VoxelUseColorBar("output/midOutput/IVHT/3-ToolPathResult-" + std::to_string(file_index) + ".obj", ToolPathResult, VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));

	for (int i = 0; i < TX; i++) {
		if (ToolPathResult[i]) 
			is_machining[i] = true;
		if(!remain_material[i])
			is_machining[i] = false;
		int x = i / NY / NZ;
		int y = (i / NZ) % NY;
		int z = i % NZ;
		if (x == 1 || x == NX - 2 || y == 1 || y == NY - 2 || z == 1 || z == NZ - 2)
			is_machining[i] = false;
	}

	vbm::DistanceFunction(NX, NY, NZ, is_machining, DistanceField);

	for (int i = 0; i < TX; i++)
		DistanceField[i] = -DistanceField[i];

	extractIosSurface(DistanceField, "output/midOutput/IVHT/3.1-ToolPathVolume-" + to_string(file_index) + ".obj", -1.0, 0.5);

	int iter_count = 0;
	vector<bool> kernal(TX, false);
	copy(object_A.begin(), object_A.end(), kernal.begin());

	vector<int> inner_surface(TX, 0);  
	vector<int> outer_surface(TX, 0);  
	FloodFill(object_B, ToolPathResult, kernal, inner_surface);

	Int2VoxelUseColorBar("output/midOutput/IVHT/5-inner surface-" + std::to_string(file_index) + ".obj", inner_surface, VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));
	
	ToolPathExtraction(remain_material, inner_surface, sample_ori, IDX);
	IDX++;
	cout << "extracted tool path" << endl;

	vector<bool> kernal_B(TX, false);
	copy(object_B.begin(), object_B.end(), kernal_B.begin());

	while (iter_count < 2) {
		fill(ToolPathResult.begin(), ToolPathResult.end(), 0);
		ToolPathSimulation(inner_surface, ToolPathResult, sample_ori, non_reachability_map);


		ToolPathExtraction(remain_material, ToolPathResult, sample_ori, IDX);
		IDX++;
		cout << "extracted tool path" << endl;

		copy(object_B.begin(), object_B.end(), kernal_B.begin());
		fill(outer_surface.begin(), outer_surface.end(), 0);
		FloodFill(object_A, ToolPathResult, kernal_B, outer_surface);

		ToolPathExtraction(remain_material, outer_surface, sample_ori, IDX);
		IDX++;
		cout << "extracted tool path" << endl;

		fill(ToolPathResult.begin(), ToolPathResult.end(), 0);

		ToolPathSimulation(outer_surface, ToolPathResult, sample_ori, non_reachability_map);

		ToolPathExtraction(remain_material, ToolPathResult, sample_ori, IDX);
		IDX++;
		cout << "extracted tool path" << endl;

		copy(object_A.begin(), object_A.end(), kernal.begin());
		fill(inner_surface.begin(), inner_surface.end(), 0);
		FloodFill(object_B, ToolPathResult, kernal, inner_surface);

		ToolPathExtraction(remain_material, inner_surface, sample_ori, IDX);
		IDX++;
		cout << "extracted tool path" << endl;

		iter_count++;
	}


	copy(object_A.begin(), object_A.end(), kernal.begin());
	fill(inner_surface.begin(), inner_surface.end(), 0);
	FloodFill(object_B, ToolPathResult, remain_material, kernal, inner_surface);

	Bool2Voxel("output/midOutput/IVHT/6-kernal-" + to_string(file_index) + ".obj", kernal, VoxelColor(0.9, 0.6, 0.8), VoxelSize(ContainerSize / std::max(NX, std::max(NY, NZ))));



	vector<bool>  tmp_kernal(TX, 0);
	copy(kernal.begin(), kernal.end(), tmp_kernal.begin());
	FloodFill(object_B, ToolPathResult, tmp_kernal, inner_surface);

	fill(is_machining.begin(), is_machining.end(), false);
	for (int i = 0; i < TX; i++) {
		if (inner_surface[i])
			is_machining[i] = true;
	}
	fill(DistanceField.begin(), DistanceField.end(), 0);
	vbm::DistanceFunction(NX, NY, NZ, is_machining, DistanceField);
	for (int i = 0; i < TX; i++)
	{
		if(tmp_kernal[i])
			DistanceField[i] = -DistanceField[i];
	}
	
	extractIosSurface(DistanceField, "output/midOutput/IVHT/0-BESTSurface-" + to_string(file_index) + ".obj");

	cout<< "extracted best surface" << endl;

	vector<float> iosField(TX, -1);

	for (int i = 0; i < TX; i++) {
		if (kernal[i] && remain_material[i]) {
			iosField[i] = 1;
			for (int k = 0; k < 6; k++) {
				int nx = i / NY / NZ + sDir[k][0];
				int ny = (i / NZ) % NY + sDir[k][1];
				int nz = i % NZ + sDir[k][2];
				if (nx < 0 || nx >= NX || ny < 0 || ny >= NY || nz < 0 || nz >= NZ)
					continue;  

				int nei_idx = nx * NY * NZ + ny * NZ + nz;
				if (remain_material[nei_idx] && !object_B[nei_idx]) {
					iosField[nei_idx] = 1;
				}
			}
		}
	}



	extractIosSurface(iosField, "output/midOutput/IVHT/8-maching result-" + to_string(file_index) + ".obj");

	fill(iosField.begin(), iosField.end(), -1);
	for (int i = 0; i < TX; i++) {
		if ((kernal[i] || inner_surface[i]) && remain_material[i])
		{
			iosField[i] = 1;
			remain_material[i] = false;  

			for (int k = 0; k < 6; k++) {
				int nx = i / NY / NZ + sDir[k][0];
				int ny = (i / NZ) % NY + sDir[k][1];
				int nz = i % NZ + sDir[k][2];
				if (nx < 0 || nx >= NX || ny < 0 || ny >= NY || nz < 0 || nz >= NZ)
					continue;  

				int nei_idx = nx * NY * NZ + ny * NZ + nz;
				if (remain_material[nei_idx] && !object_B[nei_idx]) {
					iosField[nei_idx] = 1;
					remain_material[i] = false;  
				}
			}
		}
		if (ToolPathResult[i])
			remain_material[i] = false;  
	}

	extractIosSurface(iosField, "output/midOutput/IVHT/7-carving hull-" + to_string(file_index) + ".obj");



	fill(iosField.begin(), iosField.end(), -1);
	for (int i = 0; i < TX; i++)
		if (remain_material[i])
			iosField[i] = 1;

	extractIosSurface(iosField, "output/midOutput/IVHT/9-remain material-" + to_string(file_index) + ".obj");

	vector<int> ToolPathDirection(CutterNum, 0);
	for (int i = 0; i < TX; i++) {
		if (ToolPathResult[i])
			ToolPathDirection[ToolPathResult[i] - 1]++;
	}

	ofstream ofs("output/midOutput/IVHT/10-ToolPathDirection-" + to_string(file_index) + ".txt");
	for (int i = 0; i < CutterNum; i++) {
		if (ToolPathDirection[i] > 0)
			ofs << sample_ori[i][0] << " " << sample_ori[i][1] << " " << sample_ori[i][2] << endl;
	}
	ofs.close();

	
	for (int ix = 0; ix < NX; ix++)
		for (int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++) {

				if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1 || iz == 0 || iz == NZ - 1)
				{
					int idx = ix * NY * NZ + iy * NZ + iz;
					ToolPathResult[idx] = 0;

				}
			}

	ToolPathExtraction(ToolPathResult, sample_ori, file_index);
}






void PackingInfo::ComputeDisassembleOrder()
{
	cout << "ComputeDisassembleOrder..." << endl;

	disassembleOrder.clear();  
	int global_idx = 0;
	vector<bool> RemainingMaterial(TX, true);  

	for (int ix = 0; ix < NX; ix++)
		for (int iy = 0; iy < NY; iy++)
			for (int iz = 0; iz < NZ; iz++) {

				if (ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1 || iz == 0 || iz == NZ - 1)
					RemainingMaterial[ix * NY * NZ + iy * NZ + iz] = false;
			}
	
	if (ForFabrication) {
		for (int ix = 0; ix < NX; ix++)
			for (int iy = 0; iy < NY; iy++)
				for (int iz = 0; iz < NZ; iz++) {
					if (iz >= 47)
						RemainingMaterial[ix * NY * NZ + iy * NZ + iz] = false;

					if(ix == 0 || ix == NX - 1 || iy == 0 || iy == NY - 1 || iz == 0 || iz == 47)
						RemainingMaterial[ix * NY * NZ + iy * NZ + iz] = false;
				}
	}

	shared_ptr<TreeNode> current = root;
	if (root->data.size() == 1) {
		shared_ptr<TreeNode> child = make_shared<TreeNode>();
		child->data.push_back(root->data[0]);

		current->children.push_back(child);
	}

	while (!current->children.empty()) {  
		SAMPLE_ON_BALL sampleOri;
		sampleOri.num_ori_sample = CutterNum;
		sampleOri.EnumerationOritation();
		vector<Vector> sample_ori;
		for (int c = 0; c < CutterNum; c++)
		{
			Vector tmp_ori = sampleOri.sample_points[c];
			tmp_ori = -tmp_ori;
			sample_ori.push_back(tmp_ori);
		}

		vector<int> RemainingObjects;

		for (int i = 0; i < current->data.size(); i++)
			RemainingObjects.push_back(current->data[i]);

		if (CombinationMap.find(RemainingObjects) == CombinationMap.end()) {
			cout << endl << "The combination has NOT been calculated!" << std::endl;
			vector<vector<bool>> NFV_c(CutterNum, vector<bool>(NX * NY * NZ, false));

			vector<int> bestMatchKey;
			int maxMatchLength = 0;

			for (const auto& entry : CombinationMap) {
				const vector<int>& currentKey = entry.first;
				if (currentKey.size() >= RemainingObjects.size())
					continue;

				if (isSubsequence(RemainingObjects, currentKey))
				{
					int matchLength = currentKey.size();
					if (matchLength > maxMatchLength) {
						maxMatchLength = matchLength;
						bestMatchKey = currentKey;
					}
				}
			}
			auto key = CombinationMap[bestMatchKey];

			for (int i = 0; i < CutterNum; i++)
				copy(NFV_objects_comb[key.first][i].begin(), NFV_objects_comb[key.first][i].end(), NFV_c[i].begin());

			vector<int> difference;
			DiffSequence(RemainingObjects, bestMatchKey, difference);

			for (int i = 0; i < difference.size(); i++) {
				int objid = difference[i];

				if (!positions[objid].packed)
					cout << "Combination Error: The object is not packed!" << std::endl;
				if (positions[objid].OriId == -1)
					cout << "Combination Error: The object has not been placed!" << std::endl;

				vector<int> tmp;
				tmp.push_back(objid);

				key = CombinationMap[tmp];
				for (int j = 0; j < CutterNum; j++) {
					for (int k = 0; k < NFV_objects_comb[key.first][j].size(); k++) {
						if (NFV_objects_comb[key.first][j][k])
							NFV_c[j][k] = true;
					}
				}
			}

			NFV_objects_comb.push_back(NFV_c);
			CombinationMap[RemainingObjects] = make_pair(ConbinationNum, -1);   
			ConbinationNum++;
		}

		auto key = CombinationMap[RemainingObjects];
		int comb_num = key.first;

		vector<int> visibility(TX, 0);

		ComputeVis(NFV_objects_comb[comb_num], visibility);

		vector<bool> vis_copy(TX, false);
		for (int m = 0; m < TX; m++)
		{
			if (visibility[m] == 0)
				vis_copy[m] = true;
		}

		vector<int> toDisassemble;

		for (int k = 0; k < current->children.size(); k++) {
			if (current->children[k]->data.size() == 1)  
				toDisassemble.push_back(current->children[k]->data[0]);
		}

		sort(toDisassemble.begin(), toDisassemble.end(), [this](int a, int b) {
			return a > b;
			});

		cout << "The order of disassembly is: " << endl;
		for (int i = 0; i < toDisassemble.size(); i++) {
			cout << toDisassemble[i] << " ";
		}
		cout << endl;

		for (int tD = 0; tD < toDisassemble.size(); tD++) {

			cout << endl << endl << "Current object to disassemble: " << toDisassemble[tD] << "   " << global_idx << endl;

			if (CombinationMap.find(RemainingObjects) == CombinationMap.end()) {
				cout << endl << "The combination has NOT been calculated!" << std::endl;
				vector<vector<bool>> NFV_c(CutterNum, vector<bool>(NX * NY * NZ, false));

				vector<int> bestMatchKey;
				int maxMatchLength = 0;

				for (const auto& entry : CombinationMap) {
					const vector<int>& currentKey = entry.first;
					if (currentKey.size() >= RemainingObjects.size())
						continue;

					if (isSubsequence(RemainingObjects, currentKey)) {
						int matchLength = currentKey.size();
						if (matchLength > maxMatchLength) {
							maxMatchLength = matchLength;
							bestMatchKey = currentKey;
						}
					}
				}
				auto key = CombinationMap[bestMatchKey];

				for (int i = 0; i < CutterNum; i++) {
					copy(NFV_objects_comb[key.first][i].begin(), NFV_objects_comb[key.first][i].end(), NFV_c[i].begin());
				}

				vector<int> difference;
				DiffSequence(RemainingObjects, bestMatchKey, difference);

				for (int i = 0; i < difference.size(); i++) {
					int objid = difference[i];

					if (!positions[objid].packed)
						cout << "Combination Error: The object is not packed!" << std::endl;
					if (positions[objid].OriId == -1)
						cout << "Combination Error: The object has not been placed!" << std::endl;

					vector<int> tmp;
					tmp.push_back(objid);

					key = CombinationMap[tmp];
					for (int j = 0; j < CutterNum; j++) {
						for (int k = 0; k < NFV_objects_comb[key.first][j].size(); k++) {
							if (NFV_objects_comb[key.first][j][k])
								NFV_c[j][k] = true;
						}
					}
				}

				NFV_objects_comb.push_back(NFV_c);
				CombinationMap[RemainingObjects] = make_pair(ConbinationNum, -1);   
				ConbinationNum++;
			}

			key = CombinationMap[RemainingObjects];
			comb_num = key.first;

			fill(visibility.begin(), visibility.end(), 0);

			ComputeVis(NFV_objects_comb[comb_num], visibility);

			fill(vis_copy.begin(), vis_copy.end(), false);
			for (int m = 0; m < TX; m++)
			{
				if (visibility[m] == 0)
					vis_copy[m] = true;
			}
			int objid = toDisassemble[tD];
			int idx_center = positions[objid].center[0] * NY * NZ + positions[objid].center[1] * NZ + positions[objid].center[2];
			std::stack<int> dfs_stack;
			dfs_stack.push(idx_center);

			vector<int> visited(TX, 0);

			if (idx_center < 0 || idx_center >= TX)
				std::cout << "4Init stack Error: The start point is out of range!" << std::endl;
			if (visibility[idx_center] != 0)
			{
				std::cout << "4Init stack Error: The start point is Unvisible!" << std::endl;
			}

			const int sDir[26][3] = {
			{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
			{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
			{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
			{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
			{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
			{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };


			vector<int> internal_visual_hull, obstacles;
			visited[idx_center] = 1;  	
			internal_visual_hull.push_back(idx_center);

			while (!dfs_stack.empty()) {
				int idx = dfs_stack.top();
				dfs_stack.pop(); 

				int ix, iy, iz;
				ix = idx / (NY * NZ);
				iy = (idx / NZ) % NY;
				iz = idx % NZ;
				for (int j = 0; j < 26; j++)
				{
					int jx = ix + sDir[j][0];
					int jy = iy + sDir[j][1];
					int jz = iz + sDir[j][2];
					if (jx < 0 || jx >= NX || jy < 0 || jy >= NY || jz < 0 || jz >= NZ)
						continue;

					int jidx = jx * NY * NZ + jy * NZ + jz;
					if (visibility[jidx] == 0 && !visited[jidx]) {
						visited[jidx] = 1;  
						dfs_stack.push(jidx);
						internal_visual_hull.push_back(jidx);
					}
				}
			}

			vector<bool> to_remove(TX, false);  
			for (int ivh = 0; ivh < internal_visual_hull.size(); ivh++)
				to_remove[internal_visual_hull[ivh]] = true;
			vector<bool> to_remain(TX, false);  

			for (int i = 0; i < TX; i++) {
				if (visibility[i] == 0 && !to_remove[i])
				{ 
					to_remain[i] = true;
					obstacles.push_back(i);
				}
			}

			if (obstacles.size() != 0) {

				VoronoiBasedInitSolverDetailed(to_remove, to_remain, NFV_objects_comb[comb_num], sample_ori, RemainingMaterial, global_idx);

				cout << "middle test" << endl;
			}
			else {
				cout << "Remaining Material:" << endl;
				vector<int> Obstacle(TX, 0), inner_surface(TX, 0);
				for (int i = 0; i < TX; i++) {
					if (!RemainingMaterial[i])
						Obstacle[i] = 1;
				}

				FloodFill(to_remain, Obstacle, to_remove, inner_surface);

				vector<float> iosField(TX, -1);
				for (int i = 0; i < TX; i++) {
					if (to_remove[i] && RemainingMaterial[i])
					{
						iosField[i] = 1;
					}
				}
				extractIosSurface(iosField, "output/midOutput/IVHT/7-carving hull-" + to_string(global_idx) + ".obj");
			}

			cout << "Disassemble Success!" << endl;
			cout << endl;
			disassembleOrder.push_back(toDisassemble[tD]);
			for (int j = 0; j < RemainingObjects.size(); j++) {
				if (RemainingObjects[j] == toDisassemble[tD]) {
					RemainingObjects.erase(RemainingObjects.begin() + j);
					break;
				}
			}
			toDisassemble.erase(toDisassemble.begin() + tD);
			tD--;

			global_idx++;
		}


		int next_node = current->children.size() - 1;
		shared_ptr<TreeNode> next = current->children[next_node];
		if (next->data.size() > 1)
			current = next;  
		else
			break;  
	}
}

void PackingInfo::ExportDOT() {
	string dot_file = "output/midOutput/ivht.dot";

	std::ofstream file(dot_file);
	file << "digraph G {\n";
	file << "    rankdir=TB;\n";
	ExportIVHTNode(root, file);
	file << "}\n";
	file.close();
}

void PackingInfo::ExportIVHTNode(const shared_ptr<TreeNode>& node, std::ofstream& file) {
	if (!node) return;

	int nodeId = node->id;


	file << "    \"" << nodeId << "\" [label=<<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLPADDING=\"4\">\n";
	file << "        <TR>\n";
	file << "            <!-- image -->\n";
	file << "            <TD><IMG SRC=\"" << node->imagePath << "\" SCALE=\"TRUE\"/></TD>\n";
	file << "        </TR>\n";
	file << "        <TR>\n";
	file << "            <!-- data -->\n";
	file << "            <TD> Objects: "; for (int i = 0; i < node->data.size(); i++) file << node->data[i] << " "; file << "</TD>\n";
	file << "        </TR>\n";
	file << "    </TABLE>> shape=plaintext];\n";

	for (const auto& child : node->children) {
		file << "\"" << nodeId << "\" -> \"" << child->id << "\";\n";
		ExportIVHTNode(child, file);
	}
}

bool PackingInfo::CheckDisassembleWithIVHT(const int oriId, const int x, const int y, const int z, vector<int> connected_objects, int depth) {
	int curDepth;
#pragma omp critical  
	{
		curDepth = currentDepth;
	}

	if (depth > curDepth + 2)
		return false;


	clock_t start_tot_check_time = clock();

	clock_t start, end;

	vector<int> visibility(TX, 0);
	vector<int> visited(TX, 0);

	bool isExist = true;
#pragma omp critical    
	{
		if (CombinationMap.find(connected_objects) == CombinationMap.end())   
			isExist = false;
	}

	if (!isExist)   
	{
		start = clock();

		vector<vector<bool>> NFV_c(CutterNum, vector<bool>(NX * NY * NZ, false));
		vector<vector<bool>> NFV_tmp(CutterNum, vector<bool>(NX * NY * NZ, false));

		vector<int> bestMatchKey;
		int maxMatchLength = 0;

		for (const auto& entry : CombinationMap) {
			const vector<int>& currentKey = entry.first;
			if (currentKey.size() >= connected_objects.size())
				continue;

			if (isSubsequence(connected_objects, currentKey))
			{
				int matchLength = currentKey.size();
				if (matchLength > maxMatchLength) {
					maxMatchLength = matchLength;
					bestMatchKey = currentKey;
				}
			}
		}
		auto key = CombinationMap[bestMatchKey];

#pragma omp critical  
		{
			for (int i = 0; i < CutterNum; i++)
				copy(NFV_objects_comb[key.first][i].begin(), NFV_objects_comb[key.first][i].end(), NFV_c[i].begin());
		}

		vector<int> difference;
		DiffSequence(connected_objects, bestMatchKey, difference);


		for (int i = 0; i < difference.size(); i++) {
			int objid = difference[i];

			if (!positions[objid].packed)
				std::cout << "Combination Error: The object is not packed!" << std::endl;

			vector<int> tmp;
			tmp.push_back(objid);

			key = CombinationMap[tmp];

#pragma omp critical  
			{
				for (int j = 0; j < CutterNum; j++) {
					copy(NFV_objects_comb[key.first][j].begin(), NFV_objects_comb[key.first][j].end(), NFV_tmp[j].begin());
				}
			}

			for (int j = 0; j < CutterNum; j++) {
				for (int k = 0; k < NFV_tmp[j].size(); k++) {
					if (NFV_tmp[j][k])
						NFV_c[j][k] = true;
				}
			}
		}

		end = clock();
		timeRecord.compute_NFV_Recursion += (double)(end - start) / CLOCKS_PER_SEC;


#pragma omp critical
		{
			if (CombinationMap.find(connected_objects) == CombinationMap.end())
			{
				NFV_objects_comb.push_back(NFV_c);
				CombinationMap[connected_objects] = make_pair(ConbinationNum, -1);   
				ConbinationNum++;
			}
		}
	}

	auto key = CombinationMap[connected_objects];
	int comb_num = key.first;

	if (!isSlient)
	{
		std::cout << "objects: [";
		for (int i = 0; i < connected_objects.size(); i++)
			std::cout << connected_objects[i] << ", ";
		std::cout << "], compoent_num: " << key.second << " ;  ";
	}


	start = clock();
	fill(visibility.begin(), visibility.end(), 0);
	fill(visited.begin(), visited.end(), 0);

	ComputeVis(oriId, comb_num, x, y, z, visibility);
	end = clock();
	if (unique_id == 0) 
		timeRecord.compute_visibility += (double)(end - start) / CLOCKS_PER_SEC;

	start = clock();
	const int sDir[26][3] = {
		{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
		{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
		{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
		{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
		{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
		{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };

	std::stack<int> dfs_stack;
	std::stack<int> init_stack;

	int component_num = 0;

	int b_center = -1;

	for (int i = 0; i < TX; i++) {
		if ((*objects)[visitedNum][oriId][i]) {
			int start_idx = 0;
			int ix, iy, iz;
			ix = i / (NY * NZ);
			iy = i / NZ % NY;
			iz = i % NZ;
			start_idx = (x + ix) % NX * NY * NZ + (y + iy) % NY * NZ + (z + iz) % NZ;
			b_center = start_idx;

			if (start_idx < 0 || start_idx >= NX * NY * NZ)
				std::cout << "1Init stack Error: The start point is out of range!" << std::endl;
			if (visited[start_idx] != 0)
				std::cout << "1Init stack Error: The start point is Unvisible!" << std::endl;

			init_stack.push(start_idx);
			break;
		}
	}

	while (!init_stack.empty()) {
		int idx_ = init_stack.top();
		init_stack.pop();

		if (visited[idx_] != 0)
			continue;

		dfs_stack.push(idx_);
		visited[idx_] = ++component_num;


		while (!dfs_stack.empty()) {
			int idx = dfs_stack.top();
			dfs_stack.pop();

			int ix, iy, iz;
			ix = idx / (NY * NZ);
			iy = (idx / NZ) % NY;
			iz = idx % NZ;
			for (int j = 0; j < 26; j++)
			{
				int jx = ix + sDir[j][0];
				int jy = iy + sDir[j][1];
				int jz = iz + sDir[j][2];
				if (jx < 0 || jx >= NX || jy < 0 || jy >= NY || jz < 0 || jz >= NZ)
					continue;

				int jidx = jx * NY * NZ + jy * NZ + jz;
				if (visibility[jidx] == 0 && visited[jidx] == 0) {
					visited[jidx] = component_num;
					dfs_stack.push(jidx);
				}
			}
		}
	}

	end = clock();
	if (unique_id == 0) 
		timeRecord.compute_connectedness += (double)(end - start) / CLOCKS_PER_SEC;


	if (!isSlient) {
		std::cout << "objects: [";
		for (int i = 0; i < connected_objects.size(); i++)
			std::cout << connected_objects[i] << ", ";
		std::cout << "Cur_Obj], compoent_num: " << component_num << " ;" << std::endl;
	}


	int current_component = visited[b_center];

	vector<int> obj_connected;

	for (int i = 0; i < visitedNum; i++)
	{
		if (!positions[i].packed) continue;

		int idx = positions[i].center[0] * NY * NZ + positions[i].center[1] * NZ + positions[i].center[2];
		if (visited[idx] != 0 && visited[idx] == current_component)
			obj_connected.push_back(i);
	}

	clock_t end_tot_check_time = clock();

	if (unique_id != 0) {
		timeRecord.Recursion_time += (double)(end_tot_check_time - start_tot_check_time) / CLOCKS_PER_SEC;
	}


	if (obj_connected.size() == 0)   
	{
#pragma omp critical  
		{
			if (depth > currentDepth)
				currentDepth = depth;
		}
		return true;
	}
	else {
		if (obj_connected.size() == connected_objects.size())   
		{
			return false;
		}
		else   
		{
			timeRecord.Recursion_num++;
			unique_id++;
			bool res = CheckDisassembleWithIVHT(oriId, x, y, z, obj_connected, depth + 1);
			return res;
		}
	}

	return false;
}


bool PackingInfo::CheckDisassembleWithIVHTLoop(const int oriId, const int x, const int y, const int z, vector<int> connected_objects, int depth){
	int curDepth;

#pragma omp critical  
	{
		curDepth = currentDepth;
	}

	vector<vector<pair<int, int>>> disassembled_objects;  
	vector<vector<int>> visibility_layers;   
	vector<vector<int>> visited_layers;   
	vector<int> component_nums;   

	bool is_disassemble = false;   

	while (true) {

		if (depth > curDepth + 2)
			return false;

		clock_t start_tot_check_time = clock();

		clock_t start, end;

		vector<int> visibility(TX, 0);
		vector<int> visited(TX, 0);

		bool isExist = true;
#pragma omp critical    
		{
			if (CombinationMap.find(connected_objects) == CombinationMap.end())   
				isExist = false;
		}

		if (!isExist)   
		{
			start = clock();

			vector<vector<bool>> NFV_c(CutterNum, vector<bool>(NX * NY * NZ, false));
			vector<vector<bool>> NFV_tmp(CutterNum, vector<bool>(NX * NY * NZ, false));

			vector<int> bestMatchKey;
			int maxMatchLength = 0;

			for (const auto& entry : CombinationMap) {
				const vector<int>& currentKey = entry.first;
				if (currentKey.size() >= connected_objects.size())
					continue;

				if (isSubsequence(connected_objects, currentKey))
				{
					int matchLength = currentKey.size();
					if (matchLength > maxMatchLength) {
						maxMatchLength = matchLength;
						bestMatchKey = currentKey;
					}
				}
			}
			auto key = CombinationMap[bestMatchKey];

#pragma omp critical  
			{
				for (int i = 0; i < CutterNum; i++)
					copy(NFV_objects_comb[key.first][i].begin(), NFV_objects_comb[key.first][i].end(), NFV_c[i].begin());
			}

			vector<int> difference;
			DiffSequence(connected_objects, bestMatchKey, difference);


			for (int i = 0; i < difference.size(); i++) {
				int objid = difference[i];

				if (!positions[objid].packed)
					std::cout << "Combination Error: The object is not packed!" << std::endl;

				vector<int> tmp;
				tmp.push_back(objid);

				key = CombinationMap[tmp];

#pragma omp critical  
				{
					for (int j = 0; j < CutterNum; j++) {
						copy(NFV_objects_comb[key.first][j].begin(), NFV_objects_comb[key.first][j].end(), NFV_tmp[j].begin());
					}
				}

				for (int j = 0; j < CutterNum; j++) {
					for (int k = 0; k < NFV_tmp[j].size(); k++) {
						if (NFV_tmp[j][k])
							NFV_c[j][k] = true;
					}
				}
			}

			end = clock();
			timeRecord.compute_NFV_Recursion += (double)(end - start) / CLOCKS_PER_SEC;


#pragma omp critical
			{
				if (CombinationMap.find(connected_objects) == CombinationMap.end())
				{
					NFV_objects_comb.push_back(NFV_c);
					CombinationMap[connected_objects] = make_pair(ConbinationNum, -1);   
					ConbinationNum++;
				}
			}
		}

		auto key = CombinationMap[connected_objects];
		int comb_num = key.first;

		if (!isSlient)
		{
			std::cout << "objects: [";
			for (int i = 0; i < connected_objects.size(); i++)
				std::cout << connected_objects[i] << ", ";
			std::cout << "], compoent_num: " << key.second << " ;  ";
		}


		start = clock();
		fill(visibility.begin(), visibility.end(), 0);
		fill(visited.begin(), visited.end(), 0);

		ComputeVis(oriId, comb_num, x, y, z, visibility);
		end = clock();
		if (unique_id == 0) 
			timeRecord.compute_visibility += (double)(end - start) / CLOCKS_PER_SEC;

		start = clock();
		const int sDir[26][3] = {
			{1,  0,  0},  {-1,  0,  0},  {0,  1,  0}, { 0, -1,  0},  {0,  0,  1},  {0,  0, -1},     
			{1,  1,  0},  {-1,  1,  0},  {1, -1,  0}, {-1, -1,  0},         
			{1,  0,  1},  {-1,  0,  1},  {1,  0, -1}, {-1,  0, -1},         
			{0,  1,  1},  { 0, -1,  1},  {0,  1, -1}, { 0, -1, -1},         
			{1,  1,  1},  {-1,  1,  1},  {1, -1,  1}, {-1, -1,  1},              
			{1,  1, -1},  {-1,  1, -1},  {1, -1, -1}, {-1, -1, -1} };

		std::stack<int> dfs_stack;
		std::stack<int> init_stack;

		int component_num = 0;

		int b_center = -1;

		for (int i = 0; i < TX; i++) {
			if ((*objects)[visitedNum][oriId][i]) {
				int start_idx = 0;
				int ix, iy, iz;
				ix = i / (NY * NZ);
				iy = i / NZ % NY;
				iz = i % NZ;
				start_idx = (x + ix) % NX * NY * NZ + (y + iy) % NY * NZ + (z + iz) % NZ;
				b_center = start_idx;

				if (start_idx < 0 || start_idx >= NX * NY * NZ)
					std::cout << "1Init stack Error: The start point is out of range!" << std::endl;
				if (visited[start_idx] != 0)
					std::cout << "1Init stack Error: The start point is Unvisible!" << std::endl;

				init_stack.push(start_idx);
				break;
			}
		}

		for (int i = 0; i < connected_objects.size(); i++) {
			int x = positions[connected_objects[i]].center[0];
			int y = positions[connected_objects[i]].center[1];
			int z = positions[connected_objects[i]].center[2];

			init_stack.push(x * NY * NZ + y * NZ + z);
			if (x * NY * NZ + y * NZ + z < 0 || x * NY * NZ + y * NZ + z >= NX * NY * NZ)
				std::cout << "2Init stack Error: The start point is out of range!" << std::endl;
			if (visibility[x * NY * NZ + y * NZ + z] != 0)
				std::cout << "2Init stack Error: The start point is Unvisible!" << std::endl;
		}

		while (!init_stack.empty()) {
			int idx_ = init_stack.top();
			init_stack.pop();

			if (visited[idx_] != 0)
				continue;

			dfs_stack.push(idx_);
			visited[idx_] = ++component_num;


			while (!dfs_stack.empty()) {
				int idx = dfs_stack.top();
				dfs_stack.pop();

				int ix, iy, iz;
				ix = idx / (NY * NZ);
				iy = (idx / NZ) % NY;
				iz = idx % NZ;
				for (int j = 0; j < 26; j++)
				{
					int jx = ix + sDir[j][0];
					int jy = iy + sDir[j][1];
					int jz = iz + sDir[j][2];
					if (jx < 0 || jx >= NX || jy < 0 || jy >= NY || jz < 0 || jz >= NZ)
						continue;

					int jidx = jx * NY * NZ + jy * NZ + jz;
					if (visibility[jidx] == 0 && visited[jidx] == 0) {
						visited[jidx] = component_num;
						dfs_stack.push(jidx);
					}
				}
			}
		}

		end = clock();
		if (unique_id == 0) 
			timeRecord.compute_connectedness += (double)(end - start) / CLOCKS_PER_SEC;


		if (!isSlient) {
			std::cout << "objects: [";
			for (int i = 0; i < connected_objects.size(); i++)
				std::cout << connected_objects[i] << ", ";
			std::cout << "Cur_Obj], compoent_num: " << component_num << " ;" << std::endl;
		}


		visibility_layers.push_back(visibility);
		visited_layers.push_back(visited);
		component_nums.push_back(component_num);

		int current_component = visited[b_center];

		vector<int> obj_connected;

		for (int i = 0; i < visitedNum; i++)
		{
			if (!positions[i].packed) continue;

			int idx = positions[i].center[0] * NY * NZ + positions[i].center[1] * NZ + positions[i].center[2];
			if (visited[idx] != 0 && visited[idx] == current_component)
				obj_connected.push_back(i);
		}

		clock_t end_tot_check_time = clock();

		if (unique_id != 0) {
			timeRecord.Recursion_time += (double)(end_tot_check_time - start_tot_check_time) / CLOCKS_PER_SEC;
		}


		if (obj_connected.size() == 0)   
		{
#pragma omp critical  
			{
				if (depth > currentDepth)
					currentDepth = depth;
			}
			is_disassemble = true;
			vector<vector<int>> components(component_num);
			for (int i = 0; i < visitedNum; i++) {
				if (!positions[i].packed) continue;

				int idx = positions[i].center[0] * NY * NZ + positions[i].center[1] * NZ + positions[i].center[2];
				if (visited[idx] != 0)
					components[visited[idx] - 1].push_back(i);
			}

			vector<pair<int, int>> disassembledobject;
			for (int i = 0; i < component_num; i++) {
				if (components[i].empty())
					continue;
				if (components[i].size() == 1)
					disassembledobject.push_back(make_pair(components[i][0], i));
			}
			disassembledobject.push_back(make_pair(visitedNum, current_component - 1));

			disassembled_objects.push_back(disassembledobject);
			break;
		}
		else {
			if (obj_connected.size() == connected_objects.size())   
			{
				return false;
			}
			else   
			{
				timeRecord.Recursion_num++;
				unique_id++;

				depth++;
#if 0  
				connected_objects.clear();
				for (int i = 0; i < obj_connected.size(); i++)
					connected_objects.push_back(obj_connected[i]);
#else  
				vector<vector<int>> components(component_num);
				for (int i = 0; i < visitedNum; i++) {
					if (!positions[i].packed) continue;

					int idx = positions[i].center[0] * NY * NZ + positions[i].center[1] * NZ + positions[i].center[2];
					if (visited[idx] != 0)
						components[visited[idx] - 1].push_back(i);
				}

				vector<int> UnremovedObjects;
				vector<pair<int, int>> disassembledobject;
				for (int i = 0; i < component_num; i++) {
					if (components[i].empty())
						continue;
					if (components[i].size() == 1 && i != current_component - 1)
					{
						disassembledobject.push_back(make_pair(components[i][0], i));
						continue;  
					}

					for (int j = 0; j < components[i].size(); j++)
						UnremovedObjects.push_back(components[i][j]);
				}
				disassembled_objects.push_back(disassembledobject);

				if (!UnremovedObjects.empty()) {
					if (UnremovedObjects.size() == connected_objects.size())
						return false;  

					sort(UnremovedObjects.begin(), UnremovedObjects.end());
				}

				connected_objects.clear();
				for (int i = 0; i < UnremovedObjects.size(); i++)
					connected_objects.push_back(UnremovedObjects[i]);
#endif
			}
		}
	}

	if (is_disassemble) {
		SAMPLE_ON_BALL sampleOri;
		sampleOri.num_ori_sample = CutterNum;
		sampleOri.EnumerationOritation();
		vector<Vector> sample_ori;
		for (int c = 0; c < CutterNum; c++)
		{
			Vector tmp_ori = sampleOri.sample_points[c];
			tmp_ori = -tmp_ori;
			sample_ori.push_back(tmp_ori);
		}

		for (int i = 0; i < visibility_layers.size(); i++) {
			vector<pair<int, int>> disassembledobject;
			disassembledobject = disassembled_objects[i];

			int component_num = component_nums[i];
			vector<vector<int>> internal_visual_hulls(component_num);
			for (int j = 0; j < TX; j++) {
				if(visited_layers[i][j] != 0){
					int component_id = visited_layers[i][j] - 1;
					internal_visual_hulls[component_id].push_back(j);
				}
			}

			while (!disassembledobject.empty()) {
				int remain_num = disassembledobject.size();
				
				for (int j = 0; j < disassembledobject.size();) {
					int obj_id = disassembledobject[j].first;
					int component_id = disassembledobject[j].second;
					vector<int> internal_visual_hull = internal_visual_hulls[component_id];
					bool isvalid = true;
					for (int c = 0; c < CutterNum; c++) {  
						isvalid = true;
						for (int ivh = 0; ivh < internal_visual_hull.size(); ivh++)
						{
							int idx = internal_visual_hull[ivh];
							Vector direction = sample_ori[c];
							int ix, iy, iz;
							ix = idx / (NY * NZ);
							iy = (idx / NZ) % NY;
							iz = idx % NZ;

							int ox = ix, oy = iy, oz = iz;

							double dir_x = direction[0];
							double dir_y = direction[1];
							double dir_z = direction[2];

							int step_x = (dir_x > 0) ? 1 : (dir_x < 0) ? -1 : 0;
							int step_y = (dir_y > 0) ? 1 : (dir_y < 0) ? -1 : 0;
							int step_z = (dir_z > 0) ? 1 : (dir_z < 0) ? -1 : 0;

							double t_max_x = (step_x == 0)
								? std::numeric_limits<double>::infinity()
								: (step_x > 0
									? 1.0 / dir_x
									: 0.0);

							double t_max_y = (step_y == 0)
								? std::numeric_limits<double>::infinity()
								: (step_y > 0
									? 1.0 / dir_y
									: 0.0);

							double t_max_z = (step_z == 0)
								? std::numeric_limits<double>::infinity()
								: (step_z > 0
									? 1.0 / dir_z
									: 0.0);

							double t_delta_x = (dir_x == 0)
								? std::numeric_limits<double>::infinity()
								: 1.0 / std::abs(dir_x);

							double t_delta_y = (dir_y == 0)
								? std::numeric_limits<double>::infinity()
								: 1.0 / std::abs(dir_y);

							double t_delta_z = (dir_z == 0)
								? std::numeric_limits<double>::infinity()
								: 1.0 / std::abs(dir_z);

							bool is_intersect = false;
							while (ox >= 0 && ox < NX && oy >= 0 && oy < NY && oz >= 0 && oz < NZ) {
								int idx = ox * NY * NZ + oy * NZ + oz;

								if (visibility_layers[i][idx] == 0 && visited_layers[i][idx] != component_id + 1) 
								{
									is_intersect = true;
									break;
								}

								if (t_max_x < t_max_y && t_max_x < t_max_z) {
									ox += step_x;
									t_max_x += t_delta_x;
								}
								else if (t_max_y < t_max_z) {
									oy += step_y;
									t_max_y += t_delta_y;
								}
								else {
									oz += step_z;
									t_max_z += t_delta_z;
								}
							}

							if (is_intersect) {
								isvalid = false;
								break;
							}
						}
						if (isvalid) {
							break;
						}
					} 

					if (isvalid) {
						disassembledobject.erase(disassembledobject.begin() + j);
						for (int ivh = 0; ivh < internal_visual_hull.size(); ivh++)
						{
							int idx = internal_visual_hull[ivh];
							visited_layers[i][idx] = 0;
						}
					}
					else {
						j++;
					}
				}
				if (remain_num == disassembledobject.size()) 
					return false;  
			}
		}
	}
	return true;
}



