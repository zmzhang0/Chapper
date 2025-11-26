#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <regex>
#include "CGALib.h"

using VEctor = std::array<double, 3>;
using TRiangle = std::array<int, 3>;

class Slicer_2
{
public:
	std::vector<VEctor> positions;
	std::vector<TRiangle> triangles;
	std::vector<bool> jud_plane;

	Point origin;
	Vector normal;

	void clear()
	{
		positions.clear();
		triangles.clear();
	}


	bool load_off(std::string filename)
	{
		clear();
		std::ifstream file(filename);
		if (!file)
			return false;
		std::string line;
		file >> line;
		int Vnum = 0, Fnum = 0, Enum = 0;
		file >> Vnum >> Fnum >> Enum;

		for (int i = 0; i < Vnum; i++)
		{
			VEctor p;

			file >> p[0] >> p[1] >> p[2];
			positions.push_back(p);
		}
		for (int i = 0; i < Fnum; i++) {
			TRiangle t;

			int p_num = 3;

			file >> p_num >> t[0] >> t[1] >> t[2];

			triangles.push_back(t);
		}
		file.clear();
		file.close();

		return true;
	}

	bool load(std::string filename)
	{
		clear();
		std::ifstream file(filename);
		if (!file)
			return false;
		std::string line;
		while (std::getline(file, line))
		{
			std::istringstream iss(line);
			std::string prefix;
			iss >> prefix;

			if (prefix == "v")
			{
				VEctor p;
				iss >> p[0] >> p[1] >> p[2];
				positions.push_back(p);
			}

			else if (prefix == "f")
			{
				TRiangle t;
				iss >> t[0] >> t[1] >> t[2];
				for (int i = 0; i < 3; i++)
					t[i]--;
				triangles.push_back(t);
			}

			else
			{
			}
		}

		file.clear();
		file.close();

		return true;
	}

	bool save(std::string filename)
	{
		std::ofstream file(filename);
		if (!file)
			return false;
		for (VEctor p : positions)
		{
			file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		}
		for (TRiangle t : triangles)
		{
			for (int i = 0; i < 3; i++)
				t[i]++;
			file << "f " << t[0] << " " << t[1] << " " << t[2] << std::endl;
		}
		file.clear();
		file.close();
		return true;
	}

	bool save_off(std::string filename)
	{
		std::ofstream file(filename);
		if (!file)
			return false;
		file << "OFF" << std::endl;
		file << positions.size() << " " << triangles.size() << " " << 0 << std::endl;
		file << std::endl;

		for (VEctor p : positions)
		{
			file << p[0] << " " << p[1] << " " << p[2] << std::endl;
		}
		for (TRiangle t : triangles)
		{
			int num_of_p = 3;
			file << num_of_p << " " << t[0] << " " << t[1] << " " << t[2] << std::endl;
		}
		file.clear();
		file.close();

		return true;
	}

	void cut()
	{
		intersections.clear();
		if (origin[0] == 0 && origin[1] == 0 && origin[2] == 0)
			return;
		for (size_t tid = 0; tid < triangles.size(); tid += !do_triangle(tid))
			;
	}

private:
	static constexpr double precision = 0.0000001;  

	std::map<std::pair<int, int>, int> intersections;
	double get_lambda(VEctor p, VEctor q)
	{
		double num = (origin[0] - q[0]) * normal[0] + (origin[1] - q[1]) * normal[1] + (origin[2] - q[2]) * normal[2];

		double den = (p[0] - q[0]) * normal[0] + (p[1] - q[1]) * normal[1] + (p[2] - q[2]) * normal[2];
		if (abs(num) < 0.0000000001 || abs(den) < 0.0000000001)                   
			return 0.0;
		else
			return num / den;
	}

	int get_intersection(int i, int j)
	{
		if (i > j)
			std::swap(i, j);

		VEctor p = positions[i];
		VEctor q = positions[j];
		double lambda = get_lambda(p, q);
		if (!std::isfinite(lambda) || lambda < precision || lambda > 1 - precision)
			return -1;

		auto key = std::make_pair(i, j);
		auto found = intersections.find(key);
		if (found != intersections.end())
			return found->second;

		int m = positions.size();
		positions.push_back(
			{ lambda * p[0] + (1 - lambda) * q[0],
			 lambda * p[1] + (1 - lambda) * q[1],
			 lambda * p[2] + (1 - lambda) * q[2] });
		intersections[key] = m;
		return m;
	}

	bool do_triangle(int tid)
	{
		for (int n = 0; n < 3; n++)
		{
			int i = triangles[tid][n];
			int j = triangles[tid][(n + 1) % 3];
			int k = triangles[tid][(n + 2) % 3];

			int m = get_intersection(j, k);
			if (m != -1)
			{
				triangles.push_back({ i, j, m });
				triangles[tid] = { i, m, k };
				jud_plane.push_back(true);
				return true;
			}
		}
		return false;
	}

};