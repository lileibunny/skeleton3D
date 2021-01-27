#pragma once
#include "Vec.h"
#include <vector>
#include <queue>
using namespace std;

class Solid_sphere
{

private:
	queue<int> bfs_que;

public:
	Point3f center;
	float radius;
	vector<int> sample_point_indexs;
	vector<int> outer_boundary;
	bool* vist;
	bool debug_flag;


	Solid_sphere(Point3f center, int center_index)
	{
		this->center = center;
		radius = 0;
		sample_point_indexs.clear();
		outer_boundary.push_back(center_index);
		vist = NULL;
		debug_flag = false;
	}
	void free_memery()
	{
		if (vist != NULL)
		{
			delete[] vist;
		}
		vist = NULL;
		sample_point_indexs.clear();
		sample_point_indexs.~vector();

		outer_boundary.clear();
		outer_boundary.~vector();
	}

	void expand(float delta_radius, vector<Point3f>& sample_space, vector< vector<int> >& sample_neighbors)
	{

		if (debug_flag)
		{
			cout << "expand in" << endl;
		}

		//init
		if (radius==0)
		{
			vist = new bool[sample_space.size()];
			memset(vist, 0, sample_space.size());
		}
		radius += delta_radius;
		vector<int> new_outer_boundary;
		

		//BFS
		int sample_space_size = sample_space.size();
	
	
		for (int i = 0; i < outer_boundary.size(); i++)
		{
			float dist_to_center = (center - sample_space[outer_boundary[i]]).length();

			if (dist_to_center <= radius)
			{
				sample_point_indexs.push_back(outer_boundary[i]);
				bfs_que.push(outer_boundary[i]);
			}
			else
			{
				new_outer_boundary.push_back(outer_boundary[i]);
			}

			vist[outer_boundary[i]] = true;
		}

		while (!bfs_que.empty()) //BFS mainloop
		{
			
			int u = bfs_que.front();
			bfs_que.pop();
			for (int i = 0; i < sample_neighbors[u].size(); i++)
			{
				int v = sample_neighbors[u][i];

				if (vist[v])
				{
					continue;
				}

				float dist_to_center = (center - sample_space[v]).length();
				if (dist_to_center <= radius)
				{
					sample_point_indexs.push_back(v);
					bfs_que.push(v);
				}
				else
				{
					new_outer_boundary.push_back(v);
				}

				vist[v] = true;
			}
		}

		outer_boundary = new_outer_boundary;
		new_outer_boundary.clear();
		new_outer_boundary.~vector();

		if (debug_flag)
		{
			cout << "expand out" << endl;
		}
	}

};