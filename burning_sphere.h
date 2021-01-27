#pragma once
#include "Vec.h"
#include <vector>
#include <set>
using namespace std;

class BurningSphere
{

private:
	vector<Point3f>& sample_space;
	vector< vector<int> >& sample_neighbors;

public:
	int center_index;
	vector<int> inner_points_index;
	vector<int> outer_boundary;

	float* burning_time;
	//vector<int> pre;
	int* pre;
	//vector<vector<int>> path_tree;
	vector<int> son_num;
	float max_burning_time;
	
	bool debug_flag;

	BurningSphere(const BurningSphere& it) : sample_space(it.sample_space), sample_neighbors(it.sample_neighbors)
	{
		center_index = it.center_index;
		inner_points_index = it.inner_points_index;
		outer_boundary = it.outer_boundary;
		burning_time = it.burning_time;
		max_burning_time = it.max_burning_time;

		
	}

	BurningSphere(int center_index, vector<Point3f>& tsample_space, vector< vector<int> >& tsample_neighbors, float* init_burning_time): sample_space(tsample_space), sample_neighbors(tsample_neighbors)
	{
		burning_time = init_burning_time;		

		this->center_index = center_index;
		burning_time[center_index]=-1;
		max_burning_time = 0;
		inner_points_index.push_back(center_index);

		//pre.clear();
		pre = new int[sample_space.size()];
		memset(pre, -1, sample_space.size()*sizeof(int));
		//path_tree.clear();
		//for (int i = 0; i < sample_space.size(); i++)
		//{
		//	pre.push_back(-1);
		//	son_num.push_back(0);

		//	//vector<int> son_set;
		//	//path_tree.push_back(son_set);

		//}

		int u = center_index;
		for (int i = 0; i < sample_neighbors[u].size(); i++)
		{
			int v = sample_neighbors[u][i];
			outer_boundary.push_back(v);

			float dist = (sample_space[u] - sample_space[v]).length();
			burning_time[v] = dist;

			pre[v] = u;
			//path_tree[u].push_back(v);
		}

		debug_flag = false;
	}

	void free_memery()
	{
		for (int i = 0; i < inner_points_index.size(); i++)
		{
			burning_time[inner_points_index[i]] = 0.0;
		}
		for (int i = 0; i < outer_boundary.size(); i++)
		{
			burning_time[outer_boundary[i]] = 0.0;
		}

		inner_points_index.~vector();
		outer_boundary.~vector();
		//pre.~vector();
		son_num.~vector();

		delete[] pre;
	}

	void burning()
	{
		max_burning_time += 1.0;
	
		//burning
		vector<int> new_outer_boundary;
		vector<int> old_new_outer_boundary;

		for (int i = 0; i < outer_boundary.size(); i++)
		{
			if (burning_time[outer_boundary[i]] <= max_burning_time)
			{
				int u = outer_boundary[i]; 
				inner_points_index.push_back(u); 

				for (int j = 0; j < sample_neighbors[u].size(); j++)
				{
					int v = sample_neighbors[u][j];
					float dist = (sample_space[u] - sample_space[v]).length();

					if (burning_time[v]==0)
					{
						new_outer_boundary.push_back(v); 

						burning_time[v] = burning_time[u] + dist;
						pre[v] = u;

					}
					else
					{
						if (burning_time[v] >= burning_time[u] + dist )
						{
							burning_time[v] = burning_time[u] + dist;
							pre[v] = u;

							
						}
					}
				}
			}
			else
			{
				new_outer_boundary.push_back(outer_boundary[i]);
			}
		}

		outer_boundary = new_outer_boundary;

	}
};