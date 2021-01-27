#include "Main.h"
#include "skeleton.h"
#include "stdio.h"
#include "pba3D.h"
#include "common.h"
#include "graph.h"
#include "sphere.h"
#include "burning_sphere.h"
#include "morphology.h"

#include<set>
#include <algorithm>
#include <ctime>
#include <queue>
using namespace std;


bool* vist;
vector<int> visited_points_index;
int visitNodeCount;


inline int find(int x, int father[]){
	return x == father[x] ? x : (father[x] = find(father[x], father));   //�ݹ���Ҽ��ϵĴ���Ԫ�أ���·��ѹ����
}

set<int> boundary_points;

bool Skeleton::is_noise_segment(int* path_pre, vector<vector<int>>& front_segments, int segment_index)
{
	bool* hash_table = new bool[medialSurface.size()];
	int crossed_path_num = 0;
	int uncrossed_path_num = 0;
	int current_segment_size = front_segments[segment_index].size();

	vector<int> all_other_path_end;
	for (int i = 0; i < front_segments.size(); i++)
	{ 
		if (i == segment_index)
		{
			continue;
		}

		for (int j = 0; j < front_segments[i].size(); j++)
		{
			all_other_path_end.push_back(front_segments[i][j]);
		}
	}

	for (int p = 0; p < front_segments[segment_index].size(); p += 1)
	{
		bool is_crossed = false;

		vector<int> main_burning_path; 
		vector<int> other_burning_path;

		int path_point = front_segments[segment_index][p]; 
		
		while (path_point != -1)
		{
			main_burning_path.push_back(path_point);  
			
			path_point = path_pre[path_point];
		} 

		drawing_paths.push_back(main_burning_path);

		memset(hash_table, 0, medialSurface.size()*sizeof(bool));
		int path_size = main_burning_path.size();

		for (int i = 0; i < path_size; i++)
		{
			hash_table[main_burning_path[i]] = 1;
		}


		for (int j = 0; j < all_other_path_end.size(); j++)
		{
			if (is_crossed)
			{
				break;
			}

			other_burning_path.clear();
			int cross_num = 0;

			path_point = all_other_path_end[j];
			
			while (path_point != -1)
			{
				other_burning_path.push_back(path_point); 
			
				if (hash_table[path_point] == 1)
				{
					cross_num++;

					if (cross_num >= 2)
					{
						is_crossed = true;
						break;
					}
				}

				path_point = path_pre[path_point];
			} 
		}

		if (is_crossed == false)
		{
			delete[] hash_table;
			return false;
		}

	}
	
	delete[] hash_table;
	return true;
}
void Skeleton::symmetry_field_generation()
{
	if (symmetry_field_generation_flag == true)
	{
		return;
	}

	//prepare data
	max_symmetry_field_value = symmetry_field_threshold;
	int denoise_factor = MAT_noise_size;
	int n = medialSurface.size();
	int* medialSurface_to_boundary = new int[n];
	symmetry_field.clear();

	//parameter setting
	const int max_graphNode_num = 10000;
	const int max_neighbor_num = 30;  //26 is enough

	int father[max_graphNode_num];
	int graph_edges[max_graphNode_num*max_neighbor_num][2];
	int n_graph_edge;

	Graph<int> propagating_boundary_graph;
	vector<int> propagating_boundary;

	vector<vector<int>> fire_front_segments;
	set<int> fire_front;
	set<int> fire_front_noises;


	float* init_burning_time = new float[n];
	memset(init_burning_time, 0, sizeof(float)*n);

	int time_count1 = 0;
	int time_count2 = 0;
	int time_tic, time_toc;

	//main loop
	int symmetry_field_time_tic = clock();
	int O_k = 0;
	for (int i = 0; i < medialSurface.size(); i++)
	{
		symmetry_field.push_back(0);
	}

	for (int i = 0; i < medialSurface.size(); i++) 
	{
		memset(medialSurface_to_boundary, -1, sizeof(int)*n);
		memset(init_burning_time, 0, sizeof(float)*n);

		BurningSphere burning_sphere(i, medialSurface, medialSurface_net, init_burning_time);

		int break_flag = 0;
		vector<vector<int>> fire_front_set;
		fire_front_set.push_back(burning_sphere.outer_boundary);

		while (1)
		{
			burning_sphere.burning();
			fire_front_set.push_back(burning_sphere.outer_boundary);

			propagating_boundary = burning_sphere.outer_boundary;
			int outer_boundary_size = propagating_boundary.size();

			//check if the outer boundary is connected by Union-Find data structure
			for (int j = 0; j < outer_boundary_size; j++)
			{
				medialSurface_to_boundary[propagating_boundary[j]] = j;
				father[j] = j;
			}

			n_graph_edge = 0;
			for (int j = 0; j < outer_boundary_size; ++j)
			{
				int point_id = propagating_boundary[j];
				int neighbor_num = medialSurface_net[point_id].size();
				for (int k = 0; k < neighbor_num; ++k)
				{
					int v = medialSurface_to_boundary[medialSurface_net[point_id][k]];
					if (v > j)
					{
						graph_edges[n_graph_edge][0] = j;
						graph_edges[n_graph_edge++][1] = v;
					}
				}

			}

			//build the Union-Find data structure
			for (int j = 0; j < n_graph_edge; j++)
			{
				int root_u = find(graph_edges[j][0], father);
				int root_v = find(graph_edges[j][1], father);
				father[root_u] = root_v;
			}

			int n_subGraph = 0;
			vector<int> segment_roots;
			for (int j = 0; j < outer_boundary_size; j++)
			{

				medialSurface_to_boundary[propagating_boundary[j]] = -1;

				if (father[j] == j)
				{
					segment_roots.push_back(j);
					n_subGraph++;  //count the number of subGraphs
				}
			}

			fire_front_segments.clear();
			for (int j = 0; j < segment_roots.size(); j++)
			{
				fire_front_segments.push_back(vector<int>());
			}

			if (n_subGraph > 1)
			{
				fire_front.clear();
				for (int j = 0; j < outer_boundary_size; ++j)
				{
					fire_front.insert(propagating_boundary[j]);
					int root = find(j, father);
					for (int k = 0; k < segment_roots.size(); k++)
					{
						if (segment_roots[k] == root)
						{
							fire_front_segments[k].push_back(propagating_boundary[j]);
						}
					}
				}

				vector<int> new_fire_front_segment;
				int front_seg_num = fire_front_segments.size();

				special_leaves = new_fire_front_segment;
				

				for (int j = 0; j < front_seg_num; j++)
				{				
					if (is_noise_segment(burning_sphere.pre, fire_front_segments, j))
					{
						n_subGraph--;

						fire_front_noises.clear();
						for (int k = 0; k < fire_front_segments[j].size(); k++)
						{
							fire_front_noises.insert(fire_front_segments[j][k]);
						}

						propagating_boundary.clear();
						for (int k = 0; k < burning_sphere.outer_boundary.size(); k++)
						{
							if (fire_front_noises.find(burning_sphere.outer_boundary[k]) == fire_front_noises.end())
							{
								propagating_boundary.push_back(burning_sphere.outer_boundary[k]);
							}
							else
							{
								init_burning_time[burning_sphere.outer_boundary[k]] = 0.0;
							}
						}

						//update the outer_boundary
						burning_sphere.outer_boundary = propagating_boundary;
					}
				}
				if (debug_flag)
				{
					cout << endl;
				}

			}

			if (n_subGraph>1)  //the graph is not connected
			{
				break_flag++;
			}
			else
			{
				break_flag = 0; //reset 
			}
			float fire_burning_time = burning_sphere.max_burning_time;

			if (break_flag >= denoise_factor)
			{
				symmetry_field[i] = max_symmetry_field_value - (fire_burning_time - denoise_factor + 1);
				break;
			}

			//reach the threshold
			if (burning_sphere.max_burning_time >= max_symmetry_field_value)
			{
				float radius = sqrt(dist->at_dist(medialSurface[i][0], medialSurface[i][1], medialSurface[i][2]));
				symmetry_field[i] = 0;
				break;
			}
		}

		fire_front_space.push_back(fire_front_set);



		////////////////free memeory///////////////////////
		burning_sphere.free_memery();

		propagating_boundary.~vector();
		fire_front_segments.~vector();

	}
	delete[] medialSurface_to_boundary;
	delete[] init_burning_time;

	int symmetry_field_time_toc = clock();
	time_count1 = symmetry_field_time_toc - symmetry_field_time_tic;
	cout << "symmetry field generation finished!  used time: " << time_count1 << " ms" << endl;

	//symmetry_field_smoothing();
	symmetry_field_generation_flag = true;
}

void Skeleton::symmetry_field_smoothing()
{
	int singularity_num;
	do
	{
		singularity_num = 0;
		for (int i = 0; i < medialSurface.size(); i++)
		{
			if (symmetry_field[i] == 0)
			{
				singularity_num++;

				float symmetry_field_mean = 0;
				int neighbor_num = 0;
				for (int j = 0; j<medialSurface_net[i].size(); j++)
				{
					int neighbor = medialSurface_net[i][j];
					symmetry_field_mean += symmetry_field[neighbor];
					neighbor_num++;
				}
				if (neighbor_num > 0)
				{
					symmetry_field_mean /= neighbor_num;
				}
				symmetry_field[i] = symmetry_field_mean;
			}
		}
	}while (singularity_num > 0);

}


void Skeleton::symmetry_field_gradient_generation()
{
	float* init_burning_time = new float[medialSurface.size()];
	memset(init_burning_time, 0, sizeof(float)*medialSurface.size());

#if 1
	float neighborhood_radius = 2.5;

	symmetry_field_gradient.clear();
	Point3f gradient_vector;
	for (int i = 0; i<medialSurface.size(); i++)
	{
		/*Solid_sphere tiny_spherical_neighborhood(medialSurface[i], i);
		tiny_spherical_neighborhood.expand(neighborhood_radius, medialSurface, medialSurface_net);*/
		BurningSphere tiny_spherical_neighborhood(i, medialSurface, medialSurface_net, init_burning_time);
		while (tiny_spherical_neighborhood.max_burning_time < neighborhood_radius)
		{
			tiny_spherical_neighborhood.burning();
		}
		int neighbor_num = tiny_spherical_neighborhood.inner_points_index.size();

		float max_value = -1;
		float min_value = 1 << 30;
		for (int j = 0; j < neighbor_num; j++)
		{
			int neighbor_id = tiny_spherical_neighborhood.inner_points_index[j];

			if (symmetry_field[neighbor_id]>max_value)
			{
				max_value = symmetry_field[neighbor_id];
			}

			if (symmetry_field[neighbor_id]<min_value)
			{
				min_value = symmetry_field[neighbor_id];
			}
		}

		Point3f gradient_vector(0, 0, 0);
		Point3f local_maximal_point(0, 0, 0);
		Point3f local_minimal_point(0, 0, 0);
		float cnt1 = 0, cnt2 = 0;
		for (int j = 0; j < neighbor_num; j++)
		{
			int neighbor_id = tiny_spherical_neighborhood.inner_points_index[j];
			Point3f neighbor_point = medialSurface[neighbor_id];
			/*
			if (symmetry_field[neighbor_id]>=symmetry_field[i])
			{
			Point3f orientation_vector = neighbor_point - medialSurface[i];
			float weight = symmetry_field[neighbor_id] - symmetry_field[i]; weight = 1;
			gradient_vector += weight*orientation_vector;
			}*/

			if (symmetry_field[neighbor_id] == max_value)
			{
				local_maximal_point += neighbor_point;
				cnt1++;
			}
			if (symmetry_field[neighbor_id] == min_value)
			{
				local_minimal_point += neighbor_point;
				cnt2++;
			}
		}
		local_maximal_point /= cnt1;
		local_minimal_point /= cnt2;

		if (max_value < symmetry_field[i])
		{
			gradient_vector = Point3f(0, 0, 0);
		}
		else
		{
			gradient_vector = local_maximal_point - medialSurface[i];
		}
		if (gradient_vector.length() > 0)
		{
			gradient_vector /= gradient_vector.length();
		}

		symmetry_field_gradient.push_back(gradient_vector);

		tiny_spherical_neighborhood.free_memery();
	}

	vector<Point3f> smoothed_gradients;
	float smoothing_radius = 2.5;
	for (int i = 0; i < medialSurface.size(); i++)
	{
		BurningSphere tiny_spherical_neighborhood(i, medialSurface, medialSurface_net, init_burning_time);
		while (tiny_spherical_neighborhood.max_burning_time < smoothing_radius)
		{
			tiny_spherical_neighborhood.burning();
		}
		int neighbor_num = tiny_spherical_neighborhood.inner_points_index.size();

		Point3f smoothed_gradient(0, 0, 0);
		for (int j = 0; j < neighbor_num; j++)
		{
			int neighbor_id = tiny_spherical_neighborhood.inner_points_index[j];
			smoothed_gradient += symmetry_field_gradient[neighbor_id];
		}

		smoothed_gradient /= neighbor_num;
		if (smoothed_gradient.length()>0)
		{
			smoothed_gradient /= smoothed_gradient.length();
		}
		smoothed_gradients.push_back(smoothed_gradient);

		tiny_spherical_neighborhood.free_memery();
	}
	//symmetry_field_gradient = smoothed_gradients;

#else 
	symmetry_field_gradient.clear();
	for (int i = 0; i < medialSurface.size(); i++)
	{
		BurningSphere tiny_spherical_neighborhood(i, medialSurface, medialSurface_net, init_burning_time);

		int neighbor_num = tiny_spherical_neighborhood.outer_boundary.size();
		Point3f sum(0, 0, 0);
		for (int j = 0; j < neighbor_num; j++)
		{
			int neighbor_index = tiny_spherical_neighborhood.outer_boundary[j];
			Point3f neighbor_dir = medialSurface[neighbor_index] - medialSurface[i];
			int dx = neighbor_dir[0];
			int dy = neighbor_dir[1];
			int dz = neighbor_dir[2];

			float sobel_x;
			if (dx != 0)
			{
				if ((dy || dz) == 0)
				{
					sobel_x = dx * 6;
				}
				else if ((dy && dz) == 1)
				{
					sobel_x = dx * 1;
				}
				else
				{
					sobel_x = dx * 3;
				}

			}
			else
			{
				sobel_x = 0;
			}

			float sobel_y;
			if (dy != 0)
			{
				if ((dx || dz) == 0)
				{
					sobel_y = dy * 6;
				}
				else if ((dx && dz) == 1)
				{
					sobel_y = dy * 1;
				}
				else
				{
					sobel_y = dy * 3;
				}

			}
			else
			{
				sobel_y = 0;
			}

			float sobel_z;
			if (dz != 0)
			{
				if ((dy || dx) == 0)
				{
					sobel_z = dz * 6;
				}
				else if ((dy && dx) == 1)
				{
					sobel_z = dz * 1;
				}
				else
				{
					sobel_z = dz * 3;
				}

			}
			else
			{
				sobel_z = 0;
			}

			float function_value = symmetry_field[neighbor_index];
			sum[0] += sobel_x * function_value;
			sum[1] += sobel_y * function_value;
			sum[2] += sobel_z * function_value;
		}
		if (sum.length() > 0)
		{
			sum /= sum.length();
		}

		symmetry_field_gradient.push_back(sum);
		tiny_spherical_neighborhood.free_memery();
	}
#endif

	delete[] init_burning_time;

}


void Skeleton::dfs_in_ridges(int root)
{
	vist[root] = true;
	visitNodeCount++;
	visited_points_index.push_back(root);

	for (int i = 0; i < symmetry_field_ridges.size(); i++)
	{
		if (i == root || vist[i] == true) continue;

		int dx = medialSurface[symmetry_field_ridges[root]][0] - medialSurface[symmetry_field_ridges[i]][0];
		int dy = medialSurface[symmetry_field_ridges[root]][1] - medialSurface[symmetry_field_ridges[i]][1];
		int dz = medialSurface[symmetry_field_ridges[root]][2] - medialSurface[symmetry_field_ridges[i]][2];

		if (abs(dx) <= 1 && abs(dy) <= 1 && abs(dz) <= 1)
		{
			dfs_in_ridges(i);
		}
	}
}

void  Skeleton::symmetry_field_ridges_detection()
{
	vector<float> div_values;
	float* init_burning_time = new float[medialSurface.size()];
	memset(init_burning_time, 0, sizeof(float)*medialSurface.size());

	float neighborhood_radius = 3.5;
	float neighbor_size_mean = 0;
	for (int i = 0; i < medialSurface.size(); i++)
	{
		BurningSphere tiny_spherical_neighborhood(i, medialSurface, medialSurface_net, init_burning_time);
		while (tiny_spherical_neighborhood.max_burning_time < neighborhood_radius)
		{
			tiny_spherical_neighborhood.burning();
		}
		neighbor_size_mean += tiny_spherical_neighborhood.inner_points_index.size() - 1;
		tiny_spherical_neighborhood.free_memery();
	}
	neighbor_size_mean /= medialSurface.size();
	neighbor_size_mean *= 0.0;

	symmetry_field_ridges.clear();
	for (int i = 0; i < medialSurface.size(); i++)
	{
		if (symmetry_field[i] <= 0)
		{
			continue;
		}

		bool is_minimum_point = true;
		float smaller_value_rate = 0;
		float larger_value_rate = 0;

		//Solid_sphere tiny_spherical_neighborhood(medialSurface[i], i);
		//tiny_spherical_neighborhood.expand(neighborhood_radius, medialSurface, medialSurface_net);

		BurningSphere tiny_spherical_neighborhood(i, medialSurface, medialSurface_net, init_burning_time);
		while (tiny_spherical_neighborhood.max_burning_time < neighborhood_radius)
		{
			tiny_spherical_neighborhood.burning();
		}

		int neighbor_num = tiny_spherical_neighborhood.inner_points_index.size();
		for (int k = 0; k <neighbor_num; k++)
		{
			int neighbor_id = tiny_spherical_neighborhood.inner_points_index[k];
			if (neighbor_id == i)
			{
				continue;
			}

			if (symmetry_field[neighbor_id] < symmetry_field[i])
			{
				is_minimum_point = false;
			}

			if (symmetry_field[neighbor_id]<symmetry_field[i])
			{
				smaller_value_rate++;
			}

			if (symmetry_field[neighbor_id] > symmetry_field[i] + 0.0001)
			{
				larger_value_rate++;
			}
		}		
		neighbor_num--;

		smaller_value_rate /= neighbor_num;
		larger_value_rate /= neighbor_num;

		if (larger_value_rate <= 0.125 /*&& neighbor_num >= neighbor_size_mean*/ /*|| (larger_value_rate >= 0.49&&larger_value_rate <= 0.51 && smaller_value_rate >= 0.49 && smaller_value_rate <= 0.51)*/)
		{
			symmetry_field_ridges.push_back(i);
		}

		tiny_spherical_neighborhood.free_memery();
	}
	cout << "ridge!!   neighbor_size_mean: " << neighbor_size_mean << endl;

	delete[] init_burning_time;

	//connect all ridges
	symmetry_field_ridges_connection();

	//build the adjacent relation of symmetry field ridge points
	symmetry_field_ridges_net.clear();
	for (int i = 0; i < symmetry_field_ridges.size(); i++)
	{
		symmetry_field_ridges_net.push_back(*(new vector<int>));
	}
	for (int i = 0; i < symmetry_field_ridges.size(); i++)
	{
		Point3f ridge_point_u = medialSurface[symmetry_field_ridges[i]];

		for (int j = i + 1; j < symmetry_field_ridges.size(); j++)
		{
			Point3f ridge_point_v = medialSurface[symmetry_field_ridges[j]];
			if ((ridge_point_u - ridge_point_v).length() < 2.0)
			{
				symmetry_field_ridges_net[i].push_back(j);
				symmetry_field_ridges_net[j].push_back(i);
			}
		}

	}

	cout << endl << "symmetry_field_ridges detection finished! the number of ridge points is " << symmetry_field_ridges.size() << endl;
}


//connect ridges by doing path tracing at each ridge point along the gradient vector 
void Skeleton::symmetry_field_ridges_connection()
{
	//return;
#if 1
	int path_tracing_time_tic = clock();
	int n = medialSurface.size();
	bool* is_ridge = new bool[n];
	memset(is_ridge, 0, n);
	for (int i = 0; i < symmetry_field_ridges.size(); i++)
	{
		is_ridge[symmetry_field_ridges[i]] = true;
	}

	//construct the graph by medial surface points and their burnning time
	vector<vector<pair<int, float>> > graph;
	float burnning_time_factor = 3;
	for (int i = 0; i < n; i++)
	{
		float burnning_time_u = max_symmetry_field_value - symmetry_field[i];
		burnning_time_u = pow(burnning_time_factor, burnning_time_u / 2);

		if (is_ridge[i])
		{
			burnning_time_u =0;
		}

		vector<pair<int, float>> neighbors;
		for (int j = 0; j < medialSurface_net[i].size(); j++)
		{
			int neighbor = medialSurface_net[i][j];
			float burnning_time_v = max_symmetry_field_value - symmetry_field[neighbor];
			burnning_time_v = pow(burnning_time_factor, burnning_time_v / 3);
			if (is_ridge[neighbor])
			{
				burnning_time_v =0;
			}

			float dist = /*burnning_time_u + */burnning_time_v;
			dist = (medialSurface[i] - medialSurface[neighbor]).length();
			neighbors.push_back(pair<int, float>(neighbor, dist));
		}

		graph.push_back(neighbors);
	}


	//connecting...

	//connect all ridges
	int m = symmetry_field_ridges.size();
	for (int i = 0; i < m; i++)
	{
		int ridge_index1 = symmetry_field_ridges[i];
		Point3f ridge_pt1 = medialSurface[symmetry_field_ridges[i]];

		int* pre = new int[n];
		long  double* shortest_path_dist = new long double[n];
		spfa(graph, ridge_index1, shortest_path_dist, pre);

		for (int j = 0; j < m; j++)
		{
			int ridge_index2 = symmetry_field_ridges[j];
			Point3f ridge_pt2 = medialSurface[symmetry_field_ridges[j]];
			float dist = (ridge_pt1 - ridge_pt2).length();

			if (dist >= 2 && dist <= 2 * sqrt(3.0))
			{
				int new_ridge_num = 0;			
				int path_point = ridge_index2;

				//cout << "ridge pt2: " << ridge_index2 << endl; continue;
				while (path_point != -1)
				{
					if (is_ridge[path_point] == false )
					{
						symmetry_field_ridges.push_back(path_point);
						new_ridge_num++;
						is_ridge[path_point] = true;

					//	cout << "new ridge point: " << path_point << endl;
					}

					path_point = pre[path_point];
				}				
			}

		}

		delete[] pre;
		delete[] shortest_path_dist;

	}

	delete[] is_ridge;

	int path_tracing_time_toc = clock();
	cout << "ridges connection time: " << path_tracing_time_toc - path_tracing_time_tic << " ms" << endl;
//	cout << "the number of ridges is:  " << symmetry_field_ridges.size() << endl;

#else

	int path_tracing_time_tic = clock();
	int n = medialSurface.size();
	bool* is_ridge = new bool[n];
	memset(is_ridge, 0, n);
	for (int i = 0; i < symmetry_field_ridges.size(); i++)
	{
		is_ridge[symmetry_field_ridges[i]] = true;
	}

	vector<int> new_ridges;
	vector<int> highest_neighbors;
	//vector<int> 



	while (1)
	{
		new_ridges.clear();
		for (int i = 0; i < symmetry_field_ridges.size(); i++)
		{
			int ridge_index = symmetry_field_ridges[i];
			float max_value = symmetry_field[ridge_index];
			float second_max_value = symmetry_field[ridge_index];

			//find the highest neighbor for path tracing 
			for (int j = 0; j < medialSurface_net[ridge_index].size(); j++)
			{
				int neighbor_index = medialSurface_net[ridge_index][j];
				if (symmetry_field[neighbor_index] > max_value)
				{
					max_value = symmetry_field[neighbor_index];
				}
			}
			/*for (int j = 0; j < medialSurface_net[ridge_index].size(); j++)
			{
				int neighbor_index = medialSurface_net[ridge_index][j];
				if (symmetry_field[neighbor_index] == max_value)
				{
					continue;
				}

				if (symmetry_field[neighbor_index] > second_max_value )
				{
					second_max_value = symmetry_field[neighbor_index];
				}
			}*/

			highest_neighbors.clear();
			for (int j = 0; j < medialSurface_net[ridge_index].size(); j++)
			{
				int neighbor_index = medialSurface_net[ridge_index][j];
				if (symmetry_field[neighbor_index] == max_value && symmetry_field[neighbor_index]> symmetry_field[ridge_index])
				{
					highest_neighbors.push_back(neighbor_index);
				}
			}

			//save highest neighbors as new ridge points
			for (int j = 0; j < highest_neighbors.size(); j++)
			{
				int highest_neighbor_index = highest_neighbors[j];
				if (is_ridge[highest_neighbor_index] == false)
				{
					new_ridges.push_back(highest_neighbor_index);
					is_ridge[highest_neighbor_index] = true;
				}
			}
		}

		//terminal condition
		if (new_ridges.empty())
		{
			break;
		}

		for (int i = 0; i < new_ridges.size(); i++)
		{
			symmetry_field_ridges.push_back(new_ridges[i]);
		}
		cout << "new ridge point num: " << new_ridges.size() << endl;
	}

	delete[] is_ridge;
	int path_tracing_time_toc = clock();
	cout << "path tracing time: " << path_tracing_time_toc - path_tracing_time_tic << " ms" << endl;

#endif
}

#include<queue>
void Skeleton::spfa(vector<vector<pair<int, float>> >& graph, int src, long double* d, int* pre)
{
	queue<int> que;
	que.push(src);
	bool* inq = new bool[medialSurface.size()];
	memset(inq, false, medialSurface.size());
	inq[src] = true;
	memset(pre, -1, medialSurface.size()*sizeof(int));
	for (int i = 0; i < medialSurface.size(); i++)
	{
		d[i] = -1;
	}
	d[src] = 0;

	while (!que.empty())
	{
		int u = que.front();
		que.pop();
		inq[u] = false;

		for (int i = 0; i < graph[u].size(); i++)
		{
			int v = graph[u][i].first;
			float dist = graph[u][i].second;

			if (d[v]==-1 || d[v]>d[u] + dist)
			{
				d[v] = d[u] + dist;
				pre[v] = u;

				if (inq[v] == false)
				{
					que.push(v);
					inq[v] = true;
				}
			}
		}
	}

	delete[] inq;
}
