#include "skeleton.h"
#include "stdio.h"
#include "graph.h"
#include "sphere.h"
#include "burning_sphere.h"

#include <algorithm>
#include <ctime>
#include <queue>
using namespace std;


Skeleton::Skeleton()
{
	
}

void Skeleton::init()
{
	symmetry_field_ridges.clear();
	skeleton_spheres_generation_flag = false;
	symmetry_field_generation_flag = false;

	skeleton_sphere_indexs.clear();
	skeletonEdges.clear();
	fire_front_space.clear();
}

void Skeleton::skeleton_spheres_generation()
{
	//****1. select fewest  skeleton spheres to cover all ridges through a greedy alogrithm of the minimum set cover problem
	skeleton_spheres.clear();
	covered_points_of_sphere.clear();
	//prepare datas 
	int n = symmetry_field_ridges.size();
	int m = n;
	//max_symmetry_field_value++;
	float shrinking_offset = 0;
	float lamda = 0;
	float shrinking_scale = 1;
	float minimal_radius = 3;

	float density_factor = 0;  //smaller value results in denser distribution of skeleton spheres
	float centeredness_factor = 100; 

	bool* is_covered = new bool[medialSurface.size()];
	memset(is_covered, 0, medialSurface.size());

	vector<int>* sphere_points = new vector<int>[medialSurface.size()];
	float* centeredness_metric = new float[medialSurface.size()];


	for (int i = 0; i < n; i++)
	{
		Point3f sphere_center = medialSurface[symmetry_field_ridges[i]];
		float dist_field_value = sqrt(dist->at_dist(sphere_center[0], sphere_center[1], sphere_center[2]));
		float symmetry_field_value=0;
		if (lamda!=0)
		{
			symmetry_field_value = max_symmetry_field_value - symmetry_field[i];
		}
		if (symmetry_field_value >= max_symmetry_field_value)
		{
			symmetry_field_value = 0;
		}

		float sphere_radius = (symmetry_field_value*lamda + dist_field_value*(1 - lamda)) * shrinking_scale;
		if (symmetry_field_value>sphere_radius)
		{
			//sphere_radius = symmetry_field_value;
		}

		sphere_radius -= shrinking_offset;
		sphere_radius = max(minimal_radius, sphere_radius);

		Point3f center_of_neighbors(0, 0, 0);
		float neighbor_num = 0;
		for (int j = 0; j < n; j++)
		{
			int ridge_index = symmetry_field_ridges[j];
			float dist = (sphere_center - medialSurface[ridge_index]).length()-0.001;
			if (dist <= sphere_radius)
			{   
				sphere_points[symmetry_field_ridges[i]].push_back(ridge_index);
				center_of_neighbors += medialSurface[ridge_index];
				neighbor_num++;
			}
		}
		if (neighbor_num>0)
		{
			center_of_neighbors /= neighbor_num;
		}

		centeredness_metric[symmetry_field_ridges[i]] = (medialSurface[symmetry_field_ridges[i]] - center_of_neighbors).length();
	}

	//main loop of the greedy algorithm
	vector<int> pending_points = symmetry_field_ridges;
	skeleton_sphere_indexs.clear();
	skeleton_sphere_centers.clear();

	vector<int> new_pending_points;

	int greedy_algorithm_time_tic = clock();
	while (pending_points.size())
	{
		float max_key_value = -(1 << 30);
		int optimal_MIS__index = -1;

		for (int i = 0; i < pending_points.size(); i++)
		{
			int covered_point_num = 0;
			int uncovered_point_num = 0;
			Point3f sphere_center = medialSurface[pending_points[i]];
			float symmetry_field_value = 0;
			if (lamda != 0)
			{
				symmetry_field_value = max_symmetry_field_value - symmetry_field[pending_points[i]];
			}
			if (symmetry_field_value >= max_symmetry_field_value)
			{
				symmetry_field_value = 0;
			}
			float dist_field_value = sqrt(dist->at_dist(sphere_center[0], sphere_center[1], sphere_center[2]));
			float sphere_radius = (symmetry_field_value*lamda + dist_field_value*(1 - lamda)) * shrinking_scale;
			if (symmetry_field_value>sphere_radius)
			{
				//sphere_radius = symmetry_field_value;
			}
			
			sphere_radius -= shrinking_offset;
			sphere_radius = max(minimal_radius, sphere_radius);
			

			for (int j = 0; j < sphere_points[pending_points[i]].size(); j++)
			{
				int ridge_index = sphere_points[pending_points[i]][j];
				
				if (is_covered[ridge_index])
				{
					covered_point_num++;
				}
				else
				{
					uncovered_point_num++;
				}
				
			}

			float current_key_value = 1 * uncovered_point_num - centeredness_factor*centeredness_metric[pending_points[i]];  // the cost function
			if (current_key_value > max_key_value)
			{
				max_key_value = current_key_value;
				optimal_MIS__index = pending_points[i];
			}
		}

		bool finish_flag = true;

		//flag current optimal sphere and carry out covering
		if (optimal_MIS__index != -1)
		{
			skeleton_sphere_indexs.push_back(optimal_MIS__index);
			skeleton_sphere_centers.push_back(medialSurface[optimal_MIS__index]);

			Point3f sphere_center = medialSurface[optimal_MIS__index];
			float symmetry_field_value = 0;
			if (lamda != 0)
			{
				symmetry_field_value = max_symmetry_field_value - symmetry_field[optimal_MIS__index];
			}
			if (symmetry_field_value >= max_symmetry_field_value)
			{
				symmetry_field_value = 0;
			}
			float dist_field_value = sqrt(dist->at_dist(sphere_center[0], sphere_center[1], sphere_center[2]));
			float sphere_radius = (symmetry_field_value*lamda + dist_field_value*(1 - lamda)) * shrinking_scale;
			if (symmetry_field_value>sphere_radius)
			{
				//sphere_radius = symmetry_field_value;
			}
			
			sphere_radius -= shrinking_offset;
			sphere_radius = max(minimal_radius, sphere_radius);


			skeleton_spheres.push_back(pair<Point3f, float>(sphere_center, sphere_radius));
			vector<int> covered_points;
			new_pending_points.clear();
			for (int i = 0; i < pending_points.size(); i++)
			{
				float dist = (sphere_center - medialSurface[pending_points[i]]).length() - density_factor;
				if (dist <= sphere_radius)
				{
					//cover it
					if (!is_covered[pending_points[i]])
					{
						is_covered[pending_points[i]] = true;
						covered_points.push_back(pending_points[i]);
					}
				}
				else// if (pending_points[i] != optimal_MIS__index)
				{
					new_pending_points.push_back(pending_points[i]);
				}

				if (!is_covered[pending_points[i]])
				{
					finish_flag = false;
				}
			}
			covered_points_of_sphere.push_back(covered_points);

		}

		if (finish_flag)
		{
			break;
		}

		//update pending sphere set
		pending_points = new_pending_points;

	}
	int greedy_algorithm_time_toc = clock();
	cout << "skeleton spheres generation finished!  used time: " << (greedy_algorithm_time_toc - greedy_algorithm_time_tic) << " ms.  "<<endl;
	cout << "skeleton spheres num: " << skeleton_spheres.size() << endl;
	skeleton_spheres_generation_flag = true;

	delete[] is_covered;


	//get connectivity information of skeleton spheres
	vector<int> point_to_skeleton_sphere;
	for (int i = 0; i < medialSurface.size(); i++)
	{
		point_to_skeleton_sphere.push_back(-1);
	}
	for (int i = 0; i < skeleton_spheres.size(); i++)
	{
		for (int j = 0; j < covered_points_of_sphere[i].size(); j++)
		{
			point_to_skeleton_sphere[covered_points_of_sphere[i][j]] = i;
		}
	}

	vector<bool> vist;
	skeleton_sphere_neighbors.clear();
	for (int i = 0; i < skeleton_spheres.size(); i++)
	{
		vist.clear();
		for (int j = 0; j < skeleton_spheres.size(); j++)
		{
			vist.push_back(false);
		}
		skeleton_sphere_neighbors.push_back(*(new vector<int>));
		for (int j = 0; j < covered_points_of_sphere[i].size(); j++)
		{
			int covered_point_index = covered_points_of_sphere[i][j];
			for (int k = 0; k < medialSurface_net[covered_point_index].size(); k++)
			{
				int neighbor_index = medialSurface_net[covered_point_index][k];
				int adjacent_sphere = point_to_skeleton_sphere[neighbor_index];
				if (adjacent_sphere != i && adjacent_sphere != -1)
				{
					if (vist[adjacent_sphere] == false)
					{
						skeleton_sphere_neighbors[i].push_back(adjacent_sphere);
						vist[adjacent_sphere] = true;
					}
				}
			}
		}
	}

}

void Skeleton::skeleton_spheres_regularization()
{
	cout << "skeleton_spheres_regularization(...)" << endl;

	int n = skeleton_sphere_indexs.size();
	float regularization_radius = 3;


	vector<Point3f> average_term;
	vector<Point3f> repulsion_force_term;
	for (int i = 0; i < n; i++)
	{
		float sphere_radius = skeleton_spheres[i].second;
		Point3f sphere_center = skeleton_spheres[i].first;

		int sphere_center_index=skeleton_sphere_indexs[i];;
		regularization_radius = sphere_radius+1;

		double weightSum = 0;
		Point3f cur_average_term(0, 0, 0);

		
		cur_average_term=Point3f(0, 0, 0);
		weightSum = 0;

		for (int j = 0; j < symmetry_field_ridges.size(); j++)
		{
			//cout << "ok" << endl;
			int ridge_index =  symmetry_field_ridges[j];
			float dist = (sphere_center - medialSurface[ridge_index]).length();
			if (dist <= regularization_radius)
			{
				float weight = 1.0;
				double mi = 4 * dist*dist / sphere_radius / sphere_radius;
				weight /= exp(mi); //weight = pow(symmetry_field[ridge_index],0);
				weight = 1;
				weightSum += weight;
				cur_average_term += weight*medialSurface[ridge_index];
			}
		}
	
	
		//for (int j = 0; j < covered_points_of_sphere[i].size(); j++)
		//{
		//	int covered_point_index = covered_points_of_sphere[i][j]; 
		//	float dist = (sphere_center - skeleton_points[covered_point_index]).length();
		//	
		//	float weight = 1.0;
		////	double mi = 4 * dist*dist / sphere_radius / sphere_radius;
		//	//weight /= exp(mi); weight = 1;
		//	weightSum += weight;
		//	cur_average_term += weight*skeleton_points[covered_point_index];
		//	
		//}


		/*for (int j = 0; j < skeleton_points.size(); j++)
		{
			float dist = (sphere_center - skeleton_points[j]).length();
			if (dist <= sphere_radius)
			{
				float weight = 1.0;
				double mi = 4 * dist*dist / sphere_radius / sphere_radius;
				weight /= exp(mi); 
				weightSum += weight;
				cur_average_term += weight*skeleton_points[j];
			}
		}
*/
		if (weightSum > 0)
		{
			cur_average_term /= weightSum;
		}
		average_term.push_back(cur_average_term);


		Point3f cur_repulsion_force(0, 0, 0);
		float weight_sum = 0;
		for (int j = 0; j < skeleton_spheres.size(); j++)
		{
			if (j == i)
			{
				continue;
			}
			/*float neighbor_sphere_radius = skeleton_spheres[skeleton_sphere_neighbors[i][j]].second;
			Point3f neighbor_sphere_center = skeleton_spheres[skeleton_sphere_neighbors[i][j]].first;*/

			float neighbor_sphere_radius = skeleton_spheres[j].second;
			Point3f neighbor_sphere_center = skeleton_spheres[j].first;


			float dist = (sphere_center - neighbor_sphere_center).length();
			float radius_sum = sphere_radius + neighbor_sphere_radius;
			float intersection_diameter = radius_sum - dist;
			float weight = intersection_diameter;
			//weight /= radius_sum;
			if (intersection_diameter>0)
			{
				Vec<3, float> repulsion_force_direction = (sphere_center - neighbor_sphere_center) / dist;
				cur_repulsion_force += weight*repulsion_force_direction;
				weight_sum += weight;
			}
			else
			{
				/*Vec<3, float> repulsion_force_direction = (neighbor_sphere_center - sphere_center) / dist;
				weight *= -2;
				cur_repulsion_force += weight*repulsion_force_direction;
				weight_sum += weight;*/
			}
		}

		
		repulsion_force_term.push_back(cur_repulsion_force);

	}


	float move_dist = 0;
	for (int i = 0; i < n; i++)
	{
		float repulsion_force_power = 0.1;
		Point3f new_center = average_term[i] + repulsion_force_power*repulsion_force_term[i];
		move_dist += (new_center - skeleton_sphere_centers[i]).length();
		skeleton_sphere_centers[i] = new_center;
		skeleton_spheres[i].first = new_center;
	}

	move_dist /= n;
	cout << "cur_average_move_dist: " << move_dist << endl;

}


void Skeleton::skeleton_generation()
{
	vector<int> point_to_skeleton_sphere;
	for (int i = 0; i < medialSurface.size(); i++)
	{
		point_to_skeleton_sphere.push_back(-1);
	}
	for (int i = 0; i < skeleton_spheres.size(); i++)
	{
		for (int j = 0; j < covered_points_of_sphere[i].size(); j++)
		{
			point_to_skeleton_sphere[covered_points_of_sphere[i][j]] = i;
		}
	}

	//construct the skeleton as a graph structure
	skeletonNodes.clear();
	for (int i = 0; i < skeleton_spheres.size(); i++)
	{
		skeletonNodes.push_back(SkeletonNode(skeleton_spheres[i].first, skeleton_spheres[i].second));
	}

	skeletonEdges.clear();
	int m = skeletonNodes.size();
	set<int>* neighbor_nodes = new set<int>[m];

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < covered_points_of_sphere[i].size(); j++)
		{
			int covered_point_index = covered_points_of_sphere[i][j];
			for (int k = 0; k < medialSurface_net[covered_point_index].size(); k++)
			{
				int neighbor_index = medialSurface_net[covered_point_index][k];
				int adjacent_sphere = point_to_skeleton_sphere[neighbor_index];
				if (adjacent_sphere != i&& adjacent_sphere != -1)
				{
					if (neighbor_nodes[i].find(adjacent_sphere) == neighbor_nodes[i].end())
					{
						skeletonNodes[i].neighbors.push_back(adjacent_sphere);
						neighbor_nodes[i].insert(adjacent_sphere);
					}
					
					skeletonEdges.push_back(pair<int, int>(i, adjacent_sphere));
				}
			}
		}
		//continue;

		for (int j = i + 1; j < m; j++)
		{
			bool is_connected = false;

			/*	for (int ii = 0; ii < skeletonNodes[i].covered_centers.size(); ii++)
			{
			for (int jj = 0; jj < skeletonNodes[j].covered_centers.size(); jj++)
			{

			float dist = (skeletonNodes[i].covered_centers[ii] - skeletonNodes[j].covered_centers[jj]).length();
			if (dist <= 1.01)
			{
			is_connected = true;
			break;
			}
			}

			if (is_connected)
			{
			break;
			}
			}*/

			float dist = (skeletonNodes[i].position - skeletonNodes[j].position).length() - 0.1;
			float radius_sum = skeletonNodes[i].radius + skeletonNodes[j].radius;
			if (dist <= radius_sum)
			{
				//is_connected = true;
			}

			if (is_connected)
			{
				skeletonEdges.push_back(pair<int, int>(i, j));
				if (neighbor_nodes[i].find(j) == neighbor_nodes[i].end())
				{
					skeletonNodes[i].neighbors.push_back(j);
					neighbor_nodes[i].insert(j);
				}

				skeletonEdges.push_back(pair<int, int>(j, i));
				if (neighbor_nodes[j].find(i) == neighbor_nodes[j].end())
				{
					skeletonNodes[j].neighbors.push_back(i);
					neighbor_nodes[j].insert(i);
				}
			}
		}
	}

	skeleton_refinement();
	cout << "skeleton generation results: skeletonNode num= " << skeletonNodes.size() << ",  skeletonEdge num= " << skeletonEdges.size() / 2 << endl;
}


void Skeleton::skeleton_refinement()
{
	//***********1. topology correction **************/
	//remove all triangles in the skeleton graph structure
	int n = skeletonNodes.size();
	int m = skeletonEdges.size();
	bool* is_obtuseAngle_edge = new bool[m];
	memset(is_obtuseAngle_edge, 0, m);
	for (int i = 0; i < m; i++)
	{
		int u = skeletonEdges[i].first;
		int v = skeletonEdges[i].second;
		float edge_length = (skeletonNodes[u].position - skeletonNodes[v].position).length();

		for (int j = 0; j < n; j++)
		{
			if (j == u || j == v)
			{
				continue;
			}

			bool adjacent_to_u = false;
			bool adjacent_to_v = false;

			for (int k = 0; k < skeletonNodes[j].neighbors.size(); k++)
			{
				if (skeletonNodes[j].neighbors[k] == u)
				{
					adjacent_to_u = true;
				}
				if (skeletonNodes[j].neighbors[k] == v)
				{
					adjacent_to_v = true;
				}
			}

			if (adjacent_to_u && adjacent_to_v)
			{
				Vec<3, float> triangle_edge_vector1 = skeletonNodes[u].position - skeletonNodes[j].position;
				Vec<3, float> triangle_edge_vector2 = skeletonNodes[v].position - skeletonNodes[j].position;

				float dot_product = triangle_edge_vector1 DOT triangle_edge_vector2;

				if (edge_length>triangle_edge_vector1.length() && edge_length>triangle_edge_vector2.length())
				{
					//delete edge(u,v) from the set of skeleton edges 
					is_obtuseAngle_edge[i] = true;

					//delete edge(u,v) from the skeleton graph structure 
					int erase_id = -1;
					for (int k = 0; k < skeletonNodes[u].neighbors.size(); k++)
					{
						int neighbor_id = skeletonNodes[u].neighbors[k];
						if (neighbor_id == v)
						{
							erase_id = k;
							break;
						}
					}

					if (erase_id != -1)
					{
						skeletonNodes[u].neighbors.erase(skeletonNodes[u].neighbors.begin() + erase_id);
					}

					/*	erase_id = -1;
					for (int k = 0; k < skeletonNodes[v].neighbors.size(); k++)
					{
					int neighbor_id = skeletonNodes[v].neighbors[k];
					if (neighbor_id == u)
					{
					erase_id = k;
					break;
					}
					}
					if (erase_id != -1)
					{
					skeletonNodes[v].neighbors.erase(skeletonNodes[v].neighbors.begin() + erase_id);
					}*/

				}
			}
		}
	}
	vector<pair<int, int>> temp_edges;
	for (int i = 0; i < m; i++)
	{
		if (is_obtuseAngle_edge[i] == false)
		{
			temp_edges.push_back(skeletonEdges[i]);
		}
	}
	skeletonEdges.clear();
	skeletonEdges = temp_edges;

	//*************2. simplification
	while (0)
	{
		float shortest_edge_length = 1 << 30;
		int shortest_length_id = -1;
		for (int i = 0; i < m; i++)
		{
			int u = skeletonEdges[i].first;
			int v = skeletonEdges[i].second;

			float edge_length = (skeletonNodes[u].position - skeletonNodes[v].position).length();
			float edge_degree = skeletonNodes[u].neighbors.size() + skeletonNodes[v].neighbors.size();

			if (edge_length <= 0.001 || skeletonNodes[u].neighbors.size() != 2 || skeletonNodes[v].neighbors.size() != 2)
			{
				continue;
			}

			if (edge_length < shortest_edge_length)
			{
				shortest_edge_length = edge_length;
				shortest_length_id = i;
			}
		}

		if (shortest_edge_length >= 5)
		{
			break;
		}

		int u = skeletonEdges[shortest_length_id].first;
		int v = skeletonEdges[shortest_length_id].second;

		Point3f middle_point = (float)0.5*(skeletonNodes[u].position + skeletonNodes[v].position);
		skeletonNodes[u].position = middle_point;
		skeletonNodes[v].position = middle_point;
	}

	//*************3.enhance the tails of skeleton branches
	for (int i = 0; i < n; i++)
	{
		if (skeletonNodes[i].neighbors.size() == 0)
		{
			int nearest_node;
			float min_dist = 1 << 30;
			for (int j = 0; j < n; j++)
			{
				if (j == i)
				{
					continue;
				}
				float dist = (skeletonNodes[i].position - skeletonNodes[j].position).length();
				if (dist < min_dist)
				{
					min_dist = dist;
					nearest_node = j;
				}
			}

			skeletonNodes[i].neighbors.push_back(nearest_node);
			skeletonEdges.push_back(pair<int, int>(i, nearest_node));
		}
	}

	for (int i = 0; i < n; i++)
	{
		if (skeletonNodes[i].neighbors.size() == 1)
		{
			Vec<3, float> skeleton_curve_direction = skeletonNodes[i].position - skeletonNodes[skeletonNodes[i].neighbors[0]].position;
			float radius1 = skeletonNodes[i].radius;
			float radius2 = skeletonNodes[skeletonNodes[i].neighbors[0]].radius;
			float dist = skeleton_curve_direction.length();
			skeletonNodes[i].position += (float)(radius1/1.5)*(skeleton_curve_direction / skeleton_curve_direction.length());
		}
	}
}


void Skeleton::find_quatrangle(vector<int> current_path, int current_node, int father_node)
{
	for (int i = 0; i < current_path.size(); i++)
	{
		if (current_path[i] == current_node)
		{
			int loop_size = current_path.size() - i;
			if (loop_size <= 4)
			{
				//find a quatrangle and remove the longest edge
			}
		}
	}
}