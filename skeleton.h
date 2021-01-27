#pragma once
#include "SkeletonGraph.h"
#include "burning_sphere.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class Skeleton
{

public:
	class CurveSkeleton
	{
	public:
		vector<Point3f> nodes;
		vector<pair<int, int>> edges;
		vector<vector<int>> neighbors;
	};

	CurveSkeleton* curSkeleton;

	struct SkeletonNode
	{
		Point3f position;
		float radius;
		vector<int> neighbors;
		vector<Point3f> covered_centers;

		SkeletonNode(Point3f position, float radius, vector<int> neighbors = *(new vector<int>), vector<Point3f> covered_centers = *(new vector<Point3f>))
		{
			this->position = position;
			this->radius = radius;
			this->neighbors = neighbors;
			this->covered_centers = covered_centers;
		}
		SkeletonNode()
		{
			this->position = Point3f(0, 0, 0);
			this->radius = 0;
		}
	};

public:
		int rows;
		int cols;
		int heis;
		int size;
		double scale;
		double constVal;
		float max_symmetry_field_value;

		Point3f* gradientField;
		vector<Point3f> symmetry_field_gradient;
		vector<int> symmetry_field_ridges;
		vector<vector<int>> symmetry_field_ridges_net;
		vector<int> burning_path;
		vector<vector<int>> burning_paths;
		vector<vector<vector<int>>> fire_front_space;
		vector<vector<int>> drawing_paths;

		void dfs(Point3f root);
		void dfs_in_ridges(int root);
		bool is_noise_segment(int* path_pre, vector<vector<int>>& front_segments, int segment_index);
		void find_quatrangle(vector<int> current_path,int current_node,int father_node);
		
    public:
		vector<Point3f> medialSurface;
		vector<float> symmetry_field;
		vector<Point3f> border;
		vector<int> skeleton_sphere_indexs;
		vector<Point3f> skeleton_sphere_centers;
		vector < pair < Point3f, float> >  skeleton_spheres;
		vector<vector<int>> covered_points_of_sphere;
		vector<BurningSphere> burning_spheres;
		vector<int> special_leaves;

		bool skeleton_spheres_generation_flag;
		bool symmetry_field_generation_flag;
		vector<pair<int, int>> skeletonEdges;
		vector<SkeletonNode> skeletonNodes;
		vector<vector<int>> skeleton_sphere_neighbors;
		vector<float> skeleton_radii;
		vector<Point3f> skeleton_points;

		vector<Point3i> vec;

		vector<int> shortest_path;

        ~Skeleton();
		void freeMemory();


		int      distance_field_ridges_detection(Point3i pnt);
		Point3f* distance_field_gradient_generation();
		
		
		void init();

		void symmetry_field_generation();
		void symmetry_field_gradient_generation();
		void symmetry_field_ridges_detection();
		void symmetry_field_ridges_connection();
		void spfa(vector<vector<pair<int, float>> >& graph, int src, long double* d, int* pre);
		void symmetry_field_smoothing();

		void skeleton_spheres_generation();
		void skeleton_spheres_generation2();
		void skeleton_spheres_regularization();
		void skeleton_generation();
		void skeleton_refinement();

};

