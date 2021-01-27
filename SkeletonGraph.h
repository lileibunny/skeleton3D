#ifndef SkeletonGraph_H
#define SkeletonGraph_H
#include <vector>
#include "Vec.h"
using namespace std;


class SkeletonGraph
{
public:
	SkeletonGraph(){}
	SkeletonGraph(int n);
	SkeletonGraph(vector<Point3f>& keyPoints,vector<float>& clusterDT, vector<int>* graph);
	~SkeletonGraph();
	int vn;
	vector<int>* graph;
	vector<Point3f> vertices;
	vector<float> DT;
	bool* onLoop;
	bool* vist;
	bool dfs(int s,int u,int nVist);
	bool hasLoop();
public:
	void getOnLoop();
};

#endif
