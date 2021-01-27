#include"SkeletonGraph.h"

SkeletonGraph::SkeletonGraph(int n)
{
	vn = n;
	graph = NULL;
	graph = new vector<int>[vn];
	vist = NULL;
	onLoop = NULL;
}
SkeletonGraph::SkeletonGraph(vector<Point3f>& keyPoints, vector<float>& clusterDT,vector<int>* graph)
{
	vertices = keyPoints;
	DT = clusterDT;
	this->graph = graph;
	vn = keyPoints.size();
	vist = NULL;
	onLoop = NULL;
}
SkeletonGraph::~SkeletonGraph()
{
	if (graph != NULL)
	{ 
		delete[] graph;
	}

	if (vist != NULL)
	{
		delete[] vist;
	}
	if (onLoop != NULL)
	{
		delete[] onLoop;
	}
}
void SkeletonGraph::getOnLoop()
{
	onLoop = new bool[vn];
	vist = new bool[vn];

	for (int i = 0; i < vn; i++)
	{
		memset(vist, 0, vn);

		if (dfs(i,i,0))
		{
			onLoop[i] = true;
		}
		else
		{
			onLoop[i] = false;
		}
	}
}

bool SkeletonGraph::dfs(int s,int u,int nVist)
{
	nVist++;
	vist[u] = true;

	int nNode = graph[u].size();
	for (int i = 0; i < nNode; i++)
	{
		int v = graph[u][i];

		if (v == s && nVist>=3)
		{
			return true;
		}

		if (vist[v] == false && dfs(s,v,nVist))
		{
			return true;
		}
	}

	return false;
}

bool SkeletonGraph::hasLoop()
{
	for (int i = 0; i < vn; i++)
	{
		memset(vist, 0, vn);

		if (dfs(i, i, 0))
		{
			return true;
		}
	}
	return false;
}