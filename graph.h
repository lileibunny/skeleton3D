#pragma once
#include "Vec.h"
#include <vector>
using namespace std;

template <class Node_Type>
class Graph
{
public:
	vector<Node_Type> nodes;
	vector< vector<int> > node_neighbors;

	Graph(){ }

	Graph(vector<Node_Type>& nodes, vector< vector<int> >& neighbors)
	{
		this->nodes = nodes;
		this->node_neighbors = neighbors;
	}
	~Graph()
	{

	}

	void dfs(int node_id, bool* vist)
	{
		//visit current node
		vist[node_id] = true;

		for (int i = 0; i < node_neighbors[node_id].size(); i++)
		{
			int neighbor_id = node_neighbors[node_id][i];

			if (vist[neighbor_id] == false)
			{
				dfs(neighbor_id, vist);
			}
		}
	}

	bool check_connectivity()
	{
		bool* vist = new bool[nodes.size()];
		memset(vist, 0, nodes.size());

		dfs(0, vist);
		for (int i = 0; i < nodes.size(); i++)
		{
			if (vist[i] == false)
			{
				delete[] vist;
				return false;
			}
		}

		delete[] vist;
		return true;
	}

};