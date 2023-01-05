#include<string>
#include<iostream>
using namespace std;

#include "timestamp.h"
#ifndef UGRAPH_H
#define UGRAPH_H

struct edge{
	int source;
	int neighbor; // adjacent node
	int w; //keeps auxiliary information
	edge(){
		source = 0;
		neighbor = 0;
		w = 0;
	};
	edge(int i, int j){
		source = 0;
		neighbor = i;
		w = j;
	};
	edge(int from, int to, int aweight){
		source = from;
		neighbor = to;
		w = aweight;
	};
};

class Ugraph{
public:
	Ugraph(int N);
	void backtrack(int s, int r){
		if(s == r) return;
		if(parents[r] == r) return;
		for(size_t i = 0; i < Adj[parents[r]].size(); i++){
			if(Adj[r][i].neighbor == parents[r]){
				Adj[r][i].w = -1;
				//cout << "changed weight at " << r << " vertex's edge " << i << " to -1" << endl;
			}
		}
		backtrack(s, parents[r]);
	}
	//assignment 13
	void dijkstra(int s);
	void printDistance();
	void printParents();
	void printPath(int u, int v);
	int lenghtShortestW(int u, int v);
	void printMST();
	int weightMST();
	void kruskal();
	int findSet(int v);
	void combine(int x, int y);
	void collectEdges(vector<edge> &edgesAll);
	void backtrackFill(int s, int r, vector<int>& apath){
		if(s == r) return;
		if(parents[r] == r) return;
		int parent = parents[r];
		for(size_t i = 0; i < Adj[parents[r]].size(); i++){
				//cout << Adj[r][i].neighbor << " ";
				if(Adj[parent][i].neighbor == r){
					//cout << Adj[r][i].neighbor << " should be pushed back. ";
					backtrackFill(s, parents[r], apath);
					Adj[parent][i].w = -1;
					//cout << "should be pushing back " << r << " vertex's edge and changing weight of vertex with edge " << i << " to -1" << endl;
					apath.push_back(Adj[parent][i].neighbor);
					
					//cout << "changed weight at " << r << " vertex's edge " << i << " to -1" << endl;
				}
			
		}

		
	}
	void bfs(int s);
	void dfs();
	bool bfsDP(int s, vector<int>& apath, int v);
	void dfsDPVisit(int u, int &t, vector<int>& apath);
	void dfsLC(int &max);
	void dfsLCVisit(int u, int &t, int &length);
	void dfsVisit(int u, int &t);
	void printGraph();
	void addEdge(int u, int v);
	void addEdge(int u, int v, int weight);
  void removeEdge(int u, int v);
	//Toy problem 1 (by Elena)
	bool isCycle(); //returns true if there is a cycle in a graph
	bool isCycleVisit(int u, int &t);//similar to dfsVisit
	//collects the nodes on the shortest path form the source s to node v.
	void collectPath(int v, vector<int> &apath);
	void collectPathNoSource(int s, int v, vector<int> &apath);
  //problems for proj 3
  bool distinctPaths(int u, int v);
  void printBridges();
  void bridgeDFS(int u, vector<bool>& visit, vector<int>& discovery, vector<int>& evl, int &t);
  void printCC();
  //implement using dfs and dfsVisit
  //whenever a current node has the same color as grey or black node
  //return false
  //color white neighbor into opposite to the current node color
  bool twoColoring();
  bool cycleCheck(int& cNum){
	for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
	}	
	int t = 0;

	for(int i = 0; i < size; i++){
		if(colors[i] == 'W'){
			//color[i] = 'G';
			bool res = cycleCheckVisit(i, t, cNum);
			if(res)
				return true;
		}//if
	}//for
	return false;

} //returns true if there is a cycle in a graph
bool cycleCheckVisit(int u, int &t, int& cNum){
	colors[u] = 'G';
	stamps[u].d = t;
	t++;
	for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W'){
				parents[v] = u;
				colors[v] = 'G';
				cNum++;
				//cout << cNum << endl;
				bool res = cycleCheckVisit(v, t, cNum);
				if(res){
					return true;
				}
			} else if(colors[v] == 'G'){
				//cout << "finished a cycle and cnum is "<< cNum << endl;
				return true;
			}
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
	return false;
}//similar to dfsVisit
  //helper functions from lecture
 // int countCC();
  //void countCCVisit(int u, int &t, int totalCC, vector<int> &cc);

  //problems for Assignment 12
	//Problem 1
	bool sameCycle(int s, int r);
	
    //Problem 2.
	int longestCycle(int s);
	
    //Problem 3
	bool twoPaths(int s, int r);
	
    //Problem 4
	bool isOnPath(int s, int r, int q);//returns true if q is on the shortest path from s to r
	int evenCycle();
	int evenCycleVisit(int u, int &t);
	void printCycle(int start, int cur);
	int countCC();
	void countCCVisit(int u, int &t, int adjI, vector<int> &goingle);
private:
	vector<vector<edge>> mst;//adjacency lists of mst
	vector< vector<edge> > Adj; //adjacency lists of the graph 
	vector<int> distance; //for BFS and DFS
	vector<int> parents; //for BFS and DFS
	vector<char> colors; //for DFS
	vector<TimeStamp> stamps; //for DFS: stamps[v].d returns discovery time of v, and stamps[v].f finishing time.
    int size;
    vector<char> redBlue;
};

#endif
