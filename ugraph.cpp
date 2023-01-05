#include<climits>
#include<algorithm>
#include<queue>
#include<math.h>
#include<vector>
#include<iostream>
#include<string>
#include<utility>
#include"priorityqueue.h"
using namespace std;

#include "ugraph.h"

bool lessThan(const edge &x, const edge &y){
	return (x.w < y.w) ? true : false;
}

bool Ugraph::distinctPaths(int u, int v){
	vector<int> apath1;
	
	vector<int> apath2;
	for(size_t i = 0; i < Adj.size(); i++){
		distance[i] = INT_MAX;
		for(size_t j = 0; j < Adj[i].size(); j++){
			Adj[i][j].w = 0;
			
		}
	}
	bfsDP(u, apath1, v);
	if(distance[v] == INT_MAX){
			return false;
	}
	apath1.push_back(u);
	backtrackFill(u, v, apath1);
		bfsDP(u, apath2, v);
		if(distance[v] == INT_MAX){
			return false;
		}
		apath2.push_back(u);
		backtrackFill(u, v, apath2);
	for (size_t i = 0; i < apath1.size();i++){
                cout << apath1[i] << " ";
        }
        apath1.clear();
        cout << endl;
        for (size_t i = 0; i < apath2.size(); i++){
                cout << apath2[i] << " ";
        }
		apath2.clear();
		cout << endl;
	for(size_t i = 0; i < Adj.size(); i++){
			for(size_t j = 0; j < Adj[i].size(); j++){
				Adj[i][j].w = 0;
			}
	}
	return true;
}
bool Ugraph::bfsDP(int s, vector<int> &apath, int v){
	for(int i = 0; i < size; i++){
		distance[i] = INT_MAX;
		parents[i] = i;
	}//for
	distance[s] = 0;
	queue<int> aq;
	aq.push(s);
	while(!aq.empty()){
		int u = aq.front();
        aq.pop();
		for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(Adj[u][i].w != -1){
				if(distance[v] == INT_MAX){
					distance[v] = distance[u] + 1;
					parents[v] = u;
					aq.push(v);
				}
			}
		}//for
	}//while
	
	/*for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
	}	
	int t = 0;

	for(int i = s; i < size; i++){
		if(colors[i] == 'W'){
			//cout << Adj[i][0].neighbor << " ";
			//color[i] = 'G';
			dfsDPVisit(i, t, apath);
		}//if
	}//for*/
	return true;
}//dfs
void Ugraph::dfsDPVisit(int u, int &t, vector<int> &apath){
	colors[u] = 'G';
	stamps[u].d = t;
	t++;

	for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W' && Adj[u][i].w != -1){
				//cout << v << " ";
				parents[v] = u;
				colors[v] = 'G';
				//apath.push_back(v);
				Adj[u][i].w = -1;
				dfsDPVisit(v, t, apath);
			}
			
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
}//dfsVisit	
void Ugraph::printBridges(){
	static int t = 0;
	vector<int> discovery(size);
	vector<bool> visit(size);
	vector<int> evl(size);
	for(int i = 0; i < size; i++){
		parents[i] = i;
		visit[i] = false;
	}
	for(int i = 0; i < size; i++){
		if(!visit[i]){
			bridgeDFS(i, visit, discovery, evl, t);
		}
	}
}
void Ugraph::bridgeDFS(int u, vector<bool>& visit, vector<int>& discovery, vector<int>& evl, int &t){
	visit[u] = true;
	discovery[u] = evl[u] = ++t;
	for(size_t i = 0; i < Adj[u].size(); i++){
		int v = Adj[u][i].neighbor;
		//cout << "v is " << v << endl;
		if(!visit[v]){
			parents[v] = u;
			//cout << "b4 dfs" << endl;
			bridgeDFS(v, visit, discovery, evl, t);
			evl[u] = min(evl[u], evl[v]);
			//cout << "got past the dfs" << endl;
			if(evl[v] > discovery[u]){
				cout << u <<" " << v << endl;
			}
		}
		else if(v != parents[u]){
			evl[u] = min(evl[u], discovery[v]);
		}
	}
}
void Ugraph::printCC(){
	int c = countCC();
	//cout << c << endl;
}
bool Ugraph::twoColoring(){
	bool result = true;
	for(size_t i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
		distance[i] = 0;
	}
	if(isCycle()){
		int cNum = 0;
		cycleCheck(cNum);
		//cout << cNum << endl;
		if((cNum % 2 != 0)){
			return false;
		}
	}
	return result;
}
Ugraph::Ugraph(int N){

	size = N;
	mst.resize(size);
	Adj.resize(size);
	distance.resize(size);
	parents.resize(size);
	colors.resize(size);
	stamps.resize(size);
  redBlue.resize(size); //twoColorint problem
}//default

void Ugraph::addEdge(int u, int v){
	Adj[u].push_back(edge(v, 0));
  Adj[v].push_back(edge(u, 0));
}
void Ugraph::addEdge(int u, int v, int weight){
	Adj[u].push_back(edge(u, v, weight));
  Adj[v].push_back(edge(v, u, weight));
}
void Ugraph::removeEdge(int u, int v){
  //find neighbor v in Adj[u]
  for(int j = 0; j < (int)Adj[u].size(); j++){
    int aneighbor = Adj[u][j].neighbor;
    if(aneighbor == v){
      int last = (int)Adj[u].size() - 1;
      Adj[u][j] = Adj[u][last];//overwrite
      Adj[u].resize(last);
      break; //found v, exit for loop
    }
  }
  int sizeV = (int)Adj[v].size();
  for(int j = 0; j < sizeV; j++){
    int aneighbor = Adj[v][j].neighbor;
    if(aneighbor == u){
      Adj[v][j] = Adj[v][sizeV - 1];
      Adj[v].resize(sizeV - 1);
      break;
    }
  }
}
void Ugraph::printGraph(){
	for(int i = 0; i < size; i++){
		for(int j = 0; j < (int)Adj[i].size(); j++){
			int v = Adj[i][j].neighbor;
			cout << v << " " ;
		}//for j
		cout << endl;
	}
}//printUgraph

void Ugraph::bfs(int s){
	for(int i = 0; i < size; i++){
		distance[i] = INT_MAX;
		parents[i] = i;
	}//for
	distance[s] = 0;
	queue<int> aq;
	aq.push(s);
	while(!aq.empty()){
		int u = aq.front();
		cout << u << " ";
        aq.pop();

		for(int i = 0; i < (int)Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(distance[v] == INT_MAX){
				distance[v] = distance[u] + 1;
				parents[v] = u;
				aq.push(v);
			}
		}//for
	}//while
	cout << endl;
}//bfs
void Ugraph::dfsLC(int& max){
	for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
	}	
	int t = 0;
	int length = 1;
	for(int i = 0; i < size; i++){
		if(colors[i] == 'W'){
			//color[i] = 'G';
			dfsLCVisit(i, t, length);
		}//if
	}//for
	if(length > max) max = length;
}//dfs
void Ugraph::dfsLCVisit(int u, int &t, int &length){
	colors[u] = 'G';
	length = distance[u];
	stamps[u].d = t;
	t++;
	for(int i = 0; i < (int)Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W'){
				distance[v] = distance[u] + 1;
				parents[v] = u;
				colors[v] = 'G';
				length = distance[v];
				dfsLCVisit(v, t, length);
			} else if(colors[v] == 'G' && parents[v] == v){
				distance[v] = distance[u] + 1;
				length = distance[v];
				return;
			} else if(colors[v] == 'G'){

			}
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
}//dfsVisit
void Ugraph::dfs(){
	for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
	}	
	int t = 0;

	for(int i = 0; i < size; i++){
		if(colors[i] == 'W'){
			cout << i << " ";
			//color[i] = 'G';
			dfsVisit(i, t);
		}//if
	}//for
	cout << endl;
}//dfs
int Ugraph::countCC(){
	for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
	}	
	int t = 0;
	int totalCC = 0;
	vector<vector<int>> cc;
	vector<int> idv(size);
	int id = 0;
	for(int i = 0; i < size; i++){
		if(colors[i] == 'W'){
			//color[i] = 'G';
			totalCC++;
			vector<int> goingle;
			//goingle.push_back(i);
			//cout << "new goingle created with start value " << goingle[0] << endl;
			idv[i] = id;
			countCCVisit(i, t, id, idv);
			id = id + 1;
			//cc.push_back(goingle);
			
		}//if
	}//for
	for(int i = 0; i < id; i++){
		//sort(cc[i].begin(), cc[i].end());
		for(int j = 0; j < idv.size(); j++){
			if(idv[j] == i){
				cout << j << " ";
			}
			
		}
		cout << endl;
	}
	return totalCC;
}//cc
void Ugraph::countCCVisit(int u, int &t, int adjI, vector<int> &goingle){
	colors[u] = 'G';
	stamps[u].d = t;
	t++;
	for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W'){
				goingle[v] = adjI;
				//cout << "found node " << v << " to be pushed back on " << adjI << endl;
				//goingle.push_back(v);
				//cout << goingle[0] << endl; 
				parents[v] = u;
				colors[v] = 'G';
				countCCVisit(v, t, adjI, goingle);
			}
			
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
}//ccv	
void Ugraph::dfsVisit(int u, int &t){
	colors[u] = 'G';
	stamps[u].d = t;
	t++;

	for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W'){
				cout << v << " ";
				parents[v] = u;
				colors[v] = 'G';
				dfsVisit(v, t);
			}
			
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
}//dfsVisit	
void Ugraph::dijkstra(int s){
	for(int i = 0; i < size; i++){
		distance[i] = INT_MAX;
		parents[i] = i;
	}//for
	distance[s] = 0;
	Item* array = new Item[size];
	for(int i = 0; i < size; i++){
		array[i] = Item(i, distance[i]);
	}
	PriorityQueue aq(array, size);
	//pqueued aq
	//queue<int> aq;
	while(aq.getSize() > 0){
		Item anitem = aq.getMin();
		int u = anitem.key;
        aq.pop();
		if(distance[u] != INT_MAX){
		for(int i = 0; i < (int)Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if (aq.isKey(v)){
				int weight = Adj[u][i].w;
				int curDist = distance[u] + weight;
				if(distance[v] > curDist){
					distance[v] = curDist;
					parents[v] = u;
					aq.updatePriority(v, distance[v]);
				}
			}
		}//for
		}
	}//while
	//cout << endl;
}//dijkstra
bool Ugraph::isCycle(){
	for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
	}	
	int t = 0;

	for(int i = 0; i < size; i++){
		if(colors[i] == 'W'){
			//color[i] = 'G';
			bool res = isCycleVisit(i, t);
			if(res)
				return true;
		}//if
	}//for
	return false;

} //returns true if there is a cycle in a graph
bool Ugraph::isCycleVisit(int u, int &t){
	colors[u] = 'G';
	stamps[u].d = t;
	t++;
	for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W'){
				parents[v] = u;
				colors[v] = 'G';
				bool res = isCycleVisit(v, t);
				if(res){
					return true;
				}
			} else if(colors[v] == 'G'){
				//cycle has been found
				return true;
			}
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
	return false;
}//similar to dfsVisit
//collects the nodes on the shortest path form the source s to node v.
void Ugraph::collectPath(int v, vector<int> &apath){
	if(distance[v] == INT_MAX){
		return;
	}
	if(parents[v] == v){
		apath.push_back(v);
		return;
	}
	collectPath(parents[v], apath);
	apath.push_back(v);
}
void Ugraph::collectPathNoSource(int s, int v, vector<int> &apath){ //working
	bfs(s);
	if(distance[v] == INT_MAX){
		return;
	}
	if(colors[v] == 'm'){
		return;
	}
	if(parents[v] == v){
		apath.push_back(v);
		return;
	}
	//colors[v] = 'm';
	collectPathNoSource(s, parents[v], apath);
	apath.push_back(v);

}

int Ugraph::evenCycle(){
	for(int i = 0; i < size; i++){
		parents[i] = i;
		colors[i] = 'W';
		distance[i] = -1;
	}
	int t = 0;
	int res = -1;
	for(int i = 0; i < size; i++){
		if(colors[i] == 'W'){
			distance[i] = 0;
			int res = evenCycleVisit(i, t);
			if(res > 0) break;
		}
	}
	return res;
}
int Ugraph::evenCycleVisit(int u, int &t){
	colors[u] = 'G';
	stamps[u].d = t;
	t++;

	for(size_t i = 0; i < Adj[u].size(); i++){
			int v = Adj[u][i].neighbor;
			if(colors[v] == 'W'){
				parents[v] = u;
				colors[v] = 'G';
				distance[v] = distance[u] + 1;
				int res = evenCycleVisit(v, t);
				if(res > 0) return res;
			} else if(colors[v] == 'G'){
				int length = distance[u] - distance[v] + 1;
				if(length % 2 == 0){
					printCycle(v, u);
					return length; 
				}
			}
			
	}//for
	colors[u] = 'B';
	stamps[u].f = t;
	t++;
}//dfsVisit	
void Ugraph::printCycle(int start, int cur){
	if(cur==start){
		cout << cur << " ";
		return;
	}
	printCycle(start, parents[cur]);
	cout << cur << " ";
	
}
bool Ugraph::sameCycle(int s, int r){
	for(size_t i = 0; i < Adj.size(); i++){
		for(size_t j = 0; j < Adj[i].size(); j++){
			Adj[i][j].w = 0;
		}
	}
	bfs(s);
	if(!isCycle()){
		return false;}
	backtrack(s, r);
	bfs(r);
	if(distance[s] != INT_MAX){
		return true;
	}
	return false;
}
int Ugraph::longestCycle(int s){
	int maxLength = -1;
	if(!isCycle()) return 0;
	dfsLC(maxLength);
	return maxLength;
}
bool Ugraph::twoPaths(int s, int r){
	for(size_t i = 0; i < Adj.size(); i++){
		distance[i] = INT_MAX;
		for(size_t j = 0; j < Adj[i].size(); j++){
			Adj[i][j].w = 0;
			
		}
	}
	bfs(s);
	if(isCycle()){
		backtrack(s, r);
		for(size_t i = 0; i < Adj.size(); i++){
			for(size_t j = 0; j < Adj[i].size(); j++){
				Adj[i][j].w = 0;
			}
		}
		for(size_t i = 0; i < distance.size(); i++){
				distance[i] = INT_MAX;
				parents[i] = i;
		}
		bfs(s);
		
		if(distance[s] != INT_MAX){
			return false;
		} else{
			return true;
		}
	}
	return false;
}
bool Ugraph::isOnPath(int s, int r, int q){
	vector<int> apath;
	collectPathNoSource(s, r, apath);
	for(size_t i = 0; i < apath.size(); i++){
		if(apath[i] == q) return true;
	}
	return false;
}
//Assignment 13 other than dijkstra
void Ugraph::printDistance(){
	for(int i = 0; i < size; i++){
		cout << distance[i] << " ";
	}
	cout << endl;
}
void Ugraph::printParents(){
	for(int i = 0; i < size; i++){
			cout << parents[i] << " ";
		}
		cout << endl;
}
int Ugraph::lenghtShortestW(int u, int v){
	dijkstra(u);
	printPath(u, v);
	cout << endl;
	return distance[v];
}
void Ugraph::kruskal(){
	mst.resize(size);
	//distance array stores ranks of vertices
	for(int v = 0; v < size; v++){
		parents[v] = v;
		distance[v] = 0; //rank of v is 0
	}
	vector<edge> edgesAll;
	collectEdges(edgesAll);
	std::sort(edgesAll.begin(), edgesAll.end(), lessThan);
	int count = 0; //total edges added so far to MST
	for(auto& anedge : edgesAll){
		int u = anedge.source;
		int v = anedge.neighbor;
		int w = anedge.w;
		int represU = findSet(u);
		int represV = findSet(v);
		if(represU != represV){
			mst[u].push_back(edge(u, v, w));
			mst[v].push_back(edge(v, u, w));
			combine(represU, represV);
			count++;
			if(count == (size-1)){break;}
		}
	}

}
void Ugraph::collectEdges(vector<edge> & edgesAll){
	for(int u = 0; u < size; u++){
		for(int j = 0; j < (int)Adj[u].size(); j++){
			int v = Adj[u][j].neighbor;
			if(u < v){
				edgesAll.push_back(edge(u, v, Adj[u][j].w));
			}
		}//for j
	}//for u
}
int Ugraph::findSet(int v){
	if(v != parents[v]){
		parents[v] = findSet(parents[v]);
	}
	return parents[v];
}
void Ugraph::combine(int x, int y){
	if(distance[x] > distance[y]){
		parents[y] = x;
	}
	else{
		parents[x] = y;
		if(distance[x] == distance[y]){
			distance[y]++;
		}
	}
}
void Ugraph::printPath(int u, int v){
	if(u == v){
		cout << u << " ";
		return;
	}
	printPath(u, parents[v]);
	cout << v << " ";
}
void Ugraph::printMST(){
	for(int u = 0; u < size; u++){
		for(size_t j = 0; j < mst[u].size(); j++){
			cout << mst[u][j].neighbor << " ";
		}
		cout << endl;
	}
}
int Ugraph::weightMST(){
	int sum = 0;
	for(int u = 0; u < size; u++){
		for(size_t j = 0; j < mst[u].size(); j++){
			sum+= mst[u][j].w;
		}
	}
	sum = sum/2;
	return sum;
}