#include<iostream>
#include<string>
#include<vector>
#include<queue>
#include<utility>
#include<cmath>
using namespace std;

#include "ugraph.h"


//Welcome to my UGraph! Enjoy :)
int main(){

	int N;
	cin >> N;
	Ugraph ag(N);

	string command;
	while(cin >> command){
		if(command == "addEdge"){
			int u, v, w;
			cin >> u >> v >> w;
			ag.addEdge(u, v, w);
		}else if(command == "removeEdge"){
			int u, v;
			cin >> u >> v;
			ag.removeEdge(u, v);
		}else if(command == "bfs"){
			int s;
			cin >> s;
			cout << "BFS on " << s << endl;
			ag.bfs(s);
		}else if(command == "dfs"){
			cout << "DFS" << endl;
			ag.dfs();
		}else if(command == "kruskal"){
			cout << "Running the Kruskal's algorithm." << endl;
			ag.kruskal();
			cout << "MST is:" << endl;
			ag.printMST();
			cout << "Weight of MST is: " ;
			int total = ag.weightMST();
			cout << total << endl;
		}else if(command == "dijkstra"){
			int s;
			cin >> s;
			cout << "Running the Dijkstra's algorithm on "<< s << endl;
			ag.dijkstra(s);
		}else if(command == "printPath"){
			int u, v;
			cin >> u >> v;
			cout << "Path from "<< u << " to " << v << " is:" << endl;
			ag.printPath(u, v);
			cout << endl;
		}else if(command == "lenghtShortestW"){
			int u,v;
			cin >> u >> v;
			cout << "Calculating length of the shortest path from " << u << " to " << v << endl;
			int res = ag.lenghtShortestW(u, v);
			cout << "The length is " << res << endl;

		
		}else if(command == "printDistance"){
			cout << "Distances:" << endl;
			ag.printDistance();
		}else if(command == "printParents"){
			cout << "Parents:" << endl;
			ag.printParents();
		}else if(command == "printGraph"){
			cout << "Graph is:" << endl;
			ag.printGraph();
		}else if(command == "printMST"){
			cout << "MST is:" << endl;
			ag.printMST();
		}else if(command == "weightMST"){
			int res = ag.weightMST();
			cout << "Total weight of MST is: " << res << endl;
			
		}else{
			cerr << "ERROR: the command " << command << " is not found." << endl;
		}
	}//while
  return 0;
}
