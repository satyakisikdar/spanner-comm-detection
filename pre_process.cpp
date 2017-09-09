// preprocesses the edge list to ensure that it's rehashed and weighted using Jaccard coefficient
// the rehashed weighted edgelist and the inverse maps are stored in the same directory

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>

using namespace std;
using namespace boost;


typedef adjacency_list< listS,vecS, undirectedS > Graph;
typedef graph_traits < Graph> :: vertex_descriptor vertex_descriptor;
typedef graph_traits < Graph> :: edge_descriptor edge_descriptor;
typedef Graph :: vertex_iterator vertex_iterator;
typedef Graph :: edge_iterator edge_iterator;
typedef std::pair <int ,int> Edge;
typedef graph_traits < Graph> :: adjacency_iterator my_adjacency_iterator;


void print_set(std::unordered_set<int> uset)
{
	cout << endl;
	for (auto item: uset)
		cout << item << " ";
	cout << endl;
}

std::unordered_set<int> union_(std::unordered_set<int> uset1, std::unordered_set<int> uset2)
{
	std::unordered_set<int> new_uset;
	if (uset1.size() > uset2.size())
	{
		new_uset = uset1;
		for (auto thing: uset2)
			new_uset.insert(thing);
	}
	
	else
	{
		new_uset = uset2;
		for (auto thing: uset1)
			new_uset.insert(thing);
	}

	return new_uset;
}

std::unordered_set<int> intersection(std::unordered_set<int> uset1, std::unordered_set<int> uset2)
{
	std::unordered_set<int> new_uset;

	if (uset1.size() > uset2.size())
	{
		for (auto node: uset2)
		{
			if (uset1.find(node) != uset1.end())
				new_uset.insert(node);
		}
	}
	
	else
	{
		for (auto node: uset1)
		{
			if (uset2.find(node) != uset2.end())
				new_uset.insert(node);
		}
	}
		
	return new_uset;
}


void preprocesses(string input_filename)
{
	// rehashes the edge list at input_filename 
	Graph G;
	my_adjacency_iterator start, end;
	edge_iterator e_start, e_end;
	edge_descriptor e;
	std::unordered_map<int,int> ordered_labels;
	int u, v;
	int count = 0;
	float w;

	ifstream fin(input_filename);
	

	while (fin >> u >> v)
	{
		add_edge(u, v, G);

		if (ordered_labels.find(u) == ordered_labels.end())
		{
			ordered_labels[u] = count;
			count += 1;
		}

		if (ordered_labels.find(v) == ordered_labels.end())
		{
			ordered_labels[v] = count;
			count += 1;
		}
	}
	fin.close();

	string out_filename1 = "rehashed_weighted_" + input_filename;
	string out_filename2 = "rehashed_" + input_filename;

	ofstream fout1;
	fout1.open(out_filename1);
	ofstream fout2;
	fout2.open(out_filename2);

	for (tie(e_start, e_end) = edges(G); e_start != e_end; ++ e_start)
	{
		std::unordered_set<int> u_neighbors, v_neighbors;
		e = *e_start;
		u = source(e, G);
		v = target(e, G);
		
		for (tie(start, end) = adjacent_vertices(u, G); start != end; ++ start)
		{
			int neighbor = *start;
			u_neighbors.insert(neighbor); 
		}
		
		for (tie(start, end) = adjacent_vertices(v, G); start != end; ++ start)
		{
			int neighbor = *start;
			v_neighbors.insert(neighbor); 
		}


		float num = intersection(u_neighbors, v_neighbors).size();
		float denom = union_(u_neighbors, v_neighbors).size() - 2;

		if (denom == 0)
			w = 0;
		else
			w = num / denom;

		fout1 << ordered_labels[u] << " " << ordered_labels[v] << " " << w << endl;
		fout2 << ordered_labels[u] << " " << ordered_labels[v] << endl;
	}
	fout1.close();
	fout2.close();

	cout << "\nWeighted edgelist is at " << out_filename1 << endl;
	cout << "\nRehashed unweighted edgelist is at " << out_filename2 << endl;
	

	out_filename2 = input_filename + "_inv_maps";
	fout2.open(out_filename2);

	for (auto item: ordered_labels)
		fout2 << item.first << " " << item.second << endl;

	fout2.close();
	cout << "\nThe inverse mapping is at " << out_filename2 << endl;
}

int main(int argc, char const *argv[])
{
	if (argc < 2)
	{
		cout << "\nEnter filename of the unweighted edge list to preprocesses\n";
		return 0;
	}

	preprocesses(argv[1]);
	return 0;
}