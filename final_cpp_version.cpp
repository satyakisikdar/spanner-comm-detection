#include <iostream>
#include <cstdio>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
// #include <boost/unordered_set.hpp>
#include <ctime>
#include <stack>
#include <queue>
#include <fstream>
// #include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/functional/hash.hpp>
#include <random>
#include <unordered_set>
#include <unordered_map>
// using namespace std;
using namespace boost;

typedef std::vector < std::vector <int> > vvi;

typedef adjacency_list< listS, vecS, undirectedS, no_property, property<edge_weight_t, float> > Graph;
// typedef adjacency_list< listS,vecS, undirectedS > Graph;
typedef graph_traits < Graph> ::vertex_descriptor vertex_descriptor;
typedef graph_traits < Graph> ::edge_descriptor edge_descriptor;
typedef Graph::vertex_iterator vertex_iterator;
typedef Graph::edge_iterator edge_iterator;
typedef std::pair <int, int> Edge;
typedef graph_traits < Graph> ::adjacency_iterator my_adjacency_iterator;
property_map<Graph, edge_weight_t>::type weight;


std::unordered_set<int> union_(std::unordered_set<int> &set1, std::unordered_set<int> &set2)
{
	// returns the union of two unordered sets of integers
	std::unordered_set<int> new_set = set1;

	for(auto thing : set2)
		new_set.insert(thing);
	return new_set;
}

std::vector< std::vector <int> > sample(std::vector <int> &things, int &k)
{
	// Returns the set of sampled nodes needed for the whole process of spanner construction as vector of integer vectors 
	std::vector< std::vector<int> > master_v;
	std::vector<int> v = things;

	int i = 1, n = v.size(), s;
	float f;

	std::vector<int> counts;

	for(i = 1; i < k; i ++)
	{
		f = i / (float) k;
		s = ceil(pow(n, 1 - f));
		counts.push_back(s);
	}

	master_v.push_back(v);

	std::default_random_engine engine(std::random_device {}());
	std::shuffle(master_v[0].begin(), master_v[0].end(), engine);

	for(i = 1; i < k; i ++)
	{
		std::vector<int> temp;
		for(int j = 0; j < counts[i - 1]; j ++)
			temp.push_back(master_v[i - 1][j]);
		std::shuffle(temp.begin(), temp.end(), engine);
		master_v.push_back(temp);
	}
	return master_v;
}


tuple<Graph, Graph, Graph> read_graphs(std::string filename)
{
	// Reads the graphs from the edgelist. At first, the graphs G, G_comm and G_prime are all the same
	std::ifstream fp(filename);
	Graph G, G_comm, G_prime;
	int u, v;
	float w;

	while(fp >> u >> v >> w)
	{
		add_edge(u, v, w, G);
		add_edge(u, v, w, G_comm);
		add_edge(u, v, w, G_prime);
	}
	fp.close();
	return make_tuple(G, G_comm, G_prime);
}

Graph make_spanner(std::string filename, int k)
{
	// constructs the (2k - 1)-spanner using Baswana and Sen's randomized algorithm 
	Graph G, G_comm, G_prime;
	Graph G_spanner;

	tie(G, G_comm, G_prime) = read_graphs(filename);

	std::unordered_set <std::pair<int, int>, hash< std::pair<int, int> > > old_Ei;
	std::unordered_set <std::pair<int, int>, hash< std::pair<int, int> > > new_Ei;

	std::unordered_set <int> v_prime_flag;
	std::unordered_map <int, int> membership_oldc;
	std::unordered_map <int, int> membership_newc;
	std::unordered_map <int, int> neighborhood;

	std::unordered_map <int, std::unordered_set <int> > old_c;

	std::pair <vertex_iterator, vertex_iterator> vp;

	std::vector <int> all_nodes;

	for(vp = vertices(G); vp.first != vp.second; ++vp.first) // iterating over all the nodes
	{
		auto node = *vp.first;
		old_c[node].insert(node);
		v_prime_flag.insert(node);
		membership_oldc[node] = node;
		membership_newc[node] = -1;
		all_nodes.push_back(node);
	}

	my_adjacency_iterator start, end;

	std::vector < std::vector<int> > R_vector; //this stores all the randomly sampled cluster heads for all iterations

	R_vector = sample(all_nodes, k); // R_vector stores the set of sampled nodes required at each iteration

	for(int i = 1; i < k; i ++)
	{
		/*****1. Forming a sample of clusters **********************/

		std::unordered_set<int> R_i(R_vector[i].begin(), R_vector[i].end()); // R_vector[i] gives the randomly picked cluster heads for ith iteration
		std::unordered_map <int, std::unordered_set <int> > new_c;

		for(auto item : old_c)
		{
			int v = item.first;

			if(R_i.find(v) != R_i.end())
				new_c[v] = old_c[v];
		}

		std::unordered_set <int> sampled_nodes;

		for(auto v : R_i)
			sampled_nodes = union_(sampled_nodes, old_c[v]);

		std::unordered_set <int> unsampled_nodes;


		for(vp = vertices(G_prime); vp.first != vp.second; ++vp.first)
		{
			int node = *(vp.first);
			if(sampled_nodes.find(node) == sampled_nodes.end())
				unsampled_nodes.insert(node);
		}

		neighborhood.clear();
		for(auto node : unsampled_nodes)
			neighborhood[node] = -1;

		if(i == 1)
			for(auto v : R_i)
				membership_newc[v] = v;

		else
		{
			std::pair <vertex_iterator, vertex_iterator> vp;
			for(vp = vertices(G_prime); vp.first != vp.second; ++ vp.first)
			{
				auto v = *(vp.first);
				if(membership_newc[v] != -1)
					if(R_i.find(membership_newc[v]) == R_i.end())
						membership_newc[v] = -1;
			}
		}

		new_Ei = old_Ei;
		int u, v;

		for(auto edge : old_Ei)
		{
			tie(u, v) = edge;
			if(R_i.find(membership_newc[u]) == R_i.end() || R_i.find(membership_newc[v]) == R_i.end())
			{
				std::pair<int, int> e = std::make_pair(u, v);
				new_Ei.erase(e);
			}
		}

		/***********************end of step 1**************************/

		/********2. finding nearest neighboring sampled cluster*********/

		for(int v : unsampled_nodes)
		{
			if(v_prime_flag.find(v) == v_prime_flag.end()) // if v is not in G_prime, we don't consider it
				continue;

			std::unordered_map <int, std::pair<float, int> > min_table;
			edge_descriptor e;
			vertex_descriptor src, dest;

			for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
			{
				int neighbor = *start;
				float wt;
				src = vertex(v, G);
				dest = vertex(neighbor, G);

				std::tie(e, std::ignore) = edge(src, dest, G);

				wt = get(weight, e);

				if(membership_newc[neighbor] != -1) // checking if the neighbor is sampled
				{
					if(min_table.find(membership_newc[neighbor]) == min_table.end())
						min_table[membership_newc[neighbor]] = std::make_pair(wt, neighbor);


					else if(wt < min_table[membership_newc[neighbor]].first)
					{
						min_table[membership_newc[neighbor]].first = wt;
						min_table[membership_newc[neighbor]].second = neighbor;
					}
				}
			}

			float min_ = INT_MAX;
			int nearest_neighbor = -1;

			for(auto item : min_table)
			{
				if(item.second.first  < min_)
				{
					min_ = item.second.first;
					nearest_neighbor = item.second.second;
				}
			}

			neighborhood[v] = nearest_neighbor;
		}

		/************end of step 2 ********************************/


		/* **********3. adding edges to spanner ******************/

		for(int v : unsampled_nodes)
		{
			if(v_prime_flag.find(v) == v_prime_flag.end())
				continue;

			if(neighborhood[v] == -1) //v is not adjacent to any sampled nodes
			{
				std::unordered_map <int, std::pair<float, int> > min_table;

				edge_descriptor e;
				vertex_descriptor src, dest;
				std::vector <std::pair<int, int> > edges_to_be_removed;
				my_adjacency_iterator start, end;


				for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
				{
					int neighbor = *start;
					float wt;

					src = vertex(v, G);
					dest = vertex(neighbor, G);

					std::tie(e, std::ignore) = edge(src, dest, G);

					wt = get(weight, e);

					if(membership_oldc[neighbor] != -1)
					{
						if(min_table.find(neighbor) == min_table.end())
							min_table[neighbor] = std::make_pair(wt, neighbor);

						if(wt <= min_table[neighbor].first)
						{
							min_table[neighbor].first = wt;
							min_table[neighbor].second = neighbor;

							edges_to_be_removed.push_back(std::make_pair(v, neighbor));
						}
					}
				}

				for(auto edge : edges_to_be_removed)
					remove_edge(edge.first, edge.second, G_prime);

				for(auto item : min_table)
				{
					add_edge(v, item.second.second, item.second.first, G_spanner);
					remove_edge(v, item.second.second, G_comm);
				}
			}


			else
			{
				my_adjacency_iterator start, end;
				float wt;

				edge_descriptor e;
				std::tie(e, std::ignore) = edge(v, neighborhood[v], G);
				wt = get(weight, e);

				add_edge(v, neighborhood[v], wt, G_spanner);
				remove_edge(v, neighborhood[v], G_comm);

				new_Ei.insert(std::make_pair(v, neighborhood[v]));

				std::vector <int> neighbors;
				for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
				{
					int neighbor = *start;
					neighbors.push_back(neighbor);
				}

				for(auto neighbor : neighbors)
				{
					if(membership_newc[neighbor] == membership_newc[neighborhood[v]])
						remove_edge(v, neighbor, G_prime);
				}

				std::unordered_map <int, std::pair<float, int> > min_table;
				edge_descriptor e1;
				vertex_descriptor src, dest;

				for(auto neighbor : neighbors)
				{
					float wt1, wt2;
					src = vertex(v, G);
					dest = vertex(neighbor, G);

					std::tie(e1, std::ignore) = edge(src, dest, G);
					wt1 = get(weight, e1);

					dest = vertex(neighborhood[v], G);
					std::tie(e1, std::ignore) = edge(src, dest, G);
					wt2 = get(weight, e1);

					if(wt1 < wt2)
					{
						if(min_table.find(membership_oldc[neighbor]) == min_table.end())
							min_table[membership_oldc[neighbor]] = std::make_pair(wt1, neighbor);

						else if(wt1 < min_table[membership_oldc[neighbor]].first)
						{
							min_table[membership_oldc[neighbor]].first = wt2;
							min_table[membership_oldc[neighbor]].second = neighbor;
							remove_edge(v, neighbor, G_prime);
						}
					}
				}

				for(auto item : min_table)
				{
					add_edge(v, item.second.second, item.second.first, G_spanner);
					remove_edge(v, item.second.second, G_comm);
					remove_edge(v, item.second.second, G_prime);
				}


			}

		}


		/************end of step 3*****************************************/


		/******4. removing intra cluster edges ***********************/

		std::vector< std::pair <int, int> > intracluster_edges;

		for(auto item : new_Ei)
		{
			int u = item.first;
			int v = item.second;

			if(membership_newc[u] != -1)
			{
				membership_newc[v] = membership_newc[u];
				new_c[membership_newc[u]].insert(v);
			}

			else if(membership_newc[v] != -1)
			{
				membership_newc[u] = membership_newc[v];
				new_c[membership_newc[v]].insert(u);
			}
		}

		edge_iterator e_start, e_end;

		for(std::tie(e_start, e_end) = edges(G_prime); e_start != e_end; ++ e_start)
		{
			int u = source(*e_start, G_prime), v = target(*e_start, G_prime);

			if(membership_newc[u] == membership_newc[v])
				intracluster_edges.push_back(std::make_pair(u, v));
		}



		for(auto edge : intracluster_edges)
		{
			// // std:: cout << "\n(" << edge.first << ", " << edge.second << ")";
			remove_edge(edge.first, edge.second, G_prime);
		}

		/*************end of step 4*************************************************/

		/**** updations needed at the end of iteration******************************/

		v_prime_flag.clear();

		int count = 0;
		// // std:: cout << "E_prime: \n";
		for(std::tie(e_start, e_end) = edges(G_prime); e_start != e_end; ++ e_start)
		{
			int u = source(*e_start, G_prime), v = target(*e_start, G_prime);
			v_prime_flag.insert(u);
			v_prime_flag.insert(v);
			count ++;
		}

		for(auto item : new_Ei)
		{
			v_prime_flag.insert(item.first);
			v_prime_flag.insert(item.second);
		}

		old_c = new_c;
		membership_oldc = membership_newc;

		/***********The end of iteration i *****************************/
	}

	/******* Phase 2 **************************/

	vertex_iterator v_start, v_end;

	for(std::tie(v_start, v_end) = vertices(G_prime); v_start != v_end; ++ v_start)
	{
		int v = *v_start;
		if(v_prime_flag.find(v) == v_prime_flag.end())
			continue;

		std::unordered_map <int, std::pair<float, int> > min_table;
		edge_descriptor e;
		vertex_descriptor src, dest;

		std::vector <std::pair<int, int> > edges_to_be_removed;


		for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
		{
			int neighbor = *start;
			float wt;
			src = vertex(v, G);
			dest = vertex(neighbor, G);

			std::tie(e, std::ignore) = edge(src, dest, G);

			wt = get(weight, e);

			if(membership_newc[neighbor] == -1)
				continue;

			if(min_table.find(membership_newc[neighbor]) == min_table.end())
				min_table[membership_newc[neighbor]] = std::make_pair(wt, neighbor);

			if(wt <= min_table[membership_newc[neighbor]].first)
			{
				min_table[membership_newc[neighbor]].first = wt;
				min_table[membership_newc[neighbor]].second = neighbor;

				edges_to_be_removed.push_back(std::make_pair(v, neighbor));
			}

			if(min_table.find(membership_newc[neighbor]) == min_table.end())
				edges_to_be_removed.push_back(std::make_pair(v, neighbor));
		}

		for(auto edge : edges_to_be_removed)
			remove_edge(edge.first, edge.second, G_prime);

		for(auto item : min_table)
		{
			add_edge(v, item.second.second, item.second.first, G_spanner);
			remove_edge(v, item.second.second, G_comm);
		}
	}
	return G_comm;
}

float THRESHOLD = 0.66;


std::stack <int> Broker_Stack;
std::queue <int> Community_Queue;
std::unordered_set <int> isolated;
std::unordered_set <int> influenced;
std::unordered_map<int, int> labels;
std::unordered_set <int> traversed;
std::unordered_set <int> brokers;
std::unordered_set <int> marked;
std::unordered_map<int, float> scores;


float INS_score(Graph &g, int node)
{
	int count = 0;
	int total_count = 0;

	my_adjacency_iterator start, end;
	vertex_descriptor v = vertex(node, g);


	for(std::tie(start, end) = adjacent_vertices(v, g); start != end; ++ start)
	{
		auto neighbor = *start;
		if(influenced.find(neighbor) != influenced.end())
			count += 1;
		total_count += 1;
	}

	return count * 1.0 / total_count;
}


void INS(Graph &g, int starting_node)
{
	vertex_iterator v_st, v_end;
	influenced.insert(starting_node);
	traversed.insert(starting_node);
	Broker_Stack.push(starting_node);
	brokers.insert(starting_node);
	scores[starting_node] = 0;
	labels[starting_node] = starting_node;


	my_adjacency_iterator start, end;

	while(Broker_Stack.size() + Community_Queue.size() > 0)
	{
		int node;
		if(! Community_Queue.empty())
		{
			node = Community_Queue.front();
			Community_Queue.pop();
		}
		else
		{
			node = Broker_Stack.top();
			Broker_Stack.pop();
		}

		isolated.erase(node);

		vertex_descriptor v = vertex(node, g);

		for(std::tie(start, end) = adjacent_vertices(v, g); start != end; ++ start)
			influenced.insert(*start);

		for(std::tie(start, end) = adjacent_vertices(v, g); start != end; ++ start)
		{
			int neighbor = *start;

			if(traversed.find(neighbor) != traversed.end())
				continue;
			else
			{
				traversed.insert(neighbor);
			}

			auto score = INS_score(g, neighbor);
			// std:: cout << "node: " << neighbor << " score: " << score << std:: endl;
			scores[neighbor] = score;

			if(score < THRESHOLD) //broker 
			{
				labels[neighbor] = neighbor;
				Broker_Stack.push(neighbor);
				brokers.insert(neighbor);
			}

			else
			{
				if(marked.find(neighbor) == marked.end())
				{
					labels[neighbor] = labels[node];
					marked.insert(neighbor);
				}

				isolated.erase(neighbor);
				if(score != 1)
					Community_Queue.push(neighbor);
			}
		}
	}
}


void run_INS(Graph &G, std::string filename)
{
	for(unsigned int i = 0; i < num_vertices(G); i ++)
		isolated.insert(i);

	std::vector<int> component(num_vertices(G));
	int num = connected_components(G, &component[0]);

	vvi comp_lst(num);
	for(unsigned int i = 0; i != component.size(); ++i)
		comp_lst[component[i]].push_back(i);

	std::vector <int> starting_nodes;
	int min_deg = INT_MAX;
	int starting_node, deg;

	for(int i = 0; i < num; i ++)
	{
		starting_node = -1;
		min_deg = INT_MAX;
		for(auto node : comp_lst[i])
		{
			vertex_descriptor v = vertex(node, G);
			deg = degree(v, G);
			if(deg == 1)
			{
				starting_node = node;
				min_deg = deg;
				break;
			}
			else if(deg < min_deg)
			{
				starting_node = node;
				min_deg = deg;
			}
		}
		starting_nodes.push_back(starting_node);
	}

	for(auto start_node : starting_nodes)
		INS(G, start_node);

	std::string st = filename;
	std::ofstream fout(st + "_cover.part");

	for(auto it = labels.begin(); it != labels.end(); ++ it)
	{
		fout << it->first << " " << it->second << std::endl;
	}
	
	for(auto thing : isolated)
	{
		fout << thing << " " << thing << std::endl;
	}
	fout.close();

}




int main(int argc, char const *argv[])
{
	if(argc < 3)
	{
		std::cout << "Enter filename and k";
		return 0;
	}
	Graph G_comm;
	G_comm = make_spanner(argv[1], std::stoi(argv[2]));
	
	run_INS(G_comm, argv[1]);
	return 0;
}
