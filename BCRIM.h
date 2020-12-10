// bcrim header file 
#pragma once
#ifndef BCRIM_H_
#define BCRIM_H_

#include<bits/stdc++.h>
using namespace std;

struct Neighbour {
	int node;
	double p;
};

class Bcrim {
private:
	int num_comm; //community
	int num_node; //nodes
	int num_seed; //seeds
	int round; //monte carlo simulation
	int *hub; // stores if a node is hub node
	int *node; //store node's name
	int *is_seed; // if a node is included in the seed set
	int *visited; // visited nodes
	int *comm; //node i in community; comm[i] c;
	int *comm_update; //if community is to be updated
	double *threshold; //activation threshold of i
	double *weight; //weighted fraction of neighbours of i = wieght[i].
	double *weight2;  //current  weighted fraction of i's neighbour
	double *IS; // the current influence of each comm
	double **MIS; // the marginal influence spread of each community of a new seed node
	vector< int > seed_set; // seed set
	vector< int > *comm_seed; // the chosen seed nodes in a comm
	vector< int > *comm_set; // the nodes in a comm
	vector< int > *H; // the extended nodes in a comm
	vector< Neighbour > *nbr; // network cut by commuunities
	vector< Neighbour > *G; // the whole network (hub)
	vector< Neighbour > *GG; // the whole network
public:
	Bcrim();
	void initialize(int t, int k, int n, int c, string model);
	void load(string path, string model);
	double IC(int t, int cm);
	double LT(int t, int cm);
	double CIC(int u);
	double CLT(int u);
	void Extend_seedsIC(int u);
	void Extend_seedsLT(int u);
	void output_to_file(string path, string model, int k, double timecost, double influencespread);
	void output_is(string path, string model, int k, double influencespread, double hubratio);
	void influence_maximization(string path, int n, int k, int t, int c, string model);
	double RandCasIC(int k, int t); // influence spread on network
	double RandCasLT(int k, int t); // influence spread on network
	void clr();
};

#endif // !BCRIM_H_
