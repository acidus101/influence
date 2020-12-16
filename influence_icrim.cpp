// include header file for bcrim class
// TC = recordes timestamp when the last seed node that will affect this comm
// TN = records timestamp when node's MIS is computed
#include "ICRIM.h"

std::default_random_engine generatorI((unsigned int)(time(NULL)));
std::uniform_real_distribution<double> distI(0.0, 1.0);

// constructor for icrim
Icrim::Icrim(){
    // intialize private attributes
	hub = NULL;
	comm = NULL;
	node = NULL;
	visited = NULL;
	is_seed = NULL;
	threshold = NULL;
	comm_update = NULL;
	weight = NULL;
	weight2 = NULL;
	IS = NULL;
	MIS = NULL;
	TN = NULL;
	TC = NULL;
	comm_seed = NULL;
	comm_set = NULL;
	H = NULL;
	nbr = NULL;
	G = NULL;
	GG = NULL;
}

// allocating space to class attributes
void Icrim::initialize(int t, int k, int n, int c, string model){
	seed_set.clear();
	num_comm = c; // the number of communities
	round = t;
	num_seed = k;
	num_node = n;

	hub = new int[num_node + 2];
	node = new int[num_node + 2];
	is_seed = new int[num_node + 2];
	visited = new int[num_node + 2];
	comm = new int[num_node + 2];
	if (model.compare("LT") == 0) {
		threshold = new double[num_node + 2];
		weight = new double[num_node + 2];
		weight2 = new double[num_node + 2];
	}

	comm_update = new int[num_comm + 2];
	IS = new double[num_comm + 2];
	MIS = new double[num_comm + 2];
	TN = new int[num_node + 2];
	TC = new int[num_comm + 2];
	comm_seed = new vector< int >[num_comm + 2];
	comm_set = new vector< int >[num_comm + 2];
	H = new vector< int >[num_comm + 2];

	nbr = new vector< NeighbourIMP >[num_node + 2];
	G = new vector< NeighbourIMP >[num_node + 2];
	GG = new vector< NeighbourIMP >[num_node + 2];

	for (int i = 0; i < num_node + 2; i++) {
		comm[i] = -1;
		hub[i] = 0;
		node[i] = -1;
		is_seed[i] = 0;
		TN[i] = 0;
	}
	for (int i = 0; i < num_comm + 2; i++) {
		comm_seed[i].clear();
		comm_set[i].clear();
		IS[i] = 0;
		TC[i] = 0;
	}

}

// loading graph into communities
void Icrim::load(string path, string model){
    cout << "graph loading..." << endl;
    string node_path = path + "node_comm.txt";// define path for community node(node_comm.txt)
    string hub_path = path + "hub.txt";// define path for hub node(hub.txt)
    string edge_path = path + model + "-edges_pp.txt";// define path for edges("model"-edeges_pp.txt)

    // open file node_comm bfs through it and update the name of nodes
    //  in a commuity(node_comm structure(node-community))
    int u, v;
	char c;
	ifstream node_file;
	node_file.open(node_path.c_str());
	if (node_file) {
		while (!node_file.eof()) {
			node_file >> u >> v;
			if (u == -1) break;
			comm[u] = v;
			comm_set[v].push_back(u);
			node[u] = u;
			u = -1;
		}
	}
	else {
		cout << "Error opening file - " + node_path << endl;
		exit(1);
	}
	node_file.close();

    // open file hub.txt bfs through it and update whether a node 
    // is a hub node
	ifstream hub_file;
	hub_file.open(hub_path.c_str());
	if (hub_file) {
		while (!hub_file.eof()) {
			hub_file >> u >> v;
			if (u == -1) break;
			hub[u] = v;
			u = -1;
		}
	}
	else {
		cout << "Error opening file - " + hub_path << endl;
		exit(1);
	}
	hub_file.close();
    
    // open file lt-egdes_pp.txt bfs through it if community node add 
    // to the neighbour of current node(u) else just update the whole network
    double p;
	ifstream edge_file;
	edge_file.open(edge_path.c_str());
	if (edge_file) {
		while (!edge_file.eof()) {
			edge_file >> u >> v >> p;
			if (u == -1) break;
			NeighbourIMP nb;
			nb.node = v;
			nb.p = p;
			GG[u].push_back(nb);
			if (comm[u] == comm[v]) {
				nbr[u].push_back(nb);
				G[u].push_back(nb);
			}
			else if (hub[u]) {
				G[u].push_back(nb);
			}
			u = -1;
		}
	}
	else {
		cout << "Error opening file - " + edge_path << endl;
		exit(1);
	}
	edge_file.close();
}

void Icrim::Extend_seedsIC(int u){
    // initialize set h to null set and using u as starting node and push
    // set s into h and follow the algorithm further to extends seeds
    queue<int> queue;
    while(!queue.empty()){
        queue.pop();
    }

    visited[u] = 0;
	for (int i = 0; i < (int)G[u].size(); i++) {
		int node = G[u][i].node;
		visited[node] = 0;
	}
	for (int i = 0; i < num_comm; i++) {
		H[i].clear();
		comm_update[i] = 0; // initially, each community is not influenced
		for (int j = 0; j < (int)comm_seed[i].size(); j++) {
			int v = comm_seed[i][j];
			visited[v] = 1;
			H[i].push_back(v);
			for (int jj = 0; jj < (int)G[v].size(); jj++) {
				int node = G[v][jj].node; // neighbour
				if (!is_seed[node]) visited[node] = 0;
			}
			queue.push(v);
		}
	}

    // the one step diffusion algorithm
    // the remporary queue contains set s
    while(!queue.empty()){
        // pop each node
		int cur_node = queue.front();
		queue.pop();
        for(int j = 0; j < (int)G[cur_node].size(); j++){// for each neighbour node
			int node = G[cur_node][j].node; // neighbour
			double p = G[cur_node][j].p; 
			if ((comm[node] != comm[cur_node]) && (!visited[node])) {
				double puw = distI(generatorI);// puw a random number(define)
				if (puw <= p) { // activate successfully (extended)
					visited[node] = 1;
					H[comm[node]].push_back(node);
				}
			}            
        }

    }
    // rest of the extendseets part
	if (!visited[u]) {
		H[comm[u]].push_back(u);
		comm_update[comm[u]] = 1;
	} else{
		comm_update[comm[u]] = 0;
    }
    // if u is in the affected community push it into set h
    if(hub[u]){
        for(int i = 0; i < (int)G[u].size(); i++){// for each neighbour node
            int node = G[u][i].node;
            double p = G[u][i].p;
            if((comm[node] != comm[u]) && (!visited[node])){
                double puw = distI(generatorI);// puw a random number [0,1]
                if(puw <= p){
                    visited[node] = 1;
                    H[comm[node]].push_back(node);
                    comm_update[comm[node]] = 1;                
                }
            }
        }
    }
    
}

void Icrim::Extend_seedsLT(int u){
    for(int i = 0; i < num_node; i++){
        threshold[i] = distI(generatorI);// initialize threshold to a random number
        weight[i] = 0.0;// intialize the added weight to zero
        visited[i] = 0;
    }
    // initialize set h to null set and using u as starting node and push
    // set s into h and follow the algorithm further to extends seeds
    queue<int> queue;
    while(!queue.empty()){
        queue.pop();
    }

	visited[u] = 2;
	weight2[u] = 0.0;
	for (int i = 0; i < (int)G[u].size(); i++) {
		int node = G[u][i].node;
		visited[node] = 2; // mark u's neighbour as 2
		weight2[node] = 0.0;
	}
	for (int i = 0; i < num_node; i++) {
		H[i].clear();
		comm_update[i] = 0; // initially, each community is not influenced
		for (int j = 0; j < (int)comm_seed[i].size(); j++) {
			int v = comm_seed[i][j];
			visited[v] = 1;
			H[i].push_back(v);
			queue.push(v);
		}
	}

    // the one step diffusion algorithm
    // the remporary queue contains set s
    while(!queue.empty()){
        // pop each node
        int cur_node = queue.front();
		queue.pop();
        for(int j = 0; j < (int)G[cur_node].size(); j++){// for each neighbour node
			int node = G[cur_node][j].node; // neighbour
			double p = G[cur_node][j].p;
			if ((comm[node] != comm[cur_node]) && (visited[node] != 1)) {
                weight[node] += p;
                if(weight[node] >= threshold[node]){
                    visited[node] = 1;
                    H[comm[node]].push_back(node);
                }
			}
        }

    }

    // rest of the extendseets part
	if (weight[u] + weight2[u] >= threshold[u]) visited[u] = 1;
	if (visited[u] != 1) {
		visited[u] = 1;
		H[comm[u]].push_back(u);
		comm_update[comm[u]] = 1;
	}else{
		comm_update[comm[u]] = 0;
    }

    // if u is in the affected community push it into set h
    if(hub[u]){
        for(int i = 0; i < (int)G[u].size(); i++){// for each neighbour node
            int node = G[u][i].node;
            double p = G[u][i].p;
            if((comm[node] != comm[u]) && (visited[node] == 2)){
                weight[node] += p;
                if((weight[node] + weight2[node] - p < threshold[node]) && (weight[node] + weight2[node] >= threshold[node])){
                    visited[node] = 1;
                    H[comm[node]].push_back(node);
                    comm_update[comm[node]] = 1;
                }else {
                    visited[node] = 0;
                }
            }else if (visited[node] != 1){
                visited[node] = 0;
            }

        }
    }

}

// the diffusion models
double Icrim::IC(int t, int cm){
    // push all the nodes into a queue mark the node visited
	double influence_spread = 0.0;
	queue< int > queue;
	while (!queue.empty()) queue.pop();
	for (int i = 0; i < (int)comm_set[cm].size(); i++) {
		visited[comm_set[cm][i]] = 0;
	}
    // push all the nodes into a queue mark the node visited
	for (int i = 0; i < (int)H[cm].size(); i++) {
		queue.push(H[cm][i]);
		visited[H[cm][i]] = 1;
	}
	influence_spread += (int)H[cm].size();
	
    while(!queue.empty()){
        // pop the front node as current
        int cur_node = queue.front();
		queue.pop();
        
        for(int j = 0; j < nbr[cur_node].size(); j++){// for each neighbour v node of u
            // extract name and probability from the neighbour
            int node = nbr[cur_node][j].node; // neighbour
			double p = nbr[cur_node][j].p;
            if(!visited[node]){// does not belong to the current active node set
                // find a random number pw [0,1]
                double pw = distI(generatorI);
                if(pw <= p){
                    // activate the node pushing into the active nodes queue 
                    queue.push(node);
 					visited[node] = 1;
					influence_spread += 1;                   
                }
            }
        }
    }
    return influence_spread;
}

double Icrim::LT(int t, int cm){
    // push all the nodes into a queue mark the node visited
    double influence_spread = 0.0;
	queue< int > queue;
	while (!queue.empty()){
        queue.pop();
    }
    //push starting nodes into queue
    for (int i = 0; i < (int)H[cm].size(); i++) {
		queue.push(H[cm][i]);
	}
	influence_spread += (int)H[cm].size();

    while(!queue.empty()){ //bfs
        // pop the front node as current
        int cur_node = queue.front();
		queue.pop();

        for(int j = 0; j < nbr[cur_node].size(); j++){// for each neighbour v node of u
             // extract name and probaability from the neighbour
            int node = nbr[cur_node][j].node; // neighbour
			double p = nbr[cur_node][j].p;
            if(!visited[node]){// does not belong to the current active node set
                weight[node] += p; // sum of activated neighbours 
                if(weight[node] >= threshold[node]){
                    // activate the node pushing into the active nodes queue 
                    queue.push(node);
                    visited[node] = 1;
                    influence_spread += 1;
                }
            }
        }
    }
    return influence_spread;
}

// for calculating coomunity based marginal influence spread
// (the edge weight is probability of the node getting activated
// and gets activated is random p is <= the original puv)
double Icrim::CIC(int u){
    // initialize marginal influence spread
    double gain = 0.0;
    
    for (int i = 0; i < num_comm; i++){
		MIS[i] = 0.0;
    }

    for(int turn = 1; turn <= round; turn++){// do t monte carlo simulations
        // calls extend seeds algo
        Extend_seedsIC(u);

        for(int j = 0; j < num_comm; j++){// for each community
            if(comm_update[j]){ // if community is afftected
                //invokes ic model to calculate influence spread of each affected comm
                double influence_spread = IC(1, j);
                // finds change by subtracting the current influence spread of each
                double delta = influence_spread - IS[j];
                // community from the influence spread obtained
                // increase influence spread(gain) of affected community by change
                gain += delta;
                MIS[j] += delta;
            }
        }
    }
    for (int i = 0; i < num_comm; i++){
		MIS[i] /= round;
    }
    return gain/round;
}

// for calculating coomunity based marginal influence spread
// (the edge weight is the influence u exerts on v 
// and gets activated on passing the threshold)
double Icrim::CLT(int u){
    // initialize marginal influence spread
    double gain = 0.0;
    
    for (int i = 0; i < num_comm; i++){
		MIS[i] = 0.0;
    }
    for(int turn = 1; turn <= round; turn++){// do t monte carlo simulations
        // calls extend seeds algo
        Extend_seedsLT(u);
        
        for(int j = 0; j < num_comm; j++){// for each community
            if(comm_update[j]){ // if community is afftected
                //invokes ic model to calculate influence spread of each affected comm
                double influence_spread = LT(1, j);
                // finds change by subtracting the current influence spread of each
                double delta = influence_spread - IS[j];
                // community from the influence spread obtained
                // increase influence spread(gain) of affected community by change
				gain += delta;
				MIS[j] += delta;
            }
        }
    }
    for (int i = 0; i < num_comm; i++){
		MIS[i] /= round;
    }
    return gain/round;
}

void Icrim::output_to_file(string path, string model, int k, double timecost, double influencespread){
	path = path + model + "-ICRIM-results.txt";
	ofstream ofile;
	ofile.open(path, ios::out | ios::app);
	if (ofile.good()) cout << "success" << endl;
	else cout << "fail" << endl;

	ofile << setiosflags(ios::fixed);
	ofile << path << endl;
	ofile << "seed size: " << k << endl;
	ofile << "time taken: " << setprecision(6) << timecost << endl;
    ofile << "community influence: " << influencespread << endl;
	for (int i = 0; i < (int)seed_set.size(); i++)
		ofile << seed_set[i] << " ";
	ofile << endl;
	ofile << "------------------------------------------------------" << endl;
	ofile.close();
}

void Icrim::output_is(string path, string model, int k, double influencespread){
	path = path + model + "-ICRIM-results.txt";
	ofstream ofile;
	ofile.open(path, ios::out | ios::app);
	if (ofile.good()) cout << "success" << endl;
	else cout << "fail" << endl;

	ofile << setiosflags(ios::fixed);
	ofile << "k: " << k << setprecision(8) << "        influence spread: " << influencespread << endl;
	ofile << "*******************************************************" << endl;
	if (k == num_seed)
		ofile << endl;
	ofile.close();
}

int Icrim::check_update(int u, int tou) {
    int flag = 0;// f=1 indicates that updating is needed
    if(tou < TC[comm[u]]){
        flag = 1;
    }else if(hub[u]){
        for(int i = 0; i < (int)G[u].size(); i++){
            // for each neighbournode of i
			int v = G[u][i].node;
            if(tou < TC[comm[v]]){
                flag = 1;
                break;
            } 
        }
    }
    return flag;
}

void Icrim::influence_maximization(string path, int n, int k, int t, int c, string model){
    // check if correct diffusion model
    if ((model.compare("IC") != 0) and (model.compare("LT") != 0)) {
		cout << "The Diffusion Model Is Wrong. Please Check" << endl;
		return;
	}

    // call initialize to allocatespace
    initialize(t, k, n, c, model);

    // (load)divide graph into communities using necd
	load(path, model);
	
    clock_t stt, edd, edd2;
	stt = clock();
    
    // set totalinfluence = 0
	double totalinfluence = 0.0;
    while(!pq.empty()){// make priority queue null
        pq.pop();
    }
    
    for(int i = 0; i < num_node; i++){// for each node

        //calculating the community based marginal influence spread
		int vnode = node[i];
		if (vnode != -1) {
			pnode pd;
			pd.node = vnode;
			pd.ms = 0.0;
			if (model.compare("IC") == 0) {
					pd.ms = CIC(vnode);
			}
			else if (model.compare("LT") == 0) {
				pd.ms = CLT(vnode);
			}
			for (int j = 0; j < num_comm; j++){
				pd.CMIS[j] = MIS[j];
            }
            // push node, marginal spread into priority queue
			pq.push(pd);
		}

    }


	int i = 0;
	int flag = 1;
    while(i < num_seed){
        // remove front tuple from priority queue
		pnode pnext = pq.top();
		pq.pop();
		int u = pnext.node;
		int tou = TN[u];

        if(check_update(u, tou)){
			if (model.compare("IC") == 0) {
				pnext.ms = CIC(u);
			}
			else if (model.compare("LT") == 0) {
				pnext.ms = CLT(u);
			}
			for (int j = 0; j < num_comm; j++){
				pnext.CMIS[j] = MIS[j];
            }
			pq.push(pnext);
			TN[u] = i;
        }else{
            // a new seed is found
            seed_set.push_back(u);// add seed
            totalinfluence += pnext.ms;// update totalinfluence
			is_seed[u] = 1;
			comm_seed[comm[u]].push_back(u);
			i++;
            // update influence spread of each community
			for (int jj = 0; jj < num_comm; jj++) {
				IS[jj] += pnext.CMIS[jj];
			}
            TC[comm[u]] = i;// update tc(community u is affected)

            if(hub[u]){
                for(int j = 0; j < (int)G[u].size(); j++){// for each neighbour
                    int v = G[u][j].node;
                    TC[comm[v]] = i; // update tc for this community
                                    // hub seed node's neighbnour comms are influenced
                }
            }
        }

		// record the time cost of finding every five seed nodes
		double dt = 0;
		if ((flag == 1) && ((i == 1) || (i % 5 == 0))) {
			edd2 = clock();
			dt = (double)(edd2 - stt) / CLOCKS_PER_SEC;
			output_to_file(path, model, i, dt, totalinfluence);
			flag = 0;
		}
		if (i > 1 && i % 5 != 0) flag = 1;
		if (((dt > 50000) && (i < num_seed / 2)) || (dt >= 100000))
		{
			num_seed = i;
			break;
		}
    }
	edd = clock();
	double dt = (double)(edd - stt) / CLOCKS_PER_SEC;

	output_to_file(path, model, num_seed, dt, totalinfluence);


	for (int i = 1; i <= num_seed; i++) {
		if ((i == 1) || (i % 5 == 0) || (i == num_seed)) {
			double ifs = 0.0;
			if (model.compare("IC") == 0){
                ifs = RandCasIC(i, 100);
            }else if (model.compare("LT") == 0){
                ifs = RandCasLT(i, 100);
            }
			output_is(path, model, i, ifs);
			cout << i << " IS: " << ifs << endl;
		}
	}

	clr();
}


double Icrim::RandCasIC(int k, int t)
{
	double totalinfluencespread = 0.0;

	for (int r = 1; r <= t; r++) {
		queue< int > queue;
		while (!queue.empty()) queue.pop();
		for (int i = 0; i < num_node; i++) {
			visited[i] = 0;
		}

		for (int i = 0; i < k; i++) {
			visited[seed_set[i]] = 1;
			queue.push(seed_set[i]);
		}
		totalinfluencespread += k;

		while (!queue.empty()) { //BFS
			int cur_node = queue.front();
			queue.pop();

			for (int j = 0; j < G[cur_node].size(); j++) {
				int node = G[cur_node][j].node; // neighbour
				double p = G[cur_node][j].p;
				if (!visited[node]) {
					double pw = distI(generatorI);
					if (pw <= p) {//activate successfully
						queue.push(node);
						visited[node] = 1;
						totalinfluencespread += 1;
					}
				}
			}
		}
	}

	return totalinfluencespread / t;
}

double Icrim::RandCasLT(int k, int t)
{
	double totalinfluencespread = 0.0;

	for (int r = 1; r <= t; r++) {
		queue< int > queue;
		while (!queue.empty()) queue.pop();
		for (int i = 0; i < num_node; i++) {
			visited[i] = 0;
			weight[i] = 0.0;
			threshold[i] = distI(generatorI);
		}

		for (int i = 0; i < k; i++) {
			visited[seed_set[i]] = 1;
			queue.push(seed_set[i]);
		}
		totalinfluencespread += k;

		while (!queue.empty()) { //BFS
			int cur_node = queue.front();
			queue.pop();

			for (int j = 0; j < GG[cur_node].size(); j++) {
				int node = GG[cur_node][j].node; // neighbour
				double p = GG[cur_node][j].p;
				if (!visited[node]) {
					weight[node] += p;
					if (weight[node] >= threshold[node]) {//activate successfully
						queue.push(node);
						visited[node] = 1;
						totalinfluencespread += 1;
					}
				}
			}
		}
	}

	return totalinfluencespread / t;
}

void Icrim::clr(){
    // to deallocate memory used
	seed_set.clear();
	delete[] hub;
	delete[] node;
	delete[] is_seed;
	delete[] visited;
	delete[] comm;
	delete[] comm_update;
	delete[] threshold;
	delete[] weight;
	delete[] weight2;
	delete[] IS;
	delete[] MIS;
	delete[] TN;
	delete[] TC;
	for (int i = 0; i < num_node + 2; i++) {
		nbr[i].clear();
		G[i].clear();
	}
	for (int i = 0; i < num_comm + 2; i++) {
		comm_seed[i].clear();
		comm_set[i].clear();
		H[i].clear();
	}
	delete[] nbr;
	delete[] G;
	delete[] GG;
	delete[] comm_seed;
	delete[] comm_set;
	delete[] H;
}
