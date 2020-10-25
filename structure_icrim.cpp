// include header file for bcrim class
// TC = recordes timestamp when the last seed node that will affect this comm
// TN = records timestamp when node's MIS is computed
#include "ICRIM.h"

std::default_random_engine generator((unsigned int)(time(NULL)));
std::uniform_real_distribution<double> dist(0.0, 1.0);

// constructor for icrim
icrim::icrim(){
    // intialize private attributes
}

// allocating space to class attributes
void icrim::initialize(){

}

// loading graph into communities
void icrim::load(string path, string M){
    // define path for community node(node_comm.txt)
    // define path for hub node(hub.txt)
    // define path for edges("model"-edeges_pp.txt)

    // open file node_comm bfs through it and update the name of nodes in a commuity(node_comm structure(node-community))

    // open file hub.txt bfs through it and update whether a node is a hub node

    // open file lt-egdes_pp.txt bfs through it if community node add to the neighbour of current node(u) else just update the whole network

}

void icrim::ExtendseedsIC(){
    // initialize set h to null set and using u as starting node and push
    // set s into h and follow the algorithm further to extends seeds
    queue<int> que;
    while(!que.empty()){
        que.pop();
    }

    // the one step diffusion algorithm
    // the remporary queue contains set s
    while(!que.empty()){
        // pop each node
        
        for(;;){// for each neighbour node
            double puw; // puw a random number(define)
            if(puw <= p){
                // include this in the set h
            }
        }

    }
    // rest of the extendseets part

    // if u is in the affected community push it into set h
    for(;;){// for each neighbour node
        if(/*does not belong to community*/){
            double puw;// puw a random number [0,1]
            if(puw <= p){
                // add this to set h
            }
        }
    }
}

void icrim::ExtendseedsLT(){
    for(;;){
        // initialize threshold to a random number
        // intialize the added weight to zero
    }
    // initialize set h to null set and using u as starting node and push
    // set s into h and follow the algorithm further to extends seeds
    queue<int> que;
    while(!que.empty()){
        que.pop();
    }

    // the one step diffusion algorithm
    // the remporary queue contains set s
    while(!que.empty()){
        // pop each node
        for(;;){// for each neighbour node
            weight[w] += p
            if(weight[w] >= threshold){
                // include this in the set h
            }
        }

    }

    // rest of the extendseets part

    // if u is in the affected community push it into set h
    for(;;){// for each neighbour node
        if(/*does not belong to community*/){
            wieght[w] += p;
            if(weight[w] >= threshold){
                // add this to set h
            }
        }
    }
}

// the diffusion models
double icrim::IC(){
    // push all the nodes into a queue mark the node visited

    while(!queue.empty()){
        // pop the front node as current

        for(;;){// for each neighbour v node of u
        // extract name and probaability from the neighbour
            if(!visited){// does not belong to the current active node set
                // find a random number pw [0,1]
                if(pw <= p){
                    // activate the node pushing into the active nodes queue 
                    // and mark visited
                    // increase the influence spread by 1
                }
            }
        }
    }
    return influence_spread;
}

double icrim::LT(){
    // push all the nodes into a queue mark the node visited

    while(!queue.empty()){
        // pop the front node as current

        for(;;){// for each neighbour v node of u
        // extract name and probaability from the neighbour
            if(!visited){// does not belong to the current active node set
                weight += p // sum of activated neighbours 
                if(weight >= threshold){
                    // activate the node pushing into the active nodes queue 
                    // and mark visited
                    // increase the influence spread by 1
                }
            }
        }
    }
    return influence_spread;
}

// for calculating coomunity based marginal influence spread
// (the edge weight is probability of the node getting activated
// and gets activated is random p is <= the original puv)
double icrim::CIC(){
    // initialize marginal influence spread
    double gain = 0.0;
    
    for(;;){// do t monte carlo simulations
        // calls extend seeds algo
        ExtendseedsIC()
        
        for(;;){// for each community
            if(){ // if community is afftected
                //invokes ic model to calculate influence spread of each affected comm
                
                // finds change by subtracting the current influence spread of each
                // community from the influence spread obtained
                
                // increase influence spread(gain) of affected community by change
            }
        }
    }
    return gain/t;
}

// for calculating coomunity based marginal influence spread
// (the edge weight is the influence u exerts on v 
// and gets activated on passing the threshold)
double icrim::CLT(){
    // initialize marginal influence spread
    double gain = 0.0;
    
    for(;;){// do t monte carlo simulations
        // calls extend seeds algo
        ExtendseedsLT()
        
        for(;;){// for each community
            if(){ // if community is afftected
                //invokes ic model to calculate influence spread of each affected comm
                
                // finds change by subtracting the current influence spread of each
                // community from the influence spread obtained
                
                // increase influence spread(gain) of affected community by change
            }
        }
    }
    return gain/t;
}

void icrim::output_to_file(){
    // write the output to a new file
}

void icrim::check_update() {
    int flag = 0;// f=1 indicates that updating is needed
    if(tou < toucom){
        f = 1;
    }else if(hub[u]){
        for(i;;){
            // for each neighbournode of i
            if(tau < taucom[v]){
                f = 1;
                break;
            } 
        }
    }
    return f;
}

void icrim::influence_maximization(string path, int n, int k, int t, int c, string M){
    // check if correct diffusion model
    // call initialize to allocatespace
    // (load)divide graph into communities using necd

    // set totalinfluence = 0

    while(!pq.empty()){// make priority queue null
        pq.pop();
    }
    
    for(;;){// for each node

        //calculating the community based marginal influence spread
        // check model(IC OR LT)

        // for ic call CIC

        // for lt call CLT

        // push node, marginal spread into priority queue

    }
    
    // for each community set tc to zero
    for(;;){

    }

    while(i < k){
        // remove front tuple from priority queue

        if(check_update(/*for tuple*/)){
            // check model(IC OR LT)

            // for ic call CIC

            // for lt call CLT
        }else{
            // a new seed is found
            // add seed
            // update totalinfluence
            // update influence spread of each community

            // update tc

            if(hub[u]){
                for(;;){// for each neighbour
                    // update tc for this community
                    // i++
                }
            }
        }
    }

}

void icrim::clr(){
    // to deallocate memory used
}

