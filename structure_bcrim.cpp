// include header file for bcrim class

// constructor for bcrim
bcrim::bcrim(){
    // intialize private attributes
}

// allocating space to class attributes
void bcrim::initialize(){

}

// loading graph into communities
void bcrim::load(string path, string M){
    // define path for community node(node_comm.txt)
    // define path for hub node(hub.txt)
    // define path for edges("model"-edeges_pp.txt)

    // open file node_comm bfs through it and update the name of nodes in a commuity(node_comm structure(node-community))

    // open file hub.txt bfs through it and update whether a node is a hub node

    // open file lt-egdes_pp.txt bfs through it if community node add to the neighbour of current node(u) else just update the whole network

}

void bcrim::ExtendseedsIC(){
    // initialize set h to null set

    // the one step diffusion algorithm

    // rest of the extendseets part
}

void bcrim::ExtendseedsLT(){
    for(;;){
        // initialize threshold to a random number
        // intialize the added weight to zero
    }
    // initialize set h to null set

    // the one step diffusion algorithm

    // rest of the extendseets part
}

// the diffusion models
double bcrim::IC(){
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

double bcrim::LT(){
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
double bcrim::CIC(){
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
double bcrim::CLT(){
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

void bcrim::output_to_file(){
    // write the output to a new file
}

void bcrim::influence_maximization(string path, int n, int k, int t, int c, string M){
    // check if correct diffusion model
    // call initialize to allocatespace
    // (load)divide graph into communities using necd

    // set totalinfluence = 0
    
    for(;;){// for each seedsize(k)
        
        for(;;){// for each node belonging to the set of nodes
            // if the node is not icluded in seed set(name present in array)

            //calculating the community based marginal influence spread
            // check model(IC OR LT)

            // for ic call CIC

            // for lt call CLT

        }
    }
    
    // push new seet node into seedset array
    // increment totalinfluence
    // mark it as seed
}

void bcrim::clear(){
    // to deallocate memory used
}

