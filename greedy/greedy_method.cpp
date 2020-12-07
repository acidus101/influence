#include "KEMPG.h"
// constructor
kempg::kempg(){

}

// allocating space to class attributes
void kempg::initialize(){

}
// the diffusion models
double kempg::IC(){
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

double kempg::LT(){
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

kempg::influence(){

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
 
}

kempg::clr(){
    // to deallocate memory used
}