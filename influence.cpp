// n= number of nodes
// m= number of edges
// k= number of seed nodes
// t= number of monter carlo simulations
// c= number of communities

// this is entry point for the application

#include "BCRIM.h"
#include "ICRIM.h"

void bcrim(string path, int nodenum, int k, int t, int c, string model){
    // intialize new member of class bcrim
    Bcrim *app = new Bcrim();

    // call influence maximization function of the new instance
    app->influence_maximization(path, nodenum, k, t, c, model);
    
    //dellocate the pointer clearing the memory
    delete app;
	app = NULL;
}

void icrim(string path, int nodenum, int k, int t, int c, string model){
    // intialize new member of class bcrim
    Icrim *app = new Icrim();

    // call influence maximization function of the new instance
    app->influence_maximization(path, nodenum, k, t, c, model);
    
    //dellocate the pointer clearing the memory
    delete app;
	app = NULL;
}

// driver code
int main(){
    string path[10];// path to data nodes
    int nodenum[10];
    int c[10];

    string model = "IC";// define model(LT OR IC) 
    int k = 30;
    
    path[0] = "data/";
    nodenum[0] = 75888;
    c[0] = 20;

    // call icrim and bcrim for the given model
    icrim(path[0], nodenum[0], k, 100, c[0], model);
	bcrim(path[0], nodenum[0], k, 100, c[0], model);
    return 0;
}