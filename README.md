## listed is decided
## structured is as the name suggests
## done is completed

# required files
    - driver code
    - bcrim file with ic and lt model
    - icrim file with ic and lt model
    - header bcrim defining the bcrim class
    - header icrim defininf the icrim class
    - driver code for greedy method
    - greedy_method file
    - greedy header file

# structure of data obtained
    - hub.txt file(hub nodes have their neighbour in different community)
        - contains names of nodes and marks them as hub nodes or not
    - node_comm.txt file()
        - contains name of the node and the community they belong to
    - M-edges_pp.txt file(M is model(IC OR LT))
        -contains the connected nodes along with their edge weight

# tasks
    - create main header file(done)
        - include header files(done)
        - bcrim function(done)
        - icrim function(done)
        - driver code(done)
    - bcrim header file(done)
        - include inbuilt header files(done)
        - add repeat check(done)
        - write neighbour nodes file(done)
        - write icrim class(done)
    - bcrim file(done)
        - bcrim class constructor(done)
        - initialize function - allocate memory to class(done)
        - load function- load graph(done)
        - ic model(done)
        - lt model (done)
        - cic algo(done)
        - clt algo(done)
        - extendseedsic algo(done)
        - extendseedslt algo (done)
        - onestepdiffusionic algo(done-listed in extendseeds algo)
        - onestepdiffusionlt algo(done- listed in extendseeds algo)
        - output_to_filfe function(done)
        - influence_maximization_bcrim(done)
        - random_number(done)
        - clear function to clear the allocated memory(done)
    - icrim header file(done)
        - include inbuilt header files(done)
        - write neighbour nodes file(done)
        - add repeat check(done)
        - write icrim class(done)
    - icrim file(done)
        - Icrim class constructor(done)
        - initialize function - allocate memory to class(done)
        - load function- load graph(done)
        - ic model(done)
        - lt model (done)
        - cic algo(done)
        - clt algo(done)
        - extendseedsic algo(done)
        - extendseedslt algo (done)
        - onestepdiffusionic algo(done-listed in extendseeds algo)
        - onestepdiffusionlt algo(done- listed in extendseeds algo)
        - output_to_filfe function(done)
        - check update function(done)
        - influence_maximization_icrim(done)
        - random_number(done)
        - clear function to clear the allocated memory(done)
    - greedy driver file(listed)
        - include header files(structured)
        - kempg function(sructured)
        - driver code(structured)
    - greedy method file(listed)
        - constructor for kempg class(structured)
        - initialize function to initialize the class attributes(listed)
        - ic model(structured)
        - lt model (structured)
        - influence function(structured)
        - clear function to deallocate(structured)
    - greedy header file(listed)
        -yet to be done