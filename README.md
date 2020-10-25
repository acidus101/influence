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
    - create main header file(structured)
        - include header files(structured)
        - bcrim function(sructured)
        - icrim function(sructured)
        - driver code(structured)
    - bcrim header file(listed)
        - include inbuilt header files(listed)
        - add repeat check(listed)
        - write neighbour nodes file(listed)
        - write icrim class(listed)
    - bcrim file(listed)
        - bcrim class constructor(structured)
        - initialize function - allocate memory to class(structured)
        - load function- load graph(structured)
        - ic model(structured)
        - lt model (structured)
        - cic algo(structured)
        - clt algo(structured)
        - extendseedsic algo(structured)
        - extendseedslt algo (structured)
        - onestepdiffusionic algo(structured-listed in extendseeds algo)
        - onestepdiffusionlt algo(structured- listed in extendseeds algo)
        - output_to_filfe function(listed)
        - influence_maximization_bcrim(structured)
        - random_number(structured)
        - clear function to clear the allocated memory(listed)
    - icrim header file(listed)
        - include inbuilt header files(listed)
        - write neighbour nodes file(listed)
        - add repeat check(listed)
        - write icrim class(listed)
    - icrim file(listed)
        - Icrim class constructor(structured)
        - initialize function - allocate memory to class(structured)
        - load function- load graph(structured)
        - ic model(structured)
        - lt model (structured)
        - cic algo(structured)
        - clt algo(structured)
        - extendseedsic algo(structured)
        - extendseedslt algo (structured)
        - onestepdiffusionic algo(structured-listed in extendseeds algo)
        - onestepdiffusionlt algo(structured- listed in extendseeds algo)
        - output_to_filfe function(listed)
        - check update function(structured)
        - influence_maximization_icrim(structured)
        - random_number(structured)
        - clear function to clear the allocated memory(listed)
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