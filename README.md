# Graphs
Data Structures

This program has the implementation of the Graph class with these main member functions:

> extract_path:
    When an algorithm like BFS or DFS (or the critical paths function 
    below) is run, the actual paths explored/discovered are encoded 
    using "predecessor" information.  Given the predecessor info, 
    this will reconstruct paths.
   
> dag_critical_paths:
    This function takes a DAG and labels each vertex with the length of the
    LONGEST input path ("critical-paths" measured by sum of edge
    weights) ending at that vertex.  It also encodes the paths 
    themselves using the predecessor values.
    
>  dag_num_paths:
     This function labels each vertex with the number of input-to-output paths 
     which include that vertex.  This function does not encode
     any particular paths.  It simply records the number of such paths.
     
>  valid_topo_order:
     This function takes a vertex ordering and determines if it is 
     indeed a valid topological ordering of the given graph.
     (If graph is not even a DAG, the function will return false).
      
>   enum_paths:
      This function takes a vertex and constructs ALL input-paths ending at that 
      vertex.  The paths are represented as strings; a vector of
      strings is populated with the paths.
