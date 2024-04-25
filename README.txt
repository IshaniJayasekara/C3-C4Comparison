2024.01.06

Network alignment algorithm - AlignMCL

1. Merge step --> Alignment graph
2. Mining --> MCL algorithm

1. Merge

    a. Union graph 
            inputs - rice network, maize network, orthologous proteins
            outputs - union graph

    b. Raw alignment graph
            inputs - union graph
            output - raw alignment graph
    
    c. Pruning
            input - raw alignment graph
            output - alignment graph

2. Mining 
    
    MCL algorithm
        input - alignment graph, eli_score, inflation, expansion
        output - clusters


1. ortho