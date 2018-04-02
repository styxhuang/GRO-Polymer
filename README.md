# GRO-Polymer
Generate polymer crosslinking network under gromacs enviroment

This process will be a cycle only in gromacs. All process will be separated into following parts:

-1 Read .gro and .itp file, and record all information according to the molecule type
-2 Output the gro and top file, that gromacs can use them to run the simulation
-3 Generate the bonds, then update the top file simulatenously


