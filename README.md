#README#

##Info##
This is a C code to compute the mincut using LP and MOSEK https://mosek.com/ 
The LP is detailed here:
https://en.wikipedia.org/wiki/Max-flow_min-cut_theorem#Linear_program_formulation

It is also an attempt to compute the correlated mincut where some edges can be correlated, that is "if edges ij and kl are correlated then ij is cut iff kl is cut". 
This is somehow a subproblem of what is define here: http://people.csail.mit.edu/stefje/papers/subcutsShort.pdf

##To compile##

Install MOSEK with the C api (license is free for academics) 
Change the Makefile replacing with your paths 
type "make"

##To execute##

type "./mclp graph.txt" for the basic mincut 
graph.txt should contain:
- the id of the source node followed by the ID of the target node on the first line,
- "source target weight" on each next lines.
The given graph.txt is an example of input.

type "./mclp_cor graph_cor.txt" for the correlated mincut 
graph_cor.txt should contain:
- the id of the source node followed by the ID of the target node on the first line,
- "n w s_1 t_1 s_2 t_2 ... s_n t_n" on each next lines (n is the number of correlated edges, w is the sum of the weights of the correlated edges, s_i t_i are the source and target for each edges".
The given graph_cor.txt is an example of input.

It will print the value (d_ij between 0 and 1) associated to each set of correlated edge and the value of the cut.

##Initial contributors##

Maximilien Danisch  
Technical consultant: Ziad Ismail  
January 2017  
http://bit.ly/maxdan94  
maximilien.danisch@telecom-paristech.fr
