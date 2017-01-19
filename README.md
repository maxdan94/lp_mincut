#README#

##Info##
This is a C code to compute the mincut using LP and MOSEK https://mosek.com/  
The LP is detailed here:
https://en.wikipedia.org/wiki/Max-flow_min-cut_theorem#Linear_program_formulation

##To compile##

Install MOSEK with the C api (license is free for academics)  
Change the Makefile replacing with your paths  
type make 

##To execute##

type "./mincut_lp graph.txt"

graph.txt should contain:
- the id of the source node followed by the ID of the target node on the first line,
- source target weight on each next lines.
graph.txt is an example of input.

It will print the value (d_ij between 0 and 1) associated to each edge and the value of the cut.

##Initial contributors##

Maximilien Danisch  
Technical consultant: Ziad Ismail  
January 2017  
http://bit.ly/maxdan94  
maximilien.danisch@telecom-paristech.fr
