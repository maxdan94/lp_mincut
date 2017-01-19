#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mosek.h"

#define EMAX 10000000 //max number of edges

typedef struct {
	unsigned s;//source node
	unsigned t;//target node
  double w;//weight
} edge;

typedef struct {
	unsigned n;//number of nodes
  unsigned e;//number of weighted edges
  unsigned emax;//maximum number of weighted edges
	edge *edges;//list of weighted edges
	unsigned s;//ID of the source node
	unsigned t;//ID of the target node
} graph;


//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
graph* readedgelist(const char* edgelist){
	FILE *file=fopen(edgelist,"r");
	unsigned e=EMAX;
	graph *g=malloc(sizeof(graph));
  g->n=0;
  g->e=0;
  g->emax=EMAX;
	g->edges=malloc(EMAX*sizeof(edge));
	fscanf(file,"%u %u\n", &(g->s), &(g->t));
	while (fscanf(file,"%u %u %le\n", &(g->edges[g->e].s), &(g->edges[g->e].t),&(g->edges[g->e].w))==3) {
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (g->e++==e) {
			e+=EMAX;
			g->edges=realloc(g->edges,e*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->edges=realloc(g->edges,g->e*sizeof(edge));
	return g;
}

/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle,
                            const char str[])
{
  printf("%s",str);
} /* printstr */

int main(int argc, const char *argv[]) {
  unsigned i;
	double res=0;
	graph *g=readedgelist(argv[1]);
	FILE* file;

  const int numvar = g->e+2*g->n;
  const int numcon = g->e+1;

  double *c = malloc(numvar*sizeof(double));
  for (i=0;i<g->e;i++){
    c[i]=g->edges[i].w;
  }
	for (i=g->e;i<numvar;i++){
		c[i]=0.;
	}

  MSKlidxt *aptrb = malloc(numcon*sizeof(MSKlidxt));
  for (i=0;i<g->e;i++){
    aptrb[i]=3*i;
  }
	aptrb[g->e]=3*g->e;
  MSKlidxt *aptre = malloc(numcon*sizeof(MSKlidxt));
  for (i=0;i<g->e;i++){
    aptre[i]=3*(i+1);
  }
	aptre[g->e]=3*g->e+2;
  MSKidxt *asub=malloc((3*g->e+2)*sizeof(MSKlidxt));
  double *aval=malloc((3*g->e+2)*sizeof(double));
  for (i=0;i<g->e;i++){
    asub[3*i]=i;
		asub[3*i+1]=g->e+g->edges[i].s;
		asub[3*i+2]=g->e+g->edges[i].t;
		aval[3*i]=1.;
		aval[3*i+1]=-1.;
		aval[3*i+2]=1.;
  }
	asub[3*i]=g->e+g->s;
	asub[3*i+1]=g->e+g->t;
	aval[3*i]=1.;
	aval[3*i+1]=-1.;

  MSKenv_t     env  = NULL;
  MSKtask_t    task = NULL;
  MSKrescodee  r;
  MSKidxt      k,l;

  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);

  if ( r==MSK_RES_OK )
  {
    /* Create the optimization task. */
    r = MSK_maketask(env,numcon,numvar,&task);

    /* Directs the log task stream to the 'printstr' function. */
    if ( r==MSK_RES_OK )
      r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

    /* Append 'numcon' empty constraints.
     The constraints will initially have no bounds. */
    if ( r == MSK_RES_OK )
      r = MSK_appendcons(task,numcon);

    /* Append 'numvar' variables.
     The variables will initially be fixed at zero (x=0). */
    if ( r == MSK_RES_OK )
      r = MSK_appendvars(task,numvar);

    for(l=0; l<numvar && r == MSK_RES_OK; ++l)
    {
      /* Set the linear term c_j in the objective.*/
      if(r == MSK_RES_OK)
        r = MSK_putcj(task,l,c[l]);

      /* Set the bounds on variable j.
       blx[j] <= x_j <= bux[j] */
      if(r == MSK_RES_OK)
        r = MSK_putvarbound(task,
                            l,           /* Index of variable.*/
                            MSK_BK_LO,      /* Bound key.*/
                            0.0,      /* Numerical value of lower bound.*/
                            +MSK_INFINITY);     /* Numerical value of upper bound.*/
    }

    /* Set the bounds on constraints.
       for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
    for(k=0; k<numcon-1 && r==MSK_RES_OK; ++k)
    {
      r = MSK_putconbound(task,k,MSK_BK_LO,0.,+MSK_INFINITY);     /* Numerical value of upper bound.*/

      /* Input row i of A */
      if(r == MSK_RES_OK)
        r = MSK_putarow(task,
                        k,                 /* Row index.*/
                        aptre[k]-aptrb[k], /* Number of non-zeros in row i.*/
                        asub+aptrb[k],     /* Pointer to column indexes of row i.*/
                        aval+aptrb[k]);    /* Pointer to values of row i.*/
    }
		if(r == MSK_RES_OK)
		r = MSK_putconbound(task,k,MSK_BK_LO,1.,+MSK_INFINITY);

		/* Input row i of A */
		if(r == MSK_RES_OK)
			r = MSK_putarow(task,
											k,                 /* Row index.*/
											aptre[k]-aptrb[k], /* Number of non-zeros in row i.*/
											asub+aptrb[k],     /* Pointer to column indexes of row i.*/
											aval+aptrb[k]);    /* Pointer to values of row i.*/


    /* Maximize objective function. */
    if (r == MSK_RES_OK)
      r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

    if ( r==MSK_RES_OK )
    {
      MSKrescodee trmcode;

      /* Run optimizer */
      r = MSK_optimizetrm(task,&trmcode);

      /* Print a summary containing information
         about the solution for debugging purposes. */
      MSK_solutionsummary (task,MSK_STREAM_LOG);

      if ( r==MSK_RES_OK )
      {
        MSKsolstae solsta;

        if (r == MSK_RES_OK)
          r = MSK_getsolsta (task,MSK_SOL_BAS,&solsta);
        switch(solsta)
        {
          case MSK_SOL_STA_OPTIMAL:
          case MSK_SOL_STA_NEAR_OPTIMAL:
          {
            double *xx = (double*) calloc(numvar,sizeof(double));
            MSK_getxx(task,MSK_SOL_BAS,xx);
            printf("Optimal primal solution\n");
            for(l=0; l<g->e; ++l){
							//if (xx[l]>1e-5){
              printf("c[%d,%d]*d[%d,%d]= %e*%e=%e\n",g->edges[l].s,g->edges[l].t,g->edges[l].s,g->edges[l].t,c[l],xx[l],xx[l]*c[l]);
							res+=xx[l]*c[l];
							//printf("objective = %le\n",res);
							//}
            }
						printf("Value of the min cut = %le\n",res);
            free(xx);
            break;
          }
          case MSK_SOL_STA_DUAL_INFEAS_CER:
          case MSK_SOL_STA_PRIM_INFEAS_CER:
          case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
          case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
            printf("Primal or dual infeasibility certificate found.\n");
            break;
          case MSK_SOL_STA_UNKNOWN:
          {
            char symname[MSK_MAX_STR_LEN];
            char desc[MSK_MAX_STR_LEN];

            /* If the solutions status is unknown, print the termination code
               indicating why the optimizer terminated prematurely. */

            MSK_getcodedesc(trmcode,
                            symname,
                            desc);

            printf("The solutuion status is unknown.\n");
            printf("The optimizer terminitated with code: %s\n",symname);
            break;
          }
          default:
            printf("Other solution status.\n");
            break;
        }
      }
    }

    if (r != MSK_RES_OK)
    {
      /* In case of an error print error code and description. */
      char symname[MSK_MAX_STR_LEN];
      char desc[MSK_MAX_STR_LEN];

      printf("An error occurred while optimizing.\n");
      MSK_getcodedesc (r,
                       symname,
                       desc);
      printf("Error %s - '%s'\n",symname,desc);
    }

    /* Delete the task and the associated data. */
    MSK_deletetask(&task);
  }

  /* Delete the environment and the associated data. */
  MSK_deleteenv(&env);

  return r;
}
