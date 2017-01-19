#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mosek.h"

#define EMAX 10000000 //max number of edges

typedef struct {
	unsigned s;//source node
	unsigned t;//target node
} Edge;

typedef struct {
	double w;//sum of weights
	unsigned s;//number of edges
	Edge *edge;//list of edges
} Coredge;

typedef struct {
	unsigned n;//number of nodes
  unsigned ce;//number of correlated edges
	unsigned e;//number of edges
  unsigned emax;//maximum number of edges
	Coredge *coredge;//list of corelated edges
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
	unsigned e=EMAX,i;
	graph *g=malloc(sizeof(graph));
  g->n=0;
  g->e=0;
	g->ce=0;
  g->emax=EMAX;
	g->coredge=malloc(EMAX*sizeof(Coredge));
	fscanf(file,"%u %u\n", &(g->s), &(g->t));
	while (fscanf(file,"%u %le ",&(g->coredge[g->e].s),&(g->coredge[g->e].w))==2){
		g->coredge[g->e].edge=malloc(g->coredge[g->e].s*sizeof(Edge));
		g->ce+=g->coredge[g->e].s;
		for (i=0;i<g->coredge[g->e].s;i++){
			fscanf(file,"%u %u ", &(g->coredge[g->e].edge[i].s), &(g->coredge[g->e].edge[i].t));
			g->n=max3(g->n,g->coredge[g->e].edge[i].s,g->coredge[g->e].edge[i].t);
		}
		if (g->e++==e) {
			e+=EMAX;
			g->coredge=realloc(g->coredge,e*sizeof(Coredge));
		}
	}
	fclose(file);
	g->n++;
	g->coredge=realloc(g->coredge,g->e*sizeof(Coredge));
	return g;
}

/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle,
                            const char str[])
{
  printf("%s",str);
} /* printstr */

int main(int argc, const char *argv[]) {
  unsigned i,j,incr;
	double res=0;

	graph *g=readedgelist(argv[1]);
	FILE* file;

  const int numvar = g->e+2*g->n;
  const int numcon = g->ce+1;

  double *c = malloc(numvar*sizeof(double));
  for (i=0;i<g->e;i++){
    c[i]=g->coredge[i].w;
  }
	for (i=g->e;i<numvar;i++){
		c[i]=0.;
	}

  MSKlidxt *aptrb = malloc(numcon*sizeof(MSKlidxt));
  for (i=0;i<g->ce;i++){
    aptrb[i]=3*i;
  }
	aptrb[g->ce]=3*g->ce;
  MSKlidxt *aptre = malloc(numcon*sizeof(MSKlidxt));
  for (i=0;i<g->ce;i++){
    aptre[i]=3*(i+1);
  }
	aptre[g->ce]=3*g->ce+2;
  MSKidxt *asub=malloc((3*g->ce+2)*sizeof(MSKlidxt));

  double *aval=malloc((3*g->ce+2)*sizeof(double));
	incr=0;
  for (i=0;i<g->e;i++){
		for (j=0;j<g->coredge[i].s;j++){
			asub[incr]=i;
			asub[incr+1]=g->e+g->coredge[i].edge[j].s;
			asub[incr+2]=g->e+g->coredge[i].edge[j].t;
			aval[incr]=1.;
			aval[incr+1]=-1.;
			aval[incr+2]=1.;
			incr+=3;
		}
  }
	asub[incr]=g->e+g->s;
	asub[incr+1]=g->e+g->t;
	aval[incr]=1.;
	aval[incr+1]=-1.;

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
			printf("outb\n");

      /* Input row i of A */
      if(r == MSK_RES_OK)
        r = MSK_putarow(task,
                        k,                 /* Row index.*/
                        aptre[k]-aptrb[k], /* Number of non-zeros in row i.*/
                        asub+aptrb[k],     /* Pointer to column indexes of row i.*/
                        aval+aptrb[k]);    /* Pointer to values of row i.*/
			printf("outa\n");

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
              printf("c[%d,%d]*d[%d,%d]= %e*%e=%e\n",g->coredge[l].edge[0].s,g->coredge[l].edge[0].t,g->coredge[l].edge[0].s,g->coredge[l].edge[0].t,c[l],xx[l],xx[l]*c[l]);
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
