#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mosek.h"

#define NNODES 10000000 //max number of nodes: will increase if needed
#define NCOMS 10000000 //max number of communities: will increase if needed
#define SCOM 10000000 //max size of community: will increase if needed
#define NCPN 10 //max number of community per nodes: will increase if needed
#define EMAX 100000000 //max number of edges in bipartit: will increase if needed

typedef struct {
	unsigned s;//size of the community
	unsigned smax;//maximum size of the community
	unsigned *nodes;//nodes in the community
} community;

community *alloccom(){
	community *com=malloc(sizeof(community));
	com->smax=SCOM;
	com->nodes=malloc(SCOM*sizeof(unsigned));
	return com;
}

typedef struct {
	unsigned n;//number of nodes
	unsigned m;//max number of nodes

	unsigned **c;//c[i]=communities containing node i.
	unsigned *s;//s[i]=number of communities containing node i.
	unsigned *smax;//smax[i]=max number of communities

	unsigned nc;//number of communities
	unsigned mc;//max number of communities
	unsigned *size;//size[i]=size of the community i

	unsigned *tmp;//tmp[i] number of nodes shared by the current community and community i
	unsigned *list;//list of the community sharing at least one node with the current community
	unsigned nlist;//length of list
} compare;

compare *alloccompare(){
	compare *comp=malloc(sizeof(compare));
	comp->n=0;
	comp->m=NNODES;
	comp->c=calloc(comp->m,sizeof(unsigned*));
	comp->s=calloc(comp->m,sizeof(unsigned));
	comp->smax=calloc(comp->m,sizeof(unsigned));

	comp->nc=0;
	comp->mc=NCOMS;
	comp->size=malloc(comp->mc*sizeof(unsigned));

	comp->tmp=calloc(comp->mc,sizeof(unsigned));
	comp->list=calloc(comp->mc,sizeof(unsigned));
	comp->nlist=0;
	return comp;
}

typedef struct {
	unsigned u;//first node
	unsigned v;//second node
  double w;//weight
} edge;

typedef struct {
	unsigned n1;//number of nodes in set 1
	unsigned n2;//number of nodes in set 2
  unsigned e;//number of weighted edges
  unsigned emax;//maximum number of weighted edges
	edge *edges;//list of weighted edges

  unsigned *cd1;//cumulative degree of nodes in set 1 cd1[0]=0, size n1+1
  unsigned *cd2;//cumulative degree of nodes in set 2
  unsigned *adj1;//IDs of adjacent edges to nodes in set 1
  unsigned *adj2;//IDs of adjacent edges to nodes in set 2

} bipartite;

bipartite *allocbip(){
	bipartite *bip=malloc(sizeof(bipartite));
	bip->n1=0;
  bip->n2=0;
  bip->e=0;
  bip->emax=EMAX;
	bip->edges=malloc(EMAX*sizeof(edge));
	return bip;
}

bool readlinecom(FILE* file,community* com){
	char c;
	com->s=0;
	while(fscanf(file,"%u%c",com->nodes+com->s,&c)==2){
		if ( ++(com->s) == com->smax) {
			com->smax+=SCOM;
			com->nodes=realloc(com->nodes,com->smax*sizeof(unsigned));
		}
		if (c=='\n') {
			return 1;
		}
	}
	return 0;
}

void com2comp(community* com,compare* comp){
	unsigned i,j,tmp;

	if (comp->nc==comp->mc){
		comp->mc+=NCOMS;
		comp->size=realloc(comp->size,comp->mc*sizeof(unsigned));
		comp->tmp=realloc(comp->tmp,comp->mc*sizeof(unsigned));
		comp->list=realloc(comp->size,comp->mc*sizeof(unsigned));
		for (i=comp->mc-NCOMS;i<comp->mc;i++){
			comp->tmp[i]=0;
		}
	}
	comp->size[comp->nc]=com->s;
	for (i=0;i<com->s;i++){
		if (com->nodes[i]>comp->m){
			tmp=com->nodes[i]-comp->m+NCOMS;
			comp->m+=tmp;
			comp->c=realloc(comp->c,comp->m*sizeof(unsigned*));
			comp->s=realloc(comp->s,comp->m*sizeof(unsigned));
			comp->smax=realloc(comp->smax,comp->m*sizeof(unsigned));
			for (j=comp->m-tmp;j<comp->m;j++){
				comp->c[j]=NULL;
				comp->s[j]=0;
				comp->smax[j]=0;
			}
		}
		if (comp->s[com->nodes[i]]==comp->smax[com->nodes[i]]){
			if (comp->s[com->nodes[i]]==0){
				comp->c[com->nodes[i]]=malloc(NCPN*sizeof(unsigned));
				comp->smax[com->nodes[i]]=NCPN;
			}
			else {
				comp->smax[com->nodes[i]]+=NCPN;
				comp->c[com->nodes[i]]=realloc(comp->c[com->nodes[i]],comp->smax[com->nodes[i]]*sizeof(unsigned));
			}
		}
		comp->c[com->nodes[i]][comp->s[com->nodes[i]]++]=comp->nc;
	}
	comp->nc++;
}

void comsim(bipartite *bip, community* com,compare *comp, double min){
	double sim;
	unsigned i,j;
  edge ed;
  unsigned c1;
  static unsigned c2=0;

	for (i=0;i<com->s;i++){
		for (j=0;j<comp->s[com->nodes[i]];j++){
			if (comp->tmp[comp->c[com->nodes[i]][j]]==0){
				comp->list[comp->nlist++]=comp->c[com->nodes[i]][j];
			}
			comp->tmp[comp->c[com->nodes[i]][j]]++;
		}
	}

	for (i=0;i<comp->nlist;i++){
    c1=comp->list[i];
		sim=(2.*(comp->tmp[c1]))/(com->s+comp->size[c1]);
		if (sim>min){
    	ed.u=c1,ed.v=c2,ed.w=sim;
    	if (bip->e==bip->emax){
      	bip->emax+=EMAX;
      	bip->edges=realloc(bip->edges,bip->emax*sizeof(edge));
    	}
			bip->edges[bip->e++]=ed;
			comp->tmp[comp->list[i]]=0;
		}
	}
	comp->nlist=0;
  c2++;
}

void mkadj(bipartite *bip){
  unsigned i;
  unsigned *deg1=calloc(bip->n1,sizeof(unsigned));
  unsigned *deg2=calloc(bip->n2,sizeof(unsigned));
  for (i=0;i<bip->e;i++){
    deg1[bip->edges[i].u]++;
    deg2[bip->edges[i].v]++;
  }
  bip->cd1=malloc((bip->n1+1)*sizeof(unsigned));
  bip->cd1[0]=0;
  for (i=0;i<bip->n1;i++){
    bip->cd1[i+1]=bip->cd1[i]+deg1[i];
    deg1[i]=0;
  }
  bip->cd2=malloc((bip->n2+1)*sizeof(unsigned));
  bip->cd2[0]=0;
  for (i=0;i<bip->n2;i++){
    bip->cd2[i+1]=bip->cd2[i]+deg2[i];
    deg2[i]=0;
  }
  bip->adj1=malloc(bip->e*sizeof(unsigned));
  bip->adj2=malloc(bip->e*sizeof(unsigned));
  for (i=0;i<bip->e;i++){
    bip->adj1[bip->cd1[bip->edges[i].u]+deg1[bip->edges[i].u]++]=i;
    bip->adj2[bip->cd2[bip->edges[i].v]+deg2[bip->edges[i].v]++]=i;
  }
  free(deg1);
  free(deg2);
}


/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle,
                            const char str[])
{
  printf("%s",str);
} /* printstr */

int main(int argc, const char *argv[]) {
  unsigned i;
	community *com=alloccom();
	compare *comp=alloccompare();
  bipartite *bip=allocbip();
	FILE* file;
	double min=atof(argv[3]);
	double res=0;

	file=fopen(argv[1],"r");
	while(readlinecom(file,com)){
		com2comp(com,comp);
    bip->n1++;
	}
	fclose(file);

	file=fopen(argv[2],"r");
	while(readlinecom(file,com)){
		comsim(bip,com,comp,min);
		bip->n2++;
	}
	fclose(file);

  mkadj(bip);

  file=fopen("resbip","w");
  fprintf(file,"%u %u %u\n",bip->n1,bip->n2,bip->e);
  for (i=0;i<bip->e;i++){
    fprintf(file,"%u %u %le %u %u\n",bip->edges[i].u,bip->edges[i].v,bip->edges[i].w,bip->adj1[i],bip->adj2[i]);
  }
  for (i=0;i<bip->n1+1;i++){
    fprintf(file,"%u\n",bip->cd1[i]);
  }
  for (i=0;i<bip->n2+1;i++){
    fprintf(file,"%u\n",bip->cd2[i]);
  }
  fclose(file);

  const int numvar = bip->e;
  const int numcon = bip->n1+bip->n2;

  double *c = malloc(numvar*sizeof(double));
  for (i=0;i<numvar;i++){
    c[i]=bip->edges[i].w;
  }
  free(bip->edges);
  MSKlidxt *aptrb = malloc(numcon*sizeof(MSKlidxt));
  for (i=0;i<bip->n1;i++){
    aptrb[i]=bip->cd1[i];
  }
  for (i=0;i<bip->n2;i++){
    aptrb[bip->n1+i]=bip->e+bip->cd2[i];
  }
  MSKlidxt *aptre = malloc(numcon*sizeof(MSKlidxt));
  for (i=0;i<bip->n1;i++){
    aptre[i]=bip->cd1[i+1];
  }
  for (i=0;i<bip->n2;i++){
    aptre[bip->n1+i]=bip->e+bip->cd2[i+1];
  }
  free(bip->cd1);
  free(bip->cd2);
  MSKidxt *asub=malloc(2*bip->e*sizeof(MSKlidxt));
  double *aval=malloc(2*bip->e*sizeof(double));
  for (i=0;i<bip->e;i++){
    asub[i]=bip->adj1[i];
    asub[bip->e+i]=bip->adj2[i];
    aval[i]=1.;
    aval[bip->e+i]=1.;
  }
  free(bip->adj1);
  free(bip->adj2);
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
    for(k=0; k<numcon && r==MSK_RES_OK; ++k)
    {
      r = MSK_putconbound(task,
                          k,           /* Index of constraint.*/
                          MSK_BK_UP,      /* Bound key.*/
                          -MSK_INFINITY,      /* Numerical value of lower bound.*/
                          1.0);     /* Numerical value of upper bound.*/

      /* Input row i of A */
      if(r == MSK_RES_OK)
        r = MSK_putarow(task,
                        k,                 /* Row index.*/
                        aptre[k]-aptrb[k], /* Number of non-zeros in row i.*/
                        asub+aptrb[k],     /* Pointer to column indexes of row i.*/
                        aval+aptrb[k]);    /* Pointer to values of row i.*/
    }

    /* Maximize objective function. */
    if (r == MSK_RES_OK)
      r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);

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
            if ( xx )
            {
              MSK_getxx(task,
                        MSK_SOL_BAS,    /* Request the basic solution. */
                        xx);

              printf("Optimal primal solution\n");
              for(l=0; l<numvar; ++l)
								if (xx[l]>1e-5){
                	printf("x[%d]: %e\n",l,xx[l]);
									res+=xx[l]*c[l];
								}
            }
            else
            {
              r = MSK_RES_ERR_SPACE;
            }
						printf("objective = %le\n",res);


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
