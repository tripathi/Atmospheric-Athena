#include "copyright.h"
/*============================================================================*/
/*! \file jeans.c
 *  \brief Problem generator for a simple spherical hydrostatic atmosphere.
 *
 * PURPOSE: Problem generator for a simple spherical hydrostatic atmosphere.
 *
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static Real PlanetPot(const Real x1, const Real x2, const Real x3);
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = (pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n;
  Real x1,x2,x3,rad;
  Real G, M, cs, Rp, Rb, rho_atm, rho_c, H;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  
  G = 1.;
  M = 1.;
  Rp = 1.;
  H = .2 ; 
  cs = sqrt(H/Gamma);
  Rb = .4;
  rho_atm = 10.;
  rho_c = rho_atm/exp(G * M / cs /cs * (1/Rp - 1/Rb));


  fprintf(stderr, "Central: %f, Atmos : %f \n", rho_c, rho_atm);

/* Initialize all cells */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	rad = sqrt(x1*x1 + x2*x2 + x3*x3);

	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;

	if (rad <= Rb)
	  {
	    /* Inner constant */
	    pGrid->U[k][j][i].d = rho_c;
	    pGrid->U[k][j][i].E = rho_c * cs * cs / Gamma_1;
	  }
	else if ( rad > Rb && rad <= Rp) 
	  {
	    /* Atmosphere */
	    pGrid->U[k][j][i].d = rho_c * exp(G * M / Gamma/ cs /cs * (1/rad - 1/Rb));    
	    pGrid->U[k][j][i].E = pGrid->U[k][j][i].d * cs * cs/ Gamma_1;
	  }
	else
	  {
	    /* Ambient gas */
	    pGrid->U[k][j][i].d = .01;
	    pGrid->U[k][j][i].E = rho_atm * cs * cs / Gamma_1;
	  }
      }
    }
  }
  /* enroll gravity of planet */
  StaticGravPot = PlanetPot;

  return;
}

/*=============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * color()   - returns first passively advected scalar s[0]
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

#if (NSCALARS > 0)
/*! \fn static Real color(const GridS *pG, const int i, const int j,const int k)
 *  \brief returns first passively advected scalar s[0] */
static Real color(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"color")==0) return color;
#endif
  return NULL;
}


VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}


void Userwork_in_loop(MeshS *pM)
{
  Real x1, x2, x3, rad;
  int is,ie,js,je,ks,ke, nl, nd, i, j, k;
  GridS *pGrid;
  Real G, M, cs, Rp, Rb, rho_atm, rho_c, H;


  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
	pGrid = pM->Domain[nl][nd].Grid;          /* ptr to Grid */

	is = pGrid->is;  ie = pGrid->ie;
	js = pGrid->js;  je = pGrid->je;
	ks = pGrid->ks;  ke = pGrid->ke;

	G = 1.;
	M = 1.;
	Rp = 1.;
	H = .2 ; 
	cs = sqrt(H/Gamma);
	Rb = .4;
	rho_atm = 10.;
	rho_c = rho_atm/exp(G * M / cs /cs * (1/Rp - 1/Rb));
  
	for (k=ks; k<=ke; k++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	      rad = sqrt(x1*x1 + x2*x2 + x3*x3);

	      if (rad <= Rb)
		{
		  /* Inner constant */
		  pGrid->U[k][j][i].d = rho_c;
		  pGrid->U[k][j][i].E = rho_c * cs * cs/ (Gamma_1);
		  pGrid->U[k][j][i].M1 = 0.0;
		  pGrid->U[k][j][i].M2 = 0.0;
		  pGrid->U[k][j][i].M3 = 0.0;
		  
		}
	      else if ( rad > Rb && rad <= (Rb + 5. * pGrid->dx1) )
		{
		  /* Atmosphere */
		  pGrid->U[k][j][i].d = rho_c * exp(G * M / Gamma/ cs /cs * (1/rad - 1/Rb));    
		  pGrid->U[k][j][i].E = pGrid->U[k][j][i].d * cs * cs / Gamma_1;
		  pGrid->U[k][j][i].M1 = 0.0;
		  pGrid->U[k][j][i].M2 = 0.0;
		  pGrid->U[k][j][i].M3 = 0.0;
		  
		}
	    }
	  }
	}
	return;
      }
    }
  }
	      
}


void Userwork_after_loop(MeshS *pM)
{
}

static Real PlanetPot(const Real x1, const Real x2, const Real x3)
{
  Real rad,phi=0.0;
  Real Ggrav=1;
  Real Rsoft = .05;
  rad = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
  
  phi = -1.0/(rad+Rsoft);
  return phi;
}
