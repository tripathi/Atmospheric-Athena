#include "copyright.h"
/*============================================================================*/
/*! \file blast.c
 *  \brief Problem generator for spherical blast wave problem.
 *
 * PURPOSE: Problem generator for spherical blast wave problem.
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.     */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static Real GM,K,C0,rp,rin,rout,rho0,Rsoft;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real x1,x2,x3;
  Real rad,rhoe,rhob,myrho,s;

  /* radius of the sphere (planet) */
  rp  = par_getd("problem","radius");
  /* numerical gravitational constant */
  GM  = par_getd_def("problem","GM",1.0);
  /* density of the ambient medium */
  rhoe= par_getd_def("problem","rhoe",1.0);
  /* density at the base of the atmosphere */
  rhob= par_getd_def("problem","rhob",100.0);

  /* integration constant determined by rho0 */
  C0 = (Gamma_1/Gamma)*pow(rhoe,-Gamma_1);

  /* coefficient K (P=Krho^gamma) determined by rhob and C0 */
  K  = GM/rp/(pow(rhob,Gamma_1)*Gamma/Gamma_1-1.0/C0);
  fprintf (stderr, "rp: %e, C0: %e, rhob: %e K: %e \n", rp, C0, rhob, K);


  if (rp < 10.0*pGrid->dx1)
    ath_error("[sphere]: rp must be larger than 10 grid cell.\n");
  if (rp > 0.25*(pDomain->RootMaxX[0]-pDomain->RootMinX[0]))
    ath_error("[sphere]: rp must be smaller than 1/4 box size.\n");

  rin  = rp - 10.0*pGrid->dx1;
  rout = rp+ 5.0*pGrid->dx1;
 /* + 5.0*pGrid->dx1; */
  Rsoft= pGrid->dx1;

  s   = 1.0/Gamma_1;
  rho0= pow(Gamma_1/Gamma*(1.0/C0+GM/K/rin),s);

/* setup uniform ambient medium with spherical over-pressured region */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(x1*x1 + x2*x2 + x3*x3);

	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;

	if (rad > rout) {
	  pGrid->U[k][j][i].d  = rho0/10000.;
	  myrho = pow(Gamma_1/Gamma*(1.0/C0+GM/K/rout),s);
	  pGrid->U[k][j][i].E = K*pow(myrho,Gamma)/Gamma_1;

	} else {
	  myrho = pow(Gamma_1/Gamma*(1.0/C0+GM/K/MAX(rad,TINY_NUMBER)),s);
	  pGrid->U[k][j][i].d  = MIN(rho0,myrho);
	  pGrid->U[k][j][i].E  = K*pow(pGrid->U[k][j][i].d,Gamma)/Gamma_1;
	}
      }
    }
  }

/* enroll gravity of planet */
  StaticGravPot = PlanetPot;

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{

  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif


void Userwork_in_loop(MeshS *pM)
{
  Real x1,x2,x3,rad,T,T0,myrho,fac,s;
  int is,ie,js,je,ks,ke, nl, nd, i, j, k;
  GridS *pGrid;

  /* fprintf(stderr,"rp: %f, rin: %f", rp, rin); */
  s=1.0/Gamma_1;

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
	pGrid = pM->Domain[nl][nd].Grid;          /* ptr to Grid */

	is = pGrid->is;  ie = pGrid->ie;
	js = pGrid->js;  je = pGrid->je;
	ks = pGrid->ks;  ke = pGrid->ke;
	
	for (k=ks; k<=ke; k++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	      rad = sqrt(x1*x1 + x2*x2 + x3*x3);

           myrho = pow(Gamma_1/Gamma*(1.0/C0+GM/K/MAX(rad,TINY_NUMBER)),s);
              myrho = MIN(myrho, rho0);
	      
	      /*Reset values within the boundary*/
	      if (rad <= rp) {
		pGrid->U[k][j][i].d  = myrho;
		pGrid->U[k][j][i].M1 = 0.0;
		pGrid->U[k][j][i].M2 = 0.0;
		pGrid->U[k][j][i].M3 = 0.0;
		pGrid->U[k][j][i].E = K*pow(myrho,Gamma)/Gamma_1; 
	      }
/*
              else if (rad < rout) {
                fac = (rad-rp)/(rout-rp);
                T   = Gamma_1*(pGrid->U[k][j][i].E
                    -0.5*(SQR(pGrid->U[k][j][i].M1)+SQR(pGrid->U[k][j][i].M2)
                         +SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d)/pGrid->U[k][j][i].d;
                T0  = K*pow(myrho,Gamma_1);
                T   = fac*T + (1.0-fac)*T0;
                pGrid->U[k][j][i].d  = fac*pGrid->U[k][j][i].d+(1.0-fac)*myrho;
		pGrid->U[k][j][i].M1 *= fac;
		pGrid->U[k][j][i].M2 *= fac;
		pGrid->U[k][j][i].M3 *= fac;
		pGrid->U[k][j][i].E = pGrid->U[k][j][i].d*T/Gamma_1
                                      -0.5*(SQR(pGrid->U[k][j][i].M1)+SQR(pGrid->U[k][j][i].M2)
                                           +SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
	      }
*/
	    }
	  }
	}
      }
    }
  }
  return;
}

void Userwork_after_loop(MeshS *pM)
{
}

/*------------------------------------------------------------------------------
 * PlanetPot:
 */

static Real PlanetPot(const Real x1, const Real x2, const Real x3)
{
  Real rad = sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  return -GM/(rad+Rsoft);
}
