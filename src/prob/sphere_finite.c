#include "copyright.h"
/*============================================================================*/
/*! \file sphere_finite.c
 *  \brief Problem generator for a spherical adiabatic atmosphere.
 *
 * PURPOSE: Problem generator for spherical adiabatic atmosphere.
 *
 * AUTHORS: X. Bai, A. Tripathi
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static Real GM,K,Cp,rp,rin,rout,rho0,Rsoft;
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
  Real rad,myrho,s;
  Real rhoe, rhop, csq;

  /* radius of the sphere (planet) */
  rp  = par_getd("problem","radius");
  /* numerical gravitational constant */
  GM  = par_getd_def("problem","GM",1.0);

  /* /\* density of the ambient medium *\/ */
  /* rhoe = par_getd_def("problem","rhoe",1.0); */

  /* density at rp */
  rhop = par_getd_def("problem","rhop",100.0);
  /* /\* temperature at rp*\/ */
  /* Tp = par_getd_def("problem","T0",100.0);    /\*NEEDS SETTING*\/ */

  if (rp < 10.0*pGrid->dx1)
    ath_error("[sphere]: rp must be larger than 10 grid cell.\n");
  if (rp > 0.25*(pDomain->RootMaxX[0]-pDomain->RootMinX[0]))
    ath_error("[sphere]: rp must be smaller than 1/4 box size.\n");

  rin  = rp - 5.0*pGrid->dx1; /*What I refer to as r0*/
  rout = rp + 15.0*pGrid->dx1; 
  Rsoft= pGrid->dx1;

  s   = 1.0/Gamma_1;


  csq = 0.083894;


  /* komu = GM / rp / Tp ;/\*NEEDS SETTING*\/ */

  
  /* coefficient K (P=Krho^gamma) determined by rhop and Tp */
  K  = pow(rhop,-Gamma_1)* csq;
  
  /*Density at inner boundary */
  rho0= pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1/rp - 1/rin),s);
  rhoe= pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1./rp),s);
  fprintf(stderr, "rho0: %e, rhoe : %e \n", rho0, rhoe);

  /* integration constant */
  Cp = pow(rho0,Gamma_1) - (Gamma_1/Gamma)*GM/K/rin;


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
	  pGrid->U[k][j][i].d  = rho0/100000.;
	  myrho = pow(Gamma_1/Gamma*GM/K/rout + Cp,s);
	  pGrid->U[k][j][i].E = K*pow(myrho,Gamma)/Gamma_1;
	} else {
	  myrho = pow(Gamma_1/Gamma*GM/K/MAX(rad,TINY_NUMBER) + Cp,s);
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
  Real x1,x2,x3,rad,myrho,s;
  int is,ie,js,je,ks,ke, nl, nd, i, j, k;
  GridS *pGrid;

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

              /* myrho = pow(Gamma_1/Gamma*(1.0/C0+GM/K/MAX(rad,TINY_NUMBER)),s); */   
	      myrho = pow(Gamma_1/Gamma*GM/K/MAX(rad,TINY_NUMBER) + Cp,s);
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
