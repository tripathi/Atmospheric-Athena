/*============================================================================*/
/*! \file ioniz_spher.c
 *  \brief Problem generator for an ionization front irradiating a spherical adiabatic atmosphere.
 *
 * PURPOSE: Problem generator for an ionization front irradiating a spherical adiabatic atmosphere.
 *
 * AUTHORS: A. Tripathi, M. Krumholz, X. Bai
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#define TII 1.0e4
#define ALPHA4 2.59e-13

static int radplanecount;
static Real GM,K,Cp,rp,rreset,rho0,Rsoft, trad, flux;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real cs, n_H, m_H;
  Real rho, pressure;

  Real x1,x2,x3; /*For setting up atmos*/
  Real rad,myrho,rhop,powindex,Ggrav; /*For setting up atmos*/
  Real rhoout, np, mp; /*For setting up atmos*/
  Real rin, rout;

  /* Set up ionizing source at box edge. 
   * The user-input parameters are n_H (initial density),
   * cs (initial sound speed in neutral gas), flux (ionizing flux,
   * in units of photons per unit time per unit area), and trad
   * (the onset time of radiation).
   */
  trad = par_getd("problem", "trad");
  m_H = par_getd("ionradiation", "m_H");
  n_H = par_getd("problem","n_H");
  cs = par_getd("problem","cs");
  flux = par_getd("problem","flux");

  radplanecount =0;
  
  /*Set up planet atmosphere.
   * User-input params are rp (planet radius),
   * mp (planet mass), and rhop (density at rp)
   */
  rp  = par_getd_def("problem","rp",1.2e10);
  mp  = par_getd_def("problem","mp",1.0e30);
  np = par_getd_def("problem","np",6.0e8);

  Ggrav = 6.67e-8;
  GM = Ggrav * mp;
  rhop = np * m_H;
  Rsoft= pGrid->dx1;

  rin = 0.5*rp;
  rreset = 0.75*rp;
    

  

  /* rin  = rp - 5.0*pGrid->dx1; /\*What I refer to as r0*\/ */
  /* rout = rp + 15.0*pGrid->dx1;  */

  powindex   = 1.0/Gamma_1;

  /* coefficient K (P=Krho^gamma) determined by rhop and cs */
  K  = pow(rhop,-Gamma_1)*cs*cs;

  /*Density at inner boundary */
  rho0= pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1.0/rp - 1.0/rin),powindex);

  /* integration constant */
  Cp = pow(rho0,Gamma_1) - (Gamma_1/Gamma)*GM/K/rin;

  /*Outer (density matching) radius for ambient gas*/
  rout = 1./(Gamma/Gamma_1/GM*K*(pow(rhop/10000, Gamma_1) - pow(rho0, Gamma_1)) + 1./rin);


  if ((rout - rreset) < 5.0*pGrid->dx1)
    ath_error("[sphere]: Insufficient separation between reset and outer radii");
  if ((rreset - rin) < 5.0*pGrid->dx1)
    ath_error("[sphere]: At least 5 cells needed for reconstruction");

  /*Density at atmosphere's edge*/
  rhoout = rhop/10000;
/* pow(Gamma_1/Gamma*GM/K/rout + Cp,powindex); */
  /* fprintf(stderr, "rhoout %e \n", rhoout/10000.); */

  /* fprintf(stderr, "K : %f, Cp: %f powindex: %f, rho_out: %f\n", K, Cp, powindex, (Gamma_1/Gamma*GM/K/rout + Cp)); */

  /* Power-law pressure and density */
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(x1*x1 + x2*x2 + x3*x3);
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;

	if (rad <= rin){
	  pGrid->U[k][j][i].d  = rho0;
	  pGrid->U[k][j][i].E  = K*pow(pGrid->U[k][j][i].d,Gamma)/Gamma_1;
	} else if (rad > rout){
	  pGrid->U[k][j][i].d  = rhoout;
	  pGrid->U[k][j][i].E  = K*pow(rhoout,Gamma)/Gamma_1;
	} else {
	  pGrid->U[k][j][i].d = pow(Gamma_1/Gamma*GM/K/MAX(rad,TINY_NUMBER) + Cp,powindex);
	  pGrid->U[k][j][i].E  = K*pow(pGrid->U[k][j][i].d,Gamma)/Gamma_1;
	}
      pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
      }
    }
  }

  /* Radiation originating from root domain edge */
  /*   if ((pDomain->Level == 0) && (pDomain->DomNumber==0)){ */
  /*     ath_pout(0,"On domain level %d, number %d: Adding radiator on root domain \n",  pDomain->Level, pDomain->DomNumber); */
  
   if (par_geti("problem","nradplanes") == radplanecount) {
      ath_error("Zero radplanes specified in input file\n");
    }
   if (par_geti("problem","nradplanes") != radplanecount) {
      radplanecount++;
      if (radplanecount != par_geti("problem","nradplanes")) {
	ath_error("More than 1 radplane specified in input file. Unable to add requested radplane\n");
      }
      if (radplanecount == 1 && trad < TINY_NUMBER) {
	add_radplane_3d(pGrid, -1, flux);
	radplanecount = 0;
      }
    }
/*   } */
/*   else {ath_pout(0,"On domain level %d, number %d: Not adding a radiator plane\n", pDomain->Level, pDomain->DomNumber);} */

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

void Userwork_in_loop(MeshS *pM)
{
  Real x1,x2,x3,rad,myrho,powindex;
  int is,ie,js,je,ks,ke, nl, nd, i, j, k;
  GridS *pGrid;

  powindex=1.0/Gamma_1;

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

	      /*Reset values within the boundary*/
	      if (rad <= rreset) { /*AT: May need to change this to a smaller r*/
		myrho = pow(Gamma_1/Gamma*GM/K/MAX(rad,TINY_NUMBER) + Cp,powindex);
		myrho = MIN(myrho, rho0);
		pGrid->U[k][j][i].d  = myrho;
		pGrid->U[k][j][i].M1 = 0.0;
		pGrid->U[k][j][i].M2 = 0.0;
		pGrid->U[k][j][i].M3 = 0.0;
		pGrid->U[k][j][i].E = K*pow(myrho,Gamma)/Gamma_1; 
		pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
	      }
	    }
	  }
	}
      }
    }
  }

  if (radplanecount == 1 && pM->time > trad) 
      add_radplane_3d(pGrid, -1, flux);
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
