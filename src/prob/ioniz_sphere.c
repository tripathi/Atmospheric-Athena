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

static Real GM,K,Cp,rp,rin,rout,rreset,rho0,Rsoft;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int radplanecount =0;
  Real cs, n_H, m_H, flux;
  Real rho, pressure;

  Real x1,x2,x3; /*For setting up atmos*/
  Real rad,myrho,rhop,powindex,Ggrav; /*For setting up atmos*/
  Real rhoe, np, mp; /*For setting up atmos*/

  /* Set up ionizing source at box edge. 
   * The user-input parameters are n_H (initial density),
   * cs (initial sound speed in neutral gas), flux (ionizing flux,
   * in units of photons per unit time per unit area).
   */
  m_H = par_getd("ionradiation", "m_H");
  n_H = par_getd("problem","n_H");
  cs = par_getd("problem","cs");
  flux = par_getd("problem","flux");

  
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
  rout = rp;
    

  if ((rout - rreset) < 5.0*pGrid->dx1)
    ath_error("[sphere]: Insufficient separation between reset and outer radii");
  if ((rreset - rin) < 5.0*pGrid->dx1)
    ath_error("[sphere]: At least 5 cells needed for reconstruction");
  

  /* rin  = rp - 5.0*pGrid->dx1; /\*What I refer to as r0*\/ */
  /* rout = rp + 15.0*pGrid->dx1;  */

  powindex   = 1.0/Gamma_1;

  /* coefficient K (P=Krho^gamma) determined by rhop and cs */
  K  = pow(rhop,-Gamma_1)*cs*cs;

  /*Density at inner boundary */
  rho0= pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1/rp - 1/rin),powindex);
  rhoe = pow(Gamma_1/Gamma*GM/K/(1.1*rp) + Cp,powindex);
  /*  rhoe= n_H * m_H;*/
/* pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1./rp),powindex); */
  /* fprintf(stderr, "Comparison with Xuening params: rho0: %e, rhop: %e,  rhoe : %e \n", rho0, rhop, rhoe); */

  /* integration constant */
  Cp = pow(rho0,Gamma_1) - (Gamma_1/Gamma)*GM/K/rin;

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

	if (rad > rout) {
	  pGrid->U[k][j][i].d = rhoe; /*AT: Should be changed to indpt athinput file I.C. */
	  myrho = pow(Gamma_1/Gamma*GM/K/rout + Cp,powindex);
	  pGrid->U[k][j][i].E = K*pow(myrho,Gamma)/Gamma_1;
	  /*	  fprintf(stderr,"myrho: %e, cs2myrho: %e, rhoout: %e, csinsq: %e csoutsq %e \n",  myrho, pow(Gamma_1/Gamma*GM/K/rout + Cp, powindex), pGrid->U[k][j][i].E * Gamma_1/myrho, pGrid->U[k][j][i].E * Gamma_1/pow(Gamma_1/Gamma*GM/K/rout + Cp, powindex), pGrid->U[k][j][i].E * Gamma_1/ pGrid->U[k][j][i].d);*/
	  /* fprintf(stderr, "rad/radout: %f, dens: %f, rhoout: %f, powindex %f energ: %f \n", rad/rout, pGrid->U[k][j][i].d, Gamma_1/Gamma*GM/K/rout + Cp, powindex, pGrid->U[k][j][i].E); */
	} else {
	  myrho = pow(Gamma_1/Gamma*GM/K/MAX(rad,TINY_NUMBER) + Cp,powindex);
	  pGrid->U[k][j][i].d  = MIN(rho0,myrho);
	  pGrid->U[k][j][i].E  = K*pow(pGrid->U[k][j][i].d,Gamma)/Gamma_1;
	}
	pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
      }
    }
  }

  /* Radiation originating from root domain edge */
  /*   if ((pDomain->Level == 0) && (pDomain->DomNumber==0)){ */
  /*     ath_pout(0,"On domain level %d, number %d: Adding radiator on root domain \n",  pDomain->Level, pDomain->DomNumber); */
  
   if (radplanecount == par_geti("problem","nradplanes")) {
      ath_error("Invalid number of radplanes specified in input file\n");
    }
   if (radplanecount != par_geti("problem","nradplanes")) {
      add_radplane_3d(pGrid, -1, flux);
      radplanecount++;
      if (radplanecount != par_geti("problem","nradplanes")) {
	ath_error("Invalid number of radplanes specified in input file\n");
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
	      }
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
