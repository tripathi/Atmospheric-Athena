/*============================================================================*/
/*! \file ioniz_spher.c
 *  \brief Problem generator for an ionization front irradiating a spherical adiabatic atmosphere.
 *
 * PURPOSE: Problem generator for an ionization front irradiating a spherical adiabatic atmosphere.
 *
 * AUTHORS: A. Tripathi, M. Krumholz, X. Bai
 *============================================================================*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#define TII 1.0e4
#define ALPHA4 2.59e-13

/* static int radplanecount; */
static Real GMstar, adist, GM,K,Cp,rp,rreset2,rho0,Rsoft, trad, flux;
static Real TidalPot(const Real x1, const Real x2, const Real x3);
static Real PlanetPot(const Real x1, const Real x2, const Real x3);
static Real print_flux(const GridS *pG, const int i, const int j, const int k);

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real cs, n_H, mu;
  Real rho, pressure;

  Real x1,x2,x3; /*For setting up atmos*/
  Real rad,myrho,rhop,powindex,Ggrav; /*For setting up atmos*/
  Real rhoout, np, mp; /*For setting up atmos*/
  Real rin, rout, rhoedge;

  int radplanecount =0;


  /* Set up ionizing source at box edge. 
   * The user-input parameters are n_H (initial density),
   * cs (initial sound speed in neutral gas), flux (ionizing flux,
   * in units of photons per unit time per unit area), and trad
   * (the onset time of radiation).
   */
  mu = par_getd("ionradiation", "mu");
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

  /*Derived planet quantities*/
  GM = Ggrav * mp;
  rhop = np * mu;
  Rsoft= 0.01*rp;

  /* For use in trying to add the Coriolis force - Not yet working */
#ifdef SHEARING_BOX
  adist = 7.48e11; /* .05AU in cm*/
  GMstar = Ggrav * 1.99e33; /*Using 1 solar mass in g*/
  Omega_0 = sqrt(GMstar / (pow(adist,3)));
  qshear = 0;
#endif

  /* Inner radii - change accordingly */
  rin = 0.5*rp; /* Radius within which the profile is constant */
  rreset2 = 0.5625*rp*rp; /* Reset radius^2 - within which the profile is fixed */

  /* rin  = rp - 5.0*pGrid->dx1;  */
  /* rout = rp + 15.0*pGrid->dx1;  */

  powindex   = 1.0/Gamma_1;

  /* coefficient K (P=Krho^gamma) determined by rhop and cs */
  K  = pow(rhop,-Gamma_1)*cs*cs;

  /*Density at inner boundary */
  rho0= pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1.0/rp - 1.0/rin),powindex);

  /* integration constant */
  Cp = pow(rho0,Gamma_1) - (Gamma_1/Gamma)*GM/K/rin;

  /*Outer (density matching) radius for ambient gas*/
  rhoedge = rhop/10;  
  rout = 1./(Gamma/Gamma_1/GM*K*(pow(rhoedge, Gamma_1) - pow(rho0, Gamma_1)) + 1./rin); 


  /* if ((rout - rreset) < 5.0*pGrid->dx1) */
  /*   ath_error("[sphere]: Insufficient separation between reset and outer radii"); */
  /* if ((rreset - rin) < 5.0*pGrid->dx1) */
  /*   ath_error("[sphere]: At least 5 cells needed for reconstruction"); */

  /*Density at atmosphere's edge*/
  rhoout = rhoedge/10000.;


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
	  pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
	} else if (rad > rout){
	  pGrid->U[k][j][i].d  = rhoout;
	  pGrid->U[k][j][i].E  = K*pow(rhoedge,Gamma)/Gamma_1;
	  pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d * 1.0e-4; /*Make the ambient background ionized*/
	} else {
	  pGrid->U[k][j][i].d = pow(Gamma_1/Gamma*GM/K/MAX(rad,TINY_NUMBER) + Cp,powindex);
	  pGrid->U[k][j][i].E  = K*pow(pGrid->U[k][j][i].d,Gamma)/Gamma_1;
	  pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
	}

      }
    }
  }

  /* Radiation originating from root domain edge */
  /*   if ((pDomain->Level == 0) && (pDomain->DomNumber==0)){ */
  /*     ath_pout(0,"On domain level %d, number %d: Adding radiator on root domain \n",  pDomain->Level, pDomain->DomNumber); */
#ifdef ION_RADIATION  
   if (par_geti("problem","nradplanes") == radplanecount) {
      ath_error("Zero radplanes specified in input file\n");
    }
   if (par_geti("problem","nradplanes") != radplanecount) {
      radplanecount++;
      trad = par_getd("problem", "trad");
      if (radplanecount != par_geti("problem","nradplanes")) {
	ath_error("More than 1 radplane specified in input file. Unable to add requested radplane\n");
      }
      if (radplanecount == 1 && trad < TINY_NUMBER) {
	fprintf(stderr,"I'm in here  because trad is %g \n", trad);
	add_radplane_3d(pGrid, -1, flux);
	radplanecount = 0;
      } else if (trad >= TINY_NUMBER){
	ath_error("Delaying ionization not currently enabled for use with SMR + MPI \n");
      }
   }
#endif
/*   } */
/*   else {ath_pout(0,"On domain level %d, number %d: Not adding a radiator plane\n", pDomain->Level, pDomain->DomNumber);} */

/* enroll gravity of planet */
  StaticGravPot = PlanetPot;
#ifdef SHEARING_BOX
  ShearingBoxPot = TidalPot;
#endif

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

  Real mp, np, cs, Ggrav, rin, powindex, rhop, mu;
  int nl, nd;
  GridS *pGrid;

  rp  = par_getd_def("problem","rp",1.2e10);
  mp  = par_getd_def("problem","mp",1.0e30);
  np = par_getd_def("problem","np",6.0e8);
  cs = par_getd("problem","cs");
  mu = par_getd("ionradiation", "mu");

  Ggrav = 6.67e-8;
  GM = Ggrav * mp;

  /*For use when trying to implement the Coriolis force (not yet working)*/
#ifdef SHEARING_BOX
 adist = 7.48e11; /* .05AU in cm*/
 GMstar = Ggrav * 1.99e33; /*Using 1 solar mass in g*/
 Omega_0 = sqrt(GMstar / (pow(adist,3)));
 qshear = 0;
#endif

  rin = 0.5*rp;
  rreset2 = 0.5625*rp*rp; /*(0.75Rp)^2*/

  powindex   = 1.0/Gamma_1;
  rhop = np * mu;
  K  = pow(rhop,-Gamma_1)*cs*cs;
  rho0= pow( pow(rhop,Gamma_1) - Gamma_1/Gamma*GM/K*(1.0/rp - 1.0/rin),powindex);
  Cp = pow(rho0,Gamma_1) - (Gamma_1/Gamma)*GM/K/rin;

  Rsoft= 0.01*rp;
#ifdef ION_RADIATION  
  if (par_geti("problem","nradplanes") == 1) {  
    flux = par_getd("problem","flux");
    (pM->radplanelist)->dir[0] = -1;
    (pM->radplanelist)->flux_i = flux;  
  }
#endif

  StaticGravPot = PlanetPot;
#ifdef SHEARING_BOX
  ShearingBoxPot = TidalPot;
#endif

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
#ifndef ION_RADIATION
  fprintf(stdout, "Ion radiation not turned on. Cannot output flux. \n");
#else
  if(strcmp(expr,"flux")==0) return print_flux;
#endif
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  /*Use this function to fix the profile within rreset, at each timestep*/

  Real x1,x2,x3,rad2,myrho,powindex;
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
	      rad2 = x1*x1 + x2*x2 + x3*x3;

	      if ((pGrid->U[k][j][i].d < 0) || isnan(pGrid->U[k][j][i].d)) fprintf(stderr, "Neg or NaN dens: %e at k:%d j:%d i:%d, rad: %f, lev:%d \n",pGrid->U[k][j][i].d, k, j, i, rad2, nl);
	      if (pGrid->U[k][j][i].E < 0 || isnan(pGrid->U[k][j][i].E)) fprintf(stderr, "Neg or NaN E: %e at k:%d j:%d i:%d, rad: %f, lev:%d \n",pGrid->U[k][j][i].E, k, j, i, rad2, nl);

	      /*Reset values within the boundary*/
	      if (rad2 <= rreset2) { /*AT: May need to change this to a smaller r*/
		myrho = pow(Gamma_1/Gamma*GM/K/MAX(sqrt(rad2),TINY_NUMBER) + Cp,powindex);
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

  /* if (pM->time > trad)  { */
  /*   add_radplane_3d(pGrid, -1, flux); */
  /*   radplanecount = -999; /\*To avoid entering this loop again*\/ */
  /* } */
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
#ifndef SHEARING_BOX
  Real adist = 7.48e11; /* .05AU in cm*/
  Real GMstar = 6.67e-8 * 1.99e33; /*Using 1 solar mass in g*/
  Real omega = sqrt(GMstar / (pow(adist,3)));
  Real radstar = sqrt(SQR(x1+adist) + SQR(x2) + SQR(x3));
  return -GM/(rad+Rsoft)-GMstar/radstar -.5*SQR(omega*radstar); 
#else 
  return -GM/(rad+Rsoft)
#endif

}

/*------------------------------------------------------------------------------
 *! \fn static Real TidalPot(const Real x1, const Real x2,const Real x3)
 * \brief tidal potential in shearing box
 */

static Real TidalPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifdef SHEARING_BOX
  Real radstar2 = SQR(x1+adist) + SQR(x2) + SQR(x3);
  phi -= GMstar/sqrt(radstar2) + .5*SQR(Omega_0)*radstar2;
#endif
  return phi;
}


#ifdef ION_RADIATION
/*----------------------------------------------------------------------------*/
/*! \fn static print_flux(const Grid *pG,const int i,const int j,
 *      const int k)
 *  \brief Returns flux
 */

static Real print_flux(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->EdgeFlux[k-pG->ks][j-pG->js][i-pG->is]);
}
#endif
