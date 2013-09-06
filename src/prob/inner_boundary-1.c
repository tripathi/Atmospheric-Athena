#include "copyright.h"
/*==============================================================================
 * FILE: inner_boundary-1.c
 *
 * PURPOSE: Problem generator for a planet with a hydrostatic atmosphere
 *  Input parameters are:
 *      problem/Mp      = planet mass
 *      problem/Rsoft   = softening length for planet's grav. potential
 *      problem/Xplanet = x-location of planet
 *      problem/Yplanet = y-location of planet
 *      problem/Zplanet = z-location of planet
 *      problem/Rbound  = radius within which the density is constant
 *      problem/Rp      = outer edge of planetary radius
 *      problem/rho_at  = density of planetary atmosphere at Rp
 *      problem/T_at    = temperature of planetary atmosphere at Rp
 *      problem/drat    = density ratio between atmos. and ambient gas
 *
 *  Based off the template for shk_cloud
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


static Real Mp,Rsoft,Xplanet,Yplanet,Zplanet;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);
const Real kb=1.38e-16;
const Real mH = 1.67e-24;
const Real Ggrav = 6.67e-8;
static Real rho_c, Hor2, c_s, Rp, Rb;


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  Real x1,x2,x3,diag;
  Real drat;
  Real pr;
  Real rho_at,T_at, p_c;

/* Read input parameters */
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rb = par_getd("problem","Rbound");
  Rp =par_getd("problem","Rp");
  rho_at = par_getd("problem","rho_at");
  T_at =par_getd("problem","T_at");
  drat = par_getd("problem","drat");
  
/* Set up parameters for central part of "planet" */
  Hor2=kb*T_at/(2*mH*Ggrav*Mp);
  rho_c= rho_at*exp(-1/Hor2*(1/Rp - 1/Rb));
		    
  c_s = sqrt(kb*T_at/(2*mH));
  p_c = pow(rho_c,Gamma)*(c_s)*(c_s);

/* Initialize the grid */
  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        diag = sqrt(x1*x1 + x2*x2 + x3*x3);

	/* Outside the planet */
	pGrid->U[k][j][i].d  = rho_at/drat;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;
	
#ifdef ADIABATIC
	pGrid->U[k][j][i].E = pow(rho_at,Gamma)*(c_s*c_s)/Gamma_1;  
#endif
#if (NSCALARS > 0)
	pGrid->U[k][j][i].s[0] = 0.0;
#endif
	
	/* Within the inner boundary*/	
	if (diag < Rb){
	  pGrid->U[k][j][i].d  = rho_c*0.99;
	  pGrid->U[k][j][i].M1 = 0.0;
	  pGrid->U[k][j][i].M2 = 0.0;
	  pGrid->U[k][j][i].M3 = 0.0;
	  pr = pow(rho_c*0.99,Gamma)*(c_s*c_s);
#ifdef ADIABATIC
	  pGrid->U[k][j][i].E = pr/Gamma_1;
#endif
	}else if (diag <= Rp) {

	  /* Planet's atmosphere*/
          pGrid->U[k][j][i].d  = rho_c*exp((1/diag-1/Rb)/Hor2); 
          pGrid->U[k][j][i].M1 = 0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
	  pr = pow(pGrid->U[k][j][i].d,Gamma)*(c_s*c_s);
#ifdef ADIABATIC
	  pGrid->U[k][j][i].E = pr/Gamma_1;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = drat;
#endif

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
 * color()   - returns first passively advected scalar s[0]
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  /* enroll planet graviy */
  StaticGravPot = PlanetPot;
  return;
}

#if (NSCALARS > 0)
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
  Real x1, x2, x3, diag;
  int is,ie,js,je,ks,ke, nl, nd, i, j, k;
  Real pr;
  GridS *pGrid;


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
	      diag = sqrt(x1*x1 + x2*x2 + x3*x3);
	      
	      /*Reset values within the boundary*/
	      if (diag <= Rb ){
		pGrid->U[k][j][i].d  = rho_c*0.99;
		pGrid->U[k][j][i].M1 = 0.0;
		pGrid->U[k][j][i].M2 = 0.0;
		pGrid->U[k][j][i].M3 = 0.0;
		pr = pow(rho_c*0.99,Gamma)*(c_s*c_s);
#ifdef ADIABATIC
		pGrid->U[k][j][i].E = pr/Gamma_1;
#endif
		/*also reset values up to 4 cells beyond Rb*/
	      } else if (diag > Rb && diag<= Rb + 4. * pGrid->dx1){
		pGrid->U[k][j][i].d  = rho_c*exp((1/diag-1/Rb)/Hor2);
		pGrid->U[k][j][i].M1 =0.0;
		pGrid->U[k][j][i].M2 = 0.0;
		pGrid->U[k][j][i].M3 = 0.0;
		pr = pow(pGrid->U[k][j][i].d,Gamma)*(c_s*c_s);
#ifdef ADIABATIC
		pGrid->U[k][j][i].E = pr/Gamma_1;
#endif
#if (NSCALARS > 0)
		pGrid->U[k][j][i].s[0] = drat;
#endif
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
  return;
}

/*------------------------------------------------------------------------------
 * PlanetPot:
 */

static Real PlanetPot(const Real x1, const Real x2, const Real x3)
{
  Real rad,phi=0.0;
  Real Ggrav=6.67e-8;
  rad = sqrt(SQR(x1-Xplanet) + SQR(x2-Yplanet) + SQR(x3-Zplanet));
  
  phi = -1.0*Mp*Ggrav/(rad+Rsoft);
  return phi;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
