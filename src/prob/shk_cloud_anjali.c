#include "copyright.h"
/*============================================================================*/
/*! \file shk_cloud.c
 *  \brief Problem generator for shock-cloud problem; planar shock impacting
 *   a single spherical cloud.
 *
 * PURPOSE: Problem generator for shock-cloud problem; planar shock impacting
 *   a single spherical cloud.  Input parameters are:
 *    - problem/Mach   = Mach number of incident shock
 *    - problem/beta   = ratio of Pgas/Pmag
 *    - problem/iprob  = integer flag to determine problem
 *
 *   The cloud radius is fixed at 1.0.  The center of the coordinate system
 *   defines the center of the cloud, and should be in the middle of the cloud.
 *   The shock is initially at x1=-2.0.  A typical grid domain should span
 *   x1 in [-3.0,7.0] , y and z in [-2.5,2.5] (see input file in /tst)
 *   Various test cases are possible:
 *   - (iprob=1): B parallel to shock normal
 *   - (iprob=2): B perpendicular to shock normal -- NOT YET IMPLEM        diag = sqrt(x1*x1 + x2*x2 + x3*x3);
ENTED
 *
 *   If the code is configured with nscalars>0, the cloud material is labeled
 *   with U[k][j][i].s[0]=1.						      
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
//#include "constants_wind.h"

/* postshock flow variables are shared with IIB function */

static Real dl,pl,ul;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * shk_cloud_iib() - fixes BCs on L-x1 (left edge) of grid to postshock flow.
 * PlanetPot()   - static gravitational potential of planet - A.T: Added 12/7
 *============================================================================*/

void shk_cloud_iib(GridS *pGrid);
static Real Mp,Rsoft,Xplanet,Yplanet,Zplanet;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);

static const Real Gcgs = 6.674e-8;   //cm^3 g^-1 s^-2
static const Real Mjup = 1.8986e30;  //g
static const Real Rjup = 6.9911e9;   //cm
static const Real kb = 1.3806/1.e16;   //erg K^-1
static const Real mh = 1.6737/1.e24;   //g 
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real jump1, jump2, jump3;
  Real x1,x2,x3,xshock;
  Real Mach;
  Real dr,pr,ur;
  Real gconst;


/* Read input parameters */

  xshock = -2.0;
  Mach = par_getd("problem","Mach");
  iprob = par_geti("problem","iprob");
  iprob = 1;
	
  Mp = par_getd_def("problem","Mplanet",1.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);

  gconst = Gcgs*Mp*Mjup/Rjup/Rjup;
  
/* Set parameters of ambient medium */

  dr = 1.0e-15 * exp(-gconst /(Iso_csound*Iso_csound)); /*Density in g cm^-3 */
#ifdef ADIABATIC
  pr = 1.0/Gamma; //!!!Need to scale appropriately
#else
  pr = 1.0; //Should change to appropriate scaling of kT
#endif
  ur = 0.0;

/* Uses Rankine Hugoniot relations for adiabatic gas to initialize problem */

//  jump1 = (Gamma + 1.0)/(Gamma_1 + 2.0/(Mach*Mach));
// jump2 = (2.0*Gamma*Mach*Mach - Gamma_1)/(Gamma + 1.0);
//  jump3 = 2.0*(1.0 - 1.0/(Mach*Mach))/(Gamma + 1.0);
  jump1 = 1.;
  jump2 = 1.;
  jump3 = 1.;


  dl = dr*jump1;
  pl = pr*jump2;
#ifdef ISOTHERMAL
  ul = ur + jump3*Mach*Iso_csound;
#else
  ul = ur + jump3*Mach*sqrt(Gamma*pr/dr);
#endif

/* Initialize the grid */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

          pGrid->U[k][j][i].d  = dr;
          pGrid->U[k][j][i].M1 = 0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = dr*1e5*1e5/Gamma_1;
	  // + 0.5*dr*(ur*ur);
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = 0.0;
#endif
        }

      }
    }
	
/* enroll gravitational potential of planet & shearing-box potential fns */
	StaticGravPot = PlanetPot;

/* Set IIB value function pointer */
	// if (pDomain->Disp[0] == 0) bvals_mhd_fun(pDomain,left_x1,shk_cloud_iib);

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
	Mp      = par_getd_def("problem","Mplanet",1.0);
	Xplanet = par_getd_def("problem","Xplanet",0.0);
	Yplanet = par_getd_def("problem","Yplanet",0.0);
	Zplanet = par_getd_def("problem","Zplanet",0.0);
	Rsoft   = par_getd_def("problem","Rsoft",0.1);
	
	/* enroll gravitational potential of planet & shearing-box potential fns */	
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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*------------------------------------------------------------------------------
 * PlanetPot:
 */
/*! \fn static Real PlanetPot(const Real x1, const Real x2, const Real x3)
 *  \brief static gravitational potential of planet */
static Real PlanetPot(const Real x1, const Real x2, const Real x3)
{
	Real rad,phi=0.0;

	/* Rather than introducing new variables, 
	   I just multiply the final potential by G*M/R in cgs
	   given M and R in Jupiter units */

	rad = sqrt(SQR(x1-Xplanet) + SQR(x2-Yplanet) + SQR(x3-Zplanet));
	phi = -Gcgs*Mp*Mjup/(rad+Rsoft)/Rjup; //Result in cm^2/s^2
	return phi;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void shk_cloud_iib(GridS *pGrid)
 *  \brief Sets boundary condition on left X boundary (iib) 
 *
 * Note quantities at this boundary are held fixed at the downstream state
 */

void shk_cloud_iib(GridS *pGrid)
{
  int i=0,j=0,k=0;
  int js,je,ks,ke;

  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][i].d  = dl;
        pGrid->U[k][j][i].M1 = ul*dl;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
        pGrid->U[k][j][i].E = pl/Gamma_1
          + 0.5*dl*(ul*ul);
#endif
#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = 0.0;
#endif
      }
    }
  }
}
