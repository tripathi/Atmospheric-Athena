#include "copyright.h"
/*==============================================================================
 * FILE: hot_jupiter.c
 *
 * PURPOSE: Problem generator for shock-cloud problem; planar shock impacting
 *   a single spherical cloud.  Input parameters are:
 *      problem/Mach   = Mach number of incident shock
 *      problem/drat   = density ratio of cloud to ambient
 *      problem/beta   = ratio of Pgas/Pmag
 *      problem/iprob  = integer flag to determine problem
 *   The cloud radius is fixed at 1.0.  The center of the coordinate system
 *   defines the center of the cloud, and should be in the middle of the cloud.
 *   The shock is initially at x1=-2.0.  A typical grid domain should span
 *   x1 in [-3.0,7.0] , y and z in [-2.5,2.5] (see input file in /tst)
 *   Various test cases are possible:
 *     (iprob=1): B parallel to shock normal
 *     (iprob=2): B perpendicular to shock normal -- NOT YET IMPLEMENTED
 *   If the code is configured with nscalars>0, the cloud material is labeled
 *   with U[k][j][i].s[0]=1.
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* postshock flow variables are shared with IIB function */

static Real Mp,Rsoft,Xplanet,Yplanet,Zplanet;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);
static Real dl,pl,ul;
const Real kb=1.38e-16;
const Real mH = 1.67e-24;
const Real Ggrav = 6.67e-8;
#ifdef MHD
static Real bxl,byl,bzl;

#endif /* MHD */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * shk_cloud_iib() - fixes BCs on L-x1 (left edge) of grid to postshock flow.
 *============================================================================*/

void shk_cloud_iib(GridS *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;

  Real x1,x2,x3,diag,rad,xshock;
  Real Mach,drat;
#ifdef MHD
  Real beta,bxr,byr,bzr;
#endif /* MHD */
  Real dr,pr,ur;
  Real rho_at,p_at,T_at,rho_c,p_c,Rp,Hor2,c_s,rho_H,Rb;
  Real prprev,tempval;
/* Read input parameters */

  Mach = par_getd("problem","Mach");
  drat = par_getd("problem","drat");
  
  rho_at = par_getd("problem","rho_at");
  T_at =par_getd("problem","T_at");
  Rp =par_getd("problem","Rp");
  Rb = par_getd("problem","Rbound");
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);

#ifdef MHD
  beta = par_getd("problem","beta");
#endif
  
/* Set paramters in ambient medium ("R-state" for shock) */

  dr = 1.0;
  pr = 1.0/Gamma;
  ur = 0.0;
#ifdef INNERB
  fprintf(stderr,"inner boundary is on!\n");
#endif

  /* set up parameters for central part of "planet" center is colder and denser?*/
  Hor2=kb*T_at/(2*mH*Ggrav*Mp);
  rho_c= rho_at*exp(-1/Hor2*(1/Rp - 1/Rb));
		    
  c_s = sqrt(kb*T_at/(2*mH));
  p_c = pow(rho_c,Gamma)*(c_s)*(c_s);
  fprintf(stderr,"cs is %e, and rho_c is %e and my cell is %e\n",c_s,rho_c, pGrid->dx1);

/* Initialize the grid */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        diag = sqrt(x1*x1 + x2*x2 + x3*x3);

/* outside the planet */
       
	pGrid->U[k][j][i].d  = rho_at/drat;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	
#ifdef ADIABATIC
	pGrid->U[k][j][i].E = pow(rho_at,Gamma)*(c_s*c_s)/Gamma_1;  
	  

#endif
#if (NSCALARS > 0)
	pGrid->U[k][j][i].s[0] = 0.0;
#endif
     
	
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
          pGrid->U[k][j][i].d  = rho_c*exp((1/diag-1/Rb)/Hor2); 
          pGrid->U[k][j][i].M1 =0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
	  prprev=pr;
	  pr = pow(pGrid->U[k][j][i].d,Gamma)*(c_s*c_s);
	  tempval = fmax(rho_at/drat*10,tempval);
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
  
  fprintf(stderr,"%e is sheet rho\n",tempval);


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
  Real rho_at,p_at,T_at,rho_c,p_c,Rp,Hor2,c_s,rho_H,Rb, pr;
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
	      
	      
	      if (diag < Rb ){
		pGrid->U[k][j][i].d  = rho_c*0.99;
		pGrid->U[k][j][i].M1 = 0.0;
		pGrid->U[k][j][i].M2 = 0.0;
		pGrid->U[k][j][i].M3 = 0.0;
		pr = pow(rho_c*0.99,Gamma)*(c_s*c_s);
#ifdef ADIABATIC
		pGrid->U[k][j][i].E = pr/Gamma_1;
#endif
		
	      } else if (diag <= Rp){
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

/*-----------------------------------------------------------------------------
 * shk_cloud_iib: sets boundary condition on left X boundary (iib)
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
#ifdef MHD
        pGrid->B1i[k][j][i] = bxl;
        pGrid->B2i[k][j][i] = byl;
        pGrid->B3i[k][j][i] = bzl;
        pGrid->U[k][j][i].B1c = bxl;
        pGrid->U[k][j][i].B2c = byl;
        pGrid->U[k][j][i].B3c = bzl;
#endif
#ifdef ADIABATIC
        pGrid->U[k][j][i].E = pl/Gamma_1
#ifdef MHD
          + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
          + 0.5*dl*(ul*ul);
#endif
#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = 0.0;
#endif
      }
    }
  }
}
