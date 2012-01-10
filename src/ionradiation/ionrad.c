#include "../copyright.h"
#define IONRAD_C
/*==============================================================================
 * FILE: ionrad.c
 *
 * PURPOSE: Contains function to control ionizing radiative transfer routines.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   ion_radtransfer_init() - sets pointer to appropriate ionizing
 *                            radiative transfer function
 *   ion_radtransfer_init_domain() - sets domain information for ionizing
 *                                   radiative transfer module
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "ionrad.h"

#ifdef ION_RADIATION
static int dim=0;

void ion_radtransfer_init_domain(DomainS *pD) {

  /*Set grid - new to Athena v4*/
  GridS *pG = (pD->Grid);

  /* Calculate the dimensionality and error check */
  dim = 0;
  if(pG->N[0] > 1) dim++;
  if(pG->N[1] > 1) dim++;
  if(pG->N[2] > 1) dim++;

  switch(dim) {
  case 1: break;
  case 2: break;
  case 3:
    ion_radtransfer_init_domain_3d(pG, pD);
    return;
  }

  ath_error("[ion_radtransfer_init_domain]: Unsupported dim. N[0]=%d, N[1]=%d, N[2]=%d\n",
	    pG->N[0],pG->N[1],pG->N[2]);
}

VGFun_t ion_radtransfer_init(DomainS *pD, int ires){

  /*Set grid - new to Athena v4*/
  GridS *pG = (pD->Grid);

  /* Calcualte the dimensionality and error check */
  dim = 0;
  if(pG->N[0] > 1) dim++;
  if(pG->N[1] > 1) dim++;
  if(pG->N[2] > 1) dim++;

  switch(dim){
  case 1: break;
  case 2: break;
  case 3:
    ion_radtransfer_init_3d(pG, pD, ires);
    return ion_radtransfer_3d;
  }

  ath_error("[ion_radtransfer_init]: Unsupported dim. N[0]=%d, N[1]=%d, N[2]=%d\n",
	    pG->N[0],pG->N[1],pG->N[2]);

  /* This is never executed, but lack of a return statement generates
     a warning on some compilers. */
  return NULL;
}

#endif /* ION_RADIATION */
