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

void ion_radtransfer_init_domain(MeshS *pM) {

  /*Set grid and domain for Athena v4*/

/*   DomainS *pD = &(pM->Domain[0][0]); /\*Temporarily using root domain. Will need FIXing with SMR/MPI*\/ */
/*   GridS *pG = pM->Domain[0][0].Grid; */

  /* Loop over all levels and domains per level */
  DomainS *pD;
  GridS *pG;
  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
	pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */
	pG = pM->Domain[nl][nd].Grid;          /* set ptr to Grid */

  /* Calculate the dimensionality and error check */
  dim = 0;
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;

  switch(dim) {
  case 1: break;
  case 2: break;
  case 3:
    ion_radtransfer_init_domain_3d(pG, pD);
    return;
  }
      }
    }
  }
  ath_error("[ion_radtransfer_init_domain]: Unsupported dim. Nx[0]=%d, Nx[1]=%d, Nx[2]=%d\n",
	    pG->Nx[0],pG->Nx[1],pG->Nx[2]);
}

VDFun_t ion_radtransfer_init(MeshS *pM, int ires){

  /*Set grid and domain for Athena v4*/
/*   DomainS *pD = &(pM->Domain[0][0]); /\*This works only for root domain. Use loop for more domains*\/ */
/*   GridS *pG = pM->Domain[0][0].Grid; */

  DomainS *pD;
  GridS *pG;
  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
	pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */
	pG = pM->Domain[nl][nd].Grid;          /* set ptr to Grid */

  /* Calcualte the dimensionality and error check */
  dim = 0;
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;

  switch(dim){
  case 1: break;
  case 2: break;
  case 3:
    ion_radtransfer_init_3d(pG, pD, ires);
    return ion_radtransfer_3d;
  }
      }
    }
  }

  ath_error("[ion_radtransfer_init]: Unsupported dim. Nx[0]=%d, Nx[1]=%d, Nx[2]=%d\n",
	    pG->Nx[0],pG->Nx[1],pG->Nx[2]);

  /* This is never executed, but lack of a return statement generates
     a warning on some compilers. */
  return NULL;
}

#endif /* ION_RADIATION */
