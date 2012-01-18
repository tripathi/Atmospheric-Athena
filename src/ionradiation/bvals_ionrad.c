#include "copyright.h"
/*============================================================================*/
/*! \file bvals_mhd.c
 *  \brief Sets boundary conditions (quantities in ghost zones) 
 *   for the fluxes on each edge of a Grid.
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) 
 *   for the fluxes on each edge of a Grid.  See comments at start of
 *   bvals_mhd.c for more details.  MPI is NOT implemented here.
 * The only BC functions implemented here are for
 * Outflow boundary conditions
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - bvals_ionrad() - calls appropriate functions to set ghost cells
 * - bvals_ionrad_init() - sets function pointers used by bvals_mhd()
 * - bvals_ionrad_fun() - enrolls a pointer to a user-defined BC function
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef ION_RADIATION

static VGFun_t ix1_radBCFun = NULL, ox1_radBCFun = NULL;
static VGFun_t ix2_radBCFun = NULL, ox2_radBCFun = NULL;
static VGFun_t ix3_radBCFun = NULL, ox3_radBCFun = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   outflow_flux_???()  - outflow BCs at boundary ???
 *============================================================================*/

static void outflow_flux_ix1(GridS *pG);
static void outflow_flux_ox1(GridS *pG);
static void outflow_flux_ix2(GridS *pG);
static void outflow_flux_ox2(GridS *pG);
static void outflow_flux_ix3(GridS *pG);
static void outflow_flux_ox3(GridS *pG);

static void ProlongateLater(GridS *pG);


/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mhd(DomainS *pD)
 *  \brief Calls appropriate functions to set ghost zones.  
 *
 *   The function
 *   pointers (*(pD->???_BCFun)) are set by bvals_init() to be either a
 *   user-defined function, or one of the functions corresponding to reflecting,
 *   periodic, or outflow.  If the left- or right-Grid ID numbers are >= 1
 *   (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void bvals_ionrad(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx[0] > 1){

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*(ix1_radBCFun))(pGrid);
      (*(ox1_radBCFun))(pGrid);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*(ix2_radBCFun))(pGrid);
      (*(ox2_radBCFun))(pGrid);
    } 


  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){


/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*(ix3_radBCFun))(pGrid);
      (*(ox3_radBCFun))(pGrid);
    } 

  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mhd_init(MeshS *pM)
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI.
 */

void bvals_ionrad_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,nl,nd,irefine;

/* Cycle through all the Domains that have active Grids on this proc */

  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
    pD = (DomainS*)&(pM->Domain[nl][nd]);  /* ptr to Domain */
    pG = pM->Domain[nl][nd].Grid;          /* ptr to Grid */
    irefine = 1;
    for (i=1;i<=nl;i++) irefine *= 2;   /* C pow fn only takes doubles !! */

/* Set function pointers for physical boundaries in x1-direction -------------*/

    if(pG->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

      if(ix1_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[0] != 0) {      
          ix1_radBCFun = ProlongateLater;
	  
/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) && 
               pM->BCFlag_ix1 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix1b but not ox1b and periodic BC not allowed\n",nl); 
	    
/* Domain is at L-edge of root Domain */
          } else {                    
            switch(pM->BCFlag_ix1){

            case 1: /* Reflecting, B_normal=0 */
              ix1_radBCFun = reflect_ix1;
	      break;
	      
            case 2: /* Outflow */
              ix1_radBCFun = outflow_ix1;
	      break;
	      
            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              ix1_radBCFun = periodic_ix1;
	      break;
	      
            case 5: /* Reflecting, B_normal!=0 */
              ix1_radBCFun = reflect_ix1;
	      break;
	      
            default:
              ath_perr(-1,"[bvals_init]:bc_ix1=%d unknown\n",pM->BCFlag_ix1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox1 boundary ----------------------------------------------------------*/

      if(ox1_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
          ox1_radBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[0] != 0) && (pM->BCFlag_ox1 == 4)) {      
            ath_error("[bvals_init]:level=%d Domain touching ox1b but not ix1b and periodic BC not allowed\n",nl); 


/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox1){

            case 1: /* Reflecting, B_normal=0 */
              ox1_radBCFun = reflect_ox1;
            break;

            case 2: /* Outflow */
              ox1_radBCFun = outflow_ox1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              ox1_radBCFun = periodic_ox1;

            break;

            case 5: /* Reflecting, B_normal!=0 */
              ox1_radBCFun = reflect_ox1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox1=%d unknown\n",pM->BCFlag_ox1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Set function pointers for physical boundaries in x2-direction -------------*/

    if(pG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

      if(ix2_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[1] != 0) {
          ix2_radBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) &&
               pM->BCFlag_ix2 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix2b but not ox2b and periodic BC not allowed\n",nl); 


/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix2){

            case 1: /* Reflecting, B_normal=0 */
              ix2_radBCFun = reflect_ix2;
            break;

            case 2: /* Outflow */
              ix2_radBCFun = outflow_ix2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              ix2_radBCFun = periodic_ix2;

            break;
  
            case 5: /* Reflecting, B_normal!=0 */
              ix2_radBCFun = reflect_ix2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix2=%d unknown\n",pM->BCFlag_ix2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox2 boundary ----------------------------------------------------------*/

      if(ox2_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
          ox2_radBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[1] != 0) && (pM->BCFlag_ox2 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox2b but not ix2b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox2){

            case 1: /* Reflecting, B_normal=0 */
              ox2_radBCFun = reflect_ox2;
            break;

            case 2: /* Outflow */
              ox2_radBCFun = outflow_ox2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              ox2_radBCFun = periodic_ox2;

            break;

            case 5: /* Reflecting, B_normal!=0 */
              ox2_radBCFun = reflect_ox2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox2=%d unknown\n",pM->BCFlag_ox2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Set function pointers for physical boundaries in x3-direction -------------*/

    if(pG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

      if(ix3_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[2] != 0) {
          ix3_radBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) &&
               pM->BCFlag_ix3 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix3b but not ox3b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix3){

            case 1: /* Reflecting, B_normal=0 */
              ix3_radBCFun = reflect_ix3;
            break;

            case 2: /* Outflow */
              ix3_radBCFun = outflow_ix3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              ix3_radBCFun = periodic_ix3;

            break;

            case 5: /* Reflecting, B_normal!=0 */
              ix3_radBCFun = reflect_ix3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix3=%d unknown\n",pM->BCFlag_ix3);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox3 boundary ----------------------------------------------------------*/

      if(ox3_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
          ox3_radBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[2] != 0) && (pM->BCFlag_ox3 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox3b but not ix3b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox3){

            case 1: /* Reflecting, B_normal=0 */
              ox3_radBCFun = reflect_ox3;
            break;

            case 2: /* Outflow */
              ox3_radBCFun = outflow_ox3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              ox3_radBCFun = periodic_ox3;

            break;

            case 5: /* Reflecting, B_normal!=0 */
              ox3_radBCFun = reflect_ox3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox3=%d unknown\n",pM->BCFlag_ox3);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Figure out largest size needed for send/receive buffers with MPI ----------*/


  }}}  /* End loop over all Domains with active Grids -----------------------*/


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mhd_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
 *  \brief Sets function ptrs for user-defined BCs.
 */

void bvals_ionrad_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    ix1_radBCFun = prob_bc;
    break;
  case right_x1:
    ox1_radBCFun = prob_bc;
    break;
  case left_x2:
    ix2_radBCFun = prob_bc;
    break;
  case right_x2:
    ox2_radBCFun = prob_bc;
    break;
  case left_x3:
    ix3_radBCFun = prob_bc;
    break;
  case right_x3:
    ox3_radBCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[bvals_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???:   where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   outflow_???
 *   periodic_???
 *   conduct_???
 *   pack_???
 *   unpack_???
 */

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix1(GridS *pGrid)
 *  \brief  REFLECTING boundary conditions, Inner x1 boundary (bc_ix1=1) */

static void reflect_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i]    =  pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][is-i].B1c= -pGrid->U[k][j][is-i].B1c;/* reflect 1-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* reflect normal component of B field, B1i not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][is] = 0.0;
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i]   = -pGrid->B1i[k][j][is+i];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox1(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (bc_ox1=1). */

static void reflect_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i]    =  pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].M1 = -pGrid->U[k][j][ie+i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][ie+i].B1c= -pGrid->U[k][j][ie+i].B1c;/* reflect 1-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* reflect normal component of B field */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] = 0.0;
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = -pGrid->B1i[k][j][ie-(i-2)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (bc_ix2=1) */

static void reflect_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][js-j][i].B2c= -pGrid->U[k][js-j][i].B2c;/* reflect 2-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* reflect normal component of B field, B2i not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B2i[k][js][i] = 0.0;
    }
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = -pGrid->B2i[k][js+j][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (bc_ox2=1) */

static void reflect_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][je+j][i].B2c= -pGrid->U[k][je+j][i].B2c;/* reflect 2-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* reflect normal component of B field */
  for (k=ks; k<=ke; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B2i[k][je+1][i] = 0.0;
    }
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = -pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (bc_ix3=1) */

static void reflect_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ks-k][j][i].B3c= -pGrid->U[ks-k][j][i].B3c;/* reflect 3-fld.*/
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* reflect normal component of B field, B3i not set at k=ks-nghost */
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B3i[ks][j][i] = 0.0;
    }
  }
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i]   = -pGrid->B3i[ks+k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (bc_ox3=1) */

static void reflect_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ke+k][j][i].B3c= -pGrid->U[ke+k][j][i].B3c;/* reflect 3-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* reflect normal component of B field */
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B3i[ke+1][j][i] = 0.0;
    }
  }
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i]   = -pGrid->B3i[ke-(k-2)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix1(GridS *pGrid)
 *  \brief OUTFLOW boundary condition, Inner x1 boundary (bc_ix1=2) */

static void outflow_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox1(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2) */

static void outflow_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix2(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2) */

static void outflow_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][js][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox2(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2) */

static void outflow_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][je][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je][i];
      }
    }
  }

/* j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void outflow_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ks][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void outflow_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ke][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke][j][i];
      }
    }
  }

/* k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4) */

static void periodic_ix1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][ie-(i-1)];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void periodic_ox1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][is+(i-1)];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=ie+1 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void periodic_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][je-(j-1)][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][je-(j-1)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void periodic_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][js+(j-1)][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* B2i is not set at j=je+1 */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][js+(j-1)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void periodic_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ke-(k-1)][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void periodic_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ks+(k-1)][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* B3i is not set at k=ke+1 */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix1(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Inner x1 boundary (bc_ix1=5) */

static void conduct_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i]    =  pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][is-i].B2c= -pGrid->U[k][j][is-i].B2c;/* reflect fld */
        pGrid->U[k][j][is-i].B3c= -pGrid->U[k][j][is-i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is+i];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i]   = -pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i]   = -pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox1(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Outer x1 boundary (bc_ox1=5) */

static void conduct_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i]    =  pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].M1 = -pGrid->U[k][j][ie+i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][ie+i].B2c= -pGrid->U[k][j][ie+i].B2c;/* reflect fld */
        pGrid->U[k][j][ie+i].B3c= -pGrid->U[k][j][ie+i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not set for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie-(i-2)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i]   = -pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i]   = -pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix2(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Inner x2 boundary (bc_ix2=5) */

static void conduct_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][js-j][i].B1c= -pGrid->U[k][js-j][i].B1c;/* reflect fld */
        pGrid->U[k][js-j][i].B3c= -pGrid->U[k][js-j][i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = -pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* B2i not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js+j][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = -pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox2(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Outer x2 boundary (bc_ox2=5) */

static void conduct_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][je+j][i].B1c= -pGrid->U[k][je+j][i].B1c;/* reflect fld */
        pGrid->U[k][je+j][i].B3c= -pGrid->U[k][je+j][i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i]   = -pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* j=je+1 is not set for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i]   = pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i]   = -pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix3(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Inner x3 boundary (bc_ix3=5) */

static void conduct_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ks-k][j][i].B1c= -pGrid->U[ks-k][j][i].B1c;/* reflect fld */
        pGrid->U[ks-k][j][i].B2c= -pGrid->U[ks-k][j][i].B2c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i]   = -pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i]   = -pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks+k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox3(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Outer x3 boundary (bc_ox3=5) */

static void conduct_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ke+k][j][i].B1c= -pGrid->U[ke+k][j][i].B1c;/* reflect fld */
        pGrid->U[ke+k][j][i].B2c= -pGrid->U[ke+k][j][i].B2c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i]   = -pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i]   = -pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* k=ke+1 is not set for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i]   = pGrid->B3i[ke-(k-2)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void ProlongateLater(GridS *pGrid)
 *  \brief PROLONGATION boundary conditions.  
 *
 *  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pGrid)
{
  return;
}

#endif /* ION_RADIATION */
