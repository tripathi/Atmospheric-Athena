#include "copyright.h"
/*============================================================================*/
/*! \file bvals_ionrad.c
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
 * - bvals_ionrad(DomainS *pD) - calls appropriate functions to set ghost cells
 * - bvals_ionrad_init(MeshS *pM) - sets function pointers used by bvals_ionrad()
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
static const int nghost_ionrad = 1;

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
/*! \fn void bvals_ionrad(DomainS *pD)
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
/*! \fn void bvals_ionrad_init(MeshS *pM)
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI.
 */

void bvals_ionrad_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,nl,nd,irefine;
  int dir = (pM->radplanelist)->dir[0];

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
        if((pD->Disp[0] != 0) || (dir != -1)) {      
          ix1_radBCFun = ProlongateLater;
	  
	} else {                    
	  /* Domain is at L-edge of root Domain */
	  if (dir == -1){
	    ix1_radBCFun = outflow_flux_ix1;
/* 	  } else if (dir ==1){ */
/* 	    ox1_radBCFun = outflow_flux_ox1; */
	  }

	}
      }
      
/*---- ox1 boundary ----------------------------------------------------------*/

      if(ox1_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */
	
	/* Domain boundary is in interior of root */
        if(((pD->Disp[0] + pD->Nx[0])/irefine == pM->Nx[0]) || (dir !=1)) {
          ox1_radBCFun = ProlongateLater;
	  
	  /* Domain is at R-edge of root Domain */
	} else{
	  if (dir == 1){
	    ox1_radBCFun = outflow_flux_ox1;
	  }
/* 	} else { */
/* 	  ox1_radBCFun = outflow_flux_ox1; */
	}
      }
    }

/* Set function pointers for physical boundaries in x2-direction -------------*/

    if(pG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

      if(ix2_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[1] != 0) || (dir != -2)) {
          ix2_radBCFun = ProlongateLater;
	  
	  /* Domain is at L-edge of root Domain */
	} else {
	  if (dir == -2){
	    ix2_radBCFun = outflow_flux_ix2;
	  }
	}
      }

/*---- ox2 boundary ----------------------------------------------------------*/

      if(ox2_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

	/* Domain boundary is in interior of root */
        if(((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) || (dir !=2)) {
          ox2_radBCFun = ProlongateLater;

	  /* Domain is at R-edge of root Domain */
	} else {
	  if (dir == 2){
	    ox2_radBCFun = outflow_flux_ox2;
	  }
	}
      }
    }

/* Set function pointers for physical boundaries in x3-direction -------------*/

    if(pG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/
      
      if(ix3_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */
	
/* Domain boundary is in interior of root */
        if((pD->Disp[2] != 0) || (dir !=-3)) {
          ix3_radBCFun = ProlongateLater;
	  
	  /* Domain is at L-edge of root Domain */
	} else {
	  if (dir == -3){
	    ix3_radBCFun = outflow_flux_ix3;
	  }
	}
      }
/*---- ox3 boundary ----------------------------------------------------------*/

      if(ox3_radBCFun == NULL){    /* BCFun ptr was not set in prob gen */

	/* Domain boundary is in interior of root */
        if(((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) || (dir !=3)) {
          ox3_radBCFun = ProlongateLater;

	  /* Domain is at R-edge of root Domain */
	} else {
	  if (dir == 3){
	    ox3_radBCFun = outflow_flux_ox3;
	  }
	}
      }
    }

  }}}  /* End loop over all Domains with active Grids -----------------------*/


  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   outflow_??? where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   periodic_???
 *   conduct_???
 *   pack_???
 *   unpack_???
 */


/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_flux_ix1(GridS *pGrid)
 *  \brief OUTFLOW boundary condition, Inner x1 boundary (bc_ix1=2) */

static void outflow_flux_ix1(GridS *pGrid)
{
  MeshS *pMesh = pGrid->Mesh;
  int je = pGrid->Nx[1];
  int ke = pGrid->Nx[2];
  
/*   ath_pout(0,"Check: %d vs %d and %d vs %d \n", je, pGrid->je - pGrid->js, ke, pGrid->ke - pGrid->ks); */

/*   int is = pGrid->is; */
/*   int js = pGrid->js, je = pGrid->je; */
/*   int ks = pGrid->ks, ke = pGrid->ke; */
  int i,j,k;
  /*Add in something to see if flux is coming from left.  If not, don't do this.*/
  for (k=0; k<=ke; k++) {
    for (j=0; j<=je; j++) {
/*       for (i=is-1; i<=is + nghost_ionrad; i++) { */
      pGrid->EdgeFlux[k][j][0] =(pMesh->radplanelist)->flux_i;
/*       } */
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_flux_ox1(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2) */

static void outflow_flux_ox1(GridS *pGrid)
{
  MeshS *pMesh = pGrid->Mesh;
  int ie = pGrid->Nx[0];
  int je = pGrid->Nx[1];
  int ke = pGrid->Nx[2];
  int i,j,k;

  for (k=0; k<=ke; k++) {
    for (j=0; j<=je; j++) {
/*       for (i=1; i<=nghost_ionrad; i++) { */
        pGrid->EdgeFlux[k][j][ie] =(pMesh->radplanelist)->flux_i;
/*       } */
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_flux_ix2(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2) */

static void outflow_flux_ix2(GridS *pGrid)
{
  MeshS *pMesh = pGrid->Mesh;
  int ie = pGrid->Nx[0];
  int ke = pGrid->Nx[2];
  int i,j,k;


  for (k=0; k<=ke; k++) {
/*     for (j=1; j<=nghost_ionrad; j++) { */
      for (i=0; i<=ie; i++) {
        pGrid->EdgeFlux[k][0][i] =(pMesh->radplanelist)->flux_i;;
      }
/*     } */
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_flux_ox2(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2) */

static void outflow_flux_ox2(GridS *pGrid)
{
  MeshS *pMesh = pGrid->Mesh;
  int ie = pGrid->Nx[0];
  int je = pGrid->Nx[1];
  int ke = pGrid->Nx[2];
  int i,j,k;

  for (k=0; k<=ke; k++) {
/*     for (j=1; j<=nghost_ionrad; j++) { */
      for (i=0; i<=ie; i++) {
        pGrid->EdgeFlux[k][je][i] =(pMesh->radplanelist)->flux_i;
/*       } */
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_flux_ix3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void outflow_flux_ix3(GridS *pGrid)
{
  MeshS *pMesh = pGrid->Mesh;
  int ie = pGrid->Nx[0];
  int je = pGrid->Nx[1];
  int i,j,k;

/*   for (k=1; k<=nghost_ionrad; k++) { */
    for (j=0; j<=je; j++) {
       for (i=0; i<=ie; i++) { 
        pGrid->EdgeFlux[0][j][i] = (pMesh->radplanelist)->flux_i;
/*       } */
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_flux_ox3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void outflow_flux_ox3(GridS *pGrid)
{
  MeshS *pMesh = pGrid->Mesh;
  int ie = pGrid->Nx[0];
  int je = pGrid->Nx[1];
  int ke = pGrid->Nx[2];
  int i,j,k;

/*   for (k=1; k<=nghost_ionrad; k++) { */
    for (j=0; j<=je; j++) {
      for (i=0; i<=ie; i++) { 
        pGrid->EdgeFlux[ke][j][i] = (pMesh->radplanelist)->flux_i;
/*       } */
    }
  }

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
