#include "../copyright.h"

/*==============================================================================
 * FILE: ionrad_smr.c
 *
 * PURPOSE: Contains functions to set values on coarse/fine grids appropriately
 *   when SMR is used
 *
 *   Use of these routines requires that --enable-ion-radiation and --enable-smr
 *   be set at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   ion_prolongate - prolongates coarse grid radiative flux into fine grid
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef STATIC_MESH_REFINEMENT
void ionrad_prolong_rcv(DomainS *pDomain, int dir, int arrsize)
{
  int npg;
  GridS *pGrid = pDomain->Grid;
  int dim;
  MPI_Status stat;
  int err;

  dim = (dir > 0) ? 2*(dir - 1): 2*fabs(dir) - 1;

    
  for (npg=0; npg<(pGrid->NPGrid); npg++)
    {
      pPO=(GridOvrlpS*)&(pGrid->PGrid[npg]);
/*       if (pPO->DomN == pDomain->Level - 1)  */
      err = MPI_Recv(pPO->ionFlx[dim], arrsize, MP_RL, pPO->ID, tag, MPI_COMM_WORLD, &stat);
    }
}

void ionrad_prolong_snd(GridS *pGrid)
{
  pCO=(GridOvrlpS*)&(pGrid->CGrid[ncg]);
  err = MPI_Send(pCO->ionFlx[dim], arrsize, MP_RL, pCO->ID, tag, MPI_COMM_WORLD);

}

void ionrad_prolongate(DomainS *pD)
{
  MeshS *pM = pD->Mesh;
  GridS *pG = pD->Grid;
  GridS *highergrid;
  DomainS *higherdom;
  int nd, ndtot, higherlevel, level;
  int dir = (pM->radplanelist)->dir[0];
  int upperi, upperj, upperk, fixed;
  SideS thisd, parentd, overlapd;
  int isDOverlap;
  int i, j, k, lr;

  level = pD->Level;
  higherlevel = level-1;
  ndtot = pM->DomainsPerLevel[higherlevel];

  /*Find location of this domain */
  /*relative to the origin in the higher level units*/
  for (i=0; i<3; i++){
    thisd.ijkl[i] = pD->Disp[i]/2;
    thisd.ijkr[i] = (pD->Disp[i] + pD->Nx[i])/2;
  }

  /*Find the parent domain for the current domain*/
  if (ndtot == 1) {
    highergrid = pM->Domain[higherlevel][ndtot-1].Grid;
    higherdom =(DomainS*)&(pM->Domain[higherlevel][ndtot-1]);

    /*If this is the first time searching for a parent domain, print to screen the information*/
    if (pM->time == 0) fprintf(stderr, "Dom Level %d No %d has parent Dom Level %d No %d \n", level, pD->DomNumber, pM->Domain[higherlevel][0].Level, pM->Domain[higherlevel][0].DomNumber); 

  } else {
    for (nd=0; nd<ndtot-1; nd++){
      for (i=0; i<3; i++){
	parentd.ijkl[i] = pM->Domain[higherlevel][nd].Disp[i];
	parentd.ijkr[i] = pM->Domain[higherlevel][nd].Disp[i] + pM->Domain[higherlevel][nd].Nx[i];
      }
      isDOverlap = checkOverlap(&thisd, &parentd, &overlapd);
      if (isDOverlap == 1){
	highergrid = pM->Domain[higherlevel][nd].Grid;
	higherdom = (DomainS*)&(pM->Domain[higherlevel][nd]);
	if (pM->time == 0) fprintf(stderr, "Dom Level %d No %d has parent Dom Level %d No %d \n", level, pD->DomNumber, pM->Domain[higherlevel][nd].Level, pM->Domain[higherlevel][nd].DomNumber); 

	/*Check that the system believes this domain has a child grid*/
	if (highergrid->NCGrid < 0) ath_error("Contradiction of child/parent grid \n");
	break;
      }
    }
  }

  if (highergrid == NULL) ath_error("No parent grid \n");

  lr = (dir < 0) ? 1 : -1;
  
  /* switch(dir) { */
  /* case -1: case 1: { */

    if (lr > 0) {
      fixed = 0;
      for (k=0; k<=pG->Nx[2]; k++) {
  	for (j=0; j<=pG->Nx[1]; j++) {
	  /* fprintf(stderr, "%d %d Flux: %e \n", k, j, pG->EdgeFlux[k][j][fixed]); */
  	  upperk =(int)( floor(k/2.) + thisd.ijkl[2] - higherdom->Disp[2]);
	  upperj = (int) (floor(j/2) + thisd.ijkl[1] - higherdom->Disp[1]);
	  /* fprintf(stderr, "This: %d Upper: %d This: %d Upper: %d \n", k, upperk, j, upperj); */
  	  pG->EdgeFlux[k][j][fixed] =highergrid->EdgeFlux[upperk][upperj][thisd.ijkl[0]-higherdom->Disp[0]];
  	}
      }
    } else {ath_error("Doesnt work yet \n");}

  /*   break; */
  /* } */
  /* case -2: case 2: { */
  /*   d2 = 0; */
  /*   if (lr > 0) { */
  /*    } else { */
  /*   } */
  /*   break; */
  /* } */

}

#ifdef MPI_PARALLEL
/* Need to send to all children upon completion of coarse grid and receive by all fine grids*/
/* Need to have call to this not from ionrad_3d.c in the finegrid case.  Instead it needs to be called even after the coarse.*/
/*So maybe like - at end of every level, do a send.  And if not root domain, do a receive*/
/*So if nparents != 0, then receive, and if nchildren !=0 then do a send*/
#endif

#endif

  /* GridOvrlpS *pCO, *pPO; */
  /* int npg, ncg; */
  /* /\* GridsDataS ***Ginfo = pD->GData; *\/ */
  
  /* fprintf(stderr, "---------\n Domain level: %d, number:%d, X-disp:%d, zones:%d, Nparents:%d, Nchild:%d \n", pD->Level, pD->DomNumber, pD->Disp[0], pD->Nx[0], pG->NPGrid,pG->NCGrid); */
  /* /\* fprintf(stderr, "ID %d \n", Ginfo->ID_Comm_world); *\/ */

  /* for (npg=0; npg<(pG->NPGrid); npg++){ */
  /*   pPO = (GridOvrlpS*)&(pG->PGrid[npg]); */
  /*   fprintf(stderr, "Parent: X:ijks %d ijke %d, ID:%d \n", pPO->ijks[0], pPO->ijke[0], pPO->ID); */
  /* } */

  /* for (ncg=0; ncg<(pG->NmyCGrid); ncg++){ */
  /*   pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);    /\* ptr to child Grid overlap *\/ */
  /*   fprintf(stderr, "Child: X:ijks %d ijke %d, ID:%d \n", pCO->ijks[0], pCO->ijke[0], pCO->ID); */
  /* } */
