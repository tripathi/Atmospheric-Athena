f#include "../copyright.h"

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
#include "ionrad.h"

#ifdef STATIC_MESH_REFINEMENT

int tag = 0;
void ionrad_prolong_rcv(GridS *pGrid, int dim)
{
  int npg;
  int i, j, k, fixed, arrsize, indexarith;
  MPI_Status stat;
  int err;
  GridOvrlpS *pPO;

  for (npg=0; npg<(pGrid->NPGrid); npg++)
    {
      pPO=(GridOvrlpS*)&(pGrid->PGrid[npg]);

/*       switch(dim) { */
/*       case 0: case 1: { */
	if (fmod(dim,2) == 0) {
	  fixed = pPO->ijks[0] - nghost;
	} else {
	  fixed = pPO->ijke[0] + 1 - nghost; /*Should this be +1 or +2? */
	}
	
	if(pPO->ionFlx[dim] != NULL) {
	  arrsize = (pPO->ijke[2] + 2 - pPO->ijks[2]) * (pPO->ijke[1] + 2 - pPO->ijks[1]);
	  /* Will need in nonblocking form? */
	  err = MPI_Recv(pPO->ionFlx[dim], arrsize, MP_RL, pPO->ID, pPO->ID, MPI_COMM_WORLD, &stat);


	  for (k=pPO->ijks[2] - nghost; k<= pPO->ijke[2]+1 - nghost; k++) {
	    for (j=pPO->ijks[1] - nghost; j<= pPO->ijke[1]+1 - nghost; j++) {
	      indexarith = (k-(pPO->ijks[2]-nghost))*(pPO->ijke[1] - pPO->ijks[1] + 2)+j-(pPO->ijks[1]-nghost);
	      pGrid->EdgeFlux[k][j][fixed] = pPO->ionFlx[dim][indexarith];
	      fprintf(stderr, "k:%d j:%d, index: %d \n", k, j, indexarith);
	    }
	  }
	  /*Will need to check indexing to see if it's +1 or +2*/
	}
/* 	break; */
/*       } */
/*       } */

/*     case 2: case 3: { */
/*       if (lr > 0) { */
/* 	fixed = pCO->ijks[1] - nghost; */
/*       } else { */
/* 	fixed = pCO->ijke[1] + 1 - nghost; */
/*       } */
      
/*       if(pCO->ionFlx[dim] != NULL) { */
/* 	  for (k=pCO->ijks[2] - nghost; k<= pCO->ijke[2]+1 - nghost; k++) { */
/* 	    for (i=pCO->ijks[0] - nghost; i<= pCO->ijke[0]+1 - nghost; i++) { */
/* 	      pCO->ionFlx[dim][(k-(pCO->ijks[2]-nghost))*(pCO->ijke[0] - pCO->ijks[0] + 2)+i-(pCO->ijks[0]-nghost)] = pGrid->EdgeFlux[k][fixed][i]; */
/* 	      fprintf(stderr, "k:%d i:%d, index: %d, nx: %d \n", k, i, (k-(pCO->ijks[2]-nghost))*(pCO->ijke[0] - pCO->ijks[0] + 2)+i-(pCO->ijks[0]-nghost), pCO->ijke[0] - pCO->ijks[0] +2); */
/* 	    } */
/* 	  } */
/* 	  arrsize = (pCO->ijke[2] + 2 - pCO->ijks[2]) * (pCO->ijke[0] + 2 - pCO->ijks[0]); */
/*       } */
/*       break; */
/*     } */

/*     case 4: case 5: { */
/*       if (lr > 0) { */
/* 	fixed = pCO->ijks[2] - nghost; */
/*       } else { */
/* 	fixed = pCO->ijke[2] + 1 - nghost; */
/*       } */
/*       if(pCO->ionFlx[dim] != NULL) { */
/* 	for (j=pCO->ijks[1] - nghost; j<= pCO->ijke[1]+1 - nghost; j++) { */
/* 	  for (i=pCO->ijks[0] - nghost; i<= pCO->ijke[0]+1 - nghost; i++) { */
/* 	    pCO->ionFlx[dim][(j-(pCO->ijks[1]-nghost))*(pCO->ijke[0] - pCO->ijks[0] + 2)+i-(pCO->ijks[0]-nghost)] = pGrid->EdgeFlux[fixed][j][i]; */
/* 	    fprintf(stderr, "j:%d i:%d, index: %d, nx: %d \n", j, i, (j-(pCO->ijks[1]-nghost))*(pCO->ijke[0] - pCO->ijks[0] + 2)+i-(pCO->ijks[0]-nghost), pCO->ijke[0] - pCO->ijks[0] +2); */
/* 	  } */
/* 	} */
/* 	arrsize = (pCO->ijke[0] + 2 - pCO->ijks[0]) * (pCO->ijke[1] + 2 - pCO->ijks[1]); */



    }


}

void ionrad_prolong_snd(GridS *pGrid, int dim, int arrsize)
{
  MPI_Request *send_rq=NULL;
  int ncg, ierr;
  GridOvrlpS *pCO, *pPO;

  for (ncg=0; ncg<(pGrid->NCGrid); ncg++){
    pCO=(GridOvrlpS*)&(pGrid->CGrid[ncg]);
    ierr = MPI_Isend(pCO->ionFlx[dim], arrsize, MP_RL, pCO->ID, pCO->ID, MPI_COMM_WORLD, &(send_rq[ncg]));
    fprintf(stderr, "I sent my data for %d \n", ncg);
  }
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
