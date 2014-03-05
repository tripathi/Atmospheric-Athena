#include "../copyright.h"

/*==============================================================================
 * FILE: ionrad_smr.c
 *
 * PURPOSE: Contains functions to set ionizing flux values on fine grids 
 *   when SMR is used
 *
 *   Use of these routines requires that --enable-ion-radiation and --enable-smr
 *   be set at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   ionrad_prolong_rcv - Receives radiative flux from coarser grid
 *   ionrad_prolong_snd - Sends radiative flux to finer grid
 *   ion_prolongate - prolongates coarse grid radiative flux into fine grid - DEPRECATED
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "ionrad.h"
#include <math.h>

#define RESET_COLOR "\e[m"
#define MAKE_BLUE "\e[34m"
#define MAKE_RED "\e[31m"


#ifdef STATIC_MESH_REFINEMENT

void ionrad_prolong_rcv(GridS *pGrid, int dim, int level, int domnumber)
{

  MeshS *pMesh = pGrid->Mesh;
  int npg;
  int i, j, k, fixed, arrsize, indexarith;
  int is, js, ks;
  int err, ierr;
  GridOvrlpS *pPO, *pCO;

  DomainS *pDomain;
  pDomain = (DomainS*)&(pMesh->Domain[level][domnumber]);  /* ptr to Domain */

#ifdef MPI_PARALLEL
  int allrcv =0;
  int rcvd = 0;
  int receive_done;
  int *succtest;
  int tag1, tag2, tag3;
  char temp[10];
  MPI_Status stat;
  MPI_Request *rcv_rq;
  int *rcvdarrsize;
  int dy, dz, coarsek;

  rcv_rq = (MPI_Request*) calloc_1d_array(pGrid->NPGrid,sizeof(MPI_Request));
  succtest = (int*) calloc_1d_array(pGrid->NPGrid,sizeof(int)) ;
  rcvdarrsize = (int*) calloc_1d_array(pGrid->NPGrid,sizeof(int)) ;
#else
  GridS *parentgrid; /*Parent grid of current grid*/
#endif


  /* fprintf (stderr, "Level: %d Grid disp : %d \n", level, pGrid->Disp[0]); */
/*Find my parent grid overlap structure(s if MPI) to receive data from it*/
  for (npg=0; npg<(pGrid->NPGrid); npg++)
    {
      pPO=(GridOvrlpS*)&(pGrid->PGrid[npg]);

      /*#ifndef MPI_PARALLEL*/
      /*      pPO->ionFlx[dim] = pCO->ionFlx[dim]*/

#ifdef MPI_PARALLEL
      if(pPO->ionFlx[dim] != NULL) {
	/* fprintf(stderr, "Beginning receive call for domain level %d \n", level); */
	
	/*Old tagging strategies*/
	/* tag1 = atoi(temp); */
	/* tag2 = (level - 1)*1000000 + pPO->ID * 10000 + level * 100 + domnumber; */
	/* fprintf(stderr, "rcvconcat: %d, powers:%d domno: %d, myid :%d \n", tag1, tag2, domnumber, myID_Comm_world); */
	tag3 = pPO->DomN + 100;

	
	/*AT 9/21/12: TO_DO: Insert case statement, so that arrsize is pulled from the correct indices*/
	
	/*Find the size of the array of flux values being transferred*/
	/* arrsize = (pPO->ijke[2] + 2 - pPO->ijks[2]) * (pPO->ijke[1] + 2 - pPO->ijks[1]); */
	dz = floor((pPO->ijke[2] - pPO->ijks[2])/2.)+2.;
	dy = floor((pPO->ijke[1] - pPO->ijks[1])/2.)+2.;

	rcvdarrsize[npg] = dy * dz;
	/*(pPO->ijke[2] - pPO->ijks[2]+3)/2. * (pPO->ijke[1] - pPO->ijks[1]+3)/2.;*/




	/* fprintf(stderr, MAKE_BLUE "ARR SIZE %d being received k: s: %d e: %d j: s: %d e: %d\n" RESET_COLOR, rcvdarrsize[npg], pPO->ijks[2], pPO->ijke[2], pPO->ijks[1], pPO->ijke[1]); */
	/* arrsize = (pCO->ijke[2] + 2 - pCO->ijks[2]) * (pCO->ijke[1] + 2 - pCO->ijks[1]); */


	/*Initiate non-blocking receive*/
	ierr = MPI_Irecv(pPO->ionFlx[dim], rcvdarrsize[npg], MP_RL, pPO->ID, tag3, pMesh->Domain[level][domnumber].Comm_Parent, &(rcv_rq[npg]));
	/* fprintf (stderr, "did i get here? %d \n", ierr); */

      } else {
	/*AT 9/26/12: Faking a successful test for grids that are not in the direction of propagation*/
	succtest[npg] = 1;
	rcvd++;
      }
    }

/*Loop over all possible communications to check which parent grids have been received from*/
  while(allrcv==0)
    {
      for (npg=0; npg<(pGrid->NPGrid); npg++)
	{
	  if (succtest[npg]==0) 
	    {
	      err = MPI_Test(&(rcv_rq[npg]), &receive_done, &stat);
	      if (receive_done) 
		{ 
		  succtest[npg]= 1; 
		  rcvd ++;
		  /* fprintf(stderr, " Data received from parent w/id: %d for tag %d. I'm on level %d \n", pGrid->PGrid[npg].ID, pGrid->PGrid[npg].DomN + 100, level); */
		  
		} 
	      else {
			/*fprintf(stderr, "Still waiting \n");*/
	      }
	    }
	}
      allrcv = (rcvd == pGrid->NPGrid) ? 1 : 0;
    }


/* Populate the cells on this (fine) grid with the values received from the coarse grid */
/*AT 2/28/13: TO_DO: Indexing needs to be fixed*/

  for (npg=0; npg<(pGrid->NPGrid); npg++)
    {
      pPO=(GridOvrlpS*)&(pGrid->PGrid[npg]);

      if(pPO->ionFlx[dim] != NULL) {
	/* fprintf(stderr, "Going to go fill edgeflux now \n"); */

	switch(dim) { 
	case 0: case 1: { 
	  if (fmod(dim,2) == 0) {
	    fixed = (pPO->ijks[0] - nghost) * 2.;
	  } else {
	    fixed = (pPO->ijke[0] + 1 - nghost) * 2.;
	  }  
  /* fprintf(stderr,"I'm here at line 150 \n"); */



	  for (k=pPO->ijks[2] - nghost; k<= pPO->ijke[2]+1 - nghost; k++) {
	    for (j=pPO->ijks[1] - nghost; j<= pPO->ijke[1]+1 - nghost; j++) {
	      coarsek = floor(k/2) + pDomain->Disp[2];

	      indexarith = floor((k - pPO->ijks[2] + nghost)/2.) * (floor((pPO->ijke[1] - pPO->ijks[1])/2.) + 2) + floor((j - pPO->ijks[1] + nghost)/2.);
	      
	      /* indexarith = floor(( (k-(pPO->ijks[2]-nghost))*(pPO->ijke[1] - pPO->ijks[1] + 2)+j-(pPO->ijks[1]-nghost))/2.); */
	      /* indexarith = (k-(pCO->ijks[2]-nghost))*(pCO->ijke[1] - pCO->ijks[1] + 2)+j-(pCO->ijks[1]-nghost); */



	      /* fprintf(stderr, "I'm here at line 155 k:%d j:%d \n", k, j); */

	      ks = k;

/* *2 - pDomain->Disp[2]; */

	      js = j;
/* *2 - pDomain->Disp[1]; */

	      /* fprintf(stderr, MAKE_RED"Now here at line 162. ks: %d, js:%d, indexarith: %d \n"RESET_COLOR, ks, js, indexarith); */
	      pGrid->EdgeFlux[ks][js][fixed] = pPO->ionFlx[dim][indexarith];
	      /* fprintf(stderr, "Now here at line 164 \n"); */
	      /* if (j < (pPO->ijke[1]+1 - nghost)) { */
	      /* fprintf(stderr, "Now here at line 166 \n"); */
	      /* 	if (k < (pPO->ijke[2]+1 - nghost)) {	 */
	      /* 	  /\*pPO is very likely the wrong place to store this since we really want hte parent's pCO*\/ */
	      /* 	  pGrid->EdgeFlux[ks+1][js+1][fixed] = pPO->ionFlx[dim][indexarith]; */
	      /* 	  pGrid->EdgeFlux[ks][js+1][fixed] = pPO->ionFlx[dim][indexarith]; */
	      /* 	  pGrid->EdgeFlux[ks+1][js][fixed] = pPO->ionFlx[dim][indexarith]; */
	      /* 	} */
		
	      /* 	/\*Handle edge cases*\/ */
	      /* 	else { */
	      /* 	  pGrid->EdgeFlux[ks][js+1][fixed] = pCO->ionFlx[dim][indexarith]; */

	      /* 	} */
	      /* 	/\*Handle edge cases*\/ */
	      /* }else { */
	      /* 	if (k < pCO->ijke[2]+1 - nghost ) { */
	      /* 	  pGrid->EdgeFlux[ks+1][js][fixed] = pCO->ionFlx[dim][indexarith]; */
		  
	      /* /\*Will need to check indexing to see if it's +1 or +2*\/ */
	      /* /\*	      fprintf(stderr, "Putting away received k:%d j:%d, index: %d \n", k, j, indexarith); *\/ */
	      /* 	} */
	      /* } */

  /* fprintf(stderr,"I'm here at line 185 \n"); */

	    }
	  }
	  break; 
	}

	case 2: case 3: {
	  if (fmod(dim,2) == 0) {
	    fixed = pPO->ijks[1] - nghost;
	  } else {
	    fixed = pPO->ijke[1] + 1 - nghost;
	  }	  
	  for (k=pPO->ijks[2] - nghost; k<= pPO->ijke[2]+1 - nghost; k++) {
	    for (i=pPO->ijks[0] - nghost; i<= pPO->ijke[0]+1 - nghost; i++) {
	      indexarith = (k-(pPO->ijks[2]-nghost))*(pPO->ijke[0] - pPO->ijks[0] + 2)+i-(pPO->ijks[0]-nghost);
	      pGrid->EdgeFlux[k][fixed][i] = pPO->ionFlx[dim][indexarith];
	    }
	  }
	  break;
	}

	case 4: case 5: {
	  if (fmod(dim,2) == 0) {
	    fixed = pPO->ijks[2] - nghost;
	  } else {
	    fixed = pPO->ijke[2] + 1 - nghost;
	  }
	  for (j=pPO->ijks[1] - nghost; j<= pPO->ijke[1]+1 - nghost; j++) {
	    for (i=pPO->ijks[0] - nghost; i<= pPO->ijke[0]+1 - nghost; i++) {
	      indexarith = (j-(pPO->ijks[1]-nghost))*(pPO->ijke[0] - pPO->ijks[0] + 2)+i-(pPO->ijks[0]-nghost);
	      pGrid->EdgeFlux[fixed][k][i] = pPO->ionFlx[dim][indexarith];
	    }
	  }
	  break;
	}
	}
      }
#else /*single processor*/

/*Let the parent grid be the grid whose domain is one level higher and has domain number matching pPO->DomN*/
      parentgrid = pMesh->Domain[level-1][pPO->DomN].Grid;
/*Find which child of my parent I am.  Assign that grid overlap structure to pCO*/
      for (i = 0; i <= parentgrid->NCGrid; i++) {
	pCO=(GridOvrlpS*)&(parentgrid->CGrid[i]);
	if (pCO->DomN == domnumber) {
	  break;
	}
      }



/*Use the parent grid's child grid overlap ionFlx array to fill this grid's edgeflux array*/


      /*AT 1/11/13: Changing the if and else conditions to be reversed of what they originally were so that dimensions match*/
      /* AT 12/23/12 TO_DO Add case statement back in when 1 direction works to account for other directions of propagation*/
      /* switch(dim) { */
      /* case 0: case 1: { */

/*Find the face where radiation is incoming from the coarse grid*/
      /*If in the +x direction*/
      if (fmod(dim,2) == 0) {
	fixed = (pPO->ijks[0] - nghost) * 2.;
      }
      /*If in the -x direction*/	
      else {
	fixed = (pPO->ijke[0] + 1 - nghost) * 2.;
      }
      
/*Loop over all cells orthogonal to the direction of propagation */
/*Loop indices range over all the cells, as indexed by the coarse grid*/
      if(pCO->ionFlx[dim] != NULL) {
	/* fprintf(stderr, "j Start: %d End: %d   k Start: %d End:%d \n", pCO->ijks[1] - nghost, pCO->ijke[1]- nghost, pCO->ijks[2]- nghost, pCO->ijke[2]- nghost) ; */

	for (k=pCO->ijks[2] - nghost; k<= pCO->ijke[2]+1 - nghost; k++) {
	  for (j=pCO->ijks[1] - nghost; j<= pCO->ijke[1]+1 - nghost; j++) {



	    /*Calculate the index on the ionFlx array corresponding to a given cell on this (fine) grid*/
	    /*Remember that ionFlx is a 1D array of values*/
	    indexarith = (k-(pCO->ijks[2]-nghost))*(pCO->ijke[1] - pCO->ijks[1] + 2)+j-(pCO->ijks[1]-nghost);

	    /*Find my location on the current (fine) grid*/
	    /*AT 12/23/12 BE CAREFUL WITH FACTORS OF 2 IN INDICES*/
	    ks = k*2 - pDomain->Disp[2];
	    js = j*2 - pDomain->Disp[1];

	    /* if (pCO->ionFlx[dim][indexarith] > .1) */
	    /*   fprintf(stderr,"Setting ks and js to be %e \n", pCO->ionFlx[dim][indexarith]); */

	    /*Assign my parent's ionizing flux to my EdgeFlux cell */
	    pGrid->EdgeFlux[ks][js][fixed] = pCO->ionFlx[dim][indexarith];

	    /*AT 01/14/13: The following is a clumsy way of linearly interpolating the values to 2 cells.*/
	    /*TO_DO: I also need to check that the if statements will hold up for grids of diff. size and that it works with MPI in the right area*/
	    
	    /*Assign the same value to the other 3 fine cells that make up the one coarse cell */
	    if (j < (pCO->ijke[1]+1 - nghost)) {
		if (k < (pCO->ijke[2]+1 - nghost)) {	
		pGrid->EdgeFlux[ks+1][js+1][fixed] = pCO->ionFlx[dim][indexarith];
		pGrid->EdgeFlux[ks][js+1][fixed] = pCO->ionFlx[dim][indexarith];
		pGrid->EdgeFlux[ks+1][js][fixed] = pCO->ionFlx[dim][indexarith];

		/* if ( pCO->ionFlx[dim][indexarith] > 1.) */
		/*   fprintf(stderr,"On level: %d setting [%d +=1][%d +=1][%d] ionflx[%d][%d]  %e \n", level, ks, js, fixed, dim, indexarith, pCO->ionFlx[dim][indexarith]); */

	      }
	      
	      /*Handle edge cases*/
	      else {
		pGrid->EdgeFlux[ks][js+1][fixed] = pCO->ionFlx[dim][indexarith];

		/* if ( pCO->ionFlx[dim][indexarith] > 1.) */
		/*   fprintf(stderr,MAKE_BLUE "1: On level: %d setting [%d][%d +=1][%d] ionflx[%d][%d]  %e \n" RESET_COLOR2, level, ks, js, fixed, dim, indexarith, pCO->ionFlx[dim][indexarith]); */

	      } 
	    }
	    /*Handle edge cases*/
	    else {
	      if (k < pCO->ijke[2]+1 - nghost ) {
		pGrid->EdgeFlux[ks+1][js][fixed] = pCO->ionFlx[dim][indexarith];

		/* if ( pCO->ionFlx[dim][indexarith] > 1.) */
		/*   fprintf(stderr,MAKE_RED"2: On level: %d setting [%d +=1][%d][%d] ionflx[%d][%d]  %e \n" RESET_COLOR2, level, ks, js, fixed, dim, indexarith, pCO->ionFlx[dim][indexarith]); */

	      }
	    }
	    
	  }
	}
      }

#endif
    }
}

void ionrad_prolong_snd(GridS *pGrid, int dim, int level, int domnumber)
{
  MeshS *pMesh = pGrid->Mesh;
  int ncg, ierr;
  GridOvrlpS *pCO;
  int fixed, arrsize;
  int i, j, k;
  int indexarith;

#ifdef MPI_PARALLEL
  MPI_Request *send_rq;
  send_rq = (MPI_Request*) calloc_1d_array(pGrid->NCGrid,sizeof(MPI_Request)) ;
  int tag1, tag2, tag3;
  char temp[10];
#endif

/* Loop over children grids to fill their buffer arrays*/
/* If not using MPI, this is the only operation that is performed*/
/* If using MPI, data will be sent to child grid in next step*/
  for (ncg=0; ncg<(pGrid->NCGrid); ncg++){
    pCO=(GridOvrlpS*)&(pGrid->CGrid[ncg]);

    /* fprintf(stderr, "Actually filling buffer \n"); */
      /* fprintf(stderr, "CHILD # %d of %d I think my child's x- start is %d and end %d. So the edgeflux corresponds to those - %d \n", ncg+1, pGrid->NCGrid, pCO->ijks[0], pCO->ijke[0], nghost); */

    switch(dim) {

    /*For the +x/-x direction*/
    case 0: case 1: {
      if (fmod(dim,2) == 0) {
	fixed = pCO->ijks[0] - nghost;
      } else {
	fixed = pCO->ijke[0] + 1 - nghost;
      }

      if(pCO->ionFlx[dim] != NULL) {
	/* fprintf(stderr, MAKE_BLUE"I am %d (Level %d) and sending data to i: %d, j: %d - %d, k: %d -%d relative to my overlap index\n " RESET_COLOR, myID_Comm_world, level, fixed, pCO->ijks[1] - nghost, pCO->ijke[1] - nghost, pCO->ijks[2] - nghost, pCO->ijke[2] - nghost); */
	for (k=pCO->ijks[2] - nghost; k<= pCO->ijke[2]+1 - nghost; k++) {
	  for (j=pCO->ijks[1] - nghost; j<= pCO->ijke[1]+1 - nghost; j++) {
	    indexarith = (k-(pCO->ijks[2]-nghost))*(pCO->ijke[1] - pCO->ijks[1] + 2)+j-(pCO->ijks[1]-nghost);

	    /*Store data in the child grid overlap structure ionFlx (which is a 1D array for the purposes of MPI communication)*/
	    pCO->ionFlx[dim][indexarith] = pGrid->EdgeFlux[k][j][fixed];
	    /* if (pGrid->EdgeFlux[k][j][fixed] > 1.) */
	    /*   fprintf(stderr,"On level: %d setting ionflx[%d][%d] to [%d][%d][%d] %e \n", level, dim, indexarith,k-2, j-2, fixed-2, pGrid->EdgeFlux[k][j][fixed]); */
	  }
	}
	/*TO_DO: Will need to check indexing to see if it's +1 or +2*/
	arrsize = (pCO->ijke[2] + 2 - pCO->ijks[2]) * (pCO->ijke[1] + 2 - pCO->ijks[1]);
	/* fprintf(stderr, MAKE_RED "Sending array size %d with 2s: %d 2e: %d [SENT] \n" RESET_COLOR, arrsize, pCO->ijks[2], pCO->ijke[2]); */
      } else{
	/* fprintf(stderr, MAKE_BLUE"I am %d (Level %d) and didn't fill a buffer. For reference, I have i: %d, j: %d - %d, k: %d -%d relative to my overlap index\n " RESET_COLOR, myID_Comm_world, level, fixed, pCO->ijks[1] - nghost, pCO->ijke[1] - nghost, pCO->ijks[2] - nghost, pCO->ijke[2] - nghost); */
      }
      break;
    }

    /*For the +y/-y direction*/
    case 2: case 3: {
      if (fmod(dim,2) == 0) {
	fixed = pCO->ijks[1] - nghost;
      } else {
	fixed = pCO->ijke[1] + 1 - nghost;
      }
      
      if(pCO->ionFlx[dim] != NULL) {
	  for (k=pCO->ijks[2] - nghost; k<= pCO->ijke[2]+1 - nghost; k++) {
	    for (i=pCO->ijks[0] - nghost; i<= pCO->ijke[0]+1 - nghost; i++) {
	      pCO->ionFlx[dim][(k-(pCO->ijks[2]-nghost))*(pCO->ijke[0] - pCO->ijks[0] + 2)+i-(pCO->ijks[0]-nghost)] = pGrid->EdgeFlux[k][fixed][i];
	    }
	  }
	  arrsize = (pCO->ijke[2] + 2 - pCO->ijks[2]) * (pCO->ijke[0] + 2 - pCO->ijks[0]);
      }
      break;
    }

    /*For the +z/-z direction*/
    case 4: case 5: {
      if (fmod(dim,2) == 0) {
	fixed = pCO->ijks[2] - nghost;
      } else {
	fixed = pCO->ijke[2] + 1 - nghost;
      }
      if(pCO->ionFlx[dim] != NULL) {
	for (j=pCO->ijks[1] - nghost; j<= pCO->ijke[1]+1 - nghost; j++) {
	  for (i=pCO->ijks[0] - nghost; i<= pCO->ijke[0]+1 - nghost; i++) {
	    pCO->ionFlx[dim][(j-(pCO->ijks[1]-nghost))*(pCO->ijke[0] - pCO->ijks[0] + 2)+i-(pCO->ijks[0]-nghost)] = pGrid->EdgeFlux[fixed][j][i];
	  }
	}
	arrsize = (pCO->ijke[0] + 2 - pCO->ijks[0]) * (pCO->ijke[1] + 2 - pCO->ijks[1]);
      }
      break;
    }
    }
     
    /*AT 9/26/12: I think the NULL check only needs to be here and needs to be modified above.*/
    /* The point of such a check is to ensure that we're not trying to communicate, when there's an overlapping grid not in the direction of propagation*/
    /*Send buffer arrays of radiation flux to children grids*/

    if(pCO->ionFlx[dim] != NULL) {
#ifdef MPI_PARALLEL
      /* 	tag1 = atoi(temp); */
      /* 	tag2 = level*1000000 + myID_Comm_world * 10000 + (level+1) * 100 + pCO->ID; */
	/* fprintf(stderr, "sndconcat: %d, powers:%d \n", tag1, tag2); */

      tag3 = domnumber + 100;

      /*Send data to child grid*/
      ierr = MPI_Isend(pCO->ionFlx[dim], arrsize, MP_RL, pCO->ID, tag3, pMesh->Domain[level][domnumber].Comm_Children, &send_rq[ncg]);
      
      /* fprintf(stderr, "Sent data to child ID %d using tag %d. I'm on level %d \n", pCO->ID, tag3, level); */
      /* fprintf(stderr, "Left x: %d, right x:%d I sent my data for child %d of %d\n", pGrid->lx1_id, pGrid->rx1_id,ncg+1, pGrid->NCGrid);*/
#endif
    }
  }
}


/*Older function for prolonging on a single processor.  Has been modified since it's initial use, so don't trust.*/
/* void ionrad_prolongate(DomainS *pD) */
/* { */
/*   MeshS *pM = pD->Mesh; */
/*   GridS *pG = pD->Grid; */
/*   GridS *highergrid; */
/*   DomainS *higherdom; */
/*   int nd, ndtot, higherlevel, level; */
/*   int dir = (pM->radplanelist)->dir[0]; */
/*   int upperi, upperj, upperk, fixed; */
/*   SideS thisd, parentd, overlapd; */
/*   int isDOverlap; */
/*   int i, j, k, lr; */

/*   level = pD->Level; */
/*   higherlevel = level-1; */
/*   ndtot = pM->DomainsPerLevel[higherlevel]; */

/*   /\*Find location of this domain *\/ */
/*   /\*relative to the origin in the higher level units*\/ */
/*   for (i=0; i<3; i++){ */
/*     thisd.ijkl[i] = pD->Disp[i]/2; */
/*     thisd.ijkr[i] = (pD->Disp[i] + pD->Nx[i])/2; */
/*   } */

/*   /\*Find the parent domain for the current domain*\/ */
/*   if (ndtot == 1) { */
/*     highergrid = pM->Domain[higherlevel][ndtot-1].Grid; */
/*     higherdom =(DomainS*)&(pM->Domain[higherlevel][ndtot-1]); */

/*     /\*If this is the first time searching for a parent domain, print to screen the information*\/ */
/*     if (pM->time == 0) fprintf(stderr, "Dom Level %d No %d has parent Dom Level %d No %d \n", level, pD->DomNumber, pM->Domain[higherlevel][0].Level, pM->Domain[higherlevel][0].DomNumber);  */

/*   } else { */
/*     for (nd=0; nd<ndtot-1; nd++){ */
/*       for (i=0; i<3; i++){ */
/* 	parentd.ijkl[i] = pM->Domain[higherlevel][nd].Disp[i]; */
/* 	parentd.ijkr[i] = pM->Domain[higherlevel][nd].Disp[i] + pM->Domain[higherlevel][nd].Nx[i]; */
/*       } */
/*       isDOverlap = checkOverlap(&thisd, &parentd, &overlapd); */
/*       if (isDOverlap == 1){ */
/* 	highergrid = pM->Domain[higherlevel][nd].Grid; */
/* 	higherdom = (DomainS*)&(pM->Domain[higherlevel][nd]); */
/* 	if (pM->time == 0) fprintf(stderr, "Dom Level %d No %d has parent Dom Level %d No %d \n", level, pD->DomNumber, pM->Domain[higherlevel][nd].Level, pM->Domain[higherlevel][nd].DomNumber);  */

/* 	/\*Check that the system believes this domain has a child grid*\/ */
/* 	if (highergrid->NCGrid < 0) ath_error("Contradiction of child/parent grid \n"); */
/* 	break; */
/*       } */
/*     } */
/*   } */

/*   if (highergrid == NULL) ath_error("No parent grid \n"); */

/*   lr = (dir < 0) ? 1 : -1; */
  
/*   /\* switch(dir) { *\/ */
/*   /\* case -1: case 1: { *\/ */

/*     if (lr > 0) { */
/*       fixed = 0; */
/*       for (k=0; k<=pG->Nx[2]; k++) { */
/*   	for (j=0; j<=pG->Nx[1]; j++) { */
/* 	  /\* fprintf(stderr, "%d %d Flux: %e \n", k, j, pG->EdgeFlux[k][j][fixed]); *\/ */
/*   	  upperk =(int)( floor(k/2.) + thisd.ijkl[2] - higherdom->Disp[2]); */
/* 	  upperj = (int) (floor(j/2) + thisd.ijkl[1] - higherdom->Disp[1]); */
/* 	  /\* fprintf(stderr, "This: %d Upper: %d This: %d Upper: %d \n", k, upperk, j, upperj); *\/ */
/*   	  pG->EdgeFlux[k][j][fixed] =highergrid->EdgeFlux[upperk][upperj][thisd.ijkl[0]-higherdom->Disp[0]]; */
/* 	  /\*AT: 12/23/12  - Trying to write pPO equivalent*\/ */
/* 	  /\* Instead of using highergrid, need to use pPO->ionFlx.  The question is >> -how to identify what is the higher grid >>-or if I'm using pPO how do I fill ionFlx, since I can see filling my pCO->ionFlx.  >> If I'm a child, can I access my parent's pCO? */
/*  *\/ */



/*   	} */
/*       } */
/*     } else {ath_error("Doesnt work yet \n");} */

/*   /\*   break; *\/ */
/*   /\* } *\/ */
/*   /\* case -2: case 2: { *\/ */
/*   /\*   d2 = 0; *\/ */
/*   /\*   if (lr > 0) { *\/ */
/*   /\*    } else { *\/ */
/*   /\*   } *\/ */
/*   /\*   break; *\/ */
/*   /\* } *\/ */

/* } */

/* #ifdef MPI_PARALLEL */
/* /\* Need to send to all children upon completion of coarse grid and receive by all fine grids*\/ */
/* /\* Need to have call to this not from ionrad_3d.c in the finegrid case.  Instead it needs to be called even after the coarse.*\/ */
/* /\*So maybe like - at end of every level, do a send.  And if not root domain, do a receive*\/ */
/* /\*So if nparents != 0, then receive, and if nchildren !=0 then do a send*\/ */
/* #endif */

/* #endif */

/*   /\* GridOvrlpS *pCO, *pPO; *\/ */
/*   /\* int npg, ncg; *\/ */
/*   /\* /\\* GridsDataS ***Ginfo = pD->GData; *\\/ *\/ */
  
/*   /\* fprintf(stderr, "---------\n Domain level: %d, number:%d, X-disp:%d, zones:%d, Nparents:%d, Nchild:%d \n", pD->Level, pD->DomNumber, pD->Disp[0], pD->Nx[0], pG->NPGrid,pG->NCGrid); *\/ */
/*   /\* /\\* fprintf(stderr, "ID %d \n", Ginfo->ID_Comm_world); *\\/ *\/ */

/*   /\* for (npg=0; npg<(pG->NPGrid); npg++){ *\/ */
/*   /\*   pPO = (GridOvrlpS*)&(pG->PGrid[npg]); *\/ */
/*   /\*   fprintf(stderr, "Parent: X:ijks %d ijke %d, ID:%d \n", pPO->ijks[0], pPO->ijke[0], pPO->ID); *\/ */
/*   /\* } *\/ */

/*   /\* for (ncg=0; ncg<(pG->NmyCGrid); ncg++){ *\/ */
/*   /\*   pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);    /\\* ptr to child Grid overlap *\\/ *\/ */
/*   /\*   fprintf(stderr, "Child: X:ijks %d ijke %d, ID:%d \n", pCO->ijks[0], pCO->ijke[0], pCO->ID); *\/ */
/*   /\* } *\/ */
