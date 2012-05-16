#include "copyright.h"

/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*! \fn void Prolongate(MeshS *pM)
 *  \brief Sets BC on fine Grid by prolongation (interpolation) of
 *     coarse Grid solution into fine grid ghost zones */
void Prolongate(MeshS *pM)
{
  GridS *pG;
  int nDim,nl,nd,ncg,dim,npg,rbufN,id,l,m,n,mend,nend;
  int i,ii,ics,ice,ips,ipe,igzs,igze;
  int j,jj,jcs,jce,jps,jpe,jgzs,jgze;
  int k,kk,kcs,kce,kps,kpe,kgzs,kgze;
  int ngz1,ngz2,ngz3;
  double *pRcv,*pSnd;
  GridOvrlpS *pCO, *pPO;
  ConsS ProlongedC[2][2][2];
#if (NSCALARS > 0)
  int ns;
#endif
#ifdef MHD
  Real3Vect BGZ[3][3][3], ProlongedF[3][3][3];
#endif

/* number of dimensions in Grid. */
  nDim=1;
  for (dim=1; dim<3; dim++) if (pM->Nx[dim]>1) nDim++;

/* Loop over all levels, starting at root level */

  for (nl=0; nl<(pM->NLevels); nl++){

/*=== Step 1. Send step ======================================================*/
/* Loop over all Domains, and send ghost zones to all child Grids. */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    pG=pM->Domain[nl][nd].Grid;

    //DONT UNDERSTAND THIS LINE
    for(i=0; i<maxND; i++) start_addrP[i] = 0;

    for (ncg=0; ncg<(pG->NCGrid); ncg++){
      pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);    /* ptr to child Grid overlap */

      //DONT UNDERSTAND
/* index send_buf with DomN of child, since could be multiple child Domains on
 * same processor.  Start address must be different for each DomN */
      pSnd = (double*)&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]); 

      for (dim=0; dim<(2*nDim); dim++){
        if (pCO->myFlx[dim] != NULL) {

/* Get coordinates ON THIS GRID of zones that overlap child Grid ghost zones */

          ics = pCO->ijks[0] - (nghost/2) - 1;
          ice = pCO->ijke[0] + (nghost/2) + 1;
          if (pG->Nx[1] > 1) {
            jcs = pCO->ijks[1] - (nghost/2) - 1;
            jce = pCO->ijke[1] + (nghost/2) + 1;
          } else {
            jcs = pCO->ijks[1];
            jce = pCO->ijke[1];
          }
          if (pG->Nx[2] > 1) {
            kcs = pCO->ijks[2] - (nghost/2) - 1;
            kce = pCO->ijke[2] + (nghost/2) + 1;
          } else {
            kcs = pCO->ijks[2];
            kce = pCO->ijke[2];
          }
          if (dim == 0) ice = pCO->ijks[0];
          if (dim == 1) ics = pCO->ijke[0];
          if (dim == 2) jce = pCO->ijks[1];
          if (dim == 3) jcs = pC O->ijke[1];
          if (dim == 4) kce = pCO->ijks[2];
          if (dim == 5) kcs = pCO->ijke[2];

/*--- Step 1a. ---------------------------------------------------------------*/
/* Load send buffer with values in zones that overlap child ghost zones */

          for (k=kcs; k<=kce; k++) {
          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
/*             *(pSnd++) = pG->U[k][j][i].d; */
/*             *(pSnd++) = pG->U[k][j][i].M1; */
/*             *(pSnd++) = pG->U[k][j][i].M2; */
/*             *(pSnd++) = pG->U[k][j][i].M3; */
/* #ifndef BAROTROPIC */
/*             *(pSnd++) = pG->U[k][j][i].E; */
/* #endif */
#ifdef MHD
            /* *(pSnd++) = pG->U[k][j][i].B1c; */
            /* *(pSnd++) = pG->U[k][j][i].B2c; */
            /* *(pSnd++) = pG->U[k][j][i].B3c; */
            *(pSnd++) = pG->B1i[k][j][i];
            *(pSnd++) = pG->B2i[k][j][i];
            *(pSnd++) = pG->B3i[k][j][i];
#endif
#if (NSCALARS > 0)
            for (ns=0; ns<NSCALARS; ns++) {
               *(pSnd++) = pG->U[k][j][i].s[ns];
            }
#endif
          }}}
        }
      }

/*--- Step 1b. ---------------------------------------------------------------*/
/* non-blocking send of data to child, using Domain number as tag. */

      start_addrP[pCO->DomN] += pG->CGrid[ncg].nWordsP;

    } /* end loop over child grids */
  }} /* end loop over Domains */

/*=== Step 2. Get step =======================================================*/
/* Loop over all Domains, get data sent by parent Grids, and prolongate solution
 * into ghost zones. */


  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    pG=pM->Domain[nl][nd].Grid;          /* set pointer to Grid */
    rbufN = (nl % 2);

    for (npg=0; npg<(pG->NPGrid); npg++){

/* If parent Grid is on this processor, data is at start of recv buffer */

      if (npg < pG->NmyPGrid) {
        pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
        pRcv = (double*)&(recv_bufP[rbufN][nd][0]);
      } else {
/* If not MPI_PARALLEL, and parent Grid not on this processor, then error */
        ath_error("[Prolong]: no Parent Grid on Domain[%d][%d]\n",nl,nd);
      }

/*=== Step 3. Set ghost zones ================================================*/
/* Loop over 6 boundaries, set ghost zones */

      for (dim=0; dim<(2*nDim); dim++){
        if (pPO->myFlx[dim] != NULL) {

/*--- Steps 3a.  Set GZ and BFld arrays --------------------------------------*/
/* Compute size of array containing ghost zones from parent Grid.  Set
 * starting and ending indices of GZ array. */

          if (dim == 0 || dim == 1) {
            ngz1 = (nghost/2) + 2;
            id = 0;
          } else {
            ngz1 = (pPO->ijke[0] - pPO->ijks[0] + 1)/2 + nghost + 2;
          }

          if (dim == 2 || dim == 3) {
            ngz2 = (nghost/2) + 2;
            id = 1;
          } else {
            ngz2 = (pPO->ijke[1] - pPO->ijks[1] + 1)/2 + nghost + 2;
          }

          if (dim == 4 || dim == 5) {
            ngz3 = (nghost/2) + 2;
            id = 2;
          } else {
            ngz3 = (pPO->ijke[2] - pPO->ijks[2] + 1)/2 + nghost + 2;
          }

          igzs = 0;
          igze = ngz1-1;
          if (pG->Nx[1] > 1) {
            jgzs = 0;
            jgze = ngz2-1;
            mend = 1;
          } else {
            ngz2 = 1;
            jgzs = 1;
            jgze = 1;
            mend = 0;
          }
          if (pG->Nx[2] > 1) {
            kgzs = 0;
            kgze = ngz3-1;
            nend = 1;
          } else {
            ngz3 = 1;
            kgzs = 1;
            kgze = 1;
            nend = 0;
          }

/* Load GZ array with values in receive buffer */

          for (k=kgzs; k<=kgze; k++) {
          for (j=jgzs; j<=jgze; j++) {
          for (i=igzs; i<=igze; i++) {
/*             GZ[id][k][j][i].d  = *(pRcv++); */
/*             GZ[id][k][j][i].M1 = *(pRcv++); */
/*             GZ[id][k][j][i].M2 = *(pRcv++); */
/*             GZ[id][k][j][i].M3 = *(pRcv++); */
/* #ifndef BAROTROPIC */
/*             GZ[id][k][j][i].E = *(pRcv++); */
/* #endif */
#ifdef MHD
            /* GZ[id][k][j][i].B1c = *(pRcv++); */
            /* GZ[id][k][j][i].B2c = *(pRcv++); */
            /* GZ[id][k][j][i].B3c = *(pRcv++); */
            BFld[id][k][j][i].x = *(pRcv++);
            BFld[id][k][j][i].y = *(pRcv++);
            BFld[id][k][j][i].z = *(pRcv++);
#endif
#if (NSCALARS > 0)
            for (ns=0; ns<NSCALARS; ns++) {
              GZ[id][k][j][i].s[ns] = *(pRcv++);
            }
#endif
          }}}

/* Set BC on GZ array in 1D; and on GZ and BFld arrays in 2D */

          if (nDim == 1) {
            for (i=igzs; i<=igze; i++) {
              GZ[id][1][0][i] = GZ[id][1][1][i];
              GZ[id][1][2][i] = GZ[id][1][1][i];
              GZ[id][0][1][i] = GZ[id][1][1][i];
              GZ[id][2][1][i] = GZ[id][1][1][i];
            }
          }

          if (nDim == 2) {
            for (j=jgzs; j<=jgze; j++) {
            for (i=igzs; i<=igze; i++) {
              GZ[id][0][j][i] = GZ[id][1][j][i];
              GZ[id][2][j][i] = GZ[id][1][j][i];
            }}
#ifdef MHD
            for (j=jgzs; j<=jgze; j++) {
            for (i=igzs; i<=igze; i++) {
              BFld[id][0][j][i] = BFld[id][1][j][i];
              BFld[id][2][j][i] = BFld[id][1][j][i];
            }}
#endif /* MHD */
          }

/*--- Steps 3b.  Prolongate cell-centered values -----------------------------*/
/* Get coordinates ON THIS GRID of ghost zones that overlap parent Grid */

          ips = pPO->ijks[0] - nghost;
          ipe = pPO->ijke[0] + nghost;
          if (pG->Nx[1] > 1) {
            jps = pPO->ijks[1] - nghost;
            jpe = pPO->ijke[1] + nghost;
          } else {
            jps = pPO->ijks[1];
            jpe = pPO->ijke[1];
          }
          if (pG->Nx[2] > 1) {
            kps = pPO->ijks[2] - nghost;
            kpe = pPO->ijke[2] + nghost;
          } else {
            kps = pPO->ijks[2];
            kpe = pPO->ijke[2];
          }
          if (dim == 0) {ipe = pPO->ijks[0] - 1;}
          if (dim == 1) {ips = pPO->ijke[0] + 1;}
          if (dim == 2) {jpe = pPO->ijks[1] - 1;}
          if (dim == 3) {jps = pPO->ijke[1] + 1;}
          if (dim == 4) {kpe = pPO->ijks[2] - 1;}
          if (dim == 5) {kps = pPO->ijke[2] + 1;}

/* Prolongate these values in ghost zones */

          for (k=kps, kk=1; k<=kpe; k+=2, kk++) {
          for (j=jps, jj=1; j<=jpe; j+=2, jj++) {
          for (i=ips, ii=1; i<=ipe; i+=2, ii++) {

            ProCon(GZ[id][kk][jj][ii-1],GZ[id][kk][jj][ii],GZ[id][kk][jj][ii+1],
                   GZ[id][kk][jj-1][ii],                   GZ[id][kk][jj+1][ii],
                   GZ[id][kk-1][jj][ii],                   GZ[id][kk+1][jj][ii],
                   ProlongedC);

/* 1D/2D/3D problem, set solution prolongated in x1 */

            for (n=0; n<=nend; n++) {
            for (m=0; m<=mend; m++) {
            for (l=0; l<=1; l++) {
/*               pG->U[k+n][j+m][i+l].d  = ProlongedC[n][m][l].d; */
/*               pG->U[k+n][j+m][i+l].M1 = ProlongedC[n][m][l].M1; */
/*               pG->U[k+n][j+m][i+l].M2 = ProlongedC[n][m][l].M2; */
/*               pG->U[k+n][j+m][i+l].M3 = ProlongedC[n][m][l].M3; */
/* #ifndef BAROTROPIC */
/*               pG->U[k+n][j+m][i+l].E  = ProlongedC[n][m][l].E; */
/* #endif */
#ifdef MHD
              /* pG->U[k+n][j+m][i+l].B1c = ProlongedC[n][m][l].B1c; */
              /* pG->U[k+n][j+m][i+l].B2c = ProlongedC[n][m][l].B2c; */
              /* pG->U[k+n][j+m][i+l].B3c = ProlongedC[n][m][l].B3c; */
#endif
#if (NSCALARS > 0)
              for (ns=0; ns<NSCALARS; ns++) 
                pG->U[k+n][j+m][i+l].s[ns] = ProlongedC[n][m][l].s[ns];
#endif
            }}}

#ifdef MHD
/*--- Steps 3c.  Prolongate face-centered B ----------------------------------*/
/* Set prolonged face-centered B fields for 1D (trivial case)  */

            if (nDim == 1) {
              for (l=0; l<=1; l++) {
                pG->B1i[k][j][i+l] = pG->U[k][j][i+l].B1c;
                pG->B2i[k][j][i+l] = pG->U[k][j][i+l].B2c;
                pG->B3i[k][j][i+l] = pG->U[k][j][i+l].B3c;
              }
            } else {
              for (n=0; n<3; n++) {
              for (m=0; m<3; m++) {
              for (l=0; l<3; l++) {
                ProlongedF[n][m][l].x = 0.0;
                ProlongedF[n][m][l].y = 0.0;
                ProlongedF[n][m][l].z = 0.0;
              }}}
            }

/* Load B-field ghost zone array with values read from Rcv buffer in 2D/3D */

            if (nDim == 2 || nDim ==3) {
              for (n=0; n<3; n++) {
              for (m=0; m<3; m++) {
              for (l=0; l<3; l++) {
                BGZ[n][m][l].x = BFld[id][kk+(n-1)][jj+(m-1)][ii+(l-1)].x;
                BGZ[n][m][l].y = BFld[id][kk+(n-1)][jj+(m-1)][ii+(l-1)].y;
                BGZ[n][m][l].z = BFld[id][kk+(n-1)][jj+(m-1)][ii+(l-1)].z;
              }}}

/* If edge of cell touches fine/coarse boundary, use fine grid fields for the
 * normal component at interface. ProFld will not overwrite these values.  If
 * the start/end of boundary is between MPI Grids (pPO->myFlx[]==NULL), then use
 * fine grid fields in corners as well. */

/* inner x1 boundary */
              if ((dim == 0) &&
                  (i == (ipe-1)) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[0][0][2].x = pG->B1i[k][j  ][i+2];
                ProlongedF[0][1][2].x = pG->B1i[k][j+1][i+2];
                ProlongedF[1][0][2].x = pG->B1i[k][j  ][i+2];
                ProlongedF[1][1][2].x = pG->B1i[k][j+1][i+2];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][0][2].x = pG->B1i[k+1][j  ][i+2];
                  ProlongedF[1][1][2].x = pG->B1i[k+1][j+1][i+2];
                }
              }

/* outer x1 boundary */
              if ((dim == 1) &&
                  (i == ips) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[0][0][0].x = pG->B1i[k][j  ][i];
                ProlongedF[0][1][0].x = pG->B1i[k][j+1][i];
                ProlongedF[1][0][0].x = pG->B1i[k][j  ][i];
                ProlongedF[1][1][0].x = pG->B1i[k][j+1][i];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][0][0].x = pG->B1i[k+1][j  ][i];
                  ProlongedF[1][1][0].x = pG->B1i[k+1][j+1][i];
                }
              }

/* inner x2 boundary */
              if ((dim == 2) &&
                  (j == (jpe-1)) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) ){
                ProlongedF[0][2][0].y = pG->B2i[k][j+2][i  ];
                ProlongedF[0][2][1].y = pG->B2i[k][j+2][i+1];
                ProlongedF[1][2][0].y = pG->B2i[k][j+2][i  ];
                ProlongedF[1][2][1].y = pG->B2i[k][j+2][i+1];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][2][0].y = pG->B2i[k+1][j+2][i  ];
                  ProlongedF[1][2][1].y = pG->B2i[k+1][j+2][i+1];
                }
              }

/* outer x2 boundary */
              if ((dim == 3) &&
                  (j == jps) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) ){
                ProlongedF[0][0][0].y = pG->B2i[k][j][i  ];
                ProlongedF[0][0][1].y = pG->B2i[k][j][i+1];
                ProlongedF[1][0][0].y = pG->B2i[k][j][i  ];
                ProlongedF[1][0][1].y = pG->B2i[k][j][i+1];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][0][0].y = pG->B2i[k+1][j][i  ];
                  ProlongedF[1][0][1].y = pG->B2i[k+1][j][i+1];
                }
              }

/* inner x3 boundary */
              if ((dim == 4) &&
                  (k == (kpe-1)) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[2][0][0].z = pG->B3i[k+2][j  ][i  ];
                ProlongedF[2][0][1].z = pG->B3i[k+2][j  ][i+1];
                ProlongedF[2][1][0].z = pG->B3i[k+2][j+1][i  ];
                ProlongedF[2][1][1].z = pG->B3i[k+2][j+1][i+1];
              }

/* outer x3 boundary */
              if ((dim == 5) && 
                  (k == kps) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[0][0][0].z = pG->B3i[k][j  ][i  ];
                ProlongedF[0][0][1].z = pG->B3i[k][j  ][i+1];
                ProlongedF[0][1][0].z = pG->B3i[k][j+1][i  ];
                ProlongedF[0][1][1].z = pG->B3i[k][j+1][i+1];
              }

              ProFld(BGZ, ProlongedF, pG->dx1, pG->dx2, pG->dx3);

              for (n=0; n<=nend; n++) {
              for (m=0; m<=mend; m++) {
              for (l=0; l<=1; l++) {
                if (dim != 1 || (i+l) != ips)
                  pG->B1i[k+n][j+m][i+l] = ProlongedF[n][m][l].x;
                if (dim != 3 || (j+m) != jps)
                  pG->B2i[k+n][j+m][i+l] = ProlongedF[n][m][l].y;
                if (dim != 5 || (k+n) != kps)
                  pG->B3i[k+n][j+m][i+l] = ProlongedF[n][m][l].z;

                pG->U[k+n][j+m][i+l].B1c = 
                  0.5*(ProlongedF[n][m][l].x + ProlongedF[n][m][l+1].x);
                pG->U[k+n][j+m][i+l].B2c = 
                  0.5*(ProlongedF[n][m][l].y + ProlongedF[n][m+1][l].y);
                pG->U[k+n][j+m][i+l].B3c = 
                  0.5*(ProlongedF[n][m][l].z + ProlongedF[n+1][m][l].z);
              }}}
            }

#endif /* MHD */
          }}}

        }
      } /* end loop over dims */

    } /* end loop over parent grids */
  }} /* end loop over Domains */

/*=== Step 4. Clear send_bufP ================================================*/
/* For child/parent Grids on same processor, copy send_bufP into recv_bufP to
 * prevent "send" in Step 1 above from over-writing data in buffer on the next
 * iteration of the loop over levels (for nl=nl+1). */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) { 
      pG=pM->Domain[nl][nd].Grid; 
      rbufN = ((nl+1) % 2);

/* For each Domain nd, the data for child grid on same processor must be in
 * first element of CGrid array */ 

      for (ncg=0; ncg<(pG->NmyCGrid); ncg++){
        pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);    /* ptr to child Grid overlap */

        for (i=0; i<pCO->nWordsP; i++) {
          recv_bufP[rbufN][pCO->DomN][i]=send_bufP[pCO->DomN][i];
        }
      }
    }
  }

  } /* end loop over levels */
}
