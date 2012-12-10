#ifndef IONRAD_PROTOTYPES_H
#define IONRAD_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions from the following files:
 *   ionrad.c
 *   ionrad_3d.c
 *   ionrad_chemistry.c
 *   ionradplane_3d.c
 *   ionradpoint_3d.c
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

#ifdef ION_RADIATION
/*----------------------------------------------------------------------------*/
/* bvals_ionrad.c  */
void bvals_ionrad_init(MeshS *pM);
void bvals_ionrad(DomainS *pDomain);

/*----------------------------------------------------------------------------*/
/* ionrad.c */
void ion_radtransfer_init_domain(MeshS *pM);
VDFun_t ion_radtransfer_init(MeshS *pM, int ires);

/*----------------------------------------------------------------------------*/
/* ionrad_3d.c */
void ion_radtransfer_3d(DomainS *pD);
void ion_radtransfer_init_3d(GridS *pG, DomainS *pD, int ires);
void ion_radtransfer_init_domain_3d(GridS *pG, DomainS *pD);
void set_coarse_time();
void clear_coarse_time();

/*----------------------------------------------------------------------------*/
/* ionrad_chemistry.c */
Real recomb_rate_coef(Real T);
Real coll_ion_rate_coef(Real T);
Real recomb_cool_rate_coef(Real T);
Real dmc_cool_rate(Real x, Real T);
Real osterbrock_cool_rate(Real T);
Real ki_cool_rate(Real T);
Real ki_heat_rate(void);
#endif /* ION_RADIATION */

#ifdef ION_RADPLANE
/*----------------------------------------------------------------------------*/
/* ionradplane_3d.c */
void add_radplane_3d(GridS *pGrid, int dir, Real flux);
void ion_radplane_init_domain_3d(GridS *pGrid, DomainS *pDomain);
#ifdef MPI_PARALLEL
void get_ph_rate_plane(Real initflux, int dir, Real ***ph_rate, GridS *pGrid, MPI_Comm Comm_Domain);
#else
void get_ph_rate_plane(Real initflux, int dir, Real ***ph_rate, GridS *pGrid);
#endif

/*----------------------------------------------------------------------------*/
/* ionrad_smr.c */
void ionrad_prolongate(DomainS *pD);
void ionrad_prolong_rcv(GridS *pGrid, int dim, int level, int domnumber);
void ionrad_prolong_snd(GridS *pGrid, int dim, int level, int domnumber);

#endif /* ION_RADPLANE */

#endif /* IONRAD_PROTOTYPES_H */
