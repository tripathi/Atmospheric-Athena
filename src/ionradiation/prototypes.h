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
/* ionrad.c */
void ion_radtransfer_init_domain(MeshS *pM);
VDFun_t ion_radtransfer_init(MeshS *pM, int ires);

/*----------------------------------------------------------------------------*/
/* ionrad_3d.c */
void ion_radtransfer_3d(GridS *pG);
void ion_radtransfer_init_3d(GridS *pG, DomainS *pD, int ires);
void ion_radtransfer_init_domain_3d(GridS *pG, DomainS *pD);

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
void get_ph_rate_plane(Real initflux, int dir, Real ***ph_rate, GridS *pGrid);
#endif /* ION_RADPLANE */

#endif /* IONRAD_PROTOTYPES_H */
