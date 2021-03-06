#include "../copyright.h"
/*==============================================================================
 * FILE: ionrad_3d.c
 *
 * PURPOSE: Contains functions to compute an ionization radiative transfer
 *   from update, using the algorithm described in Krumholz, Stone, 
 *   & Gardiner (2007).
 *
 *   Use of these routines requires that --enable-ion-radiation be set
 *   at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   ion_radtransfer_3d             - does an ionizing radiative transfer
 *                                      update
 *   ion_radtransfer_init_3d        - handles internal initialization
 *   ion_radtransfer_init_domain_3d - handles internal initialization
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "ionrad.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef ION_RADIATION
/* Global storage arrays */
static Real ***ph_rate;            /* Photoionization rate */
static Real ***edot;               /* Rate of change of energy */
static Real ***nHdot;              /* Rate of change of neutral
				      density */
static int  ***last_sign;          /* Last sign of nHdot -- keep track
				      of this to avoid oscillatory
				      overstability */
static int  ***sign_count;         /* Number of successive times the
				      sign of nHdot has flipped -- use
				      this to avoid oscillatory
				      overstability */
static Real ***e_init;             /* Total energies on entry to
				      routine */
static Real ***e_th_init;          /* Thermal energies on entry to
				      routine */
static Real ***x_init;             /* Ionization fraction on entry to
				      routine */
static Real tcoarse = 0; /*Keep track of higher domain time step*/
/* ------------------------------------------------------------
 * Photoionization routines
 * ------------------------------------------------------------
 *
 * These routines do the work of computing the photoionization and
 * chemisty update.
 *
 */

/* Routine to zero out initial photoionization rates */
void ph_rate_init(GridS *pGrid)
{
  int i,j,k;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) ph_rate[k][j][i] = 0.0;
    }
  }
}

/* Routine to floor temperatures */
void apply_temp_floor(GridS *pGrid) {
  int i,j,k;
  Real e_sp, e_thermal, ke, T, x, n_H, n_Hplus, n_e;
#ifdef MHD
  Real be;
#endif

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Compute temperature */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	x = n_e / (n_H + n_Hplus);
	ke = 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d;
#ifdef MHD
	be = 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		    pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		    pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c);
#endif
	e_thermal = pGrid->U[k][j][i].E - ke;
#ifdef MHD
	e_thermal -= be; 
#endif
	e_sp = e_thermal / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*m_H+(1.0-x)*mu)/ k_B;

	if (T < tfloor) {
	  e_sp = tfloor * k_B / ((x*0.5*m_H+(1.0-x)*mu) * Gamma_1);
	  pGrid->U[k][j][i].E = 
	    0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		   pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		   pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
	    / pGrid->U[k][j][i].d
#ifdef MHD
	    + 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
		     pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
		     pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
	    + e_sp * pGrid->U[k][j][i].d;
	}

	if ((T > tceil) && (tceil > 0)) {
	  e_sp = tceil * k_B / ((x*0.5*m_H+(1.0-x)*mu) * Gamma_1);
	  pGrid->U[k][j][i].E = 
	    0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		   pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		   pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
	    / pGrid->U[k][j][i].d
#ifdef MHD
	    + 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
		     pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
		     pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
	    + e_sp * pGrid->U[k][j][i].d;
	}

      }
    }
  }
}


/* Routine to keep d_n > floor and d_n < d. */
void apply_neutral_floor(GridS *pGrid) {
  int i,j,k;
  Real d_nlim;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;
	if (pGrid->U[k][j][i].s[0] < d_nlim) {
	  pGrid->U[k][j][i].s[0] = d_nlim;
	} else if (pGrid->U[k][j][i].s[0] > pGrid->U[k][j][i].d) {
	  pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
	}
      }
    }
  }
}


/* Routine to save energy, thermal energy, and ionization fraction
   passed to routine -- used for time step constraint. */
void save_energy_and_x(GridS *pGrid)
{
  int i, j, k;
  Real e_thermal, n_H, n_Hplus, n_e, x;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Compute thermal energy */
	e_thermal = pGrid->U[k][j][i].E - 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d
#ifdef MHD
	  - 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		   pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		   pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	  ;

	/* Compute ion fraction */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	x = n_e / (n_H + n_Hplus);

	/* Save thermal and total energy, and neutral fraction */
	e_init[k][j][i] = pGrid->U[k][j][i].E;
	e_th_init[k][j][i] = e_thermal;
	x_init[k][j][i] = n_e / (n_H + n_Hplus);

	/* Initialize last_sign and sign_count arrays */
	last_sign[k][j][i] = 0;
	sign_count[k][j][i] = 0;
      }
    }
  }
}


/* Routine to check if we have changed the total energy, thermal
   energy, or x_n as much as we are allowed. */
int check_range(DomainS *pDomain) {
  GridS *pGrid = pDomain->Grid;
  int i, j, k;
  Real e_thermal, n_H, n_Hplus, n_e, x;
  long cellcount = 0;
#ifdef MPI_PARALLEL
  int err;
  long cellcount_glob = 0;
  MPI_Comm Comm_Domain = pDomain->Comm_Domain;
#endif

  /* Check thermal energy */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Check D type condition */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	if (ph_rate[k][j][i] / (min_area * n_H) > 2.0*CION) continue;

	/* Check thermal energy condition */
	if (max_de_therm_step > 0) {
	  e_thermal = pGrid->U[k][j][i].E - 0.5 *
	    (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	     pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	     pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	    / pGrid->U[k][j][i].d
#ifdef MHD
	    - 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		     pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		     pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	    ;
	}
	if ((e_thermal / e_th_init[k][j][i] >= 1 + max_de_therm_step) ||
	    (e_th_init[k][j][i] / e_thermal >= 1 + max_de_therm_step)) {
	  cellcount++;
	  continue;
	}

	/* Check total energy condition */
	if (max_de_step > 0) {
	  if ((pGrid->U[k][j][i].E / e_init[k][j][i] >= 1 + max_de_step) ||
	      (e_init[k][j][i] / pGrid->U[k][j][i].E >= 1 + max_de_step)) {
	    cellcount++;
	    continue;
	  }
	}

	/* Check neutral fraction condition */
	if (max_dx_step > 0) {
	  n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	  n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	  x = n_e / (n_H + n_Hplus);
	  if ((x / x_init[k][j][i] >= 1 + max_dx_step) ||
	      (x_init[k][j][i] / x >= 1 + max_dx_step)) {
	    cellcount++;
	    continue;
	  }
	  /*if (x / x_init[k][j][i] >= 1 + max_dx_step){ 
	    cellcount++;
	    continue;
	  }*/
	}
      }
    }
  }

#ifdef MPI_PARALLEL
  err = MPI_Allreduce(&cellcount, &cellcount_glob, 1, MPI_LONG, 
		      MPI_SUM, Comm_Domain);
  if (err) ath_error("[check_range]: MPI_Allreduce returned error code %d\n"
		    ,err);
  cellcount = cellcount_glob;
#endif
  if (cellcount > MAXCELLCOUNT) return(1);
  else return(0);
}


#define MAXSIGNCOUNT 4
#define DAMPFACTOR 0.5
Real compute_chem_rates(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, j, k, n;
  Real n_H, n_Hplus, n_e, d_nlim;
  Real e_sp, T, x;
  Real dt_chem, dt_chem1, dt_chem2, dt_chem_min;
  /* int dt_chem_min_index; */
#ifdef MPI_PARALLEL
  int err;
  Real dt_chem_min_glob;
  MPI_Comm Comm_Domain = pDomain->Comm_Domain;
#endif

  /* Initialize chemistry time step to large time step */
  dt_chem_min = LARGE;
  /* dt_chem_min_index = -999; */

  /* Loop over cells to get timestep */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {


	/* Get species abundances */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	x = n_e / (n_H + n_Hplus);

	/* Get gas temperature in K */
	e_sp = (pGrid->U[k][j][i].E -
		0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		       pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		       pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
		/ pGrid->U[k][j][i].d
#ifdef MHD
		- 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
			 pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
			 pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
		) / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*m_H+(1.0-x)*mu)/ k_B;
	if (T < tfloor) T = tfloor;

	/* Get rate of change of neutral density */
	nHdot[k][j][i] = 
	  recomb_rate_coef(T) * time_unit * n_e * n_Hplus
	  - ph_rate[k][j][i] * n_H;
	
	/* AT 2/23/14: Removing collisional ionization*/
	/* - coll_ion_rate_coef(T) * time_unit * n_e * n_H; */

	/* Check if the sign has flipped -- oscillatory overstability
	   check */
	if (nHdot[k][j][i] < 0.0) {
	  if (last_sign[k][j][i] == 1) sign_count[k][j][i]++;
	  else if (sign_count[k][j][i] > 0) sign_count[k][j][i]--;
	  last_sign[k][j][i] = -1;
	} else if (nHdot[k][j][i] > 0.0) {
	  if (last_sign[k][j][i] == -1) sign_count[k][j][i]++;
	  else if (sign_count[k][j][i] > 0) sign_count[k][j][i]--;
	  last_sign[k][j][i] = 1;
	} else {
	  sign_count[k][j][i] = last_sign[k][j][i] = 0;
	}

	/* If sign has flipped too many times successively, this cell
	 * is probably experiencing oscillatory overstability. To
	 * combat this, decrease the update amount by a chosen factor
	 * for every repeat of the same sign past the trigger number
	 */
	for (n=MAXSIGNCOUNT; n<sign_count[k][j][i]; n++) {
	  edot[k][j][i] *= DAMPFACTOR;
	  nHdot[k][j][i] *= DAMPFACTOR;
	}

	/* Get ionization fraction floor */
	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;

	/* Compute chemistry time step for this cell and find the min */
	if (nHdot[k][j][i] == 0.0) {
	  dt_chem1 = dt_chem2 = LARGE;
	} else if (nHdot[k][j][i] > 0.0) {
	  dt_chem1 = max_dx_iter / (1+max_dx_iter) * n_e / nHdot[k][j][i];
	  dt_chem2 = max_dx_iter * n_H / nHdot[k][j][i];
	} else if (pGrid->U[k][j][i].s[0] > 1.0001*d_nlim) {
	  dt_chem1 = -max_dx_iter * n_e / nHdot[k][j][i];
	  dt_chem2 = -max_dx_iter / (1+max_dx_iter) * n_H / nHdot[k][j][i];
	} else {
	  dt_chem1 = dt_chem2 = LARGE;
	}
	dt_chem = (dt_chem1 < dt_chem2) ? dt_chem1 : dt_chem2;
	/* dt_chem_min = (dt_chem < dt_chem_min) ? dt_chem : dt_chem_min; */
	if (dt_chem < dt_chem_min){
	  /* fprintf(stderr, "Cell i: %d j %d k %d, dt_chem: %e dt_chem_min: %e \n", i, j, k, dt_chem, dt_chem_min); */
	  dt_chem_min = dt_chem;
	  /* dt_chem_min_index = i; */
	}
	    
	if (dt_chem < 0) {
	  ath_error("[compute_chem_rates]: cell %d %d %d: dt_chem = %e, d = %e, d_n = %e, T = %e, nH = %e, nH+ = %e, ne = %e, nHdot = %e\n", i,j,k, dt_chem, pGrid->U[k][j][i].d, pGrid->U[k][j][i].s[0], T, n_H, n_Hplus, n_e, nHdot[k][j][i]);
	}
      }
    }
  }

  /* fprintf(stderr, "Min dt_chem at %d with left dens: %e, dens %e, right dens: %e \n", dt_chem_min_index, pGrid->U[4][4][dt_chem_min_index].d, pGrid->U[4][4][dt_chem_min_index-1].d, pGrid->U[4][4][dt_chem_min_index+1].d); */
#ifdef MPI_PARALLEL
  /* Sync chemistry timestep across processors */
  err = MPI_Allreduce(&dt_chem_min, &dt_chem_min_glob, 1, MP_RL, 
		      MPI_MIN, Comm_Domain);
  if(err) ath_error("[compute_chem_rates]: MPI_Allreduce returned error code %d\n"
		    ,err);
  dt_chem_min = dt_chem_min_glob;
#endif

  /*fclose(fp);*/

  return(dt_chem_min);
}
#undef MAXSIGNCOUNT
#undef DAMPFACTOR


Real compute_therm_rates(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, j, k;
  Real n_H, n_Hplus, n_e, e_thermal;
  Real e_sp, T, x, e_sp_min, e_th_min, e_min, d_nlim;
  Real dt_therm, dt_therm1, dt_therm2, dt_therm_min;
  /* int dt_therm_min_index; */
#ifdef MPI_PARALLEL
  int err;
  Real dt_therm_min_glob;
  MPI_Comm Comm_Domain = pDomain->Comm_Domain;
#endif

  /* Initialize thermal time step to large value */
  dt_therm_min = LARGE;
  /* dt_therm_min_index = -999; */

  /* Loop over cells to get timestep */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Get species abundances */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	x = n_e / (n_H + n_Hplus);

	/* Get gas temperature in K */
	e_thermal = pGrid->U[k][j][i].E - 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d
#ifdef MHD
	  - 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		   pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		   pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	  ;
	e_sp = e_thermal / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*m_H+(1.0-x)*mu)/ k_B;

	/* Check temperature floor. If cell is below temperature
	   floor, skip it. We'll fix it later */
	if (T < tfloor) {
	  edot[k][j][i] = 0.0;
	  continue;
	}

	/* If we're at the ionization floor and trying to ionize
	   further, we won't update the temperature, so skip this
	   cell. */
	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;
	if ((nHdot[k][j][i] < 0) &&
	    (pGrid->U[k][j][i].s[0] < 1.0001*d_nlim)) {
	  edot[k][j][i] = 0.0;
	  continue;
	}

	/* Get rate of change of gas energy. Only use molecular cooling
	   in cells with < COOLFRAC or > 1 - COOLFRAC ionization fraction, to
	   avoid artificial over-cooling in mixed or transition cells. */
	edot[k][j][i] = ph_rate[k][j][i] * e_gamma * n_H
	  /* AT 2/23/14: Removing metal cooling*/
	  /* - osterbrock_cool_rate(T) * n_e*n_Hplus */
	  - recomb_cool_rate_coef(T) * time_unit * n_Hplus * n_e
	  /* AT 3/13/14: Adding Lya cooling (which has a negative coefficient)*/
	  + lya_cool_rate(n_H, n_e, T) * time_unit;
	/* AT 2/23/14: Removing molecular terms*/
	/* if ((n_Hplus / (n_H+n_Hplus) < COOLFRAC) ||  */
	/*     (n_Hplus / (n_H+n_Hplus) > 1.0-COOLFRAC)) { */
	/*   edot[k][j][i] += ki_heat_rate() * time_unit * n_H */
	/*     - ki_cool_rate(T) * time_unit * n_H * n_H; */
	/* } */

	/* Compute thermal time step for this cell and find the
	   min. Note that, if we're cooling, we need to take into
	   account the effect of the floor. */
	if (edot[k][j][i] == 0.0) {
	  dt_therm1 = dt_therm2 = LARGE;
	} else if (edot[k][j][i] > 0.0) {

	  /* We're heating, no need to consider floor */
	  dt_therm1 = max_de_iter * pGrid->U[k][j][i].E / edot[k][j][i];
	  dt_therm2 = max_de_therm_iter * e_thermal / edot[k][j][i];

	} else {

	  /* We're cooling. Start by computing the total and thermal
	     energy the gas would have if it were at the temperature
	     floor. */
	  e_sp_min = tfloor * k_B / ((x*0.5*m_H+(1.0-x)*mu) * Gamma_1);
	  e_th_min = e_sp_min * pGrid->U[k][j][i].d;
	  e_min = 
	    0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		   pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		   pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
	    / pGrid->U[k][j][i].d
#ifdef MHD
	    + 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
		     pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
		     pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
	    + e_th_min;

	  /* If cooling to the temperature floor would not violate the
	     constraint on the maximum allowable change in either
	     total or thermal energy, there is no constraint, and we
	     can continue to the next cell */
	  if ((e_thermal/(1.0+max_de_therm_iter) < e_th_min) &&
	      (pGrid->U[k][j][i].E/(1.0+max_de_iter) < e_min))
	    continue;

	  /* If we're here, cooling to the temperature floor would
	     violate our time step constraint. Therefore compute the
	     time step normally. */
	  dt_therm1 = -max_de_iter / (1+max_de_iter) * pGrid->U[k][j][i].E 
	    / edot[k][j][i];
	  dt_therm2 = -max_de_therm_iter / (1+max_de_therm_iter) * e_thermal /
	    edot[k][j][i];

	}

	/* Set time step to minimum */
	dt_therm = (dt_therm1 < dt_therm2) ? dt_therm1 : dt_therm2;
	/* dt_therm_min = (dt_therm < dt_therm_min) ? dt_therm : dt_therm_min; */
	if (dt_therm < dt_therm_min){
	  /* fprintf(stderr, "Cell i: %d j %d k %d, dt_therm: %e dt_therm_min: %e \n", i, j, k, dt_therm, dt_therm_min); */
	  dt_therm_min = dt_therm;
	  /* dt_therm_min_index = i; */
	}
      }
    }
  }
  /* fprintf(stderr, "Min dt_therm at %d \n", dt_therm_min_index); */
#ifdef MPI_PARALLEL
  /* Sync thermal timestep across processors */
  err = MPI_Allreduce(&dt_therm_min, &dt_therm_min_glob, 1, MP_RL, 
		      MPI_MIN, Comm_Domain);
  if(err) ath_error("[compute_therm_rates]: MPI_Allreduce returned error code %d\n"
		    ,err);
  dt_therm_min = dt_therm_min_glob;
#endif

  return(dt_therm_min);
}


void ionization_update(GridS *pGrid, Real dt)
{
  int i, j, k;
  Real d_nlim;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;

	if ((nHdot[k][j][i] > 0) ||
	    (pGrid->U[k][j][i].s[0] > 1.0001*d_nlim)) {

	  /* Update gas energy */
	  pGrid->U[k][j][i].E += edot[k][j][i] * dt;

	  /* Update neutral density */
	  pGrid->U[k][j][i].s[0] += nHdot[k][j][i] * dt * m_H;

	}
      }
    }
  }
}


Real compute_dt_hydro(DomainS *pDomain) {
  GridS *pGrid = pDomain->Grid;
  int i,j,k;
  Real di,v1,v2,v3,qsq,p,asq,cf1sq,cf2sq,cf3sq,max_dti=0.0,dt;
#ifdef MHD
  Real b1,b2,b3,bsq,tsum,tdif;
#endif /* MHD */
#ifdef MPI_PARALLEL
  MPI_Comm Comm_Domain = pDomain->Comm_Domain;
  Real dt_glob;
  int err;
#endif /* MPI_PARALLEL */

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	di = 1.0/(pGrid->U[k][j][i].d);
	v1 = pGrid->U[k][j][i].M1*di;
	v2 = pGrid->U[k][j][i].M2*di;
	v3 = pGrid->U[k][j][i].M3*di;
	qsq = v1*v1 + v2*v2 + v3*v3;

#ifdef MHD
	/* Use maximum of face-centered fields (always larger than
	   cell-centered B) */
	b1 = pGrid->U[k][j][i].B1c 
	  + fabs((double)(pGrid->B1i[k][j][i] - pGrid->U[k][j][i].B1c));
	b2 = pGrid->U[k][j][i].B2c 
	  + fabs((double)(pGrid->B2i[k][j][i] - pGrid->U[k][j][i].B2c));
	b3 = pGrid->U[k][j][i].B3c 
	  + fabs((double)(pGrid->B3i[k][j][i] - pGrid->U[k][j][i].B3c));
	bsq = b1*b1 + b2*b2 + b3*b3;
	/* compute sound speed squared */
#ifdef ADIABATIC
	p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq
			 - 0.5*bsq), TINY_NUMBER);
	asq = Gamma*p*di;
#else
	asq = Iso_csound2;
#endif /* ADIABATIC */
	/* compute fast magnetosonic speed squared in each direction */
	tsum = bsq*di + asq;
	tdif = bsq*di - asq;
	cf1sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b2*b2+b3*b3)*di));
	cf2sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b3*b3)*di));
	cf3sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b2*b2)*di));
	
#else /* MHD */
	
	/* compute sound speed squared */
#ifdef ADIABATIC
	p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq),
		TINY_NUMBER);
	asq = Gamma*p*di;
#else
	asq = Iso_csound2;
#endif /* ADIABATIC */
	/* compute fast magnetosonic speed squared in each direction */
	cf1sq = asq;
	cf2sq = asq;
	cf3sq = asq;
	
#endif /* MHD */

	/* compute maximum inverse of dt (corresponding to minimum dt) */
	if (pGrid->Nx[0] > 1)
	  max_dti = MAX(max_dti,(fabs(v1)+sqrt((double)cf1sq))/pGrid->dx1);
	if (pGrid->Nx[1] > 1)
	  max_dti = MAX(max_dti,(fabs(v2)+sqrt((double)cf2sq))/pGrid->dx2);
	if (pGrid->Nx[2] > 1)
	  max_dti = MAX(max_dti,(fabs(v3)+sqrt((double)cf3sq))/pGrid->dx3);
      }
    }
  }

  /* get timestep. */
  dt = CourNo/max_dti;

#ifdef MPI_PARALLEL
  err = MPI_Allreduce(&dt, &dt_glob, 1, MP_RL, MPI_MIN, Comm_Domain);
  if(err) ath_error("[compute_dt_hydro]: MPI_Allreduce returned error code %d\n"
		    ,err);
  dt = dt_glob;
#endif /* MPI_PARALLEL */
  return(dt);
}


/* Routine to manually set energy passed to hydro solver
 to test if the ionization fraction is working*/
void set_energy_manually(GridS *pGrid)
{
  int i, j, k;
  Real e_nonthermal, n_H, n_Hplus, n_e, x;
  Real Tcalc, mucalc;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Compute kinetic energy */
	e_nonthermal = 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d
#ifdef MHD
	  + 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		   pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		   pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	  ;

	/* Compute ion fraction */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	x = n_e / (n_H + n_Hplus);

	/* Save thermal and total energy, and neutral fraction */
	if (x <= .5) {
	  /*Neutral gas*/
	  Tcalc = 55. ;
	  mucalc = mu;
	} else {
	  /*Ionized gas*/
	  Tcalc=10000.;
	  mucalc = m_H/2.;
	}
	pGrid->U[k][j][i].E = e_nonthermal + k_B*Tcalc/Gamma_1/mucalc*pGrid->U[k][j][i].d;
      }
    }
  }
}


/* ------------------------------------------------------------
 * Initialize and store routines
 * ------------------------------------------------------------
 *
 * These are called at problem setup to allocate memory and read
 * setup parameters, or when checkpointing to dump internal data
 * needed on restart.
 *
 */

void ion_radtransfer_init_3d(GridS *pGrid, DomainS *pDomain, int ires, int sizei, int sizej, int sizek) {

  /* Read input values  */
  sigma_ph = par_getd("ionradiation", "sigma_ph");
  m_H = par_getd("ionradiation", "m_H");
  mu = par_getd("ionradiation", "mu");
  e_gamma = par_getd("ionradiation", "e_gamma");
  alpha_C = par_getd("ionradiation", "alpha_C");
  k_B = par_getd("ionradiation", "k_B");
  time_unit = par_getd("ionradiation", "time_unit");
  max_de_iter = par_getd("ionradiation", "max_de_iter");
  max_de_therm_iter = par_getd("ionradiation", "max_de_therm_iter");
  max_dx_iter = par_getd("ionradiation", "max_dx_iter");
  max_de_step = par_getd("ionradiation", "max_de_step");
  max_de_therm_step = par_getd("ionradiation", "max_de_therm_step");
  max_dx_step = par_getd("ionradiation", "max_dx_step");
  tfloor = par_getd("ionradiation", "tfloor");
  tceil = par_getd("ionradiation", "tceil");
  maxiter = par_getd("ionradiation", "maxiter");

  sizei = sizei + 2*nghost;
  sizej = sizej + 2*nghost;
  sizek = sizek + 2*nghost;

  /* Allocate memory for rate arrays */
  ph_rate = (Real***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  fprintf(stderr, "I'm level :%d, sized k: %d j:%d i:%d \n", pDomain->Level, sizek, sizej, sizei);
  edot = (Real***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  nHdot = (Real***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  last_sign = (int***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  sign_count = (int***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  e_init = (Real***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  e_th_init = (Real***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));
  x_init = (Real***) 
    calloc_3d_array(sizek, sizej, sizei, sizeof(Real));

  /* Offset pointers to account for ghost cells */
  /* ph_rate -= pGrid->ks; */
  /* edot -= pGrid->ks; */
  /* nHdot -= pGrid->ks; */
  /* last_sign -= pGrid->ks; */
  /* sign_count -= pGrid->ks; */
  /* e_init -= pGrid->ks; */
  /* e_th_init -= pGrid->ks; */
  /* x_init -= pGrid->ks; */
  /* for (k=pGrid->ks; k<=pGrid->ke; k++) { */
  /*   ph_rate[k] -= pGrid->js; */
  /*   edot[k] -= pGrid->js; */
  /*   nHdot[k] -= pGrid->js; */
  /*   last_sign[k] -= pGrid->js; */
  /*   sign_count[k] -= pGrid->js; */
  /*   e_init[k] -= pGrid->js; */
  /*   e_th_init[k] -= pGrid->js; */
  /*   x_init[k] -= pGrid->js; */
  /*   for (j=pGrid->js; j<=pGrid->je; j++) { */
  /*     ph_rate[k][j] -= pGrid->is; */
  /*     edot[k][j] -= pGrid->is; */
  /*     nHdot[k][j] -= pGrid->is; */
  /*     last_sign[k][j] -= pGrid->is; */
  /*     sign_count[k][j] -= pGrid->is; */
  /*     e_init[k][j] -= pGrid->is; */
  /*     e_th_init[k][j] -= pGrid->is; */
  /*     x_init[k][j] -= pGrid->is; */
  /*   } */
  /* } */


  return;
}


void ion_radtransfer_init_domain_3d(GridS *pGrid, DomainS *pDomain) {

  /*A Tripathi 06/01/12: CHECK to see if this is correct and/or necessary*/
  /* Store parallel grid information for internal use */

/*AT 3/5/14: Commenting out these calls to Ngrid. Local defs now in ionradplane_3d.c*/
/* #ifdef MPI_PARALLEL */
/*   pD = pDomain; */
/*   NGrid_x1 = pDomain->NGrid[0]; */
/*   NGrid_x2 = pDomain->NGrid[1]; */
/*   NGrid_x3 = pDomain->NGrid[2]; */
/* #endif */

  /* Store information specific to point sources and planes */
#ifdef ION_RADPLANE
  ion_radplane_init_domain_3d(pGrid, pDomain);
#endif
}

void set_coarse_time(){
#ifdef MPI_PARALLEL
  int err;
  Real tcoarse_glob;
  err = MPI_Allreduce(&tcoarse, &tcoarse_glob, 1, MP_RL, MPI_MAX, MPI_COMM_WORLD);
  tcoarse = tcoarse_glob;
  return;
#endif
}

void clear_coarse_time(){
#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  tcoarse = 0;
  return;
#endif
}

/* ------------------------------------------------------------
 * Main integration routine
 * ------------------------------------------------------------
 *
 * This is the driver routine for the radiation integration step.
 *
 */

void ion_radtransfer_3d(DomainS *pDomain) 
{
  MeshS *pMesh = pDomain->Mesh;
  GridS *pGrid = pDomain->Grid;
  Real dt_chem, dt_therm, dt_hydro, dt, dt_done;
  int n, niter, hydro_done;
  int nchem, ntherm;
  int finegrid, coarsetime_done;
  int dir, dim;
  finegrid = 0;
  niter = 0;

  /* fprintf(stderr,"My level is %d. My time is %f [grid] %f [mesh] \n ", pDomain->Level, pMesh->time, pGrid->time); */

  /*Set the direction of propagation to that of the first radiation source*/
  /*AT 9/24/12: Make this dir not be hardwired to dir[0] and also consistent with the get_ph_rate_plane call*/
  dir = (pMesh->radplanelist)->dir[0];
  dim = (dir < 0) ? 2*(fabs(dir) - 1): 2*dir - 1;

  /*Set the finegrid flag if on level number > 0*/
  if(pDomain->Level != 0) finegrid = 1;

#ifdef STATIC_MESH_REFINEMENT
  if (finegrid) { 
    /*If not on root domain, go to the receive function.  Send call in line 1046*/
    ionrad_prolong_rcv(pGrid, dim, pDomain->Level, pDomain->DomNumber);
  }
  else { 
    /*AT 4/3/13: This step may be redundant, given the existence of the clear_coarse_time called in the main.*/
    tcoarse = 0;
  } /*AT 11/19/12: Will need to fix placement of calls to be valid for MPI+-SMR*/
#endif

  /* Set all temperatures below the floor to the floor */
  apply_temp_floor(pGrid);

  /* Set neutral densities below floor to sensible values. This is necessary
     because the hydro update can make the neutral density negative. */
  apply_neutral_floor(pGrid);

  /* Save total and thermal energies and neutral fractions passed in
     -- we use these to compute time steps. Also initialize the last_sign
     and sign_count arrays */
  save_energy_and_x(pGrid);

  /* Begin the radiation sub-cycle */
  dt_done = 0.0;
  hydro_done = 0;
  /*Flag to check that the timestep of the root level has not been exceded by the fine grids*/
  coarsetime_done = 0; 
  nchem = ntherm = 0;

  /*Do radiation sub cycle differently depending on whether on a coarse or fine grid*/

  /*If on the coarsest level, run under the regular stopping conditions condition*/
  /*If on a finer level, run under the time condition*/
  /*This ONLY treats the root level as special.*/
  while(finegrid || !hydro_done){
    
    /* Initialize photoionization rate array */
    ph_rate_init(pGrid);

    /* Compute photoionization rate from all sources */
    /* AT 4/3/13: Again, remove the for loop if there's only going to be 1 source and adjust radplane structure*/
#ifdef ION_RADPLANE
    for (n=0; n<(pMesh->radplanelist)->nradplane; n++) 
      {
	get_ph_rate_plane((pMesh->radplanelist)->flux_i,(pMesh->radplanelist)->dir[n],ph_rate, pDomain);
      }
#endif


    /* Compute rates and time step for chemistry update */
    dt_chem = compute_chem_rates(pDomain);

    /* Compute rates and time step for thermal energy update */
    dt_therm = compute_therm_rates(pDomain);


    /* Set time step to smaller of thermal and chemical time
       steps, and record whether this is a thermal or chemical step */
    if (dt_chem < dt_therm) nchem++;
    else ntherm++;
    dt = MIN(dt_therm, dt_chem);


    /* If necessary, scale back time step to avoid exceeding hydro
       time step. */
    if (!finegrid){
      if (dt_done + dt > pGrid->dt) {
	dt = pGrid->dt - dt_done;
	hydro_done = 1;
      }
    /* If necessary and on a fine grid, scale back time step to avoid 
       exceeding coarse time step. */
    } else {
       if (dt_done + dt >tcoarse) {
	dt = tcoarse - dt_done;
	coarsetime_done = 1;
       }
    }

    /* Do an update */
    ionization_update(pGrid, dt);
    dt_done += dt;
    niter++;

    /* Set all temperatures below the floor to the floor */
    apply_temp_floor(pGrid);
    apply_neutral_floor(pGrid);

    /*Manually set thermal energy due to ionization*/
    /* set_energy_manually(pGrid); */

    /*Check stopping criteria for coarse grid */
    if (!finegrid){

      /* Check new energies and ionization fractions against initial
	 values to see if we've changed them as much as possible. If so,
	 exit loop. */
      if (check_range(pDomain)) {
	pGrid->dt = dt_done;
	/* fprintf(stderr,"In check range \n"); */
	break;
      }

      /* Have we advanced the full hydro time step? If so, exit loop. */
      if (hydro_done) {
	/* fprintf(stderr,"In hydro done \n"); */
	break;
      }
      
      /* Compute a new hydro time step based on the new temperature
	 distribution. If it's smaller than the time step we've already 
	 advanced, then exit. */
      dt_hydro = compute_dt_hydro(pDomain);
      if (dt_hydro < dt_done) {
      	/* fprintf(stderr,"dt_hydro %e dt done %e \n", dt_hydro, dt_done); */
      	pGrid->dt = dt_done;
      	break;
      }

    } else {
      /* if (check_range(pDomain)) {fprintf(stderr, "I've hit check range! on domn: %d \n", pDomain->DomNumber);} */

      /*Check time to stop fine grid */
      if (coarsetime_done) {
	pGrid->dt = dt_done;
/*     fprintf(stderr,"Exceeding coarse time limit \n"); */
	break;
      }
    }
  }
  

  if (!finegrid){

  /* Have we exceeded the maximum number of iterations? If so, return
     to hydro with the time step re-set to what we managed to do and
     return this time to the finer levels. */    
    if (niter==maxiter) {
	pGrid->dt = dt_done;
	/* fprintf(stderr,"Reached maxiter \n"); */
    }

    /*Set coarsetime to the timestep done on this root grid*/
    tcoarse = dt_done;
  }

    /*Set mesh timestep equal to grid timestep*/
    /*AT:4/15/13: Moved outside of !finegrid conditional*/
    pMesh->dt = pGrid->dt;

#ifdef STATIC_MESH_REFINEMENT
  /*Send radiation flux to finer grids that overlap*/
  ionrad_prolong_snd(pGrid, dim, pDomain->Level, pDomain->DomNumber);
#endif

  /* Write status */
  fprintf(stderr, "Radiation done in %d iterations: %d thermal, %d chemical; new dt = %e\n", niter, ntherm, nchem, pGrid->dt);

  /* Sanity check */
  if (!finegrid && (pGrid->dt < 0)) {
    ath_error("[ion_radtransfer_3d]: dt = %e, dt_chem = %e, dt_therm = %e, dt_hydro = %e, dt_done = %e\n", pGrid->dt, dt_chem, dt_therm, dt_hydro, dt_done);
  }
}

#endif /* ION_RADIATION */
