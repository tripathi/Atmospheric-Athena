/*
 * Function testsedov.c
 *
 * Problem generator for producing a linear blast wave (1-D Sedov Taylor)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real temprat, n_H, m_H, T;
  Real rho, pressure;
  Real kb = 1.3806/1.e16;

/* Set up uniform ambient medium
 * with left edge set to higher temperature
 */

  n_H = par_getd("problem","n_H");
  m_H = par_getd("problem", "m_H");
  T = par_getd("problem", "temperature");
  temprat = par_getd("problem","tempratio");

  /* Power-law pressure and density */
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
	rho = n_H * m_H;
	pGrid->U[k][j][i].d  = rho;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;
	if (i <= is+2)
	  {
	    pGrid->U[k][j][i].E  = (n_H * kb * T*temprat)/Gamma_1;
	  } else
	  {
	    pGrid->U[k][j][i].E  = n_H * kb * T/Gamma_1;
	  }
      }
    }
  }
  return;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}


void problem_write_restart(MeshS *pM, FILE *fp){
  return;
}


void problem_read_restart(MeshS *pM, FILE *fp){
  return;
}


ConsFun_t get_usr_expr(const char *expr){
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}
