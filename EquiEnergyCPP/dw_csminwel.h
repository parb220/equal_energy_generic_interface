/*
 * Copyright (C) 1996 Christopher Sims
 * Copyright (C) 2003 Karibzhanov, Waggoner, and Zha
 *
 * This free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * If you did not received a copy of the GNU General Public License
 * with this software, see <http://www.gnu.org/licenses/>.
 */

#ifndef __CSMINWEL_H__
#define __CSMINWEL_H__

// Added by DW to eliminate optpackage.h
#define INI_H_CSMINWEL   1.0e-005  //Initial value for the diagonal of inverse Hessian in the quasi-Newton search.
                                  //1.0e-05 (sometimes used for SargentWZ USinflation project I)
                                  //5.0e-04 (for monthly TVBAR)

//--- This extern variable allows an input by the user from an input data file.
extern double GRADSTPS_CSMINWEL;

void dw_csminwel(double (*fcn)(double *x, int n, double **args, int *dims),
            double *x, int n, double *H, double *gh,
            int (*grad)(double *x, int n, double *g, double **args, int *dims),
            double *fh, double crit, int *itct, int nit,
            int *fcount, int *retcodeh, double **args, int *dims);

void dw_csminwel_SetPrintFile(char *filename);
int dw_csminwel_randomseedChanged(int seednumber);

#endif
