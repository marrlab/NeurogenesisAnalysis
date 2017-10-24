#ifndef _MY__sA_sPB_3o_2B_S_dtS__T_astS__B_atS
#define _MY__sA_sPB_3o_2B_S_dtS__T_astS__B_atS

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <udata.h>
#include <math.h>
#include <mex.h>

             void x0__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(N_Vector x0, void *user_data);
             int J__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int JSparse__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(realtype t, N_Vector x,N_Vector fx, SlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int JBand__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(long int N, long int mupper, long int mlower, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int JB__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(long int NeqB,realtype t, N_Vector x, N_Vector xB, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int JSparseB__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(realtype t, N_Vector x, N_Vector xB, N_Vector fx, SlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int JBBand__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int Jv__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp);
             int JBv__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp);
             int sx__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2);
             void sx0__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(int ip, N_Vector sx0, void *user_data);
             void y__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x);
             void dydp__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny);
             void dydx__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, double *dydx, double *y, double *p, double *k, double *x);
             void sy__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx);
             double sroot__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data);
             double srootval__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data);
             double s2root__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data);
             double s2rootval__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data);
             void deltadisc__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int idisc, N_Vector x, void *user_data);
             void sdeltadisc__sA_sPB_3o_2B_S_dtS__T_astS__B_atS(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data);


#endif /* _MY__sA_sPB_3o_2B_S_dtS__T_astS__B_atS */
