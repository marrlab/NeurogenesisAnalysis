#include "_sA_sPB_3o_2B_S_atS__T_dtS__B_atS.h"
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
#include <symbolic_functions.c>
#include <udata.h>
#include <math.h>
#include <mex.h>

#define pi 3.141592653589793


 int xdot__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*35);
xdot_tmp[0] = 1.705*p[2]*(pow(p[3],2))*x_tmp[29] + 2.2*p[2]*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[2]*p[3]*x_tmp[33];
xdot_tmp[1] = p[4]*x_tmp[4] - 1.0*p[5]*x_tmp[1];
xdot_tmp[2] = p[5]*x_tmp[1] - 2.0*p[5]*x_tmp[2] + p[4]*x_tmp[4] + 2.0*p[4]*x_tmp[6];
xdot_tmp[3] = 1.1*p[2]*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[29];
xdot_tmp[4] = p[2]*x_tmp[3] - 1.0*p[4]*x_tmp[4] + 0.9*p[2]*(pow(p[3],2))*x_tmp[29] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[29];
xdot_tmp[5] = p[4]*x_tmp[8] - 1.0*p[5]*x_tmp[5] + 1.1*p[2]*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[9];
xdot_tmp[6] = p[2]*x_tmp[5] - 1.0*p[4]*x_tmp[4] - 1.0*p[4]*x_tmp[6] + p[4]*x_tmp[7] - 1.0*p[5]*x_tmp[6] + 0.9*p[2]*(pow(p[3],2))*x_tmp[9] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[9];
xdot_tmp[7] = p[2]*x_tmp[3] + p[4]*x_tmp[4] + 2.0*p[2]*x_tmp[8] - 2.0*p[4]*x_tmp[7] + 1.305*p[2]*(pow(p[3],2))*x_tmp[29] + 1.8*p[2]*(pow(p[3],2))*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[2]*p[3]*x_tmp[34];
xdot_tmp[8] = p[2]*x_tmp[0] - 1.0*p[4]*x_tmp[8] + 0.495*p[2]*(pow(p[3],2))*x_tmp[29] + 0.9*p[2]*(pow(p[3],2))*x_tmp[33] + 1.1*p[2]*(pow(p[3],2))*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[34];
xdot_tmp[9] = p[2]*x_tmp[31] - 1.0*p[5]*x_tmp[9] + p[4]*x_tmp[34] - 1.0*p[2]*(pow(p[3],2))*x_tmp[9] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[9];
xdot_tmp[10] = p[0]*x_tmp[26] - 1.0*p[1]*x_tmp[10] + 1.1*p[2]*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[30];
xdot_tmp[11] = 1.1*p[2]*(pow(p[3],2))*x_tmp[24] - 0.0002*x_tmp[11] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[24];
xdot_tmp[12] = p[4]*x_tmp[14] - 1.0*p[5]*x_tmp[12] - 0.0002*x_tmp[12];
xdot_tmp[13] = p[4]*x_tmp[25] - 1.0*p[5]*x_tmp[13] - 1.0*p[0]*x_tmp[13] + p[1]*x_tmp[31] + 0.0002*x_tmp[12];
xdot_tmp[14] = p[2]*x_tmp[11] - 1.0*p[4]*x_tmp[14] - 0.0002*x_tmp[14] + 0.9*p[2]*(pow(p[3],2))*x_tmp[24] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[24];
xdot_tmp[15] = 0.0002*x_tmp[17] - 0.0004*x_tmp[15];
xdot_tmp[16] = p[1]*x_tmp[18] - 1.0*p[0]*x_tmp[16] + 0.0002*x_tmp[15] - 0.0002*x_tmp[16] - 0.0002*x_tmp[17];
xdot_tmp[17] = -0.0002*x_tmp[17];
xdot_tmp[18] = p[0]*x_tmp[16] - 1.0*p[1]*x_tmp[18] - 0.0002*x_tmp[18];
xdot_tmp[19] = p[0]*x_tmp[20] - 2.0*p[0]*x_tmp[19] + 2.0*p[1]*x_tmp[22] + p[1]*x_tmp[28] + 0.0004*x_tmp[16] + 0.0002*x_tmp[17];
xdot_tmp[20] = p[1]*x_tmp[28] - 1.0*p[0]*x_tmp[20] + 0.0002*x_tmp[17];
xdot_tmp[21] = p[0]*x_tmp[20] + 2.0*p[0]*x_tmp[22] - 2.0*p[1]*x_tmp[21] + p[1]*x_tmp[28];
xdot_tmp[22] = p[0]*x_tmp[19] - 1.0*p[0]*x_tmp[20] - 1.0*p[0]*x_tmp[22] + p[1]*x_tmp[21] - 1.0*p[1]*x_tmp[22] - 1.0*p[1]*x_tmp[28] + 0.0002*x_tmp[18];
xdot_tmp[23] = p[2]*x_tmp[22] - 1.0*p[0]*x_tmp[23] + p[1]*x_tmp[30] + 0.0002*x_tmp[24] - 1.0*p[2]*(pow(p[3],2))*x_tmp[23] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[23];
xdot_tmp[24] = p[2]*x_tmp[18] - 0.0002*x_tmp[24] - 1.0*p[2]*(pow(p[3],2))*x_tmp[24] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[24];
xdot_tmp[25] = p[1]*x_tmp[27] - 1.0*p[0]*x_tmp[25] + p[2]*x_tmp[26] - 1.0*p[4]*x_tmp[25] + 0.0002*x_tmp[14] + 0.9*p[2]*(pow(p[3],2))*x_tmp[23] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[23];
xdot_tmp[26] = p[1]*x_tmp[10] - 1.0*p[0]*x_tmp[26] + 0.0002*x_tmp[11] + 1.1*p[2]*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[23];
xdot_tmp[27] = p[2]*x_tmp[10] + p[0]*x_tmp[25] - 1.0*p[1]*x_tmp[27] - 1.0*p[4]*x_tmp[27] + 0.9*p[2]*(pow(p[3],2))*x_tmp[30] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[30];
xdot_tmp[28] = p[0]*x_tmp[20] - 1.0*p[1]*x_tmp[28];
xdot_tmp[29] = p[2]*x_tmp[28] - 1.0*p[2]*(pow(p[3],2))*x_tmp[29] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[29];
xdot_tmp[30] = p[0]*x_tmp[23] + p[2]*x_tmp[21] - 1.0*p[1]*x_tmp[30] - 1.0*p[2]*(pow(p[3],2))*x_tmp[30] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[30];
xdot_tmp[31] = p[0]*x_tmp[13] + p[4]*x_tmp[27] - 1.0*p[1]*x_tmp[31] - 1.0*p[5]*x_tmp[31];
xdot_tmp[32] = p[2]*x_tmp[28] + 2.0*p[2]*x_tmp[30] + p[2]*(pow(p[3],2))*x_tmp[29] - 2.0*p[2]*(pow(p[3],2))*x_tmp[32] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*p[2]*x_tmp[32];
xdot_tmp[33] = p[2]*x_tmp[10] - 1.1*p[2]*(pow(p[3],2))*x_tmp[29] + 1.1*p[2]*(pow(p[3],2))*x_tmp[32] - 1.0*p[2]*(pow(p[3],2))*x_tmp[33] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[2]*p[3]*x_tmp[32];
xdot_tmp[34] = p[2]*x_tmp[27] + p[2]*x_tmp[33] - 1.0*p[4]*x_tmp[34] - 0.9*p[2]*(pow(p[3],2))*x_tmp[29] + 0.9*p[2]*(pow(p[3],2))*x_tmp[32] - 1.0*p[2]*(pow(p[3],2))*x_tmp[34] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*p[3]*x_tmp[32];

  for (ix=0; ix<35; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*315);
xBdot_tmp[0] = -1.0*p[2]*xB_tmp[8];
xBdot_tmp[1] = p[5]*xB_tmp[1] - 1.0*p[5]*xB_tmp[2];
xBdot_tmp[2] = 2.0*p[5]*xB_tmp[2];
xBdot_tmp[3] = - 1.0*p[2]*xB_tmp[4] - 1.0*p[2]*xB_tmp[7];
xBdot_tmp[4] = p[4]*xB_tmp[4] - 1.0*p[4]*xB_tmp[2] - 1.0*p[4]*xB_tmp[1] + p[4]*xB_tmp[6] - 1.0*p[4]*xB_tmp[7];
xBdot_tmp[5] = p[5]*xB_tmp[5] - 1.0*p[2]*xB_tmp[6];
xBdot_tmp[6] = xB_tmp[6]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[2];
xBdot_tmp[7] = 2.0*p[4]*xB_tmp[7] - 1.0*p[4]*xB_tmp[6];
xBdot_tmp[8] = p[4]*xB_tmp[8] - 1.0*p[4]*xB_tmp[5] - 2.0*p[2]*xB_tmp[7];
xBdot_tmp[9] = xB_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[6] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[5];
xBdot_tmp[10] = p[1]*xB_tmp[10] - 1.0*p[1]*xB_tmp[26] - 1.0*p[2]*xB_tmp[27] - 1.0*p[2]*xB_tmp[33];
xBdot_tmp[11] = 0.0002*xB_tmp[11] - 1.0*p[2]*xB_tmp[14] - 0.0002*xB_tmp[26];
xBdot_tmp[12] = (p[5] + 0.0002)*xB_tmp[12] - 0.0002*xB_tmp[13];
xBdot_tmp[13] = xB_tmp[13]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[31];
xBdot_tmp[14] = (p[4] + 0.0002)*xB_tmp[14] - 1.0*p[4]*xB_tmp[12] - 0.0002*xB_tmp[25];
xBdot_tmp[15] = 0.0004*xB_tmp[15] - 0.0002*xB_tmp[16];
xBdot_tmp[16] = (p[0] + 0.0002)*xB_tmp[16] - 1.0*p[0]*xB_tmp[18] - 0.0004*xB_tmp[19];
xBdot_tmp[17] = 0.0002*xB_tmp[16] - 0.0002*xB_tmp[15] + 0.0002*xB_tmp[17] - 0.0002*xB_tmp[19] - 0.0002*xB_tmp[20];
xBdot_tmp[18] = (p[1] + 0.0002)*xB_tmp[18] - 1.0*p[2]*xB_tmp[24] - 1.0*p[1]*xB_tmp[16] - 0.0002*xB_tmp[22];
xBdot_tmp[19] = 2.0*p[0]*xB_tmp[19] - 1.0*p[0]*xB_tmp[22];
xBdot_tmp[20] = p[0]*xB_tmp[20] - 1.0*p[0]*xB_tmp[19] - 1.0*p[0]*xB_tmp[21] + p[0]*xB_tmp[22] - 1.0*p[0]*xB_tmp[28];
xBdot_tmp[21] = 2.0*p[1]*xB_tmp[21] - 1.0*p[1]*xB_tmp[22] - 1.0*p[2]*xB_tmp[30];
xBdot_tmp[22] = xB_tmp[22]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[21] - 1.0*p[2]*xB_tmp[23] - 2.0*p[1]*xB_tmp[19];
xBdot_tmp[23] = xB_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[30] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[25] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[26];
xBdot_tmp[24] = xB_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[14] - 0.0002*xB_tmp[23] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[11];
xBdot_tmp[25] = xB_tmp[25]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[27] - 1.0*p[4]*xB_tmp[13];
xBdot_tmp[26] = p[0]*xB_tmp[26] - 1.0*p[0]*xB_tmp[10] - 1.0*p[2]*xB_tmp[25];
xBdot_tmp[27] = xB_tmp[27]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[31] - 1.0*p[2]*xB_tmp[34] - 1.0*p[1]*xB_tmp[25];
xBdot_tmp[28] = p[1]*xB_tmp[22] - 1.0*p[1]*xB_tmp[20] - 1.0*p[1]*xB_tmp[21] - 1.0*p[1]*xB_tmp[19] + p[1]*xB_tmp[28] - 1.0*p[2]*xB_tmp[29] - 1.0*p[2]*xB_tmp[32];
xBdot_tmp[29] = xB_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[7]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[0]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[4] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[3] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[32] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[8] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[33] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[34];
xBdot_tmp[30] = xB_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[23] - 2.0*p[2]*xB_tmp[32] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[27] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[10];
xBdot_tmp[31] = xB_tmp[31]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[13] - 1.0*p[2]*xB_tmp[9];
xBdot_tmp[32] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[32] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[34] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[33];
xBdot_tmp[33] = xB_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[34] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[8] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[0];
xBdot_tmp[34] = xB_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[9] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[7] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[8];
xBdot_tmp[35] = -1.0*p[2]*xB_tmp[43];
xBdot_tmp[36] = p[5]*xB_tmp[36] - 1.0*p[5]*xB_tmp[37];
xBdot_tmp[37] = 2.0*p[5]*xB_tmp[37];
xBdot_tmp[38] = - 1.0*p[2]*xB_tmp[39] - 1.0*p[2]*xB_tmp[42];
xBdot_tmp[39] = p[4]*xB_tmp[39] - 1.0*p[4]*xB_tmp[37] - 1.0*p[4]*xB_tmp[36] + p[4]*xB_tmp[41] - 1.0*p[4]*xB_tmp[42];
xBdot_tmp[40] = p[5]*xB_tmp[40] - 1.0*p[2]*xB_tmp[41];
xBdot_tmp[41] = xB_tmp[41]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[37];
xBdot_tmp[42] = 2.0*p[4]*xB_tmp[42] - 1.0*p[4]*xB_tmp[41];
xBdot_tmp[43] = p[4]*xB_tmp[43] - 1.0*p[4]*xB_tmp[40] - 2.0*p[2]*xB_tmp[42];
xBdot_tmp[44] = xB_tmp[44]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[41] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[40];
xBdot_tmp[45] = p[1]*xB_tmp[45] - 1.0*p[1]*xB_tmp[61] - 1.0*p[2]*xB_tmp[62] - 1.0*p[2]*xB_tmp[68];
xBdot_tmp[46] = 0.0002*xB_tmp[46] - 1.0*p[2]*xB_tmp[49] - 0.0002*xB_tmp[61];
xBdot_tmp[47] = (p[5] + 0.0002)*xB_tmp[47] - 0.0002*xB_tmp[48];
xBdot_tmp[48] = xB_tmp[48]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[66];
xBdot_tmp[49] = (p[4] + 0.0002)*xB_tmp[49] - 1.0*p[4]*xB_tmp[47] - 0.0002*xB_tmp[60];
xBdot_tmp[50] = 0.0004*xB_tmp[50] - 0.0002*xB_tmp[51];
xBdot_tmp[51] = (p[0] + 0.0002)*xB_tmp[51] - 1.0*p[0]*xB_tmp[53] - 0.0004*xB_tmp[54];
xBdot_tmp[52] = 0.0002*xB_tmp[51] - 0.0002*xB_tmp[50] + 0.0002*xB_tmp[52] - 0.0002*xB_tmp[54] - 0.0002*xB_tmp[55];
xBdot_tmp[53] = (p[1] + 0.0002)*xB_tmp[53] - 1.0*p[2]*xB_tmp[59] - 1.0*p[1]*xB_tmp[51] - 0.0002*xB_tmp[57];
xBdot_tmp[54] = 2.0*p[0]*xB_tmp[54] - 1.0*p[0]*xB_tmp[57];
xBdot_tmp[55] = p[0]*xB_tmp[55] - 1.0*p[0]*xB_tmp[54] - 1.0*p[0]*xB_tmp[56] + p[0]*xB_tmp[57] - 1.0*p[0]*xB_tmp[63];
xBdot_tmp[56] = 2.0*p[1]*xB_tmp[56] - 1.0*p[1]*xB_tmp[57] - 1.0*p[2]*xB_tmp[65];
xBdot_tmp[57] = xB_tmp[57]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[56] - 1.0*p[2]*xB_tmp[58] - 2.0*p[1]*xB_tmp[54];
xBdot_tmp[58] = xB_tmp[58]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[65] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[60] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[61];
xBdot_tmp[59] = xB_tmp[59]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[49] - 0.0002*xB_tmp[58] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[46];
xBdot_tmp[60] = xB_tmp[60]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[62] - 1.0*p[4]*xB_tmp[48];
xBdot_tmp[61] = p[0]*xB_tmp[61] - 1.0*p[0]*xB_tmp[45] - 1.0*p[2]*xB_tmp[60];
xBdot_tmp[62] = xB_tmp[62]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[66] - 1.0*p[2]*xB_tmp[69] - 1.0*p[1]*xB_tmp[60];
xBdot_tmp[63] = p[1]*xB_tmp[57] - 1.0*p[1]*xB_tmp[55] - 1.0*p[1]*xB_tmp[56] - 1.0*p[1]*xB_tmp[54] + p[1]*xB_tmp[63] - 1.0*p[2]*xB_tmp[64] - 1.0*p[2]*xB_tmp[67];
xBdot_tmp[64] = xB_tmp[64]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[42]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[35]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[39] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[38] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[67] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[43] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[68] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[69];
xBdot_tmp[65] = xB_tmp[65]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[58] - 2.0*p[2]*xB_tmp[67] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[62] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[45];
xBdot_tmp[66] = xB_tmp[66]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[48] - 1.0*p[2]*xB_tmp[44];
xBdot_tmp[67] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[67] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[69] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[68];
xBdot_tmp[68] = xB_tmp[68]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[69] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[43] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[35];
xBdot_tmp[69] = xB_tmp[69]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[44] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[42] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[43];
xBdot_tmp[70] = -1.0*p[2]*xB_tmp[78];
xBdot_tmp[71] = p[5]*xB_tmp[71] - 1.0*p[5]*xB_tmp[72];
xBdot_tmp[72] = 2.0*p[5]*xB_tmp[72];
xBdot_tmp[73] = - 1.0*p[2]*xB_tmp[74] - 1.0*p[2]*xB_tmp[77];
xBdot_tmp[74] = p[4]*xB_tmp[74] - 1.0*p[4]*xB_tmp[72] - 1.0*p[4]*xB_tmp[71] + p[4]*xB_tmp[76] - 1.0*p[4]*xB_tmp[77];
xBdot_tmp[75] = p[5]*xB_tmp[75] - 1.0*p[2]*xB_tmp[76];
xBdot_tmp[76] = xB_tmp[76]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[72];
xBdot_tmp[77] = 2.0*p[4]*xB_tmp[77] - 1.0*p[4]*xB_tmp[76];
xBdot_tmp[78] = p[4]*xB_tmp[78] - 1.0*p[4]*xB_tmp[75] - 2.0*p[2]*xB_tmp[77];
xBdot_tmp[79] = xB_tmp[79]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[76] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[75];
xBdot_tmp[80] = p[1]*xB_tmp[80] - 1.0*p[1]*xB_tmp[96] - 1.0*p[2]*xB_tmp[97] - 1.0*p[2]*xB_tmp[103];
xBdot_tmp[81] = 0.0002*xB_tmp[81] - 1.0*p[2]*xB_tmp[84] - 0.0002*xB_tmp[96];
xBdot_tmp[82] = (p[5] + 0.0002)*xB_tmp[82] - 0.0002*xB_tmp[83];
xBdot_tmp[83] = xB_tmp[83]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[101];
xBdot_tmp[84] = (p[4] + 0.0002)*xB_tmp[84] - 1.0*p[4]*xB_tmp[82] - 0.0002*xB_tmp[95];
xBdot_tmp[85] = 0.0004*xB_tmp[85] - 0.0002*xB_tmp[86];
xBdot_tmp[86] = (p[0] + 0.0002)*xB_tmp[86] - 1.0*p[0]*xB_tmp[88] - 0.0004*xB_tmp[89];
xBdot_tmp[87] = 0.0002*xB_tmp[86] - 0.0002*xB_tmp[85] + 0.0002*xB_tmp[87] - 0.0002*xB_tmp[89] - 0.0002*xB_tmp[90];
xBdot_tmp[88] = (p[1] + 0.0002)*xB_tmp[88] - 1.0*p[2]*xB_tmp[94] - 1.0*p[1]*xB_tmp[86] - 0.0002*xB_tmp[92];
xBdot_tmp[89] = 2.0*p[0]*xB_tmp[89] - 1.0*p[0]*xB_tmp[92];
xBdot_tmp[90] = p[0]*xB_tmp[90] - 1.0*p[0]*xB_tmp[89] - 1.0*p[0]*xB_tmp[91] + p[0]*xB_tmp[92] - 1.0*p[0]*xB_tmp[98];
xBdot_tmp[91] = 2.0*p[1]*xB_tmp[91] - 1.0*p[1]*xB_tmp[92] - 1.0*p[2]*xB_tmp[100];
xBdot_tmp[92] = xB_tmp[92]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[91] - 1.0*p[2]*xB_tmp[93] - 2.0*p[1]*xB_tmp[89];
xBdot_tmp[93] = xB_tmp[93]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[100] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[95] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[96];
xBdot_tmp[94] = xB_tmp[94]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[84] - 0.0002*xB_tmp[93] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[81];
xBdot_tmp[95] = xB_tmp[95]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[97] - 1.0*p[4]*xB_tmp[83];
xBdot_tmp[96] = p[0]*xB_tmp[96] - 1.0*p[0]*xB_tmp[80] - 1.0*p[2]*xB_tmp[95];
xBdot_tmp[97] = xB_tmp[97]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[101] - 1.0*p[2]*xB_tmp[104] - 1.0*p[1]*xB_tmp[95];
xBdot_tmp[98] = p[1]*xB_tmp[92] - 1.0*p[1]*xB_tmp[90] - 1.0*p[1]*xB_tmp[91] - 1.0*p[1]*xB_tmp[89] + p[1]*xB_tmp[98] - 1.0*p[2]*xB_tmp[99] - 1.0*p[2]*xB_tmp[102];
xBdot_tmp[99] = xB_tmp[99]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[77]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[70]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[74] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[73] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[102] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[78] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[103] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[104];
xBdot_tmp[100] = xB_tmp[100]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[93] - 2.0*p[2]*xB_tmp[102] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[97] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[80];
xBdot_tmp[101] = xB_tmp[101]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[83] - 1.0*p[2]*xB_tmp[79];
xBdot_tmp[102] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[102] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[104] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[103];
xBdot_tmp[103] = xB_tmp[103]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[104] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[78] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[70];
xBdot_tmp[104] = xB_tmp[104]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[79] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[77] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[78];
xBdot_tmp[105] = -1.0*p[2]*xB_tmp[113];
xBdot_tmp[106] = p[5]*xB_tmp[106] - 1.0*p[5]*xB_tmp[107];
xBdot_tmp[107] = 2.0*p[5]*xB_tmp[107];
xBdot_tmp[108] = - 1.0*p[2]*xB_tmp[109] - 1.0*p[2]*xB_tmp[112];
xBdot_tmp[109] = p[4]*xB_tmp[109] - 1.0*p[4]*xB_tmp[107] - 1.0*p[4]*xB_tmp[106] + p[4]*xB_tmp[111] - 1.0*p[4]*xB_tmp[112];
xBdot_tmp[110] = p[5]*xB_tmp[110] - 1.0*p[2]*xB_tmp[111];
xBdot_tmp[111] = xB_tmp[111]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[107];
xBdot_tmp[112] = 2.0*p[4]*xB_tmp[112] - 1.0*p[4]*xB_tmp[111];
xBdot_tmp[113] = p[4]*xB_tmp[113] - 1.0*p[4]*xB_tmp[110] - 2.0*p[2]*xB_tmp[112];
xBdot_tmp[114] = xB_tmp[114]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[111] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[110];
xBdot_tmp[115] = p[1]*xB_tmp[115] - 1.0*p[1]*xB_tmp[131] - 1.0*p[2]*xB_tmp[132] - 1.0*p[2]*xB_tmp[138];
xBdot_tmp[116] = 0.0002*xB_tmp[116] - 1.0*p[2]*xB_tmp[119] - 0.0002*xB_tmp[131];
xBdot_tmp[117] = (p[5] + 0.0002)*xB_tmp[117] - 0.0002*xB_tmp[118];
xBdot_tmp[118] = xB_tmp[118]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[136];
xBdot_tmp[119] = (p[4] + 0.0002)*xB_tmp[119] - 1.0*p[4]*xB_tmp[117] - 0.0002*xB_tmp[130];
xBdot_tmp[120] = 0.0004*xB_tmp[120] - 0.0002*xB_tmp[121];
xBdot_tmp[121] = (p[0] + 0.0002)*xB_tmp[121] - 1.0*p[0]*xB_tmp[123] - 0.0004*xB_tmp[124];
xBdot_tmp[122] = 0.0002*xB_tmp[121] - 0.0002*xB_tmp[120] + 0.0002*xB_tmp[122] - 0.0002*xB_tmp[124] - 0.0002*xB_tmp[125];
xBdot_tmp[123] = (p[1] + 0.0002)*xB_tmp[123] - 1.0*p[2]*xB_tmp[129] - 1.0*p[1]*xB_tmp[121] - 0.0002*xB_tmp[127];
xBdot_tmp[124] = 2.0*p[0]*xB_tmp[124] - 1.0*p[0]*xB_tmp[127];
xBdot_tmp[125] = p[0]*xB_tmp[125] - 1.0*p[0]*xB_tmp[124] - 1.0*p[0]*xB_tmp[126] + p[0]*xB_tmp[127] - 1.0*p[0]*xB_tmp[133];
xBdot_tmp[126] = 2.0*p[1]*xB_tmp[126] - 1.0*p[1]*xB_tmp[127] - 1.0*p[2]*xB_tmp[135];
xBdot_tmp[127] = xB_tmp[127]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[126] - 1.0*p[2]*xB_tmp[128] - 2.0*p[1]*xB_tmp[124];
xBdot_tmp[128] = xB_tmp[128]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[135] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[130] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[131];
xBdot_tmp[129] = xB_tmp[129]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[119] - 0.0002*xB_tmp[128] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[116];
xBdot_tmp[130] = xB_tmp[130]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[132] - 1.0*p[4]*xB_tmp[118];
xBdot_tmp[131] = p[0]*xB_tmp[131] - 1.0*p[0]*xB_tmp[115] - 1.0*p[2]*xB_tmp[130];
xBdot_tmp[132] = xB_tmp[132]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[136] - 1.0*p[2]*xB_tmp[139] - 1.0*p[1]*xB_tmp[130];
xBdot_tmp[133] = p[1]*xB_tmp[127] - 1.0*p[1]*xB_tmp[125] - 1.0*p[1]*xB_tmp[126] - 1.0*p[1]*xB_tmp[124] + p[1]*xB_tmp[133] - 1.0*p[2]*xB_tmp[134] - 1.0*p[2]*xB_tmp[137];
xBdot_tmp[134] = xB_tmp[134]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[112]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[105]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[109] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[108] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[137] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[113] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[138] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[139];
xBdot_tmp[135] = xB_tmp[135]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[128] - 2.0*p[2]*xB_tmp[137] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[132] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[115];
xBdot_tmp[136] = xB_tmp[136]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[118] - 1.0*p[2]*xB_tmp[114];
xBdot_tmp[137] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[137] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[139] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[138];
xBdot_tmp[138] = xB_tmp[138]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[139] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[113] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[105];
xBdot_tmp[139] = xB_tmp[139]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[114] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[112] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[113];
xBdot_tmp[140] = -1.0*p[2]*xB_tmp[148];
xBdot_tmp[141] = p[5]*xB_tmp[141] - 1.0*p[5]*xB_tmp[142];
xBdot_tmp[142] = 2.0*p[5]*xB_tmp[142];
xBdot_tmp[143] = - 1.0*p[2]*xB_tmp[144] - 1.0*p[2]*xB_tmp[147];
xBdot_tmp[144] = p[4]*xB_tmp[144] - 1.0*p[4]*xB_tmp[142] - 1.0*p[4]*xB_tmp[141] + p[4]*xB_tmp[146] - 1.0*p[4]*xB_tmp[147];
xBdot_tmp[145] = p[5]*xB_tmp[145] - 1.0*p[2]*xB_tmp[146];
xBdot_tmp[146] = xB_tmp[146]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[142];
xBdot_tmp[147] = 2.0*p[4]*xB_tmp[147] - 1.0*p[4]*xB_tmp[146];
xBdot_tmp[148] = p[4]*xB_tmp[148] - 1.0*p[4]*xB_tmp[145] - 2.0*p[2]*xB_tmp[147];
xBdot_tmp[149] = xB_tmp[149]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[146] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[145];
xBdot_tmp[150] = p[1]*xB_tmp[150] - 1.0*p[1]*xB_tmp[166] - 1.0*p[2]*xB_tmp[167] - 1.0*p[2]*xB_tmp[173];
xBdot_tmp[151] = 0.0002*xB_tmp[151] - 1.0*p[2]*xB_tmp[154] - 0.0002*xB_tmp[166];
xBdot_tmp[152] = (p[5] + 0.0002)*xB_tmp[152] - 0.0002*xB_tmp[153];
xBdot_tmp[153] = xB_tmp[153]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[171];
xBdot_tmp[154] = (p[4] + 0.0002)*xB_tmp[154] - 1.0*p[4]*xB_tmp[152] - 0.0002*xB_tmp[165];
xBdot_tmp[155] = 0.0004*xB_tmp[155] - 0.0002*xB_tmp[156];
xBdot_tmp[156] = (p[0] + 0.0002)*xB_tmp[156] - 1.0*p[0]*xB_tmp[158] - 0.0004*xB_tmp[159];
xBdot_tmp[157] = 0.0002*xB_tmp[156] - 0.0002*xB_tmp[155] + 0.0002*xB_tmp[157] - 0.0002*xB_tmp[159] - 0.0002*xB_tmp[160];
xBdot_tmp[158] = (p[1] + 0.0002)*xB_tmp[158] - 1.0*p[2]*xB_tmp[164] - 1.0*p[1]*xB_tmp[156] - 0.0002*xB_tmp[162];
xBdot_tmp[159] = 2.0*p[0]*xB_tmp[159] - 1.0*p[0]*xB_tmp[162];
xBdot_tmp[160] = p[0]*xB_tmp[160] - 1.0*p[0]*xB_tmp[159] - 1.0*p[0]*xB_tmp[161] + p[0]*xB_tmp[162] - 1.0*p[0]*xB_tmp[168];
xBdot_tmp[161] = 2.0*p[1]*xB_tmp[161] - 1.0*p[1]*xB_tmp[162] - 1.0*p[2]*xB_tmp[170];
xBdot_tmp[162] = xB_tmp[162]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[161] - 1.0*p[2]*xB_tmp[163] - 2.0*p[1]*xB_tmp[159];
xBdot_tmp[163] = xB_tmp[163]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[170] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[165] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[166];
xBdot_tmp[164] = xB_tmp[164]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[154] - 0.0002*xB_tmp[163] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[151];
xBdot_tmp[165] = xB_tmp[165]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[167] - 1.0*p[4]*xB_tmp[153];
xBdot_tmp[166] = p[0]*xB_tmp[166] - 1.0*p[0]*xB_tmp[150] - 1.0*p[2]*xB_tmp[165];
xBdot_tmp[167] = xB_tmp[167]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[171] - 1.0*p[2]*xB_tmp[174] - 1.0*p[1]*xB_tmp[165];
xBdot_tmp[168] = p[1]*xB_tmp[162] - 1.0*p[1]*xB_tmp[160] - 1.0*p[1]*xB_tmp[161] - 1.0*p[1]*xB_tmp[159] + p[1]*xB_tmp[168] - 1.0*p[2]*xB_tmp[169] - 1.0*p[2]*xB_tmp[172];
xBdot_tmp[169] = xB_tmp[169]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[147]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[140]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[144] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[143] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[172] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[148] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[173] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[174];
xBdot_tmp[170] = xB_tmp[170]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[163] - 2.0*p[2]*xB_tmp[172] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[167] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[150];
xBdot_tmp[171] = xB_tmp[171]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[153] - 1.0*p[2]*xB_tmp[149];
xBdot_tmp[172] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[172] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[174] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[173];
xBdot_tmp[173] = xB_tmp[173]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[174] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[148] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[140];
xBdot_tmp[174] = xB_tmp[174]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[149] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[147] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[148];
xBdot_tmp[175] = -1.0*p[2]*xB_tmp[183];
xBdot_tmp[176] = p[5]*xB_tmp[176] - 1.0*p[5]*xB_tmp[177];
xBdot_tmp[177] = 2.0*p[5]*xB_tmp[177];
xBdot_tmp[178] = - 1.0*p[2]*xB_tmp[179] - 1.0*p[2]*xB_tmp[182];
xBdot_tmp[179] = p[4]*xB_tmp[179] - 1.0*p[4]*xB_tmp[177] - 1.0*p[4]*xB_tmp[176] + p[4]*xB_tmp[181] - 1.0*p[4]*xB_tmp[182];
xBdot_tmp[180] = p[5]*xB_tmp[180] - 1.0*p[2]*xB_tmp[181];
xBdot_tmp[181] = xB_tmp[181]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[177];
xBdot_tmp[182] = 2.0*p[4]*xB_tmp[182] - 1.0*p[4]*xB_tmp[181];
xBdot_tmp[183] = p[4]*xB_tmp[183] - 1.0*p[4]*xB_tmp[180] - 2.0*p[2]*xB_tmp[182];
xBdot_tmp[184] = xB_tmp[184]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[181] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[180];
xBdot_tmp[185] = p[1]*xB_tmp[185] - 1.0*p[1]*xB_tmp[201] - 1.0*p[2]*xB_tmp[202] - 1.0*p[2]*xB_tmp[208];
xBdot_tmp[186] = 0.0002*xB_tmp[186] - 1.0*p[2]*xB_tmp[189] - 0.0002*xB_tmp[201];
xBdot_tmp[187] = (p[5] + 0.0002)*xB_tmp[187] - 0.0002*xB_tmp[188];
xBdot_tmp[188] = xB_tmp[188]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[206];
xBdot_tmp[189] = (p[4] + 0.0002)*xB_tmp[189] - 1.0*p[4]*xB_tmp[187] - 0.0002*xB_tmp[200];
xBdot_tmp[190] = 0.0004*xB_tmp[190] - 0.0002*xB_tmp[191];
xBdot_tmp[191] = (p[0] + 0.0002)*xB_tmp[191] - 1.0*p[0]*xB_tmp[193] - 0.0004*xB_tmp[194];
xBdot_tmp[192] = 0.0002*xB_tmp[191] - 0.0002*xB_tmp[190] + 0.0002*xB_tmp[192] - 0.0002*xB_tmp[194] - 0.0002*xB_tmp[195];
xBdot_tmp[193] = (p[1] + 0.0002)*xB_tmp[193] - 1.0*p[2]*xB_tmp[199] - 1.0*p[1]*xB_tmp[191] - 0.0002*xB_tmp[197];
xBdot_tmp[194] = 2.0*p[0]*xB_tmp[194] - 1.0*p[0]*xB_tmp[197];
xBdot_tmp[195] = p[0]*xB_tmp[195] - 1.0*p[0]*xB_tmp[194] - 1.0*p[0]*xB_tmp[196] + p[0]*xB_tmp[197] - 1.0*p[0]*xB_tmp[203];
xBdot_tmp[196] = 2.0*p[1]*xB_tmp[196] - 1.0*p[1]*xB_tmp[197] - 1.0*p[2]*xB_tmp[205];
xBdot_tmp[197] = xB_tmp[197]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[196] - 1.0*p[2]*xB_tmp[198] - 2.0*p[1]*xB_tmp[194];
xBdot_tmp[198] = xB_tmp[198]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[205] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[200] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[201];
xBdot_tmp[199] = xB_tmp[199]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[189] - 0.0002*xB_tmp[198] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[186];
xBdot_tmp[200] = xB_tmp[200]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[202] - 1.0*p[4]*xB_tmp[188];
xBdot_tmp[201] = p[0]*xB_tmp[201] - 1.0*p[0]*xB_tmp[185] - 1.0*p[2]*xB_tmp[200];
xBdot_tmp[202] = xB_tmp[202]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[206] - 1.0*p[2]*xB_tmp[209] - 1.0*p[1]*xB_tmp[200];
xBdot_tmp[203] = p[1]*xB_tmp[197] - 1.0*p[1]*xB_tmp[195] - 1.0*p[1]*xB_tmp[196] - 1.0*p[1]*xB_tmp[194] + p[1]*xB_tmp[203] - 1.0*p[2]*xB_tmp[204] - 1.0*p[2]*xB_tmp[207];
xBdot_tmp[204] = xB_tmp[204]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[182]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[175]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[179] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[178] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[207] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[183] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[208] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[209];
xBdot_tmp[205] = xB_tmp[205]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[198] - 2.0*p[2]*xB_tmp[207] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[202] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[185];
xBdot_tmp[206] = xB_tmp[206]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[188] - 1.0*p[2]*xB_tmp[184];
xBdot_tmp[207] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[207] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[209] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[208];
xBdot_tmp[208] = xB_tmp[208]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[209] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[183] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[175];
xBdot_tmp[209] = xB_tmp[209]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[184] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[182] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[183];
xBdot_tmp[210] = -1.0*p[2]*xB_tmp[218];
xBdot_tmp[211] = p[5]*xB_tmp[211] - 1.0*p[5]*xB_tmp[212];
xBdot_tmp[212] = 2.0*p[5]*xB_tmp[212];
xBdot_tmp[213] = - 1.0*p[2]*xB_tmp[214] - 1.0*p[2]*xB_tmp[217];
xBdot_tmp[214] = p[4]*xB_tmp[214] - 1.0*p[4]*xB_tmp[212] - 1.0*p[4]*xB_tmp[211] + p[4]*xB_tmp[216] - 1.0*p[4]*xB_tmp[217];
xBdot_tmp[215] = p[5]*xB_tmp[215] - 1.0*p[2]*xB_tmp[216];
xBdot_tmp[216] = xB_tmp[216]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[212];
xBdot_tmp[217] = 2.0*p[4]*xB_tmp[217] - 1.0*p[4]*xB_tmp[216];
xBdot_tmp[218] = p[4]*xB_tmp[218] - 1.0*p[4]*xB_tmp[215] - 2.0*p[2]*xB_tmp[217];
xBdot_tmp[219] = xB_tmp[219]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[216] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[215];
xBdot_tmp[220] = p[1]*xB_tmp[220] - 1.0*p[1]*xB_tmp[236] - 1.0*p[2]*xB_tmp[237] - 1.0*p[2]*xB_tmp[243];
xBdot_tmp[221] = 0.0002*xB_tmp[221] - 1.0*p[2]*xB_tmp[224] - 0.0002*xB_tmp[236];
xBdot_tmp[222] = (p[5] + 0.0002)*xB_tmp[222] - 0.0002*xB_tmp[223];
xBdot_tmp[223] = xB_tmp[223]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[241];
xBdot_tmp[224] = (p[4] + 0.0002)*xB_tmp[224] - 1.0*p[4]*xB_tmp[222] - 0.0002*xB_tmp[235];
xBdot_tmp[225] = 0.0004*xB_tmp[225] - 0.0002*xB_tmp[226];
xBdot_tmp[226] = (p[0] + 0.0002)*xB_tmp[226] - 1.0*p[0]*xB_tmp[228] - 0.0004*xB_tmp[229];
xBdot_tmp[227] = 0.0002*xB_tmp[226] - 0.0002*xB_tmp[225] + 0.0002*xB_tmp[227] - 0.0002*xB_tmp[229] - 0.0002*xB_tmp[230];
xBdot_tmp[228] = (p[1] + 0.0002)*xB_tmp[228] - 1.0*p[2]*xB_tmp[234] - 1.0*p[1]*xB_tmp[226] - 0.0002*xB_tmp[232];
xBdot_tmp[229] = 2.0*p[0]*xB_tmp[229] - 1.0*p[0]*xB_tmp[232];
xBdot_tmp[230] = p[0]*xB_tmp[230] - 1.0*p[0]*xB_tmp[229] - 1.0*p[0]*xB_tmp[231] + p[0]*xB_tmp[232] - 1.0*p[0]*xB_tmp[238];
xBdot_tmp[231] = 2.0*p[1]*xB_tmp[231] - 1.0*p[1]*xB_tmp[232] - 1.0*p[2]*xB_tmp[240];
xBdot_tmp[232] = xB_tmp[232]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[231] - 1.0*p[2]*xB_tmp[233] - 2.0*p[1]*xB_tmp[229];
xBdot_tmp[233] = xB_tmp[233]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[240] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[235] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[236];
xBdot_tmp[234] = xB_tmp[234]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[224] - 0.0002*xB_tmp[233] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[221];
xBdot_tmp[235] = xB_tmp[235]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[237] - 1.0*p[4]*xB_tmp[223];
xBdot_tmp[236] = p[0]*xB_tmp[236] - 1.0*p[0]*xB_tmp[220] - 1.0*p[2]*xB_tmp[235];
xBdot_tmp[237] = xB_tmp[237]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[241] - 1.0*p[2]*xB_tmp[244] - 1.0*p[1]*xB_tmp[235];
xBdot_tmp[238] = p[1]*xB_tmp[232] - 1.0*p[1]*xB_tmp[230] - 1.0*p[1]*xB_tmp[231] - 1.0*p[1]*xB_tmp[229] + p[1]*xB_tmp[238] - 1.0*p[2]*xB_tmp[239] - 1.0*p[2]*xB_tmp[242];
xBdot_tmp[239] = xB_tmp[239]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[217]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[210]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[214] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[213] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[242] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[218] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[243] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[244];
xBdot_tmp[240] = xB_tmp[240]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[233] - 2.0*p[2]*xB_tmp[242] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[237] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[220];
xBdot_tmp[241] = xB_tmp[241]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[223] - 1.0*p[2]*xB_tmp[219];
xBdot_tmp[242] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[242] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[244] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[243];
xBdot_tmp[243] = xB_tmp[243]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[244] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[218] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[210];
xBdot_tmp[244] = xB_tmp[244]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[219] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[217] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[218];
xBdot_tmp[245] = -1.0*p[2]*xB_tmp[253];
xBdot_tmp[246] = p[5]*xB_tmp[246] - 1.0*p[5]*xB_tmp[247];
xBdot_tmp[247] = 2.0*p[5]*xB_tmp[247];
xBdot_tmp[248] = - 1.0*p[2]*xB_tmp[249] - 1.0*p[2]*xB_tmp[252];
xBdot_tmp[249] = p[4]*xB_tmp[249] - 1.0*p[4]*xB_tmp[247] - 1.0*p[4]*xB_tmp[246] + p[4]*xB_tmp[251] - 1.0*p[4]*xB_tmp[252];
xBdot_tmp[250] = p[5]*xB_tmp[250] - 1.0*p[2]*xB_tmp[251];
xBdot_tmp[251] = xB_tmp[251]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[247];
xBdot_tmp[252] = 2.0*p[4]*xB_tmp[252] - 1.0*p[4]*xB_tmp[251];
xBdot_tmp[253] = p[4]*xB_tmp[253] - 1.0*p[4]*xB_tmp[250] - 2.0*p[2]*xB_tmp[252];
xBdot_tmp[254] = xB_tmp[254]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[251] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[250];
xBdot_tmp[255] = p[1]*xB_tmp[255] - 1.0*p[1]*xB_tmp[271] - 1.0*p[2]*xB_tmp[272] - 1.0*p[2]*xB_tmp[278];
xBdot_tmp[256] = 0.0002*xB_tmp[256] - 1.0*p[2]*xB_tmp[259] - 0.0002*xB_tmp[271];
xBdot_tmp[257] = (p[5] + 0.0002)*xB_tmp[257] - 0.0002*xB_tmp[258];
xBdot_tmp[258] = xB_tmp[258]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[276];
xBdot_tmp[259] = (p[4] + 0.0002)*xB_tmp[259] - 1.0*p[4]*xB_tmp[257] - 0.0002*xB_tmp[270];
xBdot_tmp[260] = 0.0004*xB_tmp[260] - 0.0002*xB_tmp[261];
xBdot_tmp[261] = (p[0] + 0.0002)*xB_tmp[261] - 1.0*p[0]*xB_tmp[263] - 0.0004*xB_tmp[264];
xBdot_tmp[262] = 0.0002*xB_tmp[261] - 0.0002*xB_tmp[260] + 0.0002*xB_tmp[262] - 0.0002*xB_tmp[264] - 0.0002*xB_tmp[265];
xBdot_tmp[263] = (p[1] + 0.0002)*xB_tmp[263] - 1.0*p[2]*xB_tmp[269] - 1.0*p[1]*xB_tmp[261] - 0.0002*xB_tmp[267];
xBdot_tmp[264] = 2.0*p[0]*xB_tmp[264] - 1.0*p[0]*xB_tmp[267];
xBdot_tmp[265] = p[0]*xB_tmp[265] - 1.0*p[0]*xB_tmp[264] - 1.0*p[0]*xB_tmp[266] + p[0]*xB_tmp[267] - 1.0*p[0]*xB_tmp[273];
xBdot_tmp[266] = 2.0*p[1]*xB_tmp[266] - 1.0*p[1]*xB_tmp[267] - 1.0*p[2]*xB_tmp[275];
xBdot_tmp[267] = xB_tmp[267]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[266] - 1.0*p[2]*xB_tmp[268] - 2.0*p[1]*xB_tmp[264];
xBdot_tmp[268] = xB_tmp[268]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[275] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[270] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[271];
xBdot_tmp[269] = xB_tmp[269]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[259] - 0.0002*xB_tmp[268] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[256];
xBdot_tmp[270] = xB_tmp[270]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[272] - 1.0*p[4]*xB_tmp[258];
xBdot_tmp[271] = p[0]*xB_tmp[271] - 1.0*p[0]*xB_tmp[255] - 1.0*p[2]*xB_tmp[270];
xBdot_tmp[272] = xB_tmp[272]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[276] - 1.0*p[2]*xB_tmp[279] - 1.0*p[1]*xB_tmp[270];
xBdot_tmp[273] = p[1]*xB_tmp[267] - 1.0*p[1]*xB_tmp[265] - 1.0*p[1]*xB_tmp[266] - 1.0*p[1]*xB_tmp[264] + p[1]*xB_tmp[273] - 1.0*p[2]*xB_tmp[274] - 1.0*p[2]*xB_tmp[277];
xBdot_tmp[274] = xB_tmp[274]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[252]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[245]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[249] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[248] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[277] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[253] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[278] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[279];
xBdot_tmp[275] = xB_tmp[275]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[268] - 2.0*p[2]*xB_tmp[277] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[272] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[255];
xBdot_tmp[276] = xB_tmp[276]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[258] - 1.0*p[2]*xB_tmp[254];
xBdot_tmp[277] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[277] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[279] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[278];
xBdot_tmp[278] = xB_tmp[278]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[279] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[253] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[245];
xBdot_tmp[279] = xB_tmp[279]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[254] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[252] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[253];
xBdot_tmp[280] = -1.0*p[2]*xB_tmp[288];
xBdot_tmp[281] = p[5]*xB_tmp[281] - 1.0*p[5]*xB_tmp[282];
xBdot_tmp[282] = 2.0*p[5]*xB_tmp[282];
xBdot_tmp[283] = - 1.0*p[2]*xB_tmp[284] - 1.0*p[2]*xB_tmp[287];
xBdot_tmp[284] = p[4]*xB_tmp[284] - 1.0*p[4]*xB_tmp[282] - 1.0*p[4]*xB_tmp[281] + p[4]*xB_tmp[286] - 1.0*p[4]*xB_tmp[287];
xBdot_tmp[285] = p[5]*xB_tmp[285] - 1.0*p[2]*xB_tmp[286];
xBdot_tmp[286] = xB_tmp[286]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[282];
xBdot_tmp[287] = 2.0*p[4]*xB_tmp[287] - 1.0*p[4]*xB_tmp[286];
xBdot_tmp[288] = p[4]*xB_tmp[288] - 1.0*p[4]*xB_tmp[285] - 2.0*p[2]*xB_tmp[287];
xBdot_tmp[289] = xB_tmp[289]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[286] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[285];
xBdot_tmp[290] = p[1]*xB_tmp[290] - 1.0*p[1]*xB_tmp[306] - 1.0*p[2]*xB_tmp[307] - 1.0*p[2]*xB_tmp[313];
xBdot_tmp[291] = 0.0002*xB_tmp[291] - 1.0*p[2]*xB_tmp[294] - 0.0002*xB_tmp[306];
xBdot_tmp[292] = (p[5] + 0.0002)*xB_tmp[292] - 0.0002*xB_tmp[293];
xBdot_tmp[293] = xB_tmp[293]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[311];
xBdot_tmp[294] = (p[4] + 0.0002)*xB_tmp[294] - 1.0*p[4]*xB_tmp[292] - 0.0002*xB_tmp[305];
xBdot_tmp[295] = 0.0004*xB_tmp[295] - 0.0002*xB_tmp[296];
xBdot_tmp[296] = (p[0] + 0.0002)*xB_tmp[296] - 1.0*p[0]*xB_tmp[298] - 0.0004*xB_tmp[299];
xBdot_tmp[297] = 0.0002*xB_tmp[296] - 0.0002*xB_tmp[295] + 0.0002*xB_tmp[297] - 0.0002*xB_tmp[299] - 0.0002*xB_tmp[300];
xBdot_tmp[298] = (p[1] + 0.0002)*xB_tmp[298] - 1.0*p[2]*xB_tmp[304] - 1.0*p[1]*xB_tmp[296] - 0.0002*xB_tmp[302];
xBdot_tmp[299] = 2.0*p[0]*xB_tmp[299] - 1.0*p[0]*xB_tmp[302];
xBdot_tmp[300] = p[0]*xB_tmp[300] - 1.0*p[0]*xB_tmp[299] - 1.0*p[0]*xB_tmp[301] + p[0]*xB_tmp[302] - 1.0*p[0]*xB_tmp[308];
xBdot_tmp[301] = 2.0*p[1]*xB_tmp[301] - 1.0*p[1]*xB_tmp[302] - 1.0*p[2]*xB_tmp[310];
xBdot_tmp[302] = xB_tmp[302]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[301] - 1.0*p[2]*xB_tmp[303] - 2.0*p[1]*xB_tmp[299];
xBdot_tmp[303] = xB_tmp[303]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[310] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[305] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[306];
xBdot_tmp[304] = xB_tmp[304]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[294] - 0.0002*xB_tmp[303] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[291];
xBdot_tmp[305] = xB_tmp[305]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[307] - 1.0*p[4]*xB_tmp[293];
xBdot_tmp[306] = p[0]*xB_tmp[306] - 1.0*p[0]*xB_tmp[290] - 1.0*p[2]*xB_tmp[305];
xBdot_tmp[307] = xB_tmp[307]*(p[1] + p[4]) - 1.0*p[4]*xB_tmp[311] - 1.0*p[2]*xB_tmp[314] - 1.0*p[1]*xB_tmp[305];
xBdot_tmp[308] = p[1]*xB_tmp[302] - 1.0*p[1]*xB_tmp[300] - 1.0*p[1]*xB_tmp[301] - 1.0*p[1]*xB_tmp[299] + p[1]*xB_tmp[308] - 1.0*p[2]*xB_tmp[309] - 1.0*p[2]*xB_tmp[312];
xBdot_tmp[309] = xB_tmp[309]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*xB_tmp[287]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*xB_tmp[280]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[284] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[283] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[312] - 0.495*p[2]*(pow(p[3],2))*xB_tmp[288] + 1.1*p[2]*(pow(p[3],2))*xB_tmp[313] + 0.9*p[2]*(pow(p[3],2))*xB_tmp[314];
xBdot_tmp[310] = xB_tmp[310]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*xB_tmp[303] - 2.0*p[2]*xB_tmp[312] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[307] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[290];
xBdot_tmp[311] = xB_tmp[311]*(p[1] + p[5]) - 1.0*p[1]*xB_tmp[293] - 1.0*p[2]*xB_tmp[289];
xBdot_tmp[312] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[312] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[314] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[313];
xBdot_tmp[313] = xB_tmp[313]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[314] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[288] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[280];
xBdot_tmp[314] = xB_tmp[314]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[289] - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[287] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[288];

  for (ixB=0; ixB<315; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
{
  int iyp;
  int ip;
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  int np = *data->np;
  int ny = *data->ny;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *qBdot_tmp = N_VGetArrayPointer(qBdot);
  memset(qBdot_tmp,0,sizeof(double)*9*np);
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
qBdot_tmp[0+ip*ny] = x_tmp[13]*xB_tmp[13] + x_tmp[16]*xB_tmp[16] - 1.0*x_tmp[16]*xB_tmp[18] - 1.0*x_tmp[26]*xB_tmp[10] + x_tmp[20]*xB_tmp[20] - 1.0*x_tmp[13]*xB_tmp[31] + x_tmp[23]*xB_tmp[23] - 1.0*x_tmp[20]*xB_tmp[28] + x_tmp[25]*xB_tmp[25] - 1.0*x_tmp[25]*xB_tmp[27] + x_tmp[26]*xB_tmp[26] - 1.0*x_tmp[23]*xB_tmp[30] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[19] + xB_tmp[22]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[21];
qBdot_tmp[1+ip*ny] = x_tmp[13]*xB_tmp[48] + x_tmp[16]*xB_tmp[51] - 1.0*x_tmp[16]*xB_tmp[53] - 1.0*x_tmp[26]*xB_tmp[45] + x_tmp[20]*xB_tmp[55] - 1.0*x_tmp[13]*xB_tmp[66] + x_tmp[23]*xB_tmp[58] - 1.0*x_tmp[20]*xB_tmp[63] + x_tmp[25]*xB_tmp[60] - 1.0*x_tmp[25]*xB_tmp[62] + x_tmp[26]*xB_tmp[61] - 1.0*x_tmp[23]*xB_tmp[65] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[54] + xB_tmp[57]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[56];
qBdot_tmp[2+ip*ny] = x_tmp[13]*xB_tmp[83] + x_tmp[16]*xB_tmp[86] - 1.0*x_tmp[16]*xB_tmp[88] - 1.0*x_tmp[26]*xB_tmp[80] + x_tmp[20]*xB_tmp[90] - 1.0*x_tmp[13]*xB_tmp[101] + x_tmp[23]*xB_tmp[93] - 1.0*x_tmp[20]*xB_tmp[98] + x_tmp[25]*xB_tmp[95] - 1.0*x_tmp[25]*xB_tmp[97] + x_tmp[26]*xB_tmp[96] - 1.0*x_tmp[23]*xB_tmp[100] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[89] + xB_tmp[92]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[91];
qBdot_tmp[3+ip*ny] = x_tmp[13]*xB_tmp[118] + x_tmp[16]*xB_tmp[121] - 1.0*x_tmp[16]*xB_tmp[123] - 1.0*x_tmp[26]*xB_tmp[115] + x_tmp[20]*xB_tmp[125] - 1.0*x_tmp[13]*xB_tmp[136] + x_tmp[23]*xB_tmp[128] - 1.0*x_tmp[20]*xB_tmp[133] + x_tmp[25]*xB_tmp[130] - 1.0*x_tmp[25]*xB_tmp[132] + x_tmp[26]*xB_tmp[131] - 1.0*x_tmp[23]*xB_tmp[135] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[124] + xB_tmp[127]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[126];
qBdot_tmp[4+ip*ny] = x_tmp[13]*xB_tmp[153] + x_tmp[16]*xB_tmp[156] - 1.0*x_tmp[16]*xB_tmp[158] - 1.0*x_tmp[26]*xB_tmp[150] + x_tmp[20]*xB_tmp[160] - 1.0*x_tmp[13]*xB_tmp[171] + x_tmp[23]*xB_tmp[163] - 1.0*x_tmp[20]*xB_tmp[168] + x_tmp[25]*xB_tmp[165] - 1.0*x_tmp[25]*xB_tmp[167] + x_tmp[26]*xB_tmp[166] - 1.0*x_tmp[23]*xB_tmp[170] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[159] + xB_tmp[162]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[161];
qBdot_tmp[5+ip*ny] = x_tmp[13]*xB_tmp[188] + x_tmp[16]*xB_tmp[191] - 1.0*x_tmp[16]*xB_tmp[193] - 1.0*x_tmp[26]*xB_tmp[185] + x_tmp[20]*xB_tmp[195] - 1.0*x_tmp[13]*xB_tmp[206] + x_tmp[23]*xB_tmp[198] - 1.0*x_tmp[20]*xB_tmp[203] + x_tmp[25]*xB_tmp[200] - 1.0*x_tmp[25]*xB_tmp[202] + x_tmp[26]*xB_tmp[201] - 1.0*x_tmp[23]*xB_tmp[205] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[194] + xB_tmp[197]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[196];
qBdot_tmp[6+ip*ny] = x_tmp[13]*xB_tmp[223] + x_tmp[16]*xB_tmp[226] - 1.0*x_tmp[16]*xB_tmp[228] - 1.0*x_tmp[26]*xB_tmp[220] + x_tmp[20]*xB_tmp[230] - 1.0*x_tmp[13]*xB_tmp[241] + x_tmp[23]*xB_tmp[233] - 1.0*x_tmp[20]*xB_tmp[238] + x_tmp[25]*xB_tmp[235] - 1.0*x_tmp[25]*xB_tmp[237] + x_tmp[26]*xB_tmp[236] - 1.0*x_tmp[23]*xB_tmp[240] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[229] + xB_tmp[232]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[231];
qBdot_tmp[7+ip*ny] = x_tmp[13]*xB_tmp[258] + x_tmp[16]*xB_tmp[261] - 1.0*x_tmp[16]*xB_tmp[263] - 1.0*x_tmp[26]*xB_tmp[255] + x_tmp[20]*xB_tmp[265] - 1.0*x_tmp[13]*xB_tmp[276] + x_tmp[23]*xB_tmp[268] - 1.0*x_tmp[20]*xB_tmp[273] + x_tmp[25]*xB_tmp[270] - 1.0*x_tmp[25]*xB_tmp[272] + x_tmp[26]*xB_tmp[271] - 1.0*x_tmp[23]*xB_tmp[275] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[264] + xB_tmp[267]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[266];
qBdot_tmp[8+ip*ny] = x_tmp[13]*xB_tmp[293] + x_tmp[16]*xB_tmp[296] - 1.0*x_tmp[16]*xB_tmp[298] - 1.0*x_tmp[26]*xB_tmp[290] + x_tmp[20]*xB_tmp[300] - 1.0*x_tmp[13]*xB_tmp[311] + x_tmp[23]*xB_tmp[303] - 1.0*x_tmp[20]*xB_tmp[308] + x_tmp[25]*xB_tmp[305] - 1.0*x_tmp[25]*xB_tmp[307] + x_tmp[26]*xB_tmp[306] - 1.0*x_tmp[23]*xB_tmp[310] + (2.0*x_tmp[19] - 1.0*x_tmp[20])*xB_tmp[299] + xB_tmp[302]*(x_tmp[20] - 1.0*x_tmp[19] + x_tmp[22]) - 1.0*(x_tmp[20] + 2.0*x_tmp[22])*xB_tmp[301];

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = x_tmp[10]*xB_tmp[10] - 1.0*x_tmp[18]*xB_tmp[16] - 1.0*x_tmp[10]*xB_tmp[26] + x_tmp[18]*xB_tmp[18] - 1.0*x_tmp[31]*xB_tmp[13] - 1.0*x_tmp[28]*xB_tmp[20] - 1.0*x_tmp[27]*xB_tmp[25] - 1.0*x_tmp[30]*xB_tmp[23] + x_tmp[27]*xB_tmp[27] + x_tmp[28]*xB_tmp[28] + x_tmp[30]*xB_tmp[30] + x_tmp[31]*xB_tmp[31] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[21] + xB_tmp[22]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[19];
qBdot_tmp[1+ip*ny] = x_tmp[10]*xB_tmp[45] - 1.0*x_tmp[18]*xB_tmp[51] - 1.0*x_tmp[10]*xB_tmp[61] + x_tmp[18]*xB_tmp[53] - 1.0*x_tmp[31]*xB_tmp[48] - 1.0*x_tmp[28]*xB_tmp[55] - 1.0*x_tmp[27]*xB_tmp[60] - 1.0*x_tmp[30]*xB_tmp[58] + x_tmp[27]*xB_tmp[62] + x_tmp[28]*xB_tmp[63] + x_tmp[30]*xB_tmp[65] + x_tmp[31]*xB_tmp[66] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[56] + xB_tmp[57]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[54];
qBdot_tmp[2+ip*ny] = x_tmp[10]*xB_tmp[80] - 1.0*x_tmp[18]*xB_tmp[86] - 1.0*x_tmp[10]*xB_tmp[96] + x_tmp[18]*xB_tmp[88] - 1.0*x_tmp[31]*xB_tmp[83] - 1.0*x_tmp[28]*xB_tmp[90] - 1.0*x_tmp[27]*xB_tmp[95] - 1.0*x_tmp[30]*xB_tmp[93] + x_tmp[27]*xB_tmp[97] + x_tmp[28]*xB_tmp[98] + x_tmp[30]*xB_tmp[100] + x_tmp[31]*xB_tmp[101] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[91] + xB_tmp[92]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[89];
qBdot_tmp[3+ip*ny] = x_tmp[10]*xB_tmp[115] - 1.0*x_tmp[18]*xB_tmp[121] - 1.0*x_tmp[10]*xB_tmp[131] + x_tmp[18]*xB_tmp[123] - 1.0*x_tmp[31]*xB_tmp[118] - 1.0*x_tmp[28]*xB_tmp[125] - 1.0*x_tmp[27]*xB_tmp[130] - 1.0*x_tmp[30]*xB_tmp[128] + x_tmp[27]*xB_tmp[132] + x_tmp[28]*xB_tmp[133] + x_tmp[30]*xB_tmp[135] + x_tmp[31]*xB_tmp[136] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[126] + xB_tmp[127]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[124];
qBdot_tmp[4+ip*ny] = x_tmp[10]*xB_tmp[150] - 1.0*x_tmp[18]*xB_tmp[156] - 1.0*x_tmp[10]*xB_tmp[166] + x_tmp[18]*xB_tmp[158] - 1.0*x_tmp[31]*xB_tmp[153] - 1.0*x_tmp[28]*xB_tmp[160] - 1.0*x_tmp[27]*xB_tmp[165] - 1.0*x_tmp[30]*xB_tmp[163] + x_tmp[27]*xB_tmp[167] + x_tmp[28]*xB_tmp[168] + x_tmp[30]*xB_tmp[170] + x_tmp[31]*xB_tmp[171] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[161] + xB_tmp[162]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[159];
qBdot_tmp[5+ip*ny] = x_tmp[10]*xB_tmp[185] - 1.0*x_tmp[18]*xB_tmp[191] - 1.0*x_tmp[10]*xB_tmp[201] + x_tmp[18]*xB_tmp[193] - 1.0*x_tmp[31]*xB_tmp[188] - 1.0*x_tmp[28]*xB_tmp[195] - 1.0*x_tmp[27]*xB_tmp[200] - 1.0*x_tmp[30]*xB_tmp[198] + x_tmp[27]*xB_tmp[202] + x_tmp[28]*xB_tmp[203] + x_tmp[30]*xB_tmp[205] + x_tmp[31]*xB_tmp[206] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[196] + xB_tmp[197]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[194];
qBdot_tmp[6+ip*ny] = x_tmp[10]*xB_tmp[220] - 1.0*x_tmp[18]*xB_tmp[226] - 1.0*x_tmp[10]*xB_tmp[236] + x_tmp[18]*xB_tmp[228] - 1.0*x_tmp[31]*xB_tmp[223] - 1.0*x_tmp[28]*xB_tmp[230] - 1.0*x_tmp[27]*xB_tmp[235] - 1.0*x_tmp[30]*xB_tmp[233] + x_tmp[27]*xB_tmp[237] + x_tmp[28]*xB_tmp[238] + x_tmp[30]*xB_tmp[240] + x_tmp[31]*xB_tmp[241] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[231] + xB_tmp[232]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[229];
qBdot_tmp[7+ip*ny] = x_tmp[10]*xB_tmp[255] - 1.0*x_tmp[18]*xB_tmp[261] - 1.0*x_tmp[10]*xB_tmp[271] + x_tmp[18]*xB_tmp[263] - 1.0*x_tmp[31]*xB_tmp[258] - 1.0*x_tmp[28]*xB_tmp[265] - 1.0*x_tmp[27]*xB_tmp[270] - 1.0*x_tmp[30]*xB_tmp[268] + x_tmp[27]*xB_tmp[272] + x_tmp[28]*xB_tmp[273] + x_tmp[30]*xB_tmp[275] + x_tmp[31]*xB_tmp[276] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[266] + xB_tmp[267]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[264];
qBdot_tmp[8+ip*ny] = x_tmp[10]*xB_tmp[290] - 1.0*x_tmp[18]*xB_tmp[296] - 1.0*x_tmp[10]*xB_tmp[306] + x_tmp[18]*xB_tmp[298] - 1.0*x_tmp[31]*xB_tmp[293] - 1.0*x_tmp[28]*xB_tmp[300] - 1.0*x_tmp[27]*xB_tmp[305] - 1.0*x_tmp[30]*xB_tmp[303] + x_tmp[27]*xB_tmp[307] + x_tmp[28]*xB_tmp[308] + x_tmp[30]*xB_tmp[310] + x_tmp[31]*xB_tmp[311] + (2.0*x_tmp[21] - 1.0*x_tmp[28])*xB_tmp[301] + xB_tmp[302]*(x_tmp[22] - 1.0*x_tmp[21] + x_tmp[28]) - 1.0*(2.0*x_tmp[22] + x_tmp[28])*xB_tmp[299];

  } break;

  case 2: {
qBdot_tmp[0+ip*ny] = - 1.0*xB_tmp[9]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[24]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[23]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[30]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[29]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[6]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[4]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[14]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[25]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[27]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[32]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[7]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[8]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[34]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[0]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[33]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[5] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[11] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[3] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[10] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[26];
qBdot_tmp[1+ip*ny] = - 1.0*xB_tmp[44]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[59]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[58]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[65]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[64]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[41]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[39]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[49]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[60]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[62]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[67]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[42]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[43]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[69]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[35]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[68]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[40] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[46] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[38] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[45] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[61];
qBdot_tmp[2+ip*ny] = - 1.0*xB_tmp[79]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[94]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[93]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[100]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[99]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[76]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[74]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[84]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[95]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[97]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[102]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[77]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[78]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[104]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[70]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[103]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[75] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[81] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[73] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[80] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[96];
qBdot_tmp[3+ip*ny] = - 1.0*xB_tmp[114]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[129]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[128]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[135]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[134]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[111]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[109]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[119]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[130]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[132]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[137]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[112]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[113]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[139]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[105]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[138]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[110] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[116] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[108] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[115] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[131];
qBdot_tmp[4+ip*ny] = - 1.0*xB_tmp[149]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[164]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[163]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[170]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[169]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[146]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[144]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[154]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[165]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[167]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[172]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[147]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[148]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[174]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[140]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[173]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[145] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[151] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[143] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[150] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[166];
qBdot_tmp[5+ip*ny] = - 1.0*xB_tmp[184]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[199]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[198]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[205]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[204]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[181]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[179]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[189]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[200]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[202]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[207]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[182]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[183]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[209]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[175]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[208]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[180] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[186] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[178] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[185] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[201];
qBdot_tmp[6+ip*ny] = - 1.0*xB_tmp[219]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[234]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[233]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[240]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[239]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[216]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[214]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[224]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[235]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[237]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[242]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[217]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[218]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[244]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[210]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[243]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[215] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[221] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[213] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[220] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[236];
qBdot_tmp[7+ip*ny] = - 1.0*xB_tmp[254]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[269]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[268]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[275]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[274]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[251]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[249]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[259]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[270]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[272]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[277]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[252]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[253]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[279]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[245]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[278]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[250] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[256] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[248] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[255] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[271];
qBdot_tmp[8+ip*ny] = - 1.0*xB_tmp[289]*((pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31]) - 1.0*xB_tmp[304]*((pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] + x_tmp[18]) - 1.0*xB_tmp[303]*((pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + x_tmp[22]) - 1.0*xB_tmp[310]*((pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21]) - 1.0*xB_tmp[309]*((pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28]) - 1.0*xB_tmp[286]*(0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9]) - 1.0*xB_tmp[284]*(0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29]) - 1.0*xB_tmp[294]*(0.9*(pow(p[3],2))*x_tmp[24] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24]) - 1.0*xB_tmp[305]*(0.9*(pow(p[3],2))*x_tmp[23] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23]) - 1.0*xB_tmp[307]*(0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30]) - 1.0*xB_tmp[312]*((pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] + x_tmp[28] + 2.0*x_tmp[30]) - 1.0*xB_tmp[287]*(1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[288]*(0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34]) - 1.0*xB_tmp[314]*((pow((p[3] - 1.0),2))*x_tmp[34] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*xB_tmp[280]*(1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33]) - xB_tmp[313]*(1.0*(pow((p[3] - 1.0),2))*x_tmp[33] - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - (pow(p[3],2))*x_tmp[33] + 1.0*x_tmp[10] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32]) - 1.0*(1.1*(pow(p[3],2))*x_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9])*xB_tmp[285] - 1.0*(1.1*(pow(p[3],2))*x_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24])*xB_tmp[291] - 1.0*(1.1*(pow(p[3],2))*x_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29])*xB_tmp[283] - 1.0*(1.1*(pow(p[3],2))*x_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30])*xB_tmp[290] - 1.0*(1.1*(pow(p[3],2))*x_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23])*xB_tmp[306];

  } break;

  case 3: {
qBdot_tmp[0+ip*ny] = xB_tmp[0]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[32]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[5] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[11] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[3] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[10] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[26] + xB_tmp[34]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[9] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[23] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[24] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[29] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[30] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[6] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[14] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[4] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[25] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[27] + xB_tmp[7]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[8]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[33]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[1+ip*ny] = xB_tmp[35]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[67]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[40] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[46] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[38] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[45] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[61] + xB_tmp[69]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[44] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[58] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[59] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[64] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[65] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[41] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[49] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[39] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[60] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[62] + xB_tmp[42]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[43]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[68]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[2+ip*ny] = xB_tmp[70]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[102]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[75] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[81] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[73] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[80] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[96] + xB_tmp[104]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[79] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[93] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[94] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[99] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[100] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[76] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[84] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[74] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[95] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[97] + xB_tmp[77]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[78]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[103]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[3+ip*ny] = xB_tmp[105]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[137]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[110] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[116] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[108] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[115] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[131] + xB_tmp[139]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[114] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[128] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[129] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[134] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[135] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[111] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[119] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[109] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[130] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[132] + xB_tmp[112]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[113]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[138]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[4+ip*ny] = xB_tmp[140]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[172]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[145] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[151] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[143] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[150] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[166] + xB_tmp[174]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[149] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[163] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[164] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[169] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[170] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[146] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[154] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[144] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[165] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[167] + xB_tmp[147]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[148]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[173]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[5+ip*ny] = xB_tmp[175]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[207]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[180] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[186] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[178] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[185] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[201] + xB_tmp[209]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[184] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[198] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[199] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[204] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[205] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[181] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[189] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[179] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[200] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[202] + xB_tmp[182]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[183]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[208]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[6+ip*ny] = xB_tmp[210]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[242]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[215] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[221] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[213] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[220] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[236] + xB_tmp[244]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[219] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[233] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[234] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[239] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[240] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[216] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[224] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[214] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[235] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[237] + xB_tmp[217]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[218]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[243]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[7+ip*ny] = xB_tmp[245]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[277]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[250] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[256] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[248] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[255] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[271] + xB_tmp[279]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[254] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[268] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[269] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[274] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[275] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[251] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[259] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[249] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[270] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[272] + xB_tmp[252]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[253]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[278]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);
qBdot_tmp[8+ip*ny] = xB_tmp[280]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] - 2.31*p[2]*p[3]*x_tmp[29] - 2.2*p[2]*p[3]*x_tmp[33]) - 1.0*xB_tmp[312]*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32]) + (1.1*(p[3] - 1.0)*p[2]*x_tmp[9] - 1.1*p[2]*p[3]*x_tmp[9])*xB_tmp[285] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[24] - 1.1*p[2]*p[3]*x_tmp[24])*xB_tmp[291] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[29])*xB_tmp[283] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[30] - 1.1*p[2]*p[3]*x_tmp[30])*xB_tmp[290] + (1.1*(p[3] - 1.0)*p[2]*x_tmp[23] - 1.1*p[2]*p[3]*x_tmp[23])*xB_tmp[306] + xB_tmp[314]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[34] + 1.8*p[2]*p[3]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[34]) - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9])*xB_tmp[289] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23])*xB_tmp[303] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24])*xB_tmp[304] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29])*xB_tmp[309] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30])*xB_tmp[310] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[9] - 0.9*p[2]*p[3]*x_tmp[9])*xB_tmp[286] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[24] - 0.9*p[2]*p[3]*x_tmp[24])*xB_tmp[294] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 0.9*p[2]*p[3]*x_tmp[29])*xB_tmp[284] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[23] - 0.9*p[2]*p[3]*x_tmp[23])*xB_tmp[305] + (0.9*(p[3] - 1.0)*p[2]*x_tmp[30] - 0.9*p[2]*p[3]*x_tmp[30])*xB_tmp[307] + xB_tmp[287]*(0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] - 1.71*p[2]*p[3]*x_tmp[29] - 1.8*p[2]*p[3]*x_tmp[34]) - 1.0*xB_tmp[288]*(0.99*p[2]*p[3]*x_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34]) + xB_tmp[313]*(1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[33] + 2.2*p[2]*p[3]*x_tmp[29] - 1.1*p[2]*p[3]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[33]);

  } break;

  case 4: {
qBdot_tmp[0+ip*ny] = x_tmp[4]*xB_tmp[4] - 1.0*x_tmp[4]*xB_tmp[1] - 1.0*x_tmp[8]*xB_tmp[5] + x_tmp[8]*xB_tmp[8] - 1.0*x_tmp[14]*xB_tmp[12] + x_tmp[14]*xB_tmp[14] - 1.0*x_tmp[25]*xB_tmp[13] - 1.0*x_tmp[34]*xB_tmp[9] + x_tmp[25]*xB_tmp[25] + x_tmp[27]*xB_tmp[27] - 1.0*x_tmp[27]*xB_tmp[31] + x_tmp[34]*xB_tmp[34] + xB_tmp[6]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[2] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[7];
qBdot_tmp[1+ip*ny] = x_tmp[4]*xB_tmp[39] - 1.0*x_tmp[4]*xB_tmp[36] - 1.0*x_tmp[8]*xB_tmp[40] + x_tmp[8]*xB_tmp[43] - 1.0*x_tmp[14]*xB_tmp[47] + x_tmp[14]*xB_tmp[49] - 1.0*x_tmp[25]*xB_tmp[48] - 1.0*x_tmp[34]*xB_tmp[44] + x_tmp[25]*xB_tmp[60] + x_tmp[27]*xB_tmp[62] - 1.0*x_tmp[27]*xB_tmp[66] + x_tmp[34]*xB_tmp[69] + xB_tmp[41]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[37] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[42];
qBdot_tmp[2+ip*ny] = x_tmp[4]*xB_tmp[74] - 1.0*x_tmp[4]*xB_tmp[71] - 1.0*x_tmp[8]*xB_tmp[75] + x_tmp[8]*xB_tmp[78] - 1.0*x_tmp[14]*xB_tmp[82] + x_tmp[14]*xB_tmp[84] - 1.0*x_tmp[25]*xB_tmp[83] - 1.0*x_tmp[34]*xB_tmp[79] + x_tmp[25]*xB_tmp[95] + x_tmp[27]*xB_tmp[97] - 1.0*x_tmp[27]*xB_tmp[101] + x_tmp[34]*xB_tmp[104] + xB_tmp[76]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[72] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[77];
qBdot_tmp[3+ip*ny] = x_tmp[4]*xB_tmp[109] - 1.0*x_tmp[4]*xB_tmp[106] - 1.0*x_tmp[8]*xB_tmp[110] + x_tmp[8]*xB_tmp[113] - 1.0*x_tmp[14]*xB_tmp[117] + x_tmp[14]*xB_tmp[119] - 1.0*x_tmp[25]*xB_tmp[118] - 1.0*x_tmp[34]*xB_tmp[114] + x_tmp[25]*xB_tmp[130] + x_tmp[27]*xB_tmp[132] - 1.0*x_tmp[27]*xB_tmp[136] + x_tmp[34]*xB_tmp[139] + xB_tmp[111]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[107] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[112];
qBdot_tmp[4+ip*ny] = x_tmp[4]*xB_tmp[144] - 1.0*x_tmp[4]*xB_tmp[141] - 1.0*x_tmp[8]*xB_tmp[145] + x_tmp[8]*xB_tmp[148] - 1.0*x_tmp[14]*xB_tmp[152] + x_tmp[14]*xB_tmp[154] - 1.0*x_tmp[25]*xB_tmp[153] - 1.0*x_tmp[34]*xB_tmp[149] + x_tmp[25]*xB_tmp[165] + x_tmp[27]*xB_tmp[167] - 1.0*x_tmp[27]*xB_tmp[171] + x_tmp[34]*xB_tmp[174] + xB_tmp[146]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[142] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[147];
qBdot_tmp[5+ip*ny] = x_tmp[4]*xB_tmp[179] - 1.0*x_tmp[4]*xB_tmp[176] - 1.0*x_tmp[8]*xB_tmp[180] + x_tmp[8]*xB_tmp[183] - 1.0*x_tmp[14]*xB_tmp[187] + x_tmp[14]*xB_tmp[189] - 1.0*x_tmp[25]*xB_tmp[188] - 1.0*x_tmp[34]*xB_tmp[184] + x_tmp[25]*xB_tmp[200] + x_tmp[27]*xB_tmp[202] - 1.0*x_tmp[27]*xB_tmp[206] + x_tmp[34]*xB_tmp[209] + xB_tmp[181]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[177] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[182];
qBdot_tmp[6+ip*ny] = x_tmp[4]*xB_tmp[214] - 1.0*x_tmp[4]*xB_tmp[211] - 1.0*x_tmp[8]*xB_tmp[215] + x_tmp[8]*xB_tmp[218] - 1.0*x_tmp[14]*xB_tmp[222] + x_tmp[14]*xB_tmp[224] - 1.0*x_tmp[25]*xB_tmp[223] - 1.0*x_tmp[34]*xB_tmp[219] + x_tmp[25]*xB_tmp[235] + x_tmp[27]*xB_tmp[237] - 1.0*x_tmp[27]*xB_tmp[241] + x_tmp[34]*xB_tmp[244] + xB_tmp[216]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[212] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[217];
qBdot_tmp[7+ip*ny] = x_tmp[4]*xB_tmp[249] - 1.0*x_tmp[4]*xB_tmp[246] - 1.0*x_tmp[8]*xB_tmp[250] + x_tmp[8]*xB_tmp[253] - 1.0*x_tmp[14]*xB_tmp[257] + x_tmp[14]*xB_tmp[259] - 1.0*x_tmp[25]*xB_tmp[258] - 1.0*x_tmp[34]*xB_tmp[254] + x_tmp[25]*xB_tmp[270] + x_tmp[27]*xB_tmp[272] - 1.0*x_tmp[27]*xB_tmp[276] + x_tmp[34]*xB_tmp[279] + xB_tmp[251]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[247] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[252];
qBdot_tmp[8+ip*ny] = x_tmp[4]*xB_tmp[284] - 1.0*x_tmp[4]*xB_tmp[281] - 1.0*x_tmp[8]*xB_tmp[285] + x_tmp[8]*xB_tmp[288] - 1.0*x_tmp[14]*xB_tmp[292] + x_tmp[14]*xB_tmp[294] - 1.0*x_tmp[25]*xB_tmp[293] - 1.0*x_tmp[34]*xB_tmp[289] + x_tmp[25]*xB_tmp[305] + x_tmp[27]*xB_tmp[307] - 1.0*x_tmp[27]*xB_tmp[311] + x_tmp[34]*xB_tmp[314] + xB_tmp[286]*(x_tmp[4] + x_tmp[6] - 1.0*x_tmp[7]) - 1.0*(x_tmp[4] + 2.0*x_tmp[6])*xB_tmp[282] - 1.0*(x_tmp[4] - 2.0*x_tmp[7])*xB_tmp[287];

  } break;

  case 5: {
qBdot_tmp[0+ip*ny] = x_tmp[1]*xB_tmp[1] + x_tmp[5]*xB_tmp[5] + x_tmp[6]*xB_tmp[6] + x_tmp[9]*xB_tmp[9] + x_tmp[12]*xB_tmp[12] + x_tmp[13]*xB_tmp[13] + x_tmp[31]*xB_tmp[31] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[2];
qBdot_tmp[1+ip*ny] = x_tmp[1]*xB_tmp[36] + x_tmp[5]*xB_tmp[40] + x_tmp[6]*xB_tmp[41] + x_tmp[9]*xB_tmp[44] + x_tmp[12]*xB_tmp[47] + x_tmp[13]*xB_tmp[48] + x_tmp[31]*xB_tmp[66] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[37];
qBdot_tmp[2+ip*ny] = x_tmp[1]*xB_tmp[71] + x_tmp[5]*xB_tmp[75] + x_tmp[6]*xB_tmp[76] + x_tmp[9]*xB_tmp[79] + x_tmp[12]*xB_tmp[82] + x_tmp[13]*xB_tmp[83] + x_tmp[31]*xB_tmp[101] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[72];
qBdot_tmp[3+ip*ny] = x_tmp[1]*xB_tmp[106] + x_tmp[5]*xB_tmp[110] + x_tmp[6]*xB_tmp[111] + x_tmp[9]*xB_tmp[114] + x_tmp[12]*xB_tmp[117] + x_tmp[13]*xB_tmp[118] + x_tmp[31]*xB_tmp[136] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[107];
qBdot_tmp[4+ip*ny] = x_tmp[1]*xB_tmp[141] + x_tmp[5]*xB_tmp[145] + x_tmp[6]*xB_tmp[146] + x_tmp[9]*xB_tmp[149] + x_tmp[12]*xB_tmp[152] + x_tmp[13]*xB_tmp[153] + x_tmp[31]*xB_tmp[171] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[142];
qBdot_tmp[5+ip*ny] = x_tmp[1]*xB_tmp[176] + x_tmp[5]*xB_tmp[180] + x_tmp[6]*xB_tmp[181] + x_tmp[9]*xB_tmp[184] + x_tmp[12]*xB_tmp[187] + x_tmp[13]*xB_tmp[188] + x_tmp[31]*xB_tmp[206] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[177];
qBdot_tmp[6+ip*ny] = x_tmp[1]*xB_tmp[211] + x_tmp[5]*xB_tmp[215] + x_tmp[6]*xB_tmp[216] + x_tmp[9]*xB_tmp[219] + x_tmp[12]*xB_tmp[222] + x_tmp[13]*xB_tmp[223] + x_tmp[31]*xB_tmp[241] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[212];
qBdot_tmp[7+ip*ny] = x_tmp[1]*xB_tmp[246] + x_tmp[5]*xB_tmp[250] + x_tmp[6]*xB_tmp[251] + x_tmp[9]*xB_tmp[254] + x_tmp[12]*xB_tmp[257] + x_tmp[13]*xB_tmp[258] + x_tmp[31]*xB_tmp[276] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[247];
qBdot_tmp[8+ip*ny] = x_tmp[1]*xB_tmp[281] + x_tmp[5]*xB_tmp[285] + x_tmp[6]*xB_tmp[286] + x_tmp[9]*xB_tmp[289] + x_tmp[12]*xB_tmp[292] + x_tmp[13]*xB_tmp[293] + x_tmp[31]*xB_tmp[311] - 1.0*(x_tmp[1] - 2.0*x_tmp[2])*xB_tmp[282];

  } break;

  }
  }

  for (iyp=0; iyp<9*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*35);
x0_tmp[0] = k[29]*k[64];
x0_tmp[1] = k[6]*k[41];
x0_tmp[2] = k[34]*k[69];
x0_tmp[3] = k[4]*k[39];
x0_tmp[4] = k[5]*k[40];
x0_tmp[5] = k[31]*k[66];
x0_tmp[6] = k[33]*k[68];
x0_tmp[7] = k[32]*k[67];
x0_tmp[8] = k[30]*k[65];
x0_tmp[9] = k[28]*k[63];
x0_tmp[10] = k[22]*k[57];
x0_tmp[11] = k[11]*k[46];
x0_tmp[12] = k[13]*k[48];
x0_tmp[13] = k[19]*k[54];
x0_tmp[14] = k[12]*k[47];
x0_tmp[15] = (k[7] - 1.0)*((pow(p[6],2)) - 1.0*p[6]) + k[7]*k[42];
x0_tmp[16] = k[8]*k[43] + (k[8] - 1.0)*p[6]*p[7];
x0_tmp[17] = k[0]*k[35] - 1.0*(k[0] - 1.0)*p[6];
x0_tmp[18] = k[9]*k[44] - 1.0*(k[9] - 1.0)*p[6]*(p[6] + p[7] - 1.0);
x0_tmp[19] = (k[14] - 1.0)*((pow(p[7],2)) - 1.0*p[7]) + k[14]*k[49];
x0_tmp[20] = k[1]*k[36] - 1.0*(k[1] - 1.0)*p[7];
x0_tmp[21] = k[20]*k[55] + (k[20] - 1.0)*(p[6] + p[7])*(p[6] + p[7] - 1.0);
x0_tmp[22] = k[15]*k[50] - 1.0*(k[15] - 1.0)*p[7]*(p[6] + p[7] - 1.0);
x0_tmp[23] = k[16]*k[51];
x0_tmp[24] = k[10]*k[45];
x0_tmp[25] = k[18]*k[53];
x0_tmp[26] = k[17]*k[52];
x0_tmp[27] = k[23]*k[58];
x0_tmp[28] = (k[2] - 1.0)*(p[6] + p[7] - 1.0) + k[2]*k[37];
x0_tmp[29] = k[3]*k[38];
x0_tmp[30] = k[21]*k[56];
x0_tmp[31] = k[24]*k[59];
x0_tmp[32] = k[25]*k[60];
x0_tmp[33] = k[26]*k[61];
x0_tmp[34] = k[27]*k[62];
  
  
  return;
}


 int Jv__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(N_Vector v, N_Vector Jv, realtype t,
  	N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *v_tmp = N_VGetArrayPointer(v);
  double *Jv_tmp = N_VGetArrayPointer(Jv);
  memset(Jv_tmp,0,sizeof(double)*35);
Jv_tmp[0] = v_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*v_tmp[33];
Jv_tmp[1] = p[4]*v_tmp[4] - 1.0*p[5]*v_tmp[1];
Jv_tmp[2] = p[5]*v_tmp[1] - 2.0*p[5]*v_tmp[2] + p[4]*v_tmp[4] + 2.0*p[4]*v_tmp[6];
Jv_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[29];
Jv_tmp[4] = p[2]*v_tmp[3] - 1.0*p[4]*v_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[29];
Jv_tmp[5] = p[4]*v_tmp[8] - 1.0*p[5]*v_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[9];
Jv_tmp[6] = p[2]*v_tmp[5] - 1.0*p[4]*v_tmp[4] + p[4]*v_tmp[7] - 1.0*v_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[9];
Jv_tmp[7] = p[2]*v_tmp[3] + p[4]*v_tmp[4] + 2.0*p[2]*v_tmp[8] - 2.0*p[4]*v_tmp[7] + v_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*v_tmp[34];
Jv_tmp[8] = p[2]*v_tmp[0] - 1.0*p[4]*v_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[34] + 0.495*p[2]*(pow(p[3],2))*v_tmp[29];
Jv_tmp[9] = p[2]*v_tmp[31] + p[4]*v_tmp[34] - 1.0*v_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[10] = p[0]*v_tmp[26] - 1.0*p[1]*v_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[30];
Jv_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[24] - 0.0002*v_tmp[11];
Jv_tmp[12] = p[4]*v_tmp[14] - 1.0*(p[5] + 0.0002)*v_tmp[12];
Jv_tmp[13] = p[4]*v_tmp[25] + p[1]*v_tmp[31] - 1.0*v_tmp[13]*(p[0] + p[5]) + 0.0002*v_tmp[12];
Jv_tmp[14] = p[2]*v_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[24] - 1.0*(p[4] + 0.0002)*v_tmp[14];
Jv_tmp[15] = 0.0002*v_tmp[17] - 0.0004*v_tmp[15];
Jv_tmp[16] = p[1]*v_tmp[18] - 1.0*(p[0] + 0.0002)*v_tmp[16] + 0.0002*v_tmp[15] - 0.0002*v_tmp[17];
Jv_tmp[17] = -0.0002*v_tmp[17];
Jv_tmp[18] = p[0]*v_tmp[16] - 1.0*(p[1] + 0.0002)*v_tmp[18];
Jv_tmp[19] = p[0]*v_tmp[20] - 2.0*p[0]*v_tmp[19] + 2.0*p[1]*v_tmp[22] + p[1]*v_tmp[28] + 0.0004*v_tmp[16] + 0.0002*v_tmp[17];
Jv_tmp[20] = p[1]*v_tmp[28] - 1.0*p[0]*v_tmp[20] + 0.0002*v_tmp[17];
Jv_tmp[21] = p[0]*v_tmp[20] + 2.0*p[0]*v_tmp[22] - 2.0*p[1]*v_tmp[21] + p[1]*v_tmp[28];
Jv_tmp[22] = p[0]*v_tmp[19] - 1.0*p[0]*v_tmp[20] + p[1]*v_tmp[21] - 1.0*p[1]*v_tmp[28] - 1.0*v_tmp[22]*(p[0] + p[1]) + 0.0002*v_tmp[18];
Jv_tmp[23] = p[2]*v_tmp[22] + p[1]*v_tmp[30] + 0.0002*v_tmp[24] - 1.0*v_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[24] = p[2]*v_tmp[18] - 1.0*v_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
Jv_tmp[25] = p[1]*v_tmp[27] + p[2]*v_tmp[26] - 1.0*v_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[23] + 0.0002*v_tmp[14];
Jv_tmp[26] = p[1]*v_tmp[10] - 1.0*p[0]*v_tmp[26] + 0.0002*v_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[23];
Jv_tmp[27] = p[2]*v_tmp[10] + p[0]*v_tmp[25] - 1.0*v_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[30];
Jv_tmp[28] = p[0]*v_tmp[20] - 1.0*p[1]*v_tmp[28];
Jv_tmp[29] = p[2]*v_tmp[28] - 1.0*v_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[30] = p[0]*v_tmp[23] + p[2]*v_tmp[21] - 1.0*v_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[31] = p[0]*v_tmp[13] + p[4]*v_tmp[27] - 1.0*v_tmp[31]*(p[1] + p[5]);
Jv_tmp[32] = p[2]*v_tmp[28] + 2.0*p[2]*v_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*v_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*v_tmp[29];
Jv_tmp[33] = p[2]*v_tmp[10] - 1.0*v_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*v_tmp[32] - 1.1*p[2]*(pow(p[3],2))*v_tmp[29];
Jv_tmp[34] = p[2]*v_tmp[27] + p[2]*v_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*v_tmp[32] - 1.0*v_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*v_tmp[29];

  for (ix=0; ix<35; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(N_Vector vB, N_Vector JvB, realtype t,
  	N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *vB_tmp = N_VGetArrayPointer(vB);
  double *JvB_tmp = N_VGetArrayPointer(JvB);
  memset(JvB_tmp,0,sizeof(double)*35);
JvB_tmp[0] = -1.0*p[2]*vB_tmp[8];
JvB_tmp[1] = p[5]*vB_tmp[1] - 1.0*p[5]*vB_tmp[2];
JvB_tmp[2] = 2.0*p[5]*vB_tmp[2];
JvB_tmp[3] = - 1.0*p[2]*vB_tmp[4] - 1.0*p[2]*vB_tmp[7];
JvB_tmp[4] = p[4]*vB_tmp[4] - 1.0*p[4]*vB_tmp[2] - 1.0*p[4]*vB_tmp[1] + p[4]*vB_tmp[6] - 1.0*p[4]*vB_tmp[7];
JvB_tmp[5] = p[5]*vB_tmp[5] - 1.0*p[2]*vB_tmp[6];
JvB_tmp[6] = vB_tmp[6]*(p[4] + p[5]) - 2.0*p[4]*vB_tmp[2];
JvB_tmp[7] = 2.0*p[4]*vB_tmp[7] - 1.0*p[4]*vB_tmp[6];
JvB_tmp[8] = p[4]*vB_tmp[8] - 1.0*p[4]*vB_tmp[5] - 2.0*p[2]*vB_tmp[7];
JvB_tmp[9] = vB_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[5] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[6];
JvB_tmp[10] = p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[26] - 1.0*p[2]*vB_tmp[27] - 1.0*p[2]*vB_tmp[33];
JvB_tmp[11] = 0.0002*vB_tmp[11] - 1.0*p[2]*vB_tmp[14] - 0.0002*vB_tmp[26];
JvB_tmp[12] = (p[5] + 0.0002)*vB_tmp[12] - 0.0002*vB_tmp[13];
JvB_tmp[13] = vB_tmp[13]*(p[0] + p[5]) - 1.0*p[0]*vB_tmp[31];
JvB_tmp[14] = (p[4] + 0.0002)*vB_tmp[14] - 1.0*p[4]*vB_tmp[12] - 0.0002*vB_tmp[25];
JvB_tmp[15] = 0.0004*vB_tmp[15] - 0.0002*vB_tmp[16];
JvB_tmp[16] = (p[0] + 0.0002)*vB_tmp[16] - 1.0*p[0]*vB_tmp[18] - 0.0004*vB_tmp[19];
JvB_tmp[17] = 0.0002*vB_tmp[16] - 0.0002*vB_tmp[15] + 0.0002*vB_tmp[17] - 0.0002*vB_tmp[19] - 0.0002*vB_tmp[20];
JvB_tmp[18] = (p[1] + 0.0002)*vB_tmp[18] - 1.0*p[2]*vB_tmp[24] - 1.0*p[1]*vB_tmp[16] - 0.0002*vB_tmp[22];
JvB_tmp[19] = 2.0*p[0]*vB_tmp[19] - 1.0*p[0]*vB_tmp[22];
JvB_tmp[20] = p[0]*vB_tmp[20] - 1.0*p[0]*vB_tmp[19] - 1.0*p[0]*vB_tmp[21] + p[0]*vB_tmp[22] - 1.0*p[0]*vB_tmp[28];
JvB_tmp[21] = 2.0*p[1]*vB_tmp[21] - 1.0*p[1]*vB_tmp[22] - 1.0*p[2]*vB_tmp[30];
JvB_tmp[22] = vB_tmp[22]*(p[0] + p[1]) - 2.0*p[0]*vB_tmp[21] - 1.0*p[2]*vB_tmp[23] - 2.0*p[1]*vB_tmp[19];
JvB_tmp[23] = vB_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[25] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[26] - 1.0*p[0]*vB_tmp[30];
JvB_tmp[24] = vB_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[14] - 0.0002*vB_tmp[23] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[11];
JvB_tmp[25] = vB_tmp[25]*(p[0] + p[4]) - 1.0*p[0]*vB_tmp[27] - 1.0*p[4]*vB_tmp[13];
JvB_tmp[26] = p[0]*vB_tmp[26] - 1.0*p[0]*vB_tmp[10] - 1.0*p[2]*vB_tmp[25];
JvB_tmp[27] = vB_tmp[27]*(p[1] + p[4]) - 1.0*p[4]*vB_tmp[31] - 1.0*p[2]*vB_tmp[34] - 1.0*p[1]*vB_tmp[25];
JvB_tmp[28] = p[1]*vB_tmp[22] - 1.0*p[1]*vB_tmp[20] - 1.0*p[1]*vB_tmp[21] - 1.0*p[1]*vB_tmp[19] + p[1]*vB_tmp[28] - 1.0*p[2]*vB_tmp[29] - 1.0*p[2]*vB_tmp[32];
JvB_tmp[29] = vB_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*vB_tmp[7]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) - 1.0*vB_tmp[0]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[4] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[3] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*vB_tmp[32] - 0.495*p[2]*(pow(p[3],2))*vB_tmp[8] + 1.1*p[2]*(pow(p[3],2))*vB_tmp[33] + 0.9*p[2]*(pow(p[3],2))*vB_tmp[34];
JvB_tmp[30] = vB_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*vB_tmp[32] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[27] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[10] - 1.0*p[1]*vB_tmp[23];
JvB_tmp[31] = vB_tmp[31]*(p[1] + p[5]) - 1.0*p[1]*vB_tmp[13] - 1.0*p[2]*vB_tmp[9];
JvB_tmp[32] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*vB_tmp[32] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[34] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[33];
JvB_tmp[33] = vB_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*vB_tmp[34] - 1.0*(0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[8] - 1.0*(2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[0];
JvB_tmp[34] = vB_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[7] - 1.0*(1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[8] - 1.0*p[4]*vB_tmp[9];

  for (ix=0; ix<35; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  memset(J->data,0,sizeof(double)*1225);
J->data[8] = p[2];
J->data[36] = -1.0*p[5];
J->data[37] = p[5];
J->data[72] = -2.0*p[5];
J->data[109] = p[2];
J->data[112] = p[2];
J->data[141] = p[4];
J->data[142] = p[4];
J->data[144] = -1.0*p[4];
J->data[146] = -1.0*p[4];
J->data[147] = p[4];
J->data[180] = -1.0*p[5];
J->data[181] = p[2];
J->data[212] = 2.0*p[4];
J->data[216] = - 1.0*p[4] - 1.0*p[5];
J->data[251] = p[4];
J->data[252] = -2.0*p[4];
J->data[285] = p[4];
J->data[287] = 2.0*p[2];
J->data[288] = -1.0*p[4];
J->data[320] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[321] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[324] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[5] - 1.0*p[2]*(pow(p[3],2));
J->data[360] = -1.0*p[1];
J->data[376] = p[1];
J->data[377] = p[2];
J->data[383] = p[2];
J->data[396] = -0.0002;
J->data[399] = p[2];
J->data[411] = 0.0002;
J->data[432] = - 1.0*p[5] - 0.0002;
J->data[433] = 0.0002;
J->data[468] = - 1.0*p[0] - 1.0*p[5];
J->data[486] = p[0];
J->data[502] = p[4];
J->data[504] = - 1.0*p[4] - 0.0002;
J->data[515] = 0.0002;
J->data[540] = -0.0004;
J->data[541] = 0.0002;
J->data[576] = - 1.0*p[0] - 0.0002;
J->data[578] = p[0];
J->data[579] = 0.0004;
J->data[610] = 0.0002;
J->data[611] = -0.0002;
J->data[612] = -0.0002;
J->data[614] = 0.0002;
J->data[615] = 0.0002;
J->data[646] = p[1];
J->data[648] = - 1.0*p[1] - 0.0002;
J->data[652] = 0.0002;
J->data[654] = p[2];
J->data[684] = -2.0*p[0];
J->data[687] = p[0];
J->data[719] = p[0];
J->data[720] = -1.0*p[0];
J->data[721] = p[0];
J->data[722] = -1.0*p[0];
J->data[728] = p[0];
J->data[756] = -2.0*p[1];
J->data[757] = p[1];
J->data[765] = p[2];
J->data[789] = 2.0*p[1];
J->data[791] = 2.0*p[0];
J->data[792] = - 1.0*p[0] - 1.0*p[1];
J->data[793] = p[2];
J->data[828] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[0] - 1.0*p[2]*(pow(p[3],2));
J->data[830] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[831] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[835] = p[0];
J->data[851] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[854] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[863] = 0.0002;
J->data[864] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2)) - 0.0002;
J->data[888] = p[4];
J->data[900] = - 1.0*p[0] - 1.0*p[4];
J->data[902] = p[0];
J->data[920] = p[0];
J->data[935] = p[2];
J->data[936] = -1.0*p[0];
J->data[970] = p[1];
J->data[972] = - 1.0*p[1] - 1.0*p[4];
J->data[976] = p[4];
J->data[979] = p[2];
J->data[999] = p[1];
J->data[1000] = p[1];
J->data[1001] = p[1];
J->data[1002] = -1.0*p[1];
J->data[1008] = -1.0*p[1];
J->data[1009] = p[2];
J->data[1012] = p[2];
J->data[1015] = 1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[1018] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[1019] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[1022] = 1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[1023] = 0.495*p[2]*(pow(p[3],2));
J->data[1044] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[1047] = p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2];
J->data[1048] = -1.1*p[2]*(pow(p[3],2));
J->data[1049] = -0.9*p[2]*(pow(p[3],2));
J->data[1060] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[1073] = p[1];
J->data[1077] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[1080] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[1] - 1.0*p[2]*(pow(p[3],2));
J->data[1082] = 2.0*p[2];
J->data[1094] = p[2];
J->data[1098] = p[1];
J->data[1116] = - 1.0*p[1] - 1.0*p[5];
J->data[1152] = 2.0*(pow((p[3] - 1.0),2))*p[2] - 2.0*p[2]*(pow(p[3],2));
J->data[1153] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[1154] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[1155] = 2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3];
J->data[1163] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[1188] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[1189] = p[2];
J->data[1197] = 1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3];
J->data[1198] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[1199] = p[4];
J->data[1224] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[4] - 1.0*p[2]*(pow(p[3],2));

  for (iJ=0; iJ<1225; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(realtype t, N_Vector x,
  	N_Vector xdot, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(J);
  J->rowvals[0] = 8;
  J->rowvals[1] = 1;
  J->rowvals[2] = 2;
  J->rowvals[3] = 2;
  J->rowvals[4] = 4;
  J->rowvals[5] = 7;
  J->rowvals[6] = 1;
  J->rowvals[7] = 2;
  J->rowvals[8] = 4;
  J->rowvals[9] = 6;
  J->rowvals[10] = 7;
  J->rowvals[11] = 5;
  J->rowvals[12] = 6;
  J->rowvals[13] = 2;
  J->rowvals[14] = 6;
  J->rowvals[15] = 6;
  J->rowvals[16] = 7;
  J->rowvals[17] = 5;
  J->rowvals[18] = 7;
  J->rowvals[19] = 8;
  J->rowvals[20] = 5;
  J->rowvals[21] = 6;
  J->rowvals[22] = 9;
  J->rowvals[23] = 10;
  J->rowvals[24] = 26;
  J->rowvals[25] = 27;
  J->rowvals[26] = 33;
  J->rowvals[27] = 11;
  J->rowvals[28] = 14;
  J->rowvals[29] = 26;
  J->rowvals[30] = 12;
  J->rowvals[31] = 13;
  J->rowvals[32] = 13;
  J->rowvals[33] = 31;
  J->rowvals[34] = 12;
  J->rowvals[35] = 14;
  J->rowvals[36] = 25;
  J->rowvals[37] = 15;
  J->rowvals[38] = 16;
  J->rowvals[39] = 16;
  J->rowvals[40] = 18;
  J->rowvals[41] = 19;
  J->rowvals[42] = 15;
  J->rowvals[43] = 16;
  J->rowvals[44] = 17;
  J->rowvals[45] = 19;
  J->rowvals[46] = 20;
  J->rowvals[47] = 16;
  J->rowvals[48] = 18;
  J->rowvals[49] = 22;
  J->rowvals[50] = 24;
  J->rowvals[51] = 19;
  J->rowvals[52] = 22;
  J->rowvals[53] = 19;
  J->rowvals[54] = 20;
  J->rowvals[55] = 21;
  J->rowvals[56] = 22;
  J->rowvals[57] = 28;
  J->rowvals[58] = 21;
  J->rowvals[59] = 22;
  J->rowvals[60] = 30;
  J->rowvals[61] = 19;
  J->rowvals[62] = 21;
  J->rowvals[63] = 22;
  J->rowvals[64] = 23;
  J->rowvals[65] = 23;
  J->rowvals[66] = 25;
  J->rowvals[67] = 26;
  J->rowvals[68] = 30;
  J->rowvals[69] = 11;
  J->rowvals[70] = 14;
  J->rowvals[71] = 23;
  J->rowvals[72] = 24;
  J->rowvals[73] = 13;
  J->rowvals[74] = 25;
  J->rowvals[75] = 27;
  J->rowvals[76] = 10;
  J->rowvals[77] = 25;
  J->rowvals[78] = 26;
  J->rowvals[79] = 25;
  J->rowvals[80] = 27;
  J->rowvals[81] = 31;
  J->rowvals[82] = 34;
  J->rowvals[83] = 19;
  J->rowvals[84] = 20;
  J->rowvals[85] = 21;
  J->rowvals[86] = 22;
  J->rowvals[87] = 28;
  J->rowvals[88] = 29;
  J->rowvals[89] = 32;
  J->rowvals[90] = 0;
  J->rowvals[91] = 3;
  J->rowvals[92] = 4;
  J->rowvals[93] = 7;
  J->rowvals[94] = 8;
  J->rowvals[95] = 29;
  J->rowvals[96] = 32;
  J->rowvals[97] = 33;
  J->rowvals[98] = 34;
  J->rowvals[99] = 10;
  J->rowvals[100] = 23;
  J->rowvals[101] = 27;
  J->rowvals[102] = 30;
  J->rowvals[103] = 32;
  J->rowvals[104] = 9;
  J->rowvals[105] = 13;
  J->rowvals[106] = 31;
  J->rowvals[107] = 32;
  J->rowvals[108] = 33;
  J->rowvals[109] = 34;
  J->rowvals[110] = 0;
  J->rowvals[111] = 8;
  J->rowvals[112] = 33;
  J->rowvals[113] = 34;
  J->rowvals[114] = 7;
  J->rowvals[115] = 8;
  J->rowvals[116] = 9;
  J->rowvals[117] = 34;
  J->colptrs[0] = 0;
  J->colptrs[1] = 1;
  J->colptrs[2] = 3;
  J->colptrs[3] = 4;
  J->colptrs[4] = 6;
  J->colptrs[5] = 11;
  J->colptrs[6] = 13;
  J->colptrs[7] = 15;
  J->colptrs[8] = 17;
  J->colptrs[9] = 20;
  J->colptrs[10] = 23;
  J->colptrs[11] = 27;
  J->colptrs[12] = 30;
  J->colptrs[13] = 32;
  J->colptrs[14] = 34;
  J->colptrs[15] = 37;
  J->colptrs[16] = 39;
  J->colptrs[17] = 42;
  J->colptrs[18] = 47;
  J->colptrs[19] = 51;
  J->colptrs[20] = 53;
  J->colptrs[21] = 58;
  J->colptrs[22] = 61;
  J->colptrs[23] = 65;
  J->colptrs[24] = 69;
  J->colptrs[25] = 73;
  J->colptrs[26] = 76;
  J->colptrs[27] = 79;
  J->colptrs[28] = 83;
  J->colptrs[29] = 90;
  J->colptrs[30] = 99;
  J->colptrs[31] = 104;
  J->colptrs[32] = 107;
  J->colptrs[33] = 110;
  J->colptrs[34] = 114;
  J->colptrs[35] = 118;
J->data[0] = p[2];
J->data[1] = -1.0*p[5];
J->data[2] = p[5];
J->data[3] = -2.0*p[5];
J->data[4] = p[2];
J->data[5] = p[2];
J->data[6] = p[4];
J->data[7] = p[4];
J->data[8] = -1.0*p[4];
J->data[9] = -1.0*p[4];
J->data[10] = p[4];
J->data[11] = -1.0*p[5];
J->data[12] = p[2];
J->data[13] = 2.0*p[4];
J->data[14] = - 1.0*p[4] - 1.0*p[5];
J->data[15] = p[4];
J->data[16] = -2.0*p[4];
J->data[17] = p[4];
J->data[18] = 2.0*p[2];
J->data[19] = -1.0*p[4];
J->data[20] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[21] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[22] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[5] - 1.0*p[2]*(pow(p[3],2));
J->data[23] = -1.0*p[1];
J->data[24] = p[1];
J->data[25] = p[2];
J->data[26] = p[2];
J->data[27] = -0.0002;
J->data[28] = p[2];
J->data[29] = 0.0002;
J->data[30] = - 1.0*p[5] - 0.0002;
J->data[31] = 0.0002;
J->data[32] = - 1.0*p[0] - 1.0*p[5];
J->data[33] = p[0];
J->data[34] = p[4];
J->data[35] = - 1.0*p[4] - 0.0002;
J->data[36] = 0.0002;
J->data[37] = -0.0004;
J->data[38] = 0.0002;
J->data[39] = - 1.0*p[0] - 0.0002;
J->data[40] = p[0];
J->data[41] = 0.0004;
J->data[42] = 0.0002;
J->data[43] = -0.0002;
J->data[44] = -0.0002;
J->data[45] = 0.0002;
J->data[46] = 0.0002;
J->data[47] = p[1];
J->data[48] = - 1.0*p[1] - 0.0002;
J->data[49] = 0.0002;
J->data[50] = p[2];
J->data[51] = -2.0*p[0];
J->data[52] = p[0];
J->data[53] = p[0];
J->data[54] = -1.0*p[0];
J->data[55] = p[0];
J->data[56] = -1.0*p[0];
J->data[57] = p[0];
J->data[58] = -2.0*p[1];
J->data[59] = p[1];
J->data[60] = p[2];
J->data[61] = 2.0*p[1];
J->data[62] = 2.0*p[0];
J->data[63] = - 1.0*p[0] - 1.0*p[1];
J->data[64] = p[2];
J->data[65] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[0] - 1.0*p[2]*(pow(p[3],2));
J->data[66] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[67] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[68] = p[0];
J->data[69] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[70] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[71] = 0.0002;
J->data[72] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2)) - 0.0002;
J->data[73] = p[4];
J->data[74] = - 1.0*p[0] - 1.0*p[4];
J->data[75] = p[0];
J->data[76] = p[0];
J->data[77] = p[2];
J->data[78] = -1.0*p[0];
J->data[79] = p[1];
J->data[80] = - 1.0*p[1] - 1.0*p[4];
J->data[81] = p[4];
J->data[82] = p[2];
J->data[83] = p[1];
J->data[84] = p[1];
J->data[85] = p[1];
J->data[86] = -1.0*p[1];
J->data[87] = -1.0*p[1];
J->data[88] = p[2];
J->data[89] = p[2];
J->data[90] = 1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[91] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[92] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[93] = 1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[94] = 0.495*p[2]*(pow(p[3],2));
J->data[95] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[96] = p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2];
J->data[97] = -1.1*p[2]*(pow(p[3],2));
J->data[98] = -0.9*p[2]*(pow(p[3],2));
J->data[99] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[100] = p[1];
J->data[101] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[102] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[1] - 1.0*p[2]*(pow(p[3],2));
J->data[103] = 2.0*p[2];
J->data[104] = p[2];
J->data[105] = p[1];
J->data[106] = - 1.0*p[1] - 1.0*p[5];
J->data[107] = 2.0*(pow((p[3] - 1.0),2))*p[2] - 2.0*p[2]*(pow(p[3],2));
J->data[108] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[109] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[110] = 2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3];
J->data[111] = 0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3];
J->data[112] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[113] = p[2];
J->data[114] = 1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3];
J->data[115] = 1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3];
J->data[116] = p[4];
J->data[117] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[4] - 1.0*p[2]*(pow(p[3],2));
  return(0);
}


 int JBBand__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(long int N, realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, DlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
JB->data[29] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.705*p[2]*(pow(p[3],2));
JB->data[33] = 2.2*(p[3] - 1.0)*p[2]*p[3] - 2.2*p[2]*(pow(p[3],2));
JB->data[36] = p[5];
JB->data[39] = -1.0*p[4];
JB->data[71] = -1.0*p[5];
JB->data[72] = 2.0*p[5];
JB->data[74] = -1.0*p[4];
JB->data[76] = -2.0*p[4];
JB->data[134] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[143] = -1.0*p[2];
JB->data[144] = p[4];
JB->data[169] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[180] = p[5];
JB->data[183] = -1.0*p[4];
JB->data[184] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[214] = p[4];
JB->data[215] = -1.0*p[2];
JB->data[216] = p[4] + p[5];
JB->data[217] = -1.0*p[4];
JB->data[219] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[248] = -1.0*p[2];
JB->data[249] = -1.0*p[4];
JB->data[252] = 2.0*p[4];
JB->data[253] = -2.0*p[2];
JB->data[274] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 1.305*p[2]*(pow(p[3],2));
JB->data[279] = 1.8*(p[3] - 1.0)*p[2]*p[3] - 1.8*p[2]*(pow(p[3],2));
JB->data[280] = -1.0*p[2];
JB->data[288] = p[4];
JB->data[309] = -0.495*p[2]*(pow(p[3],2));
JB->data[313] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[314] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[324] = p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[346] = -1.0*p[2];
JB->data[349] = -1.0*p[4];
JB->data[360] = p[1];
JB->data[376] = -1.0*p[0];
JB->data[380] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[396] = 0.0002;
JB->data[409] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[432] = p[5] + 0.0002;
JB->data[434] = -1.0*p[4];
JB->data[467] = -0.0002;
JB->data[468] = p[0] + p[5];
JB->data[480] = -1.0*p[4];
JB->data[486] = -1.0*p[1];
JB->data[501] = -1.0*p[2];
JB->data[504] = p[4] + 0.0002;
JB->data[514] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[540] = 0.0004;
JB->data[542] = -0.0002;
JB->data[575] = -0.0002;
JB->data[576] = p[0] + 0.0002;
JB->data[577] = 0.0002;
JB->data[578] = -1.0*p[1];
JB->data[612] = 0.0002;
JB->data[646] = -1.0*p[0];
JB->data[648] = p[1] + 0.0002;
JB->data[681] = -0.0004;
JB->data[682] = -0.0002;
JB->data[684] = 2.0*p[0];
JB->data[685] = -1.0*p[0];
JB->data[687] = -2.0*p[1];
JB->data[693] = -1.0*p[1];
JB->data[717] = -0.0002;
JB->data[720] = p[0];
JB->data[728] = -1.0*p[1];
JB->data[755] = -1.0*p[0];
JB->data[756] = 2.0*p[1];
JB->data[757] = -2.0*p[0];
JB->data[763] = -1.0*p[1];
JB->data[788] = -0.0002;
JB->data[789] = -1.0*p[0];
JB->data[790] = p[0];
JB->data[791] = -1.0*p[1];
JB->data[792] = p[0] + p[1];
JB->data[798] = p[1];
JB->data[827] = -1.0*p[2];
JB->data[828] = p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[829] = -0.0002;
JB->data[835] = -1.0*p[1];
JB->data[858] = -1.0*p[2];
JB->data[864] = p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002;
JB->data[889] = -0.0002;
JB->data[898] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[900] = p[0] + p[4];
JB->data[901] = -1.0*p[2];
JB->data[902] = -1.0*p[1];
JB->data[920] = -1.0*p[1];
JB->data[921] = -0.0002;
JB->data[933] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[936] = p[0];
JB->data[955] = -1.0*p[2];
JB->data[970] = -1.0*p[0];
JB->data[972] = p[1] + p[4];
JB->data[975] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[1000] = -1.0*p[0];
JB->data[1008] = p[1];
JB->data[1043] = -1.0*p[2];
JB->data[1044] = p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[1071] = -1.0*p[2];
JB->data[1073] = -1.0*p[0];
JB->data[1080] = p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[1098] = -1.0*p[0];
JB->data[1112] = -1.0*p[4];
JB->data[1116] = p[1] + p[5];
JB->data[1148] = -1.0*p[2];
JB->data[1149] = - 1.0*p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[1150] = -2.0*p[2];
JB->data[1152] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[1165] = -1.0*p[2];
JB->data[1184] = 1.1*p[2]*(pow(p[3],2));
JB->data[1187] = 1.1*(p[3] - 1.0)*p[2]*p[3] - 1.1*p[2]*(pow(p[3],2));
JB->data[1188] = p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[1217] = -1.0*p[2];
JB->data[1219] = 0.9*p[2]*(pow(p[3],2));
JB->data[1222] = 0.9*(p[3] - 1.0)*p[2]*p[3] - 0.9*p[2]*(pow(p[3],2));
JB->data[1223] = -1.0*p[2];
JB->data[1224] = p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2];

  for (iJ=0; iJ<1225; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, SlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(JB);
  JB->rowvals[0] = 29;
  JB->rowvals[1] = 33;
  JB->rowvals[2] = 1;
  JB->rowvals[3] = 4;
  JB->rowvals[4] = 1;
  JB->rowvals[5] = 2;
  JB->rowvals[6] = 4;
  JB->rowvals[7] = 6;
  JB->rowvals[8] = 29;
  JB->rowvals[9] = 3;
  JB->rowvals[10] = 4;
  JB->rowvals[11] = 29;
  JB->rowvals[12] = 5;
  JB->rowvals[13] = 8;
  JB->rowvals[14] = 9;
  JB->rowvals[15] = 4;
  JB->rowvals[16] = 5;
  JB->rowvals[17] = 6;
  JB->rowvals[18] = 7;
  JB->rowvals[19] = 9;
  JB->rowvals[20] = 3;
  JB->rowvals[21] = 4;
  JB->rowvals[22] = 7;
  JB->rowvals[23] = 8;
  JB->rowvals[24] = 29;
  JB->rowvals[25] = 34;
  JB->rowvals[26] = 0;
  JB->rowvals[27] = 8;
  JB->rowvals[28] = 29;
  JB->rowvals[29] = 33;
  JB->rowvals[30] = 34;
  JB->rowvals[31] = 9;
  JB->rowvals[32] = 31;
  JB->rowvals[33] = 34;
  JB->rowvals[34] = 10;
  JB->rowvals[35] = 26;
  JB->rowvals[36] = 30;
  JB->rowvals[37] = 11;
  JB->rowvals[38] = 24;
  JB->rowvals[39] = 12;
  JB->rowvals[40] = 14;
  JB->rowvals[41] = 12;
  JB->rowvals[42] = 13;
  JB->rowvals[43] = 25;
  JB->rowvals[44] = 31;
  JB->rowvals[45] = 11;
  JB->rowvals[46] = 14;
  JB->rowvals[47] = 24;
  JB->rowvals[48] = 15;
  JB->rowvals[49] = 17;
  JB->rowvals[50] = 15;
  JB->rowvals[51] = 16;
  JB->rowvals[52] = 17;
  JB->rowvals[53] = 18;
  JB->rowvals[54] = 17;
  JB->rowvals[55] = 16;
  JB->rowvals[56] = 18;
  JB->rowvals[57] = 16;
  JB->rowvals[58] = 17;
  JB->rowvals[59] = 19;
  JB->rowvals[60] = 20;
  JB->rowvals[61] = 22;
  JB->rowvals[62] = 28;
  JB->rowvals[63] = 17;
  JB->rowvals[64] = 20;
  JB->rowvals[65] = 28;
  JB->rowvals[66] = 20;
  JB->rowvals[67] = 21;
  JB->rowvals[68] = 22;
  JB->rowvals[69] = 28;
  JB->rowvals[70] = 18;
  JB->rowvals[71] = 19;
  JB->rowvals[72] = 20;
  JB->rowvals[73] = 21;
  JB->rowvals[74] = 22;
  JB->rowvals[75] = 28;
  JB->rowvals[76] = 22;
  JB->rowvals[77] = 23;
  JB->rowvals[78] = 24;
  JB->rowvals[79] = 30;
  JB->rowvals[80] = 18;
  JB->rowvals[81] = 24;
  JB->rowvals[82] = 14;
  JB->rowvals[83] = 23;
  JB->rowvals[84] = 25;
  JB->rowvals[85] = 26;
  JB->rowvals[86] = 27;
  JB->rowvals[87] = 10;
  JB->rowvals[88] = 11;
  JB->rowvals[89] = 23;
  JB->rowvals[90] = 26;
  JB->rowvals[91] = 10;
  JB->rowvals[92] = 25;
  JB->rowvals[93] = 27;
  JB->rowvals[94] = 30;
  JB->rowvals[95] = 20;
  JB->rowvals[96] = 28;
  JB->rowvals[97] = 28;
  JB->rowvals[98] = 29;
  JB->rowvals[99] = 21;
  JB->rowvals[100] = 23;
  JB->rowvals[101] = 30;
  JB->rowvals[102] = 13;
  JB->rowvals[103] = 27;
  JB->rowvals[104] = 31;
  JB->rowvals[105] = 28;
  JB->rowvals[106] = 29;
  JB->rowvals[107] = 30;
  JB->rowvals[108] = 32;
  JB->rowvals[109] = 10;
  JB->rowvals[110] = 29;
  JB->rowvals[111] = 32;
  JB->rowvals[112] = 33;
  JB->rowvals[113] = 27;
  JB->rowvals[114] = 29;
  JB->rowvals[115] = 32;
  JB->rowvals[116] = 33;
  JB->rowvals[117] = 34;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 2;
  JB->colptrs[2] = 4;
  JB->colptrs[3] = 8;
  JB->colptrs[4] = 9;
  JB->colptrs[5] = 12;
  JB->colptrs[6] = 15;
  JB->colptrs[7] = 20;
  JB->colptrs[8] = 26;
  JB->colptrs[9] = 31;
  JB->colptrs[10] = 34;
  JB->colptrs[11] = 37;
  JB->colptrs[12] = 39;
  JB->colptrs[13] = 41;
  JB->colptrs[14] = 45;
  JB->colptrs[15] = 48;
  JB->colptrs[16] = 50;
  JB->colptrs[17] = 54;
  JB->colptrs[18] = 55;
  JB->colptrs[19] = 57;
  JB->colptrs[20] = 63;
  JB->colptrs[21] = 66;
  JB->colptrs[22] = 70;
  JB->colptrs[23] = 76;
  JB->colptrs[24] = 80;
  JB->colptrs[25] = 82;
  JB->colptrs[26] = 87;
  JB->colptrs[27] = 91;
  JB->colptrs[28] = 95;
  JB->colptrs[29] = 97;
  JB->colptrs[30] = 99;
  JB->colptrs[31] = 102;
  JB->colptrs[32] = 105;
  JB->colptrs[33] = 109;
  JB->colptrs[34] = 113;
  JB->colptrs[35] = 118;
  return(0);
}


 int sx__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(int Ns, realtype t, N_Vector x, N_Vector xdot,
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double *sxdot_tmp = N_VGetArrayPointer(sxdot);
  memset(sxdot_tmp,0,sizeof(double)*35);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + x_tmp[26] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12] - 1.0*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17] - 1.0*x_tmp[16];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18] + x_tmp[16];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17] - 2.0*x_tmp[19] + x_tmp[20];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17] - 1.0*x_tmp[20];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28] + x_tmp[20] + 2.0*x_tmp[22];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18] + x_tmp[19] - 1.0*x_tmp[20] - 1.0*x_tmp[22];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] - 1.0*x_tmp[23] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14] - 1.0*x_tmp[25];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] - 1.0*x_tmp[26] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] + x_tmp[25];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28] + x_tmp[20];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] + x_tmp[23] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]) + x_tmp[13];
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29];

  } break;

  case 1: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] - 1.0*x_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12] + x_tmp[31];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17] + x_tmp[18];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18] - 1.0*x_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17] + 2.0*x_tmp[22] + x_tmp[28];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17] + x_tmp[28];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28] - 2.0*x_tmp[21] + x_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18] + x_tmp[21] - 1.0*x_tmp[22] - 1.0*x_tmp[28];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] + x_tmp[30] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14] + x_tmp[27];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] + x_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] - 1.0*x_tmp[27];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28] - 1.0*x_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] - 1.0*x_tmp[30] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]) - 1.0*x_tmp[31];
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29];

  } break;

  case 2: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + 1.705*(pow(p[3],2))*x_tmp[29] + 2.2*(pow(p[3],2))*x_tmp[33] + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[3]*x_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6];
sxdot_tmp[3] = 1.1*(pow(p[3],2))*x_tmp[29] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[29] + x_tmp[3] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + 1.1*(pow(p[3],2))*x_tmp[9] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9] + 0.9*(pow(p[3],2))*x_tmp[9] + x_tmp[5] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[9];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 1.305*(pow(p[3],2))*x_tmp[29] + 1.8*(pow(p[3],2))*x_tmp[34] + x_tmp[3] + 2.0*x_tmp[8] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[3]*x_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + 0.495*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[33] + 1.1*(pow(p[3],2))*x_tmp[34] + x_tmp[0] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[34];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] + (pow((p[3] - 1.0),2))*x_tmp[9] - 1.0*(pow(p[3],2))*x_tmp[9] + x_tmp[31] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + 1.1*(pow(p[3],2))*x_tmp[30] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[30];
sxdot_tmp[11] = 1.1*(pow(p[3],2))*x_tmp[24] - 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[24];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] + 0.9*(pow(p[3],2))*x_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14] + x_tmp[11] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[24];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + (pow((p[3] - 1.0),2))*x_tmp[23] - 1.0*(pow(p[3],2))*x_tmp[23] + 0.0002*sx_tmp[24] + x_tmp[22] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] + (pow((p[3] - 1.0),2))*x_tmp[24] - 1.0*(pow(p[3],2))*x_tmp[24] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) + x_tmp[18];
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.9*(pow(p[3],2))*x_tmp[23] + 0.0002*sx_tmp[14] + x_tmp[26] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[23];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 1.1*(pow(p[3],2))*x_tmp[23] + 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] + 0.9*(pow(p[3],2))*x_tmp[30] + x_tmp[10] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[30];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] + (pow((p[3] - 1.0),2))*x_tmp[29] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(pow(p[3],2))*x_tmp[29] + x_tmp[28];
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] + (pow((p[3] - 1.0),2))*x_tmp[30] - 1.0*(pow(p[3],2))*x_tmp[30] + x_tmp[21] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]);
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] + (pow((p[3] - 1.0),2))*x_tmp[29] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[32] + (pow(p[3],2))*x_tmp[29] - 2.0*(pow(p[3],2))*x_tmp[32] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + x_tmp[28] + 2.0*x_tmp[30] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] + (pow((p[3] - 1.0),2))*x_tmp[33] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.1*(pow(p[3],2))*x_tmp[29] + 1.1*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[33] + x_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29] - 1.1*(p[3] - 1.0)*p[3]*x_tmp[32];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (pow((p[3] - 1.0),2))*x_tmp[34] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 0.9*(pow(p[3],2))*x_tmp[29] + 0.9*(pow(p[3],2))*x_tmp[32] - 1.0*(pow(p[3],2))*x_tmp[34] + x_tmp[27] + x_tmp[33] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29] - 0.9*(p[3] - 1.0)*p[3]*x_tmp[32];

  } break;

  case 3: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[29] - 2.2*(p[3] - 1.0)*p[2]*x_tmp[33] + 2.31*p[2]*p[3]*x_tmp[29] + 2.2*p[2]*p[3]*x_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[29] + 1.1*p[2]*p[3]*x_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[29] + 0.9*p[2]*p[3]*x_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[9] + 1.1*p[2]*p[3]*x_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[9] + 0.9*p[2]*p[3]*x_tmp[9];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[29] - 1.8*(p[3] - 1.0)*p[2]*x_tmp[34] + 1.71*p[2]*p[3]*x_tmp[29] + 1.8*p[2]*p[3]*x_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[33] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[34] + 0.99*p[2]*p[3]*x_tmp[29] + 0.9*p[2]*p[3]*x_tmp[33] + 1.1*p[2]*p[3]*x_tmp[34];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[9] - 2.0*p[2]*p[3]*x_tmp[9];
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[30] + 1.1*p[2]*p[3]*x_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[24] + 1.1*p[2]*p[3]*x_tmp[24];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[24] + 0.9*p[2]*p[3]*x_tmp[24];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[23] - 2.0*p[2]*p[3]*x_tmp[23];
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) + (2.0*p[3] - 2.0)*p[2]*x_tmp[24] - 2.0*p[2]*p[3]*x_tmp[24];
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[23] + 0.9*p[2]*p[3]*x_tmp[23];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[23] + 1.1*p[2]*p[3]*x_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[30] + 0.9*p[2]*p[3]*x_tmp[30];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[29] - 2.0*p[2]*p[3]*x_tmp[29];
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[30] - 2.0*p[2]*p[3]*x_tmp[30];
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]);
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29] + (2.0*p[3] - 2.0)*p[2]*x_tmp[29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[32] + 2.0*p[2]*p[3]*x_tmp[29] - 4.0*p[2]*p[3]*x_tmp[32];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] + (2.0*p[3] - 2.0)*p[2]*x_tmp[33] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29] - 1.1*(p[3] - 1.0)*p[2]*x_tmp[32] - 2.2*p[2]*p[3]*x_tmp[29] + 1.1*p[2]*p[3]*x_tmp[32] - 2.0*p[2]*p[3]*x_tmp[33];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[34] - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29] - 0.9*(p[3] - 1.0)*p[2]*x_tmp[32] - 1.8*p[2]*p[3]*x_tmp[29] + 0.9*p[2]*p[3]*x_tmp[32] - 2.0*p[2]*p[3]*x_tmp[34];

  } break;

  case 4: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1] + x_tmp[4];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6] + x_tmp[4] + 2.0*x_tmp[6];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29] - 1.0*x_tmp[4];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + x_tmp[8] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9] - 1.0*x_tmp[4] - 1.0*x_tmp[6] + x_tmp[7];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + x_tmp[4] - 2.0*x_tmp[7];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] - 1.0*x_tmp[8] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] + x_tmp[34] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12] + x_tmp[14];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12] + x_tmp[25];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14] - 1.0*x_tmp[14];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14] - 1.0*x_tmp[25];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30] - 1.0*x_tmp[27];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]) + x_tmp[27];
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*x_tmp[34] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29];

  } break;

  case 5: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1] - 1.0*x_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6] + x_tmp[1] - 2.0*x_tmp[2];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] - 1.0*x_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9] - 1.0*x_tmp[6];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] - 1.0*x_tmp[9] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12] - 1.0*x_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12] - 1.0*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]) - 1.0*x_tmp[31];
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29];

  } break;

  case 6: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]);
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29];

  } break;

  case 7: {
sxdot_tmp[0] = sx_tmp[29]*(1.705*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3]) + (2.2*p[2]*(pow(p[3],2)) - 2.2*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33];
sxdot_tmp[1] = p[4]*sx_tmp[4] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = p[5]*sx_tmp[1] - 2.0*p[5]*sx_tmp[2] + p[4]*sx_tmp[4] + 2.0*p[4]*sx_tmp[6];
sxdot_tmp[3] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[4] = p[2]*sx_tmp[3] - 1.0*p[4]*sx_tmp[4] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[29];
sxdot_tmp[5] = p[4]*sx_tmp[8] - 1.0*p[5]*sx_tmp[5] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[6] = p[2]*sx_tmp[5] - 1.0*p[4]*sx_tmp[4] + p[4]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[4] + p[5]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[9];
sxdot_tmp[7] = p[2]*sx_tmp[3] + p[4]*sx_tmp[4] + 2.0*p[2]*sx_tmp[8] - 2.0*p[4]*sx_tmp[7] + sx_tmp[29]*(1.305*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3]) + (1.8*p[2]*(pow(p[3],2)) - 1.8*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34];
sxdot_tmp[8] = p[2]*sx_tmp[0] - 1.0*p[4]*sx_tmp[8] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[33] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[34] + 0.495*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[9] = p[2]*sx_tmp[31] + p[4]*sx_tmp[34] - 1.0*sx_tmp[9]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[10] = p[0]*sx_tmp[26] - 1.0*p[1]*sx_tmp[10] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[11] = (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 0.0002*sx_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[14] - 1.0*(p[5] + 0.0002)*sx_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[25] + p[1]*sx_tmp[31] - 1.0*sx_tmp[13]*(p[0] + p[5]) + 0.0002*sx_tmp[12];
sxdot_tmp[14] = p[2]*sx_tmp[11] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[24] - 1.0*(p[4] + 0.0002)*sx_tmp[14];
sxdot_tmp[15] = 0.0002*sx_tmp[17] - 0.0004*sx_tmp[15];
sxdot_tmp[16] = p[1]*sx_tmp[18] - 1.0*(p[0] + 0.0002)*sx_tmp[16] + 0.0002*sx_tmp[15] - 0.0002*sx_tmp[17];
sxdot_tmp[17] = -0.0002*sx_tmp[17];
sxdot_tmp[18] = p[0]*sx_tmp[16] - 1.0*(p[1] + 0.0002)*sx_tmp[18];
sxdot_tmp[19] = p[0]*sx_tmp[20] - 2.0*p[0]*sx_tmp[19] + 2.0*p[1]*sx_tmp[22] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[16] + 0.0002*sx_tmp[17];
sxdot_tmp[20] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[20] + 0.0002*sx_tmp[17];
sxdot_tmp[21] = p[0]*sx_tmp[20] + 2.0*p[0]*sx_tmp[22] - 2.0*p[1]*sx_tmp[21] + p[1]*sx_tmp[28];
sxdot_tmp[22] = p[0]*sx_tmp[19] - 1.0*p[0]*sx_tmp[20] + p[1]*sx_tmp[21] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[22]*(p[0] + p[1]) + 0.0002*sx_tmp[18];
sxdot_tmp[23] = p[2]*sx_tmp[22] + p[1]*sx_tmp[30] + 0.0002*sx_tmp[24] - 1.0*sx_tmp[23]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[24] = p[2]*sx_tmp[18] - 1.0*sx_tmp[24]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[25] = p[1]*sx_tmp[27] + p[2]*sx_tmp[26] - 1.0*sx_tmp[25]*(p[0] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23] + 0.0002*sx_tmp[14];
sxdot_tmp[26] = p[1]*sx_tmp[10] - 1.0*p[0]*sx_tmp[26] + 0.0002*sx_tmp[11] + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[23];
sxdot_tmp[27] = p[2]*sx_tmp[10] + p[0]*sx_tmp[25] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[30];
sxdot_tmp[28] = p[0]*sx_tmp[20] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[28] - 1.0*sx_tmp[29]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[30] = p[0]*sx_tmp[23] + p[2]*sx_tmp[21] - 1.0*sx_tmp[30]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[31] = p[0]*sx_tmp[13] + p[4]*sx_tmp[27] - 1.0*sx_tmp[31]*(p[1] + p[5]);
sxdot_tmp[32] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[30] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[32] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[29];
sxdot_tmp[33] = p[2]*sx_tmp[10] - 1.0*sx_tmp[33]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (1.1*p[2]*(pow(p[3],2)) - 1.1*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.1*p[2]*(pow(p[3],2))*sx_tmp[29];
sxdot_tmp[34] = p[2]*sx_tmp[27] + p[2]*sx_tmp[33] + (0.9*p[2]*(pow(p[3],2)) - 0.9*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[32] - 1.0*sx_tmp[34]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.9*p[2]*(pow(p[3],2))*sx_tmp[29];

  } break;

  }
 for (ix=0; ix<35; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*35);
  switch (ip) {
  case 6: {
sx0_tmp[15] = (2.0*p[6] - 1.0)*(k[7] - 1.0);
sx0_tmp[16] = (k[8] - 1.0)*p[7];
sx0_tmp[17] = 1.0 - 1.0*k[0];
sx0_tmp[18] = - 1.0*(k[9] - 1.0)*(p[6] + p[7] - 1.0) - 1.0*(k[9] - 1.0)*p[6];
sx0_tmp[21] = (k[20] - 1.0)*(p[6] + p[7] - 1.0) + (k[20] - 1.0)*(p[6] + p[7]);
sx0_tmp[22] = -1.0*(k[15] - 1.0)*p[7];
sx0_tmp[28] = k[2] - 1.0;

  } break;

  case 7: {
sx0_tmp[16] = (k[8] - 1.0)*p[6];
sx0_tmp[18] = -1.0*(k[9] - 1.0)*p[6];
sx0_tmp[19] = (2.0*p[7] - 1.0)*(k[14] - 1.0);
sx0_tmp[20] = 1.0 - 1.0*k[1];
sx0_tmp[21] = (k[20] - 1.0)*(p[6] + p[7] - 1.0) + (k[20] - 1.0)*(p[6] + p[7]);
sx0_tmp[22] = - 1.0*(k[15] - 1.0)*(p[6] + p[7] - 1.0) - 1.0*(k[15] - 1.0)*p[7];
sx0_tmp[28] = k[2] - 1.0;

  } break;

  }

  return;
}


void y__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*29];
y[it+nt*1] = x[it+nt*3] + x[it+nt*4];
y[it+nt*2] = x[it+nt*1];
y[it+nt*3] = x[it+nt*32];
y[it+nt*4] = x[it+nt*33] + x[it+nt*34];
y[it+nt*5] = x[it+nt*9];
y[it+nt*6] = x[it+nt*0] + x[it+nt*7] + 2.0*x[it+nt*8];
y[it+nt*7] = x[it+nt*5] + x[it+nt*6];
y[it+nt*8] = x[it+nt*2];
    
    return;
}


void dydp__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*315);
dydx[6] = 1.0;
dydx[11] = 1.0;
dydx[26] = 1.0;
dydx[28] = 1.0;
dydx[37] = 1.0;
dydx[52] = 1.0;
dydx[61] = 1.0;
dydx[69] = 1.0;
dydx[78] = 2.0;
dydx[86] = 1.0;
dydx[261] = 1.0;
dydx[291] = 1.0;
dydx[301] = 1.0;
dydx[310] = 1.0;
  
  return;
}


void sy__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(29+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(4+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(32+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(33+np*nx)] + sx[it+nt*(34+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)] + 2.0*sx[it+nt*(8+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(5+np*nx)] + sx[it+nt*(6+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(2+np*nx)];

  } break;

  }
  
  return;
}
int root__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double dr_dp;
  switch (ip) {
  }
  return(dr_dp);
}
double s2root__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double ddr_dpdp;
  switch (ip) {
  }
  return(ddr_dpdp);
}
double srootval__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double dg_dp;
  switch (ip) {
  }
  return(dg_dp);
}
double s2rootval__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double ddg_dpdp;
  switch (ip) {
  }
  return(ddg_dpdp);
}
void deltadisc__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int idisc, N_Vector x, void *user_data){
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double deltadisc[35];
  memset(deltadisc,0,sizeof(double)*35);
  for(ix = 0; ix<35;ix++){;
  x_tmp[ix] += deltadisc[ix];
  };
}
void sdeltadisc__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
  int ix;
  int ip;
  UserData data = (UserData) user_data;
  int *plist = data->plist;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp;
  int np = *data->np;
  double deltadisc[35];
  double *sdeltadisc;
  memset(deltadisc,0,sizeof(double)*35);
  sdeltadisc = mxMalloc(sizeof(double)*35*np);
  memset(sdeltadisc,0,sizeof(double)*35*np);
  for (ip=0; ip<np; ip++) {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
     switch (plist[ip]) {
     }
  }
  for(ip = 0; ip<np;ip++){
      sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
      for(ix = 0; ix<35;ix++){
      sx_tmp[ix] += sdeltadisc[plist[ip]+np*ix];
     }
  }
  for(ix = 0; ix<35;ix++){
  x_tmp[ix] += deltadisc[ix];
  };
 mxFree(sdeltadisc);
}


void dxdotdp__sA_sPB_3o_2B_S_atS__T_dtS__B_atS(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(10+ip*nx)] = x[it+nt*26];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];
dxdotdp[(16+ip*nx)] = -1.0*x[it+nt*16];
dxdotdp[(18+ip*nx)] = x[it+nt*16];
dxdotdp[(19+ip*nx)] = x[it+nt*20] - 2.0*x[it+nt*19];
dxdotdp[(20+ip*nx)] = -1.0*x[it+nt*20];
dxdotdp[(21+ip*nx)] = x[it+nt*20] + 2.0*x[it+nt*22];
dxdotdp[(22+ip*nx)] = x[it+nt*19] - 1.0*x[it+nt*20] - 1.0*x[it+nt*22];
dxdotdp[(23+ip*nx)] = -1.0*x[it+nt*23];
dxdotdp[(25+ip*nx)] = -1.0*x[it+nt*25];
dxdotdp[(26+ip*nx)] = -1.0*x[it+nt*26];
dxdotdp[(27+ip*nx)] = x[it+nt*25];
dxdotdp[(28+ip*nx)] = x[it+nt*20];
dxdotdp[(30+ip*nx)] = x[it+nt*23];
dxdotdp[(31+ip*nx)] = x[it+nt*13];

  } break;

  case 1: {
dxdotdp[(10+ip*nx)] = -1.0*x[it+nt*10];
dxdotdp[(13+ip*nx)] = x[it+nt*31];
dxdotdp[(16+ip*nx)] = x[it+nt*18];
dxdotdp[(18+ip*nx)] = -1.0*x[it+nt*18];
dxdotdp[(19+ip*nx)] = 2.0*x[it+nt*22] + x[it+nt*28];
dxdotdp[(20+ip*nx)] = x[it+nt*28];
dxdotdp[(21+ip*nx)] = x[it+nt*28] - 2.0*x[it+nt*21];
dxdotdp[(22+ip*nx)] = x[it+nt*21] - 1.0*x[it+nt*22] - 1.0*x[it+nt*28];
dxdotdp[(23+ip*nx)] = x[it+nt*30];
dxdotdp[(25+ip*nx)] = x[it+nt*27];
dxdotdp[(26+ip*nx)] = x[it+nt*10];
dxdotdp[(27+ip*nx)] = -1.0*x[it+nt*27];
dxdotdp[(28+ip*nx)] = -1.0*x[it+nt*28];
dxdotdp[(30+ip*nx)] = -1.0*x[it+nt*30];
dxdotdp[(31+ip*nx)] = -1.0*x[it+nt*31];

  } break;

  case 2: {
dxdotdp[(0+ip*nx)] = 1.705*(pow(p[3],2))*x[it+nt*29] + 2.2*(pow(p[3],2))*x[it+nt*33] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*29] - 2.2*(p[3] - 1.0)*p[3]*x[it+nt*33];
dxdotdp[(3+ip*nx)] = 1.1*(pow(p[3],2))*x[it+nt*29] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*29];
dxdotdp[(4+ip*nx)] = 0.9*(pow(p[3],2))*x[it+nt*29] + x[it+nt*3] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*29];
dxdotdp[(5+ip*nx)] = 1.1*(pow(p[3],2))*x[it+nt*9] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*9];
dxdotdp[(6+ip*nx)] = 0.9*(pow(p[3],2))*x[it+nt*9] + x[it+nt*5] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*9];
dxdotdp[(7+ip*nx)] = 1.305*(pow(p[3],2))*x[it+nt*29] + 1.8*(pow(p[3],2))*x[it+nt*34] + x[it+nt*3] + 2.0*x[it+nt*8] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*29] - 1.8*(p[3] - 1.0)*p[3]*x[it+nt*34];
dxdotdp[(8+ip*nx)] = 0.495*(pow(p[3],2))*x[it+nt*29] + 0.9*(pow(p[3],2))*x[it+nt*33] + 1.1*(pow(p[3],2))*x[it+nt*34] + x[it+nt*0] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*33] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*34];
dxdotdp[(9+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*9] - 1.0*(pow(p[3],2))*x[it+nt*9] + x[it+nt*31];
dxdotdp[(10+ip*nx)] = 1.1*(pow(p[3],2))*x[it+nt*30] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*30];
dxdotdp[(11+ip*nx)] = 1.1*(pow(p[3],2))*x[it+nt*24] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*24];
dxdotdp[(14+ip*nx)] = 0.9*(pow(p[3],2))*x[it+nt*24] + x[it+nt*11] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*24];
dxdotdp[(23+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*23] - 1.0*(pow(p[3],2))*x[it+nt*23] + x[it+nt*22];
dxdotdp[(24+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*24] - 1.0*(pow(p[3],2))*x[it+nt*24] + x[it+nt*18];
dxdotdp[(25+ip*nx)] = 0.9*(pow(p[3],2))*x[it+nt*23] + x[it+nt*26] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*23];
dxdotdp[(26+ip*nx)] = 1.1*(pow(p[3],2))*x[it+nt*23] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*23];
dxdotdp[(27+ip*nx)] = 0.9*(pow(p[3],2))*x[it+nt*30] + x[it+nt*10] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*30];
dxdotdp[(29+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*29] - 1.0*(pow(p[3],2))*x[it+nt*29] + x[it+nt*28];
dxdotdp[(30+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*30] - 1.0*(pow(p[3],2))*x[it+nt*30] + x[it+nt*21];
dxdotdp[(32+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*29] + 2.0*(pow((p[3] - 1.0),2))*x[it+nt*32] + (pow(p[3],2))*x[it+nt*29] - 2.0*(pow(p[3],2))*x[it+nt*32] + x[it+nt*28] + 2.0*x[it+nt*30];
dxdotdp[(33+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*33] - 1.1*(pow(p[3],2))*x[it+nt*29] + 1.1*(pow(p[3],2))*x[it+nt*32] - 1.0*(pow(p[3],2))*x[it+nt*33] + x[it+nt*10] - 1.1*(p[3] - 1.0)*p[3]*x[it+nt*32];
dxdotdp[(34+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*34] - 0.9*(pow(p[3],2))*x[it+nt*29] + 0.9*(pow(p[3],2))*x[it+nt*32] - 1.0*(pow(p[3],2))*x[it+nt*34] + x[it+nt*27] + x[it+nt*33] - 0.9*(p[3] - 1.0)*p[3]*x[it+nt*32];

  } break;

  case 3: {
dxdotdp[(0+ip*nx)] = 2.31*p[2]*p[3]*x[it+nt*29] - 2.2*(p[3] - 1.0)*p[2]*x[it+nt*33] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*29] + 2.2*p[2]*p[3]*x[it+nt*33];
dxdotdp[(3+ip*nx)] = 1.1*p[2]*p[3]*x[it+nt*29] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*29];
dxdotdp[(4+ip*nx)] = 0.9*p[2]*p[3]*x[it+nt*29] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*29];
dxdotdp[(5+ip*nx)] = 1.1*p[2]*p[3]*x[it+nt*9] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*9];
dxdotdp[(6+ip*nx)] = 0.9*p[2]*p[3]*x[it+nt*9] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*9];
dxdotdp[(7+ip*nx)] = 1.71*p[2]*p[3]*x[it+nt*29] - 1.8*(p[3] - 1.0)*p[2]*x[it+nt*34] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*29] + 1.8*p[2]*p[3]*x[it+nt*34];
dxdotdp[(8+ip*nx)] = 0.99*p[2]*p[3]*x[it+nt*29] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*34] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*33] + 0.9*p[2]*p[3]*x[it+nt*33] + 1.1*p[2]*p[3]*x[it+nt*34];
dxdotdp[(9+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*9] - 2.0*p[2]*p[3]*x[it+nt*9];
dxdotdp[(10+ip*nx)] = 1.1*p[2]*p[3]*x[it+nt*30] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*30];
dxdotdp[(11+ip*nx)] = 1.1*p[2]*p[3]*x[it+nt*24] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*24];
dxdotdp[(14+ip*nx)] = 0.9*p[2]*p[3]*x[it+nt*24] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*24];
dxdotdp[(23+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*23] - 2.0*p[2]*p[3]*x[it+nt*23];
dxdotdp[(24+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*24] - 2.0*p[2]*p[3]*x[it+nt*24];
dxdotdp[(25+ip*nx)] = 0.9*p[2]*p[3]*x[it+nt*23] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*23];
dxdotdp[(26+ip*nx)] = 1.1*p[2]*p[3]*x[it+nt*23] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*23];
dxdotdp[(27+ip*nx)] = 0.9*p[2]*p[3]*x[it+nt*30] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*30];
dxdotdp[(29+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*29] - 2.0*p[2]*p[3]*x[it+nt*29];
dxdotdp[(30+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*30] - 2.0*p[2]*p[3]*x[it+nt*30];
dxdotdp[(32+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*29] + 2.0*(2.0*p[3] - 2.0)*p[2]*x[it+nt*32] + 2.0*p[2]*p[3]*x[it+nt*29] - 4.0*p[2]*p[3]*x[it+nt*32];
dxdotdp[(33+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*33] - 1.1*(p[3] - 1.0)*p[2]*x[it+nt*32] - 2.2*p[2]*p[3]*x[it+nt*29] + 1.1*p[2]*p[3]*x[it+nt*32] - 2.0*p[2]*p[3]*x[it+nt*33];
dxdotdp[(34+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*34] - 0.9*(p[3] - 1.0)*p[2]*x[it+nt*32] - 1.8*p[2]*p[3]*x[it+nt*29] + 0.9*p[2]*p[3]*x[it+nt*32] - 2.0*p[2]*p[3]*x[it+nt*34];

  } break;

  case 4: {
dxdotdp[(1+ip*nx)] = x[it+nt*4];
dxdotdp[(2+ip*nx)] = x[it+nt*4] + 2.0*x[it+nt*6];
dxdotdp[(4+ip*nx)] = -1.0*x[it+nt*4];
dxdotdp[(5+ip*nx)] = x[it+nt*8];
dxdotdp[(6+ip*nx)] = x[it+nt*7] - 1.0*x[it+nt*6] - 1.0*x[it+nt*4];
dxdotdp[(7+ip*nx)] = x[it+nt*4] - 2.0*x[it+nt*7];
dxdotdp[(8+ip*nx)] = -1.0*x[it+nt*8];
dxdotdp[(9+ip*nx)] = x[it+nt*34];
dxdotdp[(12+ip*nx)] = x[it+nt*14];
dxdotdp[(13+ip*nx)] = x[it+nt*25];
dxdotdp[(14+ip*nx)] = -1.0*x[it+nt*14];
dxdotdp[(25+ip*nx)] = -1.0*x[it+nt*25];
dxdotdp[(27+ip*nx)] = -1.0*x[it+nt*27];
dxdotdp[(31+ip*nx)] = x[it+nt*27];
dxdotdp[(34+ip*nx)] = -1.0*x[it+nt*34];

  } break;

  case 5: {
dxdotdp[(1+ip*nx)] = -1.0*x[it+nt*1];
dxdotdp[(2+ip*nx)] = x[it+nt*1] - 2.0*x[it+nt*2];
dxdotdp[(5+ip*nx)] = -1.0*x[it+nt*5];
dxdotdp[(6+ip*nx)] = -1.0*x[it+nt*6];
dxdotdp[(9+ip*nx)] = -1.0*x[it+nt*9];
dxdotdp[(12+ip*nx)] = -1.0*x[it+nt*12];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];
dxdotdp[(31+ip*nx)] = -1.0*x[it+nt*31];

  } break;

  }
  }
  
  return;
}
