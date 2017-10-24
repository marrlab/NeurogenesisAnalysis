#include "_sA_sPB_3o_2B_S_atS__T_atS__B_dtS.h"
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


 int xdot__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(realtype t, N_Vector x, N_Vector xdot, void *user_data)
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
xdot_tmp[0] = 0.55*p[2]*x_tmp[29] - 0.0002*x_tmp[0] - 1.0*p[2]*(pow(p[3],2))*x_tmp[0] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[0];
xdot_tmp[1] = p[4]*x_tmp[3] - 1.0*p[5]*x_tmp[1] - 0.0002*x_tmp[1];
xdot_tmp[2] = p[1]*x_tmp[33] - 1.0*p[5]*x_tmp[2] - 1.0*p[0]*x_tmp[2] + p[4]*x_tmp[34] + 0.0002*x_tmp[1];
xdot_tmp[3] = 0.45*p[2]*x_tmp[29] - 1.0*p[4]*x_tmp[3] - 0.0002*x_tmp[3] + 2.0*p[2]*(pow(p[3],2))*x_tmp[0] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[0];
xdot_tmp[4] = 0.0002*x_tmp[6] - 0.0004*x_tmp[4];
xdot_tmp[5] = p[1]*x_tmp[7] - 1.0*p[0]*x_tmp[5] + 0.0002*x_tmp[4] - 0.0002*x_tmp[5] - 0.0002*x_tmp[6];
xdot_tmp[6] = -0.0002*x_tmp[6];
xdot_tmp[7] = p[0]*x_tmp[5] - 1.0*p[1]*x_tmp[7] - 0.0002*x_tmp[7];
xdot_tmp[8] = p[0]*x_tmp[9] - 2.0*p[0]*x_tmp[8] + 2.0*p[1]*x_tmp[11] + p[1]*x_tmp[28] + 0.0004*x_tmp[5] + 0.0002*x_tmp[6];
xdot_tmp[9] = p[1]*x_tmp[28] - 1.0*p[0]*x_tmp[9] + 0.0002*x_tmp[6];
xdot_tmp[10] = p[0]*x_tmp[9] + 2.0*p[0]*x_tmp[11] - 2.0*p[1]*x_tmp[10] + p[1]*x_tmp[28];
xdot_tmp[11] = p[0]*x_tmp[8] - 1.0*p[0]*x_tmp[9] - 1.0*p[0]*x_tmp[11] + p[1]*x_tmp[10] - 1.0*p[1]*x_tmp[11] - 1.0*p[1]*x_tmp[28] + 0.0002*x_tmp[7];
xdot_tmp[12] = p[0]*x_tmp[31] - 1.0*p[1]*x_tmp[12] + 0.55*p[2]*x_tmp[32] - 1.0*p[2]*(pow(p[3],2))*x_tmp[12] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[12];
xdot_tmp[13] = p[4]*x_tmp[19] - 1.0*p[5]*x_tmp[13];
xdot_tmp[14] = p[5]*x_tmp[13] - 2.0*p[5]*x_tmp[14] + 2.0*p[4]*x_tmp[16] + p[4]*x_tmp[19];
xdot_tmp[15] = 0.55*p[2]*x_tmp[22] - 1.0*p[5]*x_tmp[15] + p[4]*x_tmp[20] - 1.0*p[2]*(pow(p[3],2))*x_tmp[15] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[15];
xdot_tmp[16] = 0.45*p[2]*x_tmp[22] - 1.0*p[5]*x_tmp[16] - 1.0*p[4]*x_tmp[19] - 1.0*p[4]*x_tmp[16] + p[4]*x_tmp[21] + 2.0*p[2]*(pow(p[3],2))*x_tmp[15] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[15];
xdot_tmp[17] = 0.55*p[2]*x_tmp[23] + 1.1*p[2]*x_tmp[26] - 2.0*p[2]*(pow(p[3],2))*x_tmp[17] + p[2]*(pow(p[3],2))*x_tmp[18] + 2.0*(pow((p[3] - 1.0),2))*p[2]*x_tmp[17] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[18];
xdot_tmp[18] = 0.55*p[2]*x_tmp[23] - 1.0*p[2]*(pow(p[3],2))*x_tmp[18] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[18];
xdot_tmp[19] = 0.45*p[2]*x_tmp[23] - 1.0*p[4]*x_tmp[19] + 2.0*p[2]*(pow(p[3],2))*x_tmp[18] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[18];
xdot_tmp[20] = 0.55*p[2]*x_tmp[25] - 1.0*p[4]*x_tmp[20] + 0.45*p[2]*x_tmp[26] + 2.0*p[2]*(pow(p[3],2))*x_tmp[17] - 2.0*p[2]*(pow(p[3],2))*x_tmp[18] - 1.0*p[2]*(pow(p[3],2))*x_tmp[20] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[20] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[17];
xdot_tmp[21] = p[4]*x_tmp[19] + 0.45*p[2]*x_tmp[23] - 2.0*p[4]*x_tmp[21] + 0.9*p[2]*x_tmp[25] + 4.0*p[2]*(pow(p[3],2))*x_tmp[18] + 4.0*p[2]*(pow(p[3],2))*x_tmp[20] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[20];
xdot_tmp[22] = p[4]*x_tmp[25] - 1.0*p[5]*x_tmp[22] + p[2]*x_tmp[33];
xdot_tmp[23] = p[2]*x_tmp[28];
xdot_tmp[24] = p[2]*x_tmp[28] + 2.0*p[2]*x_tmp[32];
xdot_tmp[25] = 0.45*p[2]*x_tmp[24] + p[2]*x_tmp[27] - 1.0*p[4]*x_tmp[25] + 2.0*p[2]*(pow(p[3],2))*x_tmp[26] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[26];
xdot_tmp[26] = p[2]*x_tmp[12] + 0.55*p[2]*x_tmp[24] - 1.0*p[2]*(pow(p[3],2))*x_tmp[26] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[26];
xdot_tmp[27] = p[0]*x_tmp[34] - 1.0*p[4]*x_tmp[27] - 1.0*p[1]*x_tmp[27] + 0.45*p[2]*x_tmp[32] + 2.0*p[2]*(pow(p[3],2))*x_tmp[12] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[12];
xdot_tmp[28] = p[0]*x_tmp[9] - 1.0*p[1]*x_tmp[28];
xdot_tmp[29] = p[2]*x_tmp[7] - 0.0002*x_tmp[29];
xdot_tmp[30] = p[2]*x_tmp[11] - 1.0*p[0]*x_tmp[30] + p[1]*x_tmp[32] + 0.0002*x_tmp[29];
xdot_tmp[31] = p[1]*x_tmp[12] - 1.0*p[0]*x_tmp[31] + 0.55*p[2]*x_tmp[30] + 0.0002*x_tmp[0] - 1.0*p[2]*(pow(p[3],2))*x_tmp[31] + (pow((p[3] - 1.0),2))*p[2]*x_tmp[31];
xdot_tmp[32] = p[2]*x_tmp[10] + p[0]*x_tmp[30] - 1.0*p[1]*x_tmp[32];
xdot_tmp[33] = p[0]*x_tmp[2] + p[4]*x_tmp[27] - 1.0*p[1]*x_tmp[33] - 1.0*p[5]*x_tmp[33];
xdot_tmp[34] = p[1]*x_tmp[27] + 0.45*p[2]*x_tmp[30] - 1.0*p[0]*x_tmp[34] - 1.0*p[4]*x_tmp[34] + 0.0002*x_tmp[3] + 2.0*p[2]*(pow(p[3],2))*x_tmp[31] - 2.0*(p[3] - 1.0)*p[2]*p[3]*x_tmp[31];

  for (ix=0; ix<35; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
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
xBdot_tmp[0] = xB_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[31] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[3];
xBdot_tmp[1] = (p[5] + 0.0002)*xB_tmp[1] - 0.0002*xB_tmp[2];
xBdot_tmp[2] = xB_tmp[2]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[33];
xBdot_tmp[3] = (p[4] + 0.0002)*xB_tmp[3] - 1.0*p[4]*xB_tmp[1] - 0.0002*xB_tmp[34];
xBdot_tmp[4] = 0.0004*xB_tmp[4] - 0.0002*xB_tmp[5];
xBdot_tmp[5] = (p[0] + 0.0002)*xB_tmp[5] - 1.0*p[0]*xB_tmp[7] - 0.0004*xB_tmp[8];
xBdot_tmp[6] = 0.0002*xB_tmp[5] - 0.0002*xB_tmp[4] + 0.0002*xB_tmp[6] - 0.0002*xB_tmp[8] - 0.0002*xB_tmp[9];
xBdot_tmp[7] = (p[1] + 0.0002)*xB_tmp[7] - 1.0*p[2]*xB_tmp[29] - 1.0*p[1]*xB_tmp[5] - 0.0002*xB_tmp[11];
xBdot_tmp[8] = 2.0*p[0]*xB_tmp[8] - 1.0*p[0]*xB_tmp[11];
xBdot_tmp[9] = p[0]*xB_tmp[9] - 1.0*p[0]*xB_tmp[8] - 1.0*p[0]*xB_tmp[10] + p[0]*xB_tmp[11] - 1.0*p[0]*xB_tmp[28];
xBdot_tmp[10] = 2.0*p[1]*xB_tmp[10] - 1.0*p[1]*xB_tmp[11] - 1.0*p[2]*xB_tmp[32];
xBdot_tmp[11] = xB_tmp[11]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[10] - 1.0*p[2]*xB_tmp[30] - 2.0*p[1]*xB_tmp[8];
xBdot_tmp[12] = xB_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[26] - 1.0*p[1]*xB_tmp[31] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[27];
xBdot_tmp[13] = p[5]*xB_tmp[13] - 1.0*p[5]*xB_tmp[14];
xBdot_tmp[14] = 2.0*p[5]*xB_tmp[14];
xBdot_tmp[15] = xB_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[16];
xBdot_tmp[16] = xB_tmp[16]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[14];
xBdot_tmp[17] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[17] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[20];
xBdot_tmp[18] = xB_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[19] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[21] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[17] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[20];
xBdot_tmp[19] = p[4]*xB_tmp[16] - 1.0*p[4]*xB_tmp[14] - 1.0*p[4]*xB_tmp[13] + p[4]*xB_tmp[19] - 1.0*p[4]*xB_tmp[21];
xBdot_tmp[20] = xB_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[15] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[21];
xBdot_tmp[21] = 2.0*p[4]*xB_tmp[21] - 1.0*p[4]*xB_tmp[16];
xBdot_tmp[22] = p[5]*xB_tmp[22] - 0.45*p[2]*xB_tmp[16] - 0.55*p[2]*xB_tmp[15];
xBdot_tmp[23] = - 0.55*p[2]*xB_tmp[17] - 0.55*p[2]*xB_tmp[18] - 0.45*p[2]*xB_tmp[19] - 0.45*p[2]*xB_tmp[21];
xBdot_tmp[24] = - 0.45*p[2]*xB_tmp[25] - 0.55*p[2]*xB_tmp[26];
xBdot_tmp[25] = p[4]*xB_tmp[25] - 0.9*p[2]*xB_tmp[21] - 1.0*p[4]*xB_tmp[22] - 0.55*p[2]*xB_tmp[20];
xBdot_tmp[26] = xB_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[20] - 1.1*p[2]*xB_tmp[17] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[25];
xBdot_tmp[27] = xB_tmp[27]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[34] - 1.0*p[4]*xB_tmp[33] - 1.0*p[2]*xB_tmp[25];
xBdot_tmp[28] = p[1]*xB_tmp[11] - 1.0*p[1]*xB_tmp[9] - 1.0*p[1]*xB_tmp[10] - 1.0*p[1]*xB_tmp[8] - 1.0*p[2]*xB_tmp[23] - 1.0*p[2]*xB_tmp[24] + p[1]*xB_tmp[28];
xBdot_tmp[29] = 0.0002*xB_tmp[29] - 0.45*p[2]*xB_tmp[3] - 0.55*p[2]*xB_tmp[0] - 0.0002*xB_tmp[30];
xBdot_tmp[30] = p[0]*xB_tmp[30] - 1.0*p[0]*xB_tmp[32] - 0.55*p[2]*xB_tmp[31] - 0.45*p[2]*xB_tmp[34];
xBdot_tmp[31] = xB_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[12] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[34];
xBdot_tmp[32] = p[1]*xB_tmp[32] - 2.0*p[2]*xB_tmp[24] - 0.45*p[2]*xB_tmp[27] - 1.0*p[1]*xB_tmp[30] - 0.55*p[2]*xB_tmp[12];
xBdot_tmp[33] = xB_tmp[33]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[22] - 1.0*p[1]*xB_tmp[2];
xBdot_tmp[34] = xB_tmp[34]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[27] - 1.0*p[4]*xB_tmp[2];
xBdot_tmp[35] = xB_tmp[35]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[66] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[38];
xBdot_tmp[36] = (p[5] + 0.0002)*xB_tmp[36] - 0.0002*xB_tmp[37];
xBdot_tmp[37] = xB_tmp[37]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[68];
xBdot_tmp[38] = (p[4] + 0.0002)*xB_tmp[38] - 1.0*p[4]*xB_tmp[36] - 0.0002*xB_tmp[69];
xBdot_tmp[39] = 0.0004*xB_tmp[39] - 0.0002*xB_tmp[40];
xBdot_tmp[40] = (p[0] + 0.0002)*xB_tmp[40] - 1.0*p[0]*xB_tmp[42] - 0.0004*xB_tmp[43];
xBdot_tmp[41] = 0.0002*xB_tmp[40] - 0.0002*xB_tmp[39] + 0.0002*xB_tmp[41] - 0.0002*xB_tmp[43] - 0.0002*xB_tmp[44];
xBdot_tmp[42] = (p[1] + 0.0002)*xB_tmp[42] - 1.0*p[2]*xB_tmp[64] - 1.0*p[1]*xB_tmp[40] - 0.0002*xB_tmp[46];
xBdot_tmp[43] = 2.0*p[0]*xB_tmp[43] - 1.0*p[0]*xB_tmp[46];
xBdot_tmp[44] = p[0]*xB_tmp[44] - 1.0*p[0]*xB_tmp[43] - 1.0*p[0]*xB_tmp[45] + p[0]*xB_tmp[46] - 1.0*p[0]*xB_tmp[63];
xBdot_tmp[45] = 2.0*p[1]*xB_tmp[45] - 1.0*p[1]*xB_tmp[46] - 1.0*p[2]*xB_tmp[67];
xBdot_tmp[46] = xB_tmp[46]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[45] - 1.0*p[2]*xB_tmp[65] - 2.0*p[1]*xB_tmp[43];
xBdot_tmp[47] = xB_tmp[47]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[61] - 1.0*p[1]*xB_tmp[66] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[62];
xBdot_tmp[48] = p[5]*xB_tmp[48] - 1.0*p[5]*xB_tmp[49];
xBdot_tmp[49] = 2.0*p[5]*xB_tmp[49];
xBdot_tmp[50] = xB_tmp[50]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[51];
xBdot_tmp[51] = xB_tmp[51]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[49];
xBdot_tmp[52] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[52] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[55];
xBdot_tmp[53] = xB_tmp[53]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[54] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[56] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[52] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[55];
xBdot_tmp[54] = p[4]*xB_tmp[51] - 1.0*p[4]*xB_tmp[49] - 1.0*p[4]*xB_tmp[48] + p[4]*xB_tmp[54] - 1.0*p[4]*xB_tmp[56];
xBdot_tmp[55] = xB_tmp[55]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[50] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[56];
xBdot_tmp[56] = 2.0*p[4]*xB_tmp[56] - 1.0*p[4]*xB_tmp[51];
xBdot_tmp[57] = p[5]*xB_tmp[57] - 0.45*p[2]*xB_tmp[51] - 0.55*p[2]*xB_tmp[50];
xBdot_tmp[58] = - 0.55*p[2]*xB_tmp[52] - 0.55*p[2]*xB_tmp[53] - 0.45*p[2]*xB_tmp[54] - 0.45*p[2]*xB_tmp[56];
xBdot_tmp[59] = - 0.45*p[2]*xB_tmp[60] - 0.55*p[2]*xB_tmp[61];
xBdot_tmp[60] = p[4]*xB_tmp[60] - 0.9*p[2]*xB_tmp[56] - 1.0*p[4]*xB_tmp[57] - 0.55*p[2]*xB_tmp[55];
xBdot_tmp[61] = xB_tmp[61]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[55] - 1.1*p[2]*xB_tmp[52] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[60];
xBdot_tmp[62] = xB_tmp[62]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[69] - 1.0*p[4]*xB_tmp[68] - 1.0*p[2]*xB_tmp[60];
xBdot_tmp[63] = p[1]*xB_tmp[46] - 1.0*p[1]*xB_tmp[44] - 1.0*p[1]*xB_tmp[45] - 1.0*p[1]*xB_tmp[43] - 1.0*p[2]*xB_tmp[58] - 1.0*p[2]*xB_tmp[59] + p[1]*xB_tmp[63];
xBdot_tmp[64] = 0.0002*xB_tmp[64] - 0.45*p[2]*xB_tmp[38] - 0.55*p[2]*xB_tmp[35] - 0.0002*xB_tmp[65];
xBdot_tmp[65] = p[0]*xB_tmp[65] - 1.0*p[0]*xB_tmp[67] - 0.55*p[2]*xB_tmp[66] - 0.45*p[2]*xB_tmp[69];
xBdot_tmp[66] = xB_tmp[66]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[47] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[69];
xBdot_tmp[67] = p[1]*xB_tmp[67] - 2.0*p[2]*xB_tmp[59] - 0.45*p[2]*xB_tmp[62] - 1.0*p[1]*xB_tmp[65] - 0.55*p[2]*xB_tmp[47];
xBdot_tmp[68] = xB_tmp[68]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[57] - 1.0*p[1]*xB_tmp[37];
xBdot_tmp[69] = xB_tmp[69]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[62] - 1.0*p[4]*xB_tmp[37];
xBdot_tmp[70] = xB_tmp[70]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[101] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[73];
xBdot_tmp[71] = (p[5] + 0.0002)*xB_tmp[71] - 0.0002*xB_tmp[72];
xBdot_tmp[72] = xB_tmp[72]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[103];
xBdot_tmp[73] = (p[4] + 0.0002)*xB_tmp[73] - 1.0*p[4]*xB_tmp[71] - 0.0002*xB_tmp[104];
xBdot_tmp[74] = 0.0004*xB_tmp[74] - 0.0002*xB_tmp[75];
xBdot_tmp[75] = (p[0] + 0.0002)*xB_tmp[75] - 1.0*p[0]*xB_tmp[77] - 0.0004*xB_tmp[78];
xBdot_tmp[76] = 0.0002*xB_tmp[75] - 0.0002*xB_tmp[74] + 0.0002*xB_tmp[76] - 0.0002*xB_tmp[78] - 0.0002*xB_tmp[79];
xBdot_tmp[77] = (p[1] + 0.0002)*xB_tmp[77] - 1.0*p[2]*xB_tmp[99] - 1.0*p[1]*xB_tmp[75] - 0.0002*xB_tmp[81];
xBdot_tmp[78] = 2.0*p[0]*xB_tmp[78] - 1.0*p[0]*xB_tmp[81];
xBdot_tmp[79] = p[0]*xB_tmp[79] - 1.0*p[0]*xB_tmp[78] - 1.0*p[0]*xB_tmp[80] + p[0]*xB_tmp[81] - 1.0*p[0]*xB_tmp[98];
xBdot_tmp[80] = 2.0*p[1]*xB_tmp[80] - 1.0*p[1]*xB_tmp[81] - 1.0*p[2]*xB_tmp[102];
xBdot_tmp[81] = xB_tmp[81]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[80] - 1.0*p[2]*xB_tmp[100] - 2.0*p[1]*xB_tmp[78];
xBdot_tmp[82] = xB_tmp[82]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[96] - 1.0*p[1]*xB_tmp[101] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[97];
xBdot_tmp[83] = p[5]*xB_tmp[83] - 1.0*p[5]*xB_tmp[84];
xBdot_tmp[84] = 2.0*p[5]*xB_tmp[84];
xBdot_tmp[85] = xB_tmp[85]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[86];
xBdot_tmp[86] = xB_tmp[86]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[84];
xBdot_tmp[87] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[87] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[90];
xBdot_tmp[88] = xB_tmp[88]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[89] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[91] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[87] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[90];
xBdot_tmp[89] = p[4]*xB_tmp[86] - 1.0*p[4]*xB_tmp[84] - 1.0*p[4]*xB_tmp[83] + p[4]*xB_tmp[89] - 1.0*p[4]*xB_tmp[91];
xBdot_tmp[90] = xB_tmp[90]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[85] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[91];
xBdot_tmp[91] = 2.0*p[4]*xB_tmp[91] - 1.0*p[4]*xB_tmp[86];
xBdot_tmp[92] = p[5]*xB_tmp[92] - 0.45*p[2]*xB_tmp[86] - 0.55*p[2]*xB_tmp[85];
xBdot_tmp[93] = - 0.55*p[2]*xB_tmp[87] - 0.55*p[2]*xB_tmp[88] - 0.45*p[2]*xB_tmp[89] - 0.45*p[2]*xB_tmp[91];
xBdot_tmp[94] = - 0.45*p[2]*xB_tmp[95] - 0.55*p[2]*xB_tmp[96];
xBdot_tmp[95] = p[4]*xB_tmp[95] - 0.9*p[2]*xB_tmp[91] - 1.0*p[4]*xB_tmp[92] - 0.55*p[2]*xB_tmp[90];
xBdot_tmp[96] = xB_tmp[96]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[90] - 1.1*p[2]*xB_tmp[87] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[95];
xBdot_tmp[97] = xB_tmp[97]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[104] - 1.0*p[4]*xB_tmp[103] - 1.0*p[2]*xB_tmp[95];
xBdot_tmp[98] = p[1]*xB_tmp[81] - 1.0*p[1]*xB_tmp[79] - 1.0*p[1]*xB_tmp[80] - 1.0*p[1]*xB_tmp[78] - 1.0*p[2]*xB_tmp[93] - 1.0*p[2]*xB_tmp[94] + p[1]*xB_tmp[98];
xBdot_tmp[99] = 0.0002*xB_tmp[99] - 0.45*p[2]*xB_tmp[73] - 0.55*p[2]*xB_tmp[70] - 0.0002*xB_tmp[100];
xBdot_tmp[100] = p[0]*xB_tmp[100] - 1.0*p[0]*xB_tmp[102] - 0.55*p[2]*xB_tmp[101] - 0.45*p[2]*xB_tmp[104];
xBdot_tmp[101] = xB_tmp[101]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[82] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[104];
xBdot_tmp[102] = p[1]*xB_tmp[102] - 2.0*p[2]*xB_tmp[94] - 0.45*p[2]*xB_tmp[97] - 1.0*p[1]*xB_tmp[100] - 0.55*p[2]*xB_tmp[82];
xBdot_tmp[103] = xB_tmp[103]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[92] - 1.0*p[1]*xB_tmp[72];
xBdot_tmp[104] = xB_tmp[104]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[97] - 1.0*p[4]*xB_tmp[72];
xBdot_tmp[105] = xB_tmp[105]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[136] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[108];
xBdot_tmp[106] = (p[5] + 0.0002)*xB_tmp[106] - 0.0002*xB_tmp[107];
xBdot_tmp[107] = xB_tmp[107]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[138];
xBdot_tmp[108] = (p[4] + 0.0002)*xB_tmp[108] - 1.0*p[4]*xB_tmp[106] - 0.0002*xB_tmp[139];
xBdot_tmp[109] = 0.0004*xB_tmp[109] - 0.0002*xB_tmp[110];
xBdot_tmp[110] = (p[0] + 0.0002)*xB_tmp[110] - 1.0*p[0]*xB_tmp[112] - 0.0004*xB_tmp[113];
xBdot_tmp[111] = 0.0002*xB_tmp[110] - 0.0002*xB_tmp[109] + 0.0002*xB_tmp[111] - 0.0002*xB_tmp[113] - 0.0002*xB_tmp[114];
xBdot_tmp[112] = (p[1] + 0.0002)*xB_tmp[112] - 1.0*p[2]*xB_tmp[134] - 1.0*p[1]*xB_tmp[110] - 0.0002*xB_tmp[116];
xBdot_tmp[113] = 2.0*p[0]*xB_tmp[113] - 1.0*p[0]*xB_tmp[116];
xBdot_tmp[114] = p[0]*xB_tmp[114] - 1.0*p[0]*xB_tmp[113] - 1.0*p[0]*xB_tmp[115] + p[0]*xB_tmp[116] - 1.0*p[0]*xB_tmp[133];
xBdot_tmp[115] = 2.0*p[1]*xB_tmp[115] - 1.0*p[1]*xB_tmp[116] - 1.0*p[2]*xB_tmp[137];
xBdot_tmp[116] = xB_tmp[116]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[115] - 1.0*p[2]*xB_tmp[135] - 2.0*p[1]*xB_tmp[113];
xBdot_tmp[117] = xB_tmp[117]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[131] - 1.0*p[1]*xB_tmp[136] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[132];
xBdot_tmp[118] = p[5]*xB_tmp[118] - 1.0*p[5]*xB_tmp[119];
xBdot_tmp[119] = 2.0*p[5]*xB_tmp[119];
xBdot_tmp[120] = xB_tmp[120]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[121];
xBdot_tmp[121] = xB_tmp[121]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[119];
xBdot_tmp[122] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[122] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[125];
xBdot_tmp[123] = xB_tmp[123]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[124] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[126] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[122] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[125];
xBdot_tmp[124] = p[4]*xB_tmp[121] - 1.0*p[4]*xB_tmp[119] - 1.0*p[4]*xB_tmp[118] + p[4]*xB_tmp[124] - 1.0*p[4]*xB_tmp[126];
xBdot_tmp[125] = xB_tmp[125]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[120] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[126];
xBdot_tmp[126] = 2.0*p[4]*xB_tmp[126] - 1.0*p[4]*xB_tmp[121];
xBdot_tmp[127] = p[5]*xB_tmp[127] - 0.45*p[2]*xB_tmp[121] - 0.55*p[2]*xB_tmp[120];
xBdot_tmp[128] = - 0.55*p[2]*xB_tmp[122] - 0.55*p[2]*xB_tmp[123] - 0.45*p[2]*xB_tmp[124] - 0.45*p[2]*xB_tmp[126];
xBdot_tmp[129] = - 0.45*p[2]*xB_tmp[130] - 0.55*p[2]*xB_tmp[131];
xBdot_tmp[130] = p[4]*xB_tmp[130] - 0.9*p[2]*xB_tmp[126] - 1.0*p[4]*xB_tmp[127] - 0.55*p[2]*xB_tmp[125];
xBdot_tmp[131] = xB_tmp[131]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[125] - 1.1*p[2]*xB_tmp[122] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[130];
xBdot_tmp[132] = xB_tmp[132]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[139] - 1.0*p[4]*xB_tmp[138] - 1.0*p[2]*xB_tmp[130];
xBdot_tmp[133] = p[1]*xB_tmp[116] - 1.0*p[1]*xB_tmp[114] - 1.0*p[1]*xB_tmp[115] - 1.0*p[1]*xB_tmp[113] - 1.0*p[2]*xB_tmp[128] - 1.0*p[2]*xB_tmp[129] + p[1]*xB_tmp[133];
xBdot_tmp[134] = 0.0002*xB_tmp[134] - 0.45*p[2]*xB_tmp[108] - 0.55*p[2]*xB_tmp[105] - 0.0002*xB_tmp[135];
xBdot_tmp[135] = p[0]*xB_tmp[135] - 1.0*p[0]*xB_tmp[137] - 0.55*p[2]*xB_tmp[136] - 0.45*p[2]*xB_tmp[139];
xBdot_tmp[136] = xB_tmp[136]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[117] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[139];
xBdot_tmp[137] = p[1]*xB_tmp[137] - 2.0*p[2]*xB_tmp[129] - 0.45*p[2]*xB_tmp[132] - 1.0*p[1]*xB_tmp[135] - 0.55*p[2]*xB_tmp[117];
xBdot_tmp[138] = xB_tmp[138]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[127] - 1.0*p[1]*xB_tmp[107];
xBdot_tmp[139] = xB_tmp[139]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[132] - 1.0*p[4]*xB_tmp[107];
xBdot_tmp[140] = xB_tmp[140]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[171] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[143];
xBdot_tmp[141] = (p[5] + 0.0002)*xB_tmp[141] - 0.0002*xB_tmp[142];
xBdot_tmp[142] = xB_tmp[142]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[173];
xBdot_tmp[143] = (p[4] + 0.0002)*xB_tmp[143] - 1.0*p[4]*xB_tmp[141] - 0.0002*xB_tmp[174];
xBdot_tmp[144] = 0.0004*xB_tmp[144] - 0.0002*xB_tmp[145];
xBdot_tmp[145] = (p[0] + 0.0002)*xB_tmp[145] - 1.0*p[0]*xB_tmp[147] - 0.0004*xB_tmp[148];
xBdot_tmp[146] = 0.0002*xB_tmp[145] - 0.0002*xB_tmp[144] + 0.0002*xB_tmp[146] - 0.0002*xB_tmp[148] - 0.0002*xB_tmp[149];
xBdot_tmp[147] = (p[1] + 0.0002)*xB_tmp[147] - 1.0*p[2]*xB_tmp[169] - 1.0*p[1]*xB_tmp[145] - 0.0002*xB_tmp[151];
xBdot_tmp[148] = 2.0*p[0]*xB_tmp[148] - 1.0*p[0]*xB_tmp[151];
xBdot_tmp[149] = p[0]*xB_tmp[149] - 1.0*p[0]*xB_tmp[148] - 1.0*p[0]*xB_tmp[150] + p[0]*xB_tmp[151] - 1.0*p[0]*xB_tmp[168];
xBdot_tmp[150] = 2.0*p[1]*xB_tmp[150] - 1.0*p[1]*xB_tmp[151] - 1.0*p[2]*xB_tmp[172];
xBdot_tmp[151] = xB_tmp[151]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[150] - 1.0*p[2]*xB_tmp[170] - 2.0*p[1]*xB_tmp[148];
xBdot_tmp[152] = xB_tmp[152]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[166] - 1.0*p[1]*xB_tmp[171] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[167];
xBdot_tmp[153] = p[5]*xB_tmp[153] - 1.0*p[5]*xB_tmp[154];
xBdot_tmp[154] = 2.0*p[5]*xB_tmp[154];
xBdot_tmp[155] = xB_tmp[155]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[156];
xBdot_tmp[156] = xB_tmp[156]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[154];
xBdot_tmp[157] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[157] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[160];
xBdot_tmp[158] = xB_tmp[158]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[159] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[161] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[157] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[160];
xBdot_tmp[159] = p[4]*xB_tmp[156] - 1.0*p[4]*xB_tmp[154] - 1.0*p[4]*xB_tmp[153] + p[4]*xB_tmp[159] - 1.0*p[4]*xB_tmp[161];
xBdot_tmp[160] = xB_tmp[160]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[155] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[161];
xBdot_tmp[161] = 2.0*p[4]*xB_tmp[161] - 1.0*p[4]*xB_tmp[156];
xBdot_tmp[162] = p[5]*xB_tmp[162] - 0.45*p[2]*xB_tmp[156] - 0.55*p[2]*xB_tmp[155];
xBdot_tmp[163] = - 0.55*p[2]*xB_tmp[157] - 0.55*p[2]*xB_tmp[158] - 0.45*p[2]*xB_tmp[159] - 0.45*p[2]*xB_tmp[161];
xBdot_tmp[164] = - 0.45*p[2]*xB_tmp[165] - 0.55*p[2]*xB_tmp[166];
xBdot_tmp[165] = p[4]*xB_tmp[165] - 0.9*p[2]*xB_tmp[161] - 1.0*p[4]*xB_tmp[162] - 0.55*p[2]*xB_tmp[160];
xBdot_tmp[166] = xB_tmp[166]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[160] - 1.1*p[2]*xB_tmp[157] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[165];
xBdot_tmp[167] = xB_tmp[167]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[174] - 1.0*p[4]*xB_tmp[173] - 1.0*p[2]*xB_tmp[165];
xBdot_tmp[168] = p[1]*xB_tmp[151] - 1.0*p[1]*xB_tmp[149] - 1.0*p[1]*xB_tmp[150] - 1.0*p[1]*xB_tmp[148] - 1.0*p[2]*xB_tmp[163] - 1.0*p[2]*xB_tmp[164] + p[1]*xB_tmp[168];
xBdot_tmp[169] = 0.0002*xB_tmp[169] - 0.45*p[2]*xB_tmp[143] - 0.55*p[2]*xB_tmp[140] - 0.0002*xB_tmp[170];
xBdot_tmp[170] = p[0]*xB_tmp[170] - 1.0*p[0]*xB_tmp[172] - 0.55*p[2]*xB_tmp[171] - 0.45*p[2]*xB_tmp[174];
xBdot_tmp[171] = xB_tmp[171]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[152] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[174];
xBdot_tmp[172] = p[1]*xB_tmp[172] - 2.0*p[2]*xB_tmp[164] - 0.45*p[2]*xB_tmp[167] - 1.0*p[1]*xB_tmp[170] - 0.55*p[2]*xB_tmp[152];
xBdot_tmp[173] = xB_tmp[173]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[162] - 1.0*p[1]*xB_tmp[142];
xBdot_tmp[174] = xB_tmp[174]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[167] - 1.0*p[4]*xB_tmp[142];
xBdot_tmp[175] = xB_tmp[175]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[206] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[178];
xBdot_tmp[176] = (p[5] + 0.0002)*xB_tmp[176] - 0.0002*xB_tmp[177];
xBdot_tmp[177] = xB_tmp[177]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[208];
xBdot_tmp[178] = (p[4] + 0.0002)*xB_tmp[178] - 1.0*p[4]*xB_tmp[176] - 0.0002*xB_tmp[209];
xBdot_tmp[179] = 0.0004*xB_tmp[179] - 0.0002*xB_tmp[180];
xBdot_tmp[180] = (p[0] + 0.0002)*xB_tmp[180] - 1.0*p[0]*xB_tmp[182] - 0.0004*xB_tmp[183];
xBdot_tmp[181] = 0.0002*xB_tmp[180] - 0.0002*xB_tmp[179] + 0.0002*xB_tmp[181] - 0.0002*xB_tmp[183] - 0.0002*xB_tmp[184];
xBdot_tmp[182] = (p[1] + 0.0002)*xB_tmp[182] - 1.0*p[2]*xB_tmp[204] - 1.0*p[1]*xB_tmp[180] - 0.0002*xB_tmp[186];
xBdot_tmp[183] = 2.0*p[0]*xB_tmp[183] - 1.0*p[0]*xB_tmp[186];
xBdot_tmp[184] = p[0]*xB_tmp[184] - 1.0*p[0]*xB_tmp[183] - 1.0*p[0]*xB_tmp[185] + p[0]*xB_tmp[186] - 1.0*p[0]*xB_tmp[203];
xBdot_tmp[185] = 2.0*p[1]*xB_tmp[185] - 1.0*p[1]*xB_tmp[186] - 1.0*p[2]*xB_tmp[207];
xBdot_tmp[186] = xB_tmp[186]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[185] - 1.0*p[2]*xB_tmp[205] - 2.0*p[1]*xB_tmp[183];
xBdot_tmp[187] = xB_tmp[187]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[201] - 1.0*p[1]*xB_tmp[206] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[202];
xBdot_tmp[188] = p[5]*xB_tmp[188] - 1.0*p[5]*xB_tmp[189];
xBdot_tmp[189] = 2.0*p[5]*xB_tmp[189];
xBdot_tmp[190] = xB_tmp[190]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[191];
xBdot_tmp[191] = xB_tmp[191]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[189];
xBdot_tmp[192] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[192] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[195];
xBdot_tmp[193] = xB_tmp[193]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[194] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[196] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[192] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[195];
xBdot_tmp[194] = p[4]*xB_tmp[191] - 1.0*p[4]*xB_tmp[189] - 1.0*p[4]*xB_tmp[188] + p[4]*xB_tmp[194] - 1.0*p[4]*xB_tmp[196];
xBdot_tmp[195] = xB_tmp[195]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[190] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[196];
xBdot_tmp[196] = 2.0*p[4]*xB_tmp[196] - 1.0*p[4]*xB_tmp[191];
xBdot_tmp[197] = p[5]*xB_tmp[197] - 0.45*p[2]*xB_tmp[191] - 0.55*p[2]*xB_tmp[190];
xBdot_tmp[198] = - 0.55*p[2]*xB_tmp[192] - 0.55*p[2]*xB_tmp[193] - 0.45*p[2]*xB_tmp[194] - 0.45*p[2]*xB_tmp[196];
xBdot_tmp[199] = - 0.45*p[2]*xB_tmp[200] - 0.55*p[2]*xB_tmp[201];
xBdot_tmp[200] = p[4]*xB_tmp[200] - 0.9*p[2]*xB_tmp[196] - 1.0*p[4]*xB_tmp[197] - 0.55*p[2]*xB_tmp[195];
xBdot_tmp[201] = xB_tmp[201]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[195] - 1.1*p[2]*xB_tmp[192] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[200];
xBdot_tmp[202] = xB_tmp[202]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[209] - 1.0*p[4]*xB_tmp[208] - 1.0*p[2]*xB_tmp[200];
xBdot_tmp[203] = p[1]*xB_tmp[186] - 1.0*p[1]*xB_tmp[184] - 1.0*p[1]*xB_tmp[185] - 1.0*p[1]*xB_tmp[183] - 1.0*p[2]*xB_tmp[198] - 1.0*p[2]*xB_tmp[199] + p[1]*xB_tmp[203];
xBdot_tmp[204] = 0.0002*xB_tmp[204] - 0.45*p[2]*xB_tmp[178] - 0.55*p[2]*xB_tmp[175] - 0.0002*xB_tmp[205];
xBdot_tmp[205] = p[0]*xB_tmp[205] - 1.0*p[0]*xB_tmp[207] - 0.55*p[2]*xB_tmp[206] - 0.45*p[2]*xB_tmp[209];
xBdot_tmp[206] = xB_tmp[206]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[187] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[209];
xBdot_tmp[207] = p[1]*xB_tmp[207] - 2.0*p[2]*xB_tmp[199] - 0.45*p[2]*xB_tmp[202] - 1.0*p[1]*xB_tmp[205] - 0.55*p[2]*xB_tmp[187];
xBdot_tmp[208] = xB_tmp[208]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[197] - 1.0*p[1]*xB_tmp[177];
xBdot_tmp[209] = xB_tmp[209]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[202] - 1.0*p[4]*xB_tmp[177];
xBdot_tmp[210] = xB_tmp[210]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[241] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[213];
xBdot_tmp[211] = (p[5] + 0.0002)*xB_tmp[211] - 0.0002*xB_tmp[212];
xBdot_tmp[212] = xB_tmp[212]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[243];
xBdot_tmp[213] = (p[4] + 0.0002)*xB_tmp[213] - 1.0*p[4]*xB_tmp[211] - 0.0002*xB_tmp[244];
xBdot_tmp[214] = 0.0004*xB_tmp[214] - 0.0002*xB_tmp[215];
xBdot_tmp[215] = (p[0] + 0.0002)*xB_tmp[215] - 1.0*p[0]*xB_tmp[217] - 0.0004*xB_tmp[218];
xBdot_tmp[216] = 0.0002*xB_tmp[215] - 0.0002*xB_tmp[214] + 0.0002*xB_tmp[216] - 0.0002*xB_tmp[218] - 0.0002*xB_tmp[219];
xBdot_tmp[217] = (p[1] + 0.0002)*xB_tmp[217] - 1.0*p[2]*xB_tmp[239] - 1.0*p[1]*xB_tmp[215] - 0.0002*xB_tmp[221];
xBdot_tmp[218] = 2.0*p[0]*xB_tmp[218] - 1.0*p[0]*xB_tmp[221];
xBdot_tmp[219] = p[0]*xB_tmp[219] - 1.0*p[0]*xB_tmp[218] - 1.0*p[0]*xB_tmp[220] + p[0]*xB_tmp[221] - 1.0*p[0]*xB_tmp[238];
xBdot_tmp[220] = 2.0*p[1]*xB_tmp[220] - 1.0*p[1]*xB_tmp[221] - 1.0*p[2]*xB_tmp[242];
xBdot_tmp[221] = xB_tmp[221]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[220] - 1.0*p[2]*xB_tmp[240] - 2.0*p[1]*xB_tmp[218];
xBdot_tmp[222] = xB_tmp[222]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[236] - 1.0*p[1]*xB_tmp[241] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[237];
xBdot_tmp[223] = p[5]*xB_tmp[223] - 1.0*p[5]*xB_tmp[224];
xBdot_tmp[224] = 2.0*p[5]*xB_tmp[224];
xBdot_tmp[225] = xB_tmp[225]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[226];
xBdot_tmp[226] = xB_tmp[226]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[224];
xBdot_tmp[227] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[227] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[230];
xBdot_tmp[228] = xB_tmp[228]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[229] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[231] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[227] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[230];
xBdot_tmp[229] = p[4]*xB_tmp[226] - 1.0*p[4]*xB_tmp[224] - 1.0*p[4]*xB_tmp[223] + p[4]*xB_tmp[229] - 1.0*p[4]*xB_tmp[231];
xBdot_tmp[230] = xB_tmp[230]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[225] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[231];
xBdot_tmp[231] = 2.0*p[4]*xB_tmp[231] - 1.0*p[4]*xB_tmp[226];
xBdot_tmp[232] = p[5]*xB_tmp[232] - 0.45*p[2]*xB_tmp[226] - 0.55*p[2]*xB_tmp[225];
xBdot_tmp[233] = - 0.55*p[2]*xB_tmp[227] - 0.55*p[2]*xB_tmp[228] - 0.45*p[2]*xB_tmp[229] - 0.45*p[2]*xB_tmp[231];
xBdot_tmp[234] = - 0.45*p[2]*xB_tmp[235] - 0.55*p[2]*xB_tmp[236];
xBdot_tmp[235] = p[4]*xB_tmp[235] - 0.9*p[2]*xB_tmp[231] - 1.0*p[4]*xB_tmp[232] - 0.55*p[2]*xB_tmp[230];
xBdot_tmp[236] = xB_tmp[236]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[230] - 1.1*p[2]*xB_tmp[227] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[235];
xBdot_tmp[237] = xB_tmp[237]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[244] - 1.0*p[4]*xB_tmp[243] - 1.0*p[2]*xB_tmp[235];
xBdot_tmp[238] = p[1]*xB_tmp[221] - 1.0*p[1]*xB_tmp[219] - 1.0*p[1]*xB_tmp[220] - 1.0*p[1]*xB_tmp[218] - 1.0*p[2]*xB_tmp[233] - 1.0*p[2]*xB_tmp[234] + p[1]*xB_tmp[238];
xBdot_tmp[239] = 0.0002*xB_tmp[239] - 0.45*p[2]*xB_tmp[213] - 0.55*p[2]*xB_tmp[210] - 0.0002*xB_tmp[240];
xBdot_tmp[240] = p[0]*xB_tmp[240] - 1.0*p[0]*xB_tmp[242] - 0.55*p[2]*xB_tmp[241] - 0.45*p[2]*xB_tmp[244];
xBdot_tmp[241] = xB_tmp[241]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[222] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[244];
xBdot_tmp[242] = p[1]*xB_tmp[242] - 2.0*p[2]*xB_tmp[234] - 0.45*p[2]*xB_tmp[237] - 1.0*p[1]*xB_tmp[240] - 0.55*p[2]*xB_tmp[222];
xBdot_tmp[243] = xB_tmp[243]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[232] - 1.0*p[1]*xB_tmp[212];
xBdot_tmp[244] = xB_tmp[244]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[237] - 1.0*p[4]*xB_tmp[212];
xBdot_tmp[245] = xB_tmp[245]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[276] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[248];
xBdot_tmp[246] = (p[5] + 0.0002)*xB_tmp[246] - 0.0002*xB_tmp[247];
xBdot_tmp[247] = xB_tmp[247]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[278];
xBdot_tmp[248] = (p[4] + 0.0002)*xB_tmp[248] - 1.0*p[4]*xB_tmp[246] - 0.0002*xB_tmp[279];
xBdot_tmp[249] = 0.0004*xB_tmp[249] - 0.0002*xB_tmp[250];
xBdot_tmp[250] = (p[0] + 0.0002)*xB_tmp[250] - 1.0*p[0]*xB_tmp[252] - 0.0004*xB_tmp[253];
xBdot_tmp[251] = 0.0002*xB_tmp[250] - 0.0002*xB_tmp[249] + 0.0002*xB_tmp[251] - 0.0002*xB_tmp[253] - 0.0002*xB_tmp[254];
xBdot_tmp[252] = (p[1] + 0.0002)*xB_tmp[252] - 1.0*p[2]*xB_tmp[274] - 1.0*p[1]*xB_tmp[250] - 0.0002*xB_tmp[256];
xBdot_tmp[253] = 2.0*p[0]*xB_tmp[253] - 1.0*p[0]*xB_tmp[256];
xBdot_tmp[254] = p[0]*xB_tmp[254] - 1.0*p[0]*xB_tmp[253] - 1.0*p[0]*xB_tmp[255] + p[0]*xB_tmp[256] - 1.0*p[0]*xB_tmp[273];
xBdot_tmp[255] = 2.0*p[1]*xB_tmp[255] - 1.0*p[1]*xB_tmp[256] - 1.0*p[2]*xB_tmp[277];
xBdot_tmp[256] = xB_tmp[256]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[255] - 1.0*p[2]*xB_tmp[275] - 2.0*p[1]*xB_tmp[253];
xBdot_tmp[257] = xB_tmp[257]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[271] - 1.0*p[1]*xB_tmp[276] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[272];
xBdot_tmp[258] = p[5]*xB_tmp[258] - 1.0*p[5]*xB_tmp[259];
xBdot_tmp[259] = 2.0*p[5]*xB_tmp[259];
xBdot_tmp[260] = xB_tmp[260]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[261];
xBdot_tmp[261] = xB_tmp[261]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[259];
xBdot_tmp[262] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[262] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[265];
xBdot_tmp[263] = xB_tmp[263]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[264] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[266] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[262] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[265];
xBdot_tmp[264] = p[4]*xB_tmp[261] - 1.0*p[4]*xB_tmp[259] - 1.0*p[4]*xB_tmp[258] + p[4]*xB_tmp[264] - 1.0*p[4]*xB_tmp[266];
xBdot_tmp[265] = xB_tmp[265]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[260] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[266];
xBdot_tmp[266] = 2.0*p[4]*xB_tmp[266] - 1.0*p[4]*xB_tmp[261];
xBdot_tmp[267] = p[5]*xB_tmp[267] - 0.45*p[2]*xB_tmp[261] - 0.55*p[2]*xB_tmp[260];
xBdot_tmp[268] = - 0.55*p[2]*xB_tmp[262] - 0.55*p[2]*xB_tmp[263] - 0.45*p[2]*xB_tmp[264] - 0.45*p[2]*xB_tmp[266];
xBdot_tmp[269] = - 0.45*p[2]*xB_tmp[270] - 0.55*p[2]*xB_tmp[271];
xBdot_tmp[270] = p[4]*xB_tmp[270] - 0.9*p[2]*xB_tmp[266] - 1.0*p[4]*xB_tmp[267] - 0.55*p[2]*xB_tmp[265];
xBdot_tmp[271] = xB_tmp[271]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[265] - 1.1*p[2]*xB_tmp[262] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[270];
xBdot_tmp[272] = xB_tmp[272]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[279] - 1.0*p[4]*xB_tmp[278] - 1.0*p[2]*xB_tmp[270];
xBdot_tmp[273] = p[1]*xB_tmp[256] - 1.0*p[1]*xB_tmp[254] - 1.0*p[1]*xB_tmp[255] - 1.0*p[1]*xB_tmp[253] - 1.0*p[2]*xB_tmp[268] - 1.0*p[2]*xB_tmp[269] + p[1]*xB_tmp[273];
xBdot_tmp[274] = 0.0002*xB_tmp[274] - 0.45*p[2]*xB_tmp[248] - 0.55*p[2]*xB_tmp[245] - 0.0002*xB_tmp[275];
xBdot_tmp[275] = p[0]*xB_tmp[275] - 1.0*p[0]*xB_tmp[277] - 0.55*p[2]*xB_tmp[276] - 0.45*p[2]*xB_tmp[279];
xBdot_tmp[276] = xB_tmp[276]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[257] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[279];
xBdot_tmp[277] = p[1]*xB_tmp[277] - 2.0*p[2]*xB_tmp[269] - 0.45*p[2]*xB_tmp[272] - 1.0*p[1]*xB_tmp[275] - 0.55*p[2]*xB_tmp[257];
xBdot_tmp[278] = xB_tmp[278]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[267] - 1.0*p[1]*xB_tmp[247];
xBdot_tmp[279] = xB_tmp[279]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[272] - 1.0*p[4]*xB_tmp[247];
xBdot_tmp[280] = xB_tmp[280]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*xB_tmp[311] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[283];
xBdot_tmp[281] = (p[5] + 0.0002)*xB_tmp[281] - 0.0002*xB_tmp[282];
xBdot_tmp[282] = xB_tmp[282]*(p[0] + p[5]) - 1.0*p[0]*xB_tmp[313];
xBdot_tmp[283] = (p[4] + 0.0002)*xB_tmp[283] - 1.0*p[4]*xB_tmp[281] - 0.0002*xB_tmp[314];
xBdot_tmp[284] = 0.0004*xB_tmp[284] - 0.0002*xB_tmp[285];
xBdot_tmp[285] = (p[0] + 0.0002)*xB_tmp[285] - 1.0*p[0]*xB_tmp[287] - 0.0004*xB_tmp[288];
xBdot_tmp[286] = 0.0002*xB_tmp[285] - 0.0002*xB_tmp[284] + 0.0002*xB_tmp[286] - 0.0002*xB_tmp[288] - 0.0002*xB_tmp[289];
xBdot_tmp[287] = (p[1] + 0.0002)*xB_tmp[287] - 1.0*p[2]*xB_tmp[309] - 1.0*p[1]*xB_tmp[285] - 0.0002*xB_tmp[291];
xBdot_tmp[288] = 2.0*p[0]*xB_tmp[288] - 1.0*p[0]*xB_tmp[291];
xBdot_tmp[289] = p[0]*xB_tmp[289] - 1.0*p[0]*xB_tmp[288] - 1.0*p[0]*xB_tmp[290] + p[0]*xB_tmp[291] - 1.0*p[0]*xB_tmp[308];
xBdot_tmp[290] = 2.0*p[1]*xB_tmp[290] - 1.0*p[1]*xB_tmp[291] - 1.0*p[2]*xB_tmp[312];
xBdot_tmp[291] = xB_tmp[291]*(p[0] + p[1]) - 2.0*p[0]*xB_tmp[290] - 1.0*p[2]*xB_tmp[310] - 2.0*p[1]*xB_tmp[288];
xBdot_tmp[292] = xB_tmp[292]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[2]*xB_tmp[306] - 1.0*p[1]*xB_tmp[311] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[307];
xBdot_tmp[293] = p[5]*xB_tmp[293] - 1.0*p[5]*xB_tmp[294];
xBdot_tmp[294] = 2.0*p[5]*xB_tmp[294];
xBdot_tmp[295] = xB_tmp[295]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[296];
xBdot_tmp[296] = xB_tmp[296]*(p[4] + p[5]) - 2.0*p[4]*xB_tmp[294];
xBdot_tmp[297] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*xB_tmp[297] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[300];
xBdot_tmp[298] = xB_tmp[298]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[299] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[301] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*xB_tmp[297] + 2.0*p[2]*(pow(p[3],2))*xB_tmp[300];
xBdot_tmp[299] = p[4]*xB_tmp[296] - 1.0*p[4]*xB_tmp[294] - 1.0*p[4]*xB_tmp[293] + p[4]*xB_tmp[299] - 1.0*p[4]*xB_tmp[301];
xBdot_tmp[300] = xB_tmp[300]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[4]*xB_tmp[295] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[301];
xBdot_tmp[301] = 2.0*p[4]*xB_tmp[301] - 1.0*p[4]*xB_tmp[296];
xBdot_tmp[302] = p[5]*xB_tmp[302] - 0.45*p[2]*xB_tmp[296] - 0.55*p[2]*xB_tmp[295];
xBdot_tmp[303] = - 0.55*p[2]*xB_tmp[297] - 0.55*p[2]*xB_tmp[298] - 0.45*p[2]*xB_tmp[299] - 0.45*p[2]*xB_tmp[301];
xBdot_tmp[304] = - 0.45*p[2]*xB_tmp[305] - 0.55*p[2]*xB_tmp[306];
xBdot_tmp[305] = p[4]*xB_tmp[305] - 0.9*p[2]*xB_tmp[301] - 1.0*p[4]*xB_tmp[302] - 0.55*p[2]*xB_tmp[300];
xBdot_tmp[306] = xB_tmp[306]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*xB_tmp[300] - 1.1*p[2]*xB_tmp[297] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[305];
xBdot_tmp[307] = xB_tmp[307]*(p[1] + p[4]) - 1.0*p[1]*xB_tmp[314] - 1.0*p[4]*xB_tmp[313] - 1.0*p[2]*xB_tmp[305];
xBdot_tmp[308] = p[1]*xB_tmp[291] - 1.0*p[1]*xB_tmp[289] - 1.0*p[1]*xB_tmp[290] - 1.0*p[1]*xB_tmp[288] - 1.0*p[2]*xB_tmp[303] - 1.0*p[2]*xB_tmp[304] + p[1]*xB_tmp[308];
xBdot_tmp[309] = 0.0002*xB_tmp[309] - 0.45*p[2]*xB_tmp[283] - 0.55*p[2]*xB_tmp[280] - 0.0002*xB_tmp[310];
xBdot_tmp[310] = p[0]*xB_tmp[310] - 1.0*p[0]*xB_tmp[312] - 0.55*p[2]*xB_tmp[311] - 0.45*p[2]*xB_tmp[314];
xBdot_tmp[311] = xB_tmp[311]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[0]*xB_tmp[292] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*xB_tmp[314];
xBdot_tmp[312] = p[1]*xB_tmp[312] - 2.0*p[2]*xB_tmp[304] - 0.45*p[2]*xB_tmp[307] - 1.0*p[1]*xB_tmp[310] - 0.55*p[2]*xB_tmp[292];
xBdot_tmp[313] = xB_tmp[313]*(p[1] + p[5]) - 1.0*p[2]*xB_tmp[302] - 1.0*p[1]*xB_tmp[282];
xBdot_tmp[314] = xB_tmp[314]*(p[0] + p[4]) - 1.0*p[0]*xB_tmp[307] - 1.0*p[4]*xB_tmp[282];

  for (ixB=0; ixB<315; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
qBdot_tmp[0+ip*ny] = x_tmp[2]*xB_tmp[2] + x_tmp[5]*xB_tmp[5] - 1.0*x_tmp[5]*xB_tmp[7] + x_tmp[9]*xB_tmp[9] - 1.0*x_tmp[2]*xB_tmp[33] - 1.0*x_tmp[9]*xB_tmp[28] - 1.0*x_tmp[31]*xB_tmp[12] + x_tmp[30]*xB_tmp[30] - 1.0*x_tmp[34]*xB_tmp[27] - 1.0*x_tmp[30]*xB_tmp[32] + x_tmp[31]*xB_tmp[31] + x_tmp[34]*xB_tmp[34] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[8] + xB_tmp[11]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[10];
qBdot_tmp[1+ip*ny] = x_tmp[2]*xB_tmp[37] + x_tmp[5]*xB_tmp[40] - 1.0*x_tmp[5]*xB_tmp[42] + x_tmp[9]*xB_tmp[44] - 1.0*x_tmp[2]*xB_tmp[68] - 1.0*x_tmp[9]*xB_tmp[63] - 1.0*x_tmp[31]*xB_tmp[47] + x_tmp[30]*xB_tmp[65] - 1.0*x_tmp[34]*xB_tmp[62] - 1.0*x_tmp[30]*xB_tmp[67] + x_tmp[31]*xB_tmp[66] + x_tmp[34]*xB_tmp[69] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[43] + xB_tmp[46]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[45];
qBdot_tmp[2+ip*ny] = x_tmp[2]*xB_tmp[72] + x_tmp[5]*xB_tmp[75] - 1.0*x_tmp[5]*xB_tmp[77] + x_tmp[9]*xB_tmp[79] - 1.0*x_tmp[2]*xB_tmp[103] - 1.0*x_tmp[9]*xB_tmp[98] - 1.0*x_tmp[31]*xB_tmp[82] + x_tmp[30]*xB_tmp[100] - 1.0*x_tmp[34]*xB_tmp[97] - 1.0*x_tmp[30]*xB_tmp[102] + x_tmp[31]*xB_tmp[101] + x_tmp[34]*xB_tmp[104] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[78] + xB_tmp[81]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[80];
qBdot_tmp[3+ip*ny] = x_tmp[2]*xB_tmp[107] + x_tmp[5]*xB_tmp[110] - 1.0*x_tmp[5]*xB_tmp[112] + x_tmp[9]*xB_tmp[114] - 1.0*x_tmp[2]*xB_tmp[138] - 1.0*x_tmp[9]*xB_tmp[133] - 1.0*x_tmp[31]*xB_tmp[117] + x_tmp[30]*xB_tmp[135] - 1.0*x_tmp[34]*xB_tmp[132] - 1.0*x_tmp[30]*xB_tmp[137] + x_tmp[31]*xB_tmp[136] + x_tmp[34]*xB_tmp[139] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[113] + xB_tmp[116]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[115];
qBdot_tmp[4+ip*ny] = x_tmp[2]*xB_tmp[142] + x_tmp[5]*xB_tmp[145] - 1.0*x_tmp[5]*xB_tmp[147] + x_tmp[9]*xB_tmp[149] - 1.0*x_tmp[2]*xB_tmp[173] - 1.0*x_tmp[9]*xB_tmp[168] - 1.0*x_tmp[31]*xB_tmp[152] + x_tmp[30]*xB_tmp[170] - 1.0*x_tmp[34]*xB_tmp[167] - 1.0*x_tmp[30]*xB_tmp[172] + x_tmp[31]*xB_tmp[171] + x_tmp[34]*xB_tmp[174] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[148] + xB_tmp[151]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[150];
qBdot_tmp[5+ip*ny] = x_tmp[2]*xB_tmp[177] + x_tmp[5]*xB_tmp[180] - 1.0*x_tmp[5]*xB_tmp[182] + x_tmp[9]*xB_tmp[184] - 1.0*x_tmp[2]*xB_tmp[208] - 1.0*x_tmp[9]*xB_tmp[203] - 1.0*x_tmp[31]*xB_tmp[187] + x_tmp[30]*xB_tmp[205] - 1.0*x_tmp[34]*xB_tmp[202] - 1.0*x_tmp[30]*xB_tmp[207] + x_tmp[31]*xB_tmp[206] + x_tmp[34]*xB_tmp[209] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[183] + xB_tmp[186]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[185];
qBdot_tmp[6+ip*ny] = x_tmp[2]*xB_tmp[212] + x_tmp[5]*xB_tmp[215] - 1.0*x_tmp[5]*xB_tmp[217] + x_tmp[9]*xB_tmp[219] - 1.0*x_tmp[2]*xB_tmp[243] - 1.0*x_tmp[9]*xB_tmp[238] - 1.0*x_tmp[31]*xB_tmp[222] + x_tmp[30]*xB_tmp[240] - 1.0*x_tmp[34]*xB_tmp[237] - 1.0*x_tmp[30]*xB_tmp[242] + x_tmp[31]*xB_tmp[241] + x_tmp[34]*xB_tmp[244] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[218] + xB_tmp[221]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[220];
qBdot_tmp[7+ip*ny] = x_tmp[2]*xB_tmp[247] + x_tmp[5]*xB_tmp[250] - 1.0*x_tmp[5]*xB_tmp[252] + x_tmp[9]*xB_tmp[254] - 1.0*x_tmp[2]*xB_tmp[278] - 1.0*x_tmp[9]*xB_tmp[273] - 1.0*x_tmp[31]*xB_tmp[257] + x_tmp[30]*xB_tmp[275] - 1.0*x_tmp[34]*xB_tmp[272] - 1.0*x_tmp[30]*xB_tmp[277] + x_tmp[31]*xB_tmp[276] + x_tmp[34]*xB_tmp[279] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[253] + xB_tmp[256]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[255];
qBdot_tmp[8+ip*ny] = x_tmp[2]*xB_tmp[282] + x_tmp[5]*xB_tmp[285] - 1.0*x_tmp[5]*xB_tmp[287] + x_tmp[9]*xB_tmp[289] - 1.0*x_tmp[2]*xB_tmp[313] - 1.0*x_tmp[9]*xB_tmp[308] - 1.0*x_tmp[31]*xB_tmp[292] + x_tmp[30]*xB_tmp[310] - 1.0*x_tmp[34]*xB_tmp[307] - 1.0*x_tmp[30]*xB_tmp[312] + x_tmp[31]*xB_tmp[311] + x_tmp[34]*xB_tmp[314] + (2.0*x_tmp[8] - 1.0*x_tmp[9])*xB_tmp[288] + xB_tmp[291]*(x_tmp[9] - 1.0*x_tmp[8] + x_tmp[11]) - 1.0*(x_tmp[9] + 2.0*x_tmp[11])*xB_tmp[290];

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = x_tmp[7]*xB_tmp[7] - 1.0*x_tmp[7]*xB_tmp[5] + x_tmp[12]*xB_tmp[12] - 1.0*x_tmp[33]*xB_tmp[2] - 1.0*x_tmp[28]*xB_tmp[9] - 1.0*x_tmp[12]*xB_tmp[31] + x_tmp[27]*xB_tmp[27] + x_tmp[28]*xB_tmp[28] - 1.0*x_tmp[27]*xB_tmp[34] - 1.0*x_tmp[32]*xB_tmp[30] + x_tmp[32]*xB_tmp[32] + x_tmp[33]*xB_tmp[33] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[10] + xB_tmp[11]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[8];
qBdot_tmp[1+ip*ny] = x_tmp[7]*xB_tmp[42] - 1.0*x_tmp[7]*xB_tmp[40] + x_tmp[12]*xB_tmp[47] - 1.0*x_tmp[33]*xB_tmp[37] - 1.0*x_tmp[28]*xB_tmp[44] - 1.0*x_tmp[12]*xB_tmp[66] + x_tmp[27]*xB_tmp[62] + x_tmp[28]*xB_tmp[63] - 1.0*x_tmp[27]*xB_tmp[69] - 1.0*x_tmp[32]*xB_tmp[65] + x_tmp[32]*xB_tmp[67] + x_tmp[33]*xB_tmp[68] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[45] + xB_tmp[46]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[43];
qBdot_tmp[2+ip*ny] = x_tmp[7]*xB_tmp[77] - 1.0*x_tmp[7]*xB_tmp[75] + x_tmp[12]*xB_tmp[82] - 1.0*x_tmp[33]*xB_tmp[72] - 1.0*x_tmp[28]*xB_tmp[79] - 1.0*x_tmp[12]*xB_tmp[101] + x_tmp[27]*xB_tmp[97] + x_tmp[28]*xB_tmp[98] - 1.0*x_tmp[27]*xB_tmp[104] - 1.0*x_tmp[32]*xB_tmp[100] + x_tmp[32]*xB_tmp[102] + x_tmp[33]*xB_tmp[103] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[80] + xB_tmp[81]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[78];
qBdot_tmp[3+ip*ny] = x_tmp[7]*xB_tmp[112] - 1.0*x_tmp[7]*xB_tmp[110] + x_tmp[12]*xB_tmp[117] - 1.0*x_tmp[33]*xB_tmp[107] - 1.0*x_tmp[28]*xB_tmp[114] - 1.0*x_tmp[12]*xB_tmp[136] + x_tmp[27]*xB_tmp[132] + x_tmp[28]*xB_tmp[133] - 1.0*x_tmp[27]*xB_tmp[139] - 1.0*x_tmp[32]*xB_tmp[135] + x_tmp[32]*xB_tmp[137] + x_tmp[33]*xB_tmp[138] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[115] + xB_tmp[116]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[113];
qBdot_tmp[4+ip*ny] = x_tmp[7]*xB_tmp[147] - 1.0*x_tmp[7]*xB_tmp[145] + x_tmp[12]*xB_tmp[152] - 1.0*x_tmp[33]*xB_tmp[142] - 1.0*x_tmp[28]*xB_tmp[149] - 1.0*x_tmp[12]*xB_tmp[171] + x_tmp[27]*xB_tmp[167] + x_tmp[28]*xB_tmp[168] - 1.0*x_tmp[27]*xB_tmp[174] - 1.0*x_tmp[32]*xB_tmp[170] + x_tmp[32]*xB_tmp[172] + x_tmp[33]*xB_tmp[173] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[150] + xB_tmp[151]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[148];
qBdot_tmp[5+ip*ny] = x_tmp[7]*xB_tmp[182] - 1.0*x_tmp[7]*xB_tmp[180] + x_tmp[12]*xB_tmp[187] - 1.0*x_tmp[33]*xB_tmp[177] - 1.0*x_tmp[28]*xB_tmp[184] - 1.0*x_tmp[12]*xB_tmp[206] + x_tmp[27]*xB_tmp[202] + x_tmp[28]*xB_tmp[203] - 1.0*x_tmp[27]*xB_tmp[209] - 1.0*x_tmp[32]*xB_tmp[205] + x_tmp[32]*xB_tmp[207] + x_tmp[33]*xB_tmp[208] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[185] + xB_tmp[186]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[183];
qBdot_tmp[6+ip*ny] = x_tmp[7]*xB_tmp[217] - 1.0*x_tmp[7]*xB_tmp[215] + x_tmp[12]*xB_tmp[222] - 1.0*x_tmp[33]*xB_tmp[212] - 1.0*x_tmp[28]*xB_tmp[219] - 1.0*x_tmp[12]*xB_tmp[241] + x_tmp[27]*xB_tmp[237] + x_tmp[28]*xB_tmp[238] - 1.0*x_tmp[27]*xB_tmp[244] - 1.0*x_tmp[32]*xB_tmp[240] + x_tmp[32]*xB_tmp[242] + x_tmp[33]*xB_tmp[243] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[220] + xB_tmp[221]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[218];
qBdot_tmp[7+ip*ny] = x_tmp[7]*xB_tmp[252] - 1.0*x_tmp[7]*xB_tmp[250] + x_tmp[12]*xB_tmp[257] - 1.0*x_tmp[33]*xB_tmp[247] - 1.0*x_tmp[28]*xB_tmp[254] - 1.0*x_tmp[12]*xB_tmp[276] + x_tmp[27]*xB_tmp[272] + x_tmp[28]*xB_tmp[273] - 1.0*x_tmp[27]*xB_tmp[279] - 1.0*x_tmp[32]*xB_tmp[275] + x_tmp[32]*xB_tmp[277] + x_tmp[33]*xB_tmp[278] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[255] + xB_tmp[256]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[253];
qBdot_tmp[8+ip*ny] = x_tmp[7]*xB_tmp[287] - 1.0*x_tmp[7]*xB_tmp[285] + x_tmp[12]*xB_tmp[292] - 1.0*x_tmp[33]*xB_tmp[282] - 1.0*x_tmp[28]*xB_tmp[289] - 1.0*x_tmp[12]*xB_tmp[311] + x_tmp[27]*xB_tmp[307] + x_tmp[28]*xB_tmp[308] - 1.0*x_tmp[27]*xB_tmp[314] - 1.0*x_tmp[32]*xB_tmp[310] + x_tmp[32]*xB_tmp[312] + x_tmp[33]*xB_tmp[313] + (2.0*x_tmp[10] - 1.0*x_tmp[28])*xB_tmp[290] + xB_tmp[291]*(x_tmp[11] - 1.0*x_tmp[10] + x_tmp[28]) - 1.0*(2.0*x_tmp[11] + x_tmp[28])*xB_tmp[288];

  } break;

  case 2: {
qBdot_tmp[0+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[29] - 1.0*x_tmp[11]*xB_tmp[30] - 1.0*x_tmp[10]*xB_tmp[32] - 1.0*x_tmp[28]*xB_tmp[23] - 1.0*x_tmp[33]*xB_tmp[22] - 1.0*xB_tmp[20]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[25]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[26]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[3]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[16]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[19]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[27]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[34]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[21]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[0]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[15]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[12]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[18]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[31]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[24] - 1.0*xB_tmp[17]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[1+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[64] - 1.0*x_tmp[11]*xB_tmp[65] - 1.0*x_tmp[10]*xB_tmp[67] - 1.0*x_tmp[28]*xB_tmp[58] - 1.0*x_tmp[33]*xB_tmp[57] - 1.0*xB_tmp[55]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[60]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[61]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[38]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[51]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[54]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[62]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[69]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[56]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[35]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[50]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[47]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[53]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[66]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[59] - 1.0*xB_tmp[52]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[2+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[99] - 1.0*x_tmp[11]*xB_tmp[100] - 1.0*x_tmp[10]*xB_tmp[102] - 1.0*x_tmp[28]*xB_tmp[93] - 1.0*x_tmp[33]*xB_tmp[92] - 1.0*xB_tmp[90]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[95]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[96]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[73]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[86]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[89]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[97]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[104]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[91]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[70]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[85]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[82]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[88]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[101]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[94] - 1.0*xB_tmp[87]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[3+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[134] - 1.0*x_tmp[11]*xB_tmp[135] - 1.0*x_tmp[10]*xB_tmp[137] - 1.0*x_tmp[28]*xB_tmp[128] - 1.0*x_tmp[33]*xB_tmp[127] - 1.0*xB_tmp[125]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[130]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[131]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[108]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[121]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[124]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[132]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[139]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[126]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[105]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[120]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[117]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[123]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[136]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[129] - 1.0*xB_tmp[122]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[4+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[169] - 1.0*x_tmp[11]*xB_tmp[170] - 1.0*x_tmp[10]*xB_tmp[172] - 1.0*x_tmp[28]*xB_tmp[163] - 1.0*x_tmp[33]*xB_tmp[162] - 1.0*xB_tmp[160]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[165]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[166]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[143]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[156]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[159]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[167]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[174]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[161]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[140]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[155]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[152]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[158]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[171]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[164] - 1.0*xB_tmp[157]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[5+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[204] - 1.0*x_tmp[11]*xB_tmp[205] - 1.0*x_tmp[10]*xB_tmp[207] - 1.0*x_tmp[28]*xB_tmp[198] - 1.0*x_tmp[33]*xB_tmp[197] - 1.0*xB_tmp[195]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[200]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[201]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[178]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[191]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[194]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[202]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[209]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[196]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[175]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[190]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[187]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[193]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[206]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[199] - 1.0*xB_tmp[192]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[6+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[239] - 1.0*x_tmp[11]*xB_tmp[240] - 1.0*x_tmp[10]*xB_tmp[242] - 1.0*x_tmp[28]*xB_tmp[233] - 1.0*x_tmp[33]*xB_tmp[232] - 1.0*xB_tmp[230]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[235]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[236]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[213]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[226]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[229]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[237]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[244]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[231]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[210]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[225]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[222]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[228]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[241]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[234] - 1.0*xB_tmp[227]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[7+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[274] - 1.0*x_tmp[11]*xB_tmp[275] - 1.0*x_tmp[10]*xB_tmp[277] - 1.0*x_tmp[28]*xB_tmp[268] - 1.0*x_tmp[33]*xB_tmp[267] - 1.0*xB_tmp[265]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[270]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[271]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[248]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[261]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[264]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[272]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[279]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[266]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[245]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[260]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[257]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[263]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[276]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[269] - 1.0*xB_tmp[262]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);
qBdot_tmp[8+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[309] - 1.0*x_tmp[11]*xB_tmp[310] - 1.0*x_tmp[10]*xB_tmp[312] - 1.0*x_tmp[28]*xB_tmp[303] - 1.0*x_tmp[33]*xB_tmp[302] - 1.0*xB_tmp[300]*((pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17]) - 1.0*xB_tmp[305]*(2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26]) - 1.0*xB_tmp[306]*((pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24]) - 1.0*xB_tmp[283]*(2.0*(pow(p[3],2))*x_tmp[0] + 0.45*x_tmp[29] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0]) - 1.0*xB_tmp[296]*(2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15]) - 1.0*xB_tmp[299]*(2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18]) - 1.0*xB_tmp[307]*(2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12]) - 1.0*xB_tmp[314]*(2.0*(pow(p[3],2))*x_tmp[31] + 0.45*x_tmp[30] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31]) - 1.0*xB_tmp[301]*(4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20]) - 1.0*xB_tmp[280]*((pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] + 0.55*x_tmp[29]) - 1.0*xB_tmp[295]*((pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22]) - 1.0*xB_tmp[292]*((pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32]) - 1.0*xB_tmp[298]*((pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23]) - 1.0*xB_tmp[311]*((pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.55*x_tmp[30]) - 1.0*(x_tmp[28] + 2.0*x_tmp[32])*xB_tmp[304] - 1.0*xB_tmp[297]*(2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23] + 1.1*x_tmp[26]);

  } break;

  case 3: {
qBdot_tmp[0+ip*ny] = xB_tmp[21]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[17]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[3] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[16] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[27] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[19] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[25] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[34] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[0] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[12] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[15] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[18] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[26] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[31] + xB_tmp[20]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[1+ip*ny] = xB_tmp[56]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[52]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[38] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[51] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[62] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[54] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[60] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[69] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[35] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[47] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[50] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[53] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[61] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[66] + xB_tmp[55]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[2+ip*ny] = xB_tmp[91]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[87]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[73] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[86] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[97] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[89] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[95] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[104] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[70] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[82] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[85] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[88] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[96] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[101] + xB_tmp[90]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[3+ip*ny] = xB_tmp[126]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[122]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[108] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[121] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[132] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[124] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[130] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[139] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[105] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[117] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[120] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[123] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[131] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[136] + xB_tmp[125]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[4+ip*ny] = xB_tmp[161]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[157]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[143] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[156] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[167] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[159] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[165] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[174] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[140] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[152] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[155] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[158] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[166] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[171] + xB_tmp[160]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[5+ip*ny] = xB_tmp[196]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[192]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[178] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[191] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[202] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[194] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[200] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[209] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[175] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[187] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[190] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[193] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[201] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[206] + xB_tmp[195]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[6+ip*ny] = xB_tmp[231]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[227]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[213] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[226] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[237] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[229] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[235] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[244] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[210] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[222] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[225] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[228] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[236] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[241] + xB_tmp[230]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[7+ip*ny] = xB_tmp[266]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[262]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[248] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[261] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[272] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[264] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[270] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[279] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[245] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[257] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[260] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[263] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[271] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[276] + xB_tmp[265]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);
qBdot_tmp[8+ip*ny] = xB_tmp[301]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] - 6.0*p[2]*p[3]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[20]) - 1.0*xB_tmp[297]*(2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18]) + (2.0*(p[3] - 1.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[283] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[296] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[307] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[299] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[305] + (2.0*(p[3] - 1.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[314] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0])*xB_tmp[280] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12])*xB_tmp[292] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15])*xB_tmp[295] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18])*xB_tmp[298] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26])*xB_tmp[306] - 1.0*((2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31])*xB_tmp[311] + xB_tmp[300]*(2.0*(p[3] - 1.0)*p[2]*x_tmp[17] - 1.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*p[3]*x_tmp[17] + 4.0*p[2]*p[3]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[20]);

  } break;

  case 4: {
qBdot_tmp[0+ip*ny] = x_tmp[3]*xB_tmp[3] - 1.0*x_tmp[3]*xB_tmp[1] - 1.0*x_tmp[19]*xB_tmp[13] - 1.0*x_tmp[20]*xB_tmp[15] - 1.0*x_tmp[34]*xB_tmp[2] + x_tmp[19]*xB_tmp[19] + x_tmp[20]*xB_tmp[20] - 1.0*x_tmp[25]*xB_tmp[22] + x_tmp[25]*xB_tmp[25] + x_tmp[27]*xB_tmp[27] - 1.0*x_tmp[27]*xB_tmp[33] + x_tmp[34]*xB_tmp[34] + xB_tmp[16]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[14] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[21];
qBdot_tmp[1+ip*ny] = x_tmp[3]*xB_tmp[38] - 1.0*x_tmp[3]*xB_tmp[36] - 1.0*x_tmp[19]*xB_tmp[48] - 1.0*x_tmp[20]*xB_tmp[50] - 1.0*x_tmp[34]*xB_tmp[37] + x_tmp[19]*xB_tmp[54] + x_tmp[20]*xB_tmp[55] - 1.0*x_tmp[25]*xB_tmp[57] + x_tmp[25]*xB_tmp[60] + x_tmp[27]*xB_tmp[62] - 1.0*x_tmp[27]*xB_tmp[68] + x_tmp[34]*xB_tmp[69] + xB_tmp[51]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[49] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[56];
qBdot_tmp[2+ip*ny] = x_tmp[3]*xB_tmp[73] - 1.0*x_tmp[3]*xB_tmp[71] - 1.0*x_tmp[19]*xB_tmp[83] - 1.0*x_tmp[20]*xB_tmp[85] - 1.0*x_tmp[34]*xB_tmp[72] + x_tmp[19]*xB_tmp[89] + x_tmp[20]*xB_tmp[90] - 1.0*x_tmp[25]*xB_tmp[92] + x_tmp[25]*xB_tmp[95] + x_tmp[27]*xB_tmp[97] - 1.0*x_tmp[27]*xB_tmp[103] + x_tmp[34]*xB_tmp[104] + xB_tmp[86]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[84] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[91];
qBdot_tmp[3+ip*ny] = x_tmp[3]*xB_tmp[108] - 1.0*x_tmp[3]*xB_tmp[106] - 1.0*x_tmp[19]*xB_tmp[118] - 1.0*x_tmp[20]*xB_tmp[120] - 1.0*x_tmp[34]*xB_tmp[107] + x_tmp[19]*xB_tmp[124] + x_tmp[20]*xB_tmp[125] - 1.0*x_tmp[25]*xB_tmp[127] + x_tmp[25]*xB_tmp[130] + x_tmp[27]*xB_tmp[132] - 1.0*x_tmp[27]*xB_tmp[138] + x_tmp[34]*xB_tmp[139] + xB_tmp[121]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[119] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[126];
qBdot_tmp[4+ip*ny] = x_tmp[3]*xB_tmp[143] - 1.0*x_tmp[3]*xB_tmp[141] - 1.0*x_tmp[19]*xB_tmp[153] - 1.0*x_tmp[20]*xB_tmp[155] - 1.0*x_tmp[34]*xB_tmp[142] + x_tmp[19]*xB_tmp[159] + x_tmp[20]*xB_tmp[160] - 1.0*x_tmp[25]*xB_tmp[162] + x_tmp[25]*xB_tmp[165] + x_tmp[27]*xB_tmp[167] - 1.0*x_tmp[27]*xB_tmp[173] + x_tmp[34]*xB_tmp[174] + xB_tmp[156]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[154] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[161];
qBdot_tmp[5+ip*ny] = x_tmp[3]*xB_tmp[178] - 1.0*x_tmp[3]*xB_tmp[176] - 1.0*x_tmp[19]*xB_tmp[188] - 1.0*x_tmp[20]*xB_tmp[190] - 1.0*x_tmp[34]*xB_tmp[177] + x_tmp[19]*xB_tmp[194] + x_tmp[20]*xB_tmp[195] - 1.0*x_tmp[25]*xB_tmp[197] + x_tmp[25]*xB_tmp[200] + x_tmp[27]*xB_tmp[202] - 1.0*x_tmp[27]*xB_tmp[208] + x_tmp[34]*xB_tmp[209] + xB_tmp[191]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[189] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[196];
qBdot_tmp[6+ip*ny] = x_tmp[3]*xB_tmp[213] - 1.0*x_tmp[3]*xB_tmp[211] - 1.0*x_tmp[19]*xB_tmp[223] - 1.0*x_tmp[20]*xB_tmp[225] - 1.0*x_tmp[34]*xB_tmp[212] + x_tmp[19]*xB_tmp[229] + x_tmp[20]*xB_tmp[230] - 1.0*x_tmp[25]*xB_tmp[232] + x_tmp[25]*xB_tmp[235] + x_tmp[27]*xB_tmp[237] - 1.0*x_tmp[27]*xB_tmp[243] + x_tmp[34]*xB_tmp[244] + xB_tmp[226]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[224] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[231];
qBdot_tmp[7+ip*ny] = x_tmp[3]*xB_tmp[248] - 1.0*x_tmp[3]*xB_tmp[246] - 1.0*x_tmp[19]*xB_tmp[258] - 1.0*x_tmp[20]*xB_tmp[260] - 1.0*x_tmp[34]*xB_tmp[247] + x_tmp[19]*xB_tmp[264] + x_tmp[20]*xB_tmp[265] - 1.0*x_tmp[25]*xB_tmp[267] + x_tmp[25]*xB_tmp[270] + x_tmp[27]*xB_tmp[272] - 1.0*x_tmp[27]*xB_tmp[278] + x_tmp[34]*xB_tmp[279] + xB_tmp[261]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[259] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[266];
qBdot_tmp[8+ip*ny] = x_tmp[3]*xB_tmp[283] - 1.0*x_tmp[3]*xB_tmp[281] - 1.0*x_tmp[19]*xB_tmp[293] - 1.0*x_tmp[20]*xB_tmp[295] - 1.0*x_tmp[34]*xB_tmp[282] + x_tmp[19]*xB_tmp[299] + x_tmp[20]*xB_tmp[300] - 1.0*x_tmp[25]*xB_tmp[302] + x_tmp[25]*xB_tmp[305] + x_tmp[27]*xB_tmp[307] - 1.0*x_tmp[27]*xB_tmp[313] + x_tmp[34]*xB_tmp[314] + xB_tmp[296]*(x_tmp[16] + x_tmp[19] - 1.0*x_tmp[21]) - 1.0*(2.0*x_tmp[16] + x_tmp[19])*xB_tmp[294] - 1.0*(x_tmp[19] - 2.0*x_tmp[21])*xB_tmp[301];

  } break;

  case 5: {
qBdot_tmp[0+ip*ny] = x_tmp[1]*xB_tmp[1] + x_tmp[2]*xB_tmp[2] + x_tmp[13]*xB_tmp[13] + x_tmp[15]*xB_tmp[15] + x_tmp[16]*xB_tmp[16] + x_tmp[22]*xB_tmp[22] + x_tmp[33]*xB_tmp[33] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[14];
qBdot_tmp[1+ip*ny] = x_tmp[1]*xB_tmp[36] + x_tmp[2]*xB_tmp[37] + x_tmp[13]*xB_tmp[48] + x_tmp[15]*xB_tmp[50] + x_tmp[16]*xB_tmp[51] + x_tmp[22]*xB_tmp[57] + x_tmp[33]*xB_tmp[68] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[49];
qBdot_tmp[2+ip*ny] = x_tmp[1]*xB_tmp[71] + x_tmp[2]*xB_tmp[72] + x_tmp[13]*xB_tmp[83] + x_tmp[15]*xB_tmp[85] + x_tmp[16]*xB_tmp[86] + x_tmp[22]*xB_tmp[92] + x_tmp[33]*xB_tmp[103] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[84];
qBdot_tmp[3+ip*ny] = x_tmp[1]*xB_tmp[106] + x_tmp[2]*xB_tmp[107] + x_tmp[13]*xB_tmp[118] + x_tmp[15]*xB_tmp[120] + x_tmp[16]*xB_tmp[121] + x_tmp[22]*xB_tmp[127] + x_tmp[33]*xB_tmp[138] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[119];
qBdot_tmp[4+ip*ny] = x_tmp[1]*xB_tmp[141] + x_tmp[2]*xB_tmp[142] + x_tmp[13]*xB_tmp[153] + x_tmp[15]*xB_tmp[155] + x_tmp[16]*xB_tmp[156] + x_tmp[22]*xB_tmp[162] + x_tmp[33]*xB_tmp[173] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[154];
qBdot_tmp[5+ip*ny] = x_tmp[1]*xB_tmp[176] + x_tmp[2]*xB_tmp[177] + x_tmp[13]*xB_tmp[188] + x_tmp[15]*xB_tmp[190] + x_tmp[16]*xB_tmp[191] + x_tmp[22]*xB_tmp[197] + x_tmp[33]*xB_tmp[208] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[189];
qBdot_tmp[6+ip*ny] = x_tmp[1]*xB_tmp[211] + x_tmp[2]*xB_tmp[212] + x_tmp[13]*xB_tmp[223] + x_tmp[15]*xB_tmp[225] + x_tmp[16]*xB_tmp[226] + x_tmp[22]*xB_tmp[232] + x_tmp[33]*xB_tmp[243] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[224];
qBdot_tmp[7+ip*ny] = x_tmp[1]*xB_tmp[246] + x_tmp[2]*xB_tmp[247] + x_tmp[13]*xB_tmp[258] + x_tmp[15]*xB_tmp[260] + x_tmp[16]*xB_tmp[261] + x_tmp[22]*xB_tmp[267] + x_tmp[33]*xB_tmp[278] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[259];
qBdot_tmp[8+ip*ny] = x_tmp[1]*xB_tmp[281] + x_tmp[2]*xB_tmp[282] + x_tmp[13]*xB_tmp[293] + x_tmp[15]*xB_tmp[295] + x_tmp[16]*xB_tmp[296] + x_tmp[22]*xB_tmp[302] + x_tmp[33]*xB_tmp[313] - 1.0*(x_tmp[13] - 2.0*x_tmp[14])*xB_tmp[294];

  } break;

  }
  }

  for (iyp=0; iyp<9*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*35);
x0_tmp[0] = k[11]*k[46];
x0_tmp[1] = k[13]*k[48];
x0_tmp[2] = k[19]*k[54];
x0_tmp[3] = k[12]*k[47];
x0_tmp[4] = (k[7] - 1.0)*((pow(p[6],2)) - 1.0*p[6]) + k[7]*k[42];
x0_tmp[5] = k[8]*k[43] + (k[8] - 1.0)*p[6]*p[7];
x0_tmp[6] = k[0]*k[35] - 1.0*(k[0] - 1.0)*p[6];
x0_tmp[7] = k[9]*k[44] - 1.0*(k[9] - 1.0)*p[6]*(p[6] + p[7] - 1.0);
x0_tmp[8] = (k[14] - 1.0)*((pow(p[7],2)) - 1.0*p[7]) + k[14]*k[49];
x0_tmp[9] = k[1]*k[36] - 1.0*(k[1] - 1.0)*p[7];
x0_tmp[10] = k[20]*k[55] + (k[20] - 1.0)*(p[6] + p[7])*(p[6] + p[7] - 1.0);
x0_tmp[11] = k[15]*k[50] - 1.0*(k[15] - 1.0)*p[7]*(p[6] + p[7] - 1.0);
x0_tmp[12] = k[22]*k[57];
x0_tmp[13] = k[6]*k[41];
x0_tmp[14] = k[34]*k[69];
x0_tmp[15] = k[31]*k[66];
x0_tmp[16] = k[33]*k[68];
x0_tmp[17] = k[29]*k[64];
x0_tmp[18] = k[4]*k[39];
x0_tmp[19] = k[5]*k[40];
x0_tmp[20] = k[30]*k[65];
x0_tmp[21] = k[32]*k[67];
x0_tmp[22] = k[28]*k[63];
x0_tmp[23] = k[3]*k[38];
x0_tmp[24] = k[25]*k[60];
x0_tmp[25] = k[27]*k[62];
x0_tmp[26] = k[26]*k[61];
x0_tmp[27] = k[23]*k[58];
x0_tmp[28] = (k[2] - 1.0)*(p[6] + p[7] - 1.0) + k[2]*k[37];
x0_tmp[29] = k[10]*k[45];
x0_tmp[30] = k[16]*k[51];
x0_tmp[31] = k[17]*k[52];
x0_tmp[32] = k[21]*k[56];
x0_tmp[33] = k[24]*k[59];
x0_tmp[34] = k[18]*k[53];
  
  
  return;
}


 int Jv__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(N_Vector v, N_Vector Jv, realtype t,
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
Jv_tmp[0] = 0.55*p[2]*v_tmp[29] - 1.0*v_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
Jv_tmp[1] = p[4]*v_tmp[3] - 1.0*(p[5] + 0.0002)*v_tmp[1];
Jv_tmp[2] = p[1]*v_tmp[33] + p[4]*v_tmp[34] - 1.0*v_tmp[2]*(p[0] + p[5]) + 0.0002*v_tmp[1];
Jv_tmp[3] = 0.45*p[2]*v_tmp[29] - 1.0*(p[4] + 0.0002)*v_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[0];
Jv_tmp[4] = 0.0002*v_tmp[6] - 0.0004*v_tmp[4];
Jv_tmp[5] = p[1]*v_tmp[7] - 1.0*(p[0] + 0.0002)*v_tmp[5] + 0.0002*v_tmp[4] - 0.0002*v_tmp[6];
Jv_tmp[6] = -0.0002*v_tmp[6];
Jv_tmp[7] = p[0]*v_tmp[5] - 1.0*(p[1] + 0.0002)*v_tmp[7];
Jv_tmp[8] = p[0]*v_tmp[9] - 2.0*p[0]*v_tmp[8] + 2.0*p[1]*v_tmp[11] + p[1]*v_tmp[28] + 0.0004*v_tmp[5] + 0.0002*v_tmp[6];
Jv_tmp[9] = p[1]*v_tmp[28] - 1.0*p[0]*v_tmp[9] + 0.0002*v_tmp[6];
Jv_tmp[10] = p[0]*v_tmp[9] + 2.0*p[0]*v_tmp[11] - 2.0*p[1]*v_tmp[10] + p[1]*v_tmp[28];
Jv_tmp[11] = p[0]*v_tmp[8] - 1.0*p[0]*v_tmp[9] + p[1]*v_tmp[10] - 1.0*p[1]*v_tmp[28] - 1.0*v_tmp[11]*(p[0] + p[1]) + 0.0002*v_tmp[7];
Jv_tmp[12] = p[0]*v_tmp[31] + 0.55*p[2]*v_tmp[32] - 1.0*v_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[13] = p[4]*v_tmp[19] - 1.0*p[5]*v_tmp[13];
Jv_tmp[14] = p[5]*v_tmp[13] - 2.0*p[5]*v_tmp[14] + 2.0*p[4]*v_tmp[16] + p[4]*v_tmp[19];
Jv_tmp[15] = 0.55*p[2]*v_tmp[22] + p[4]*v_tmp[20] - 1.0*v_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[16] = 0.45*p[2]*v_tmp[22] - 1.0*p[4]*v_tmp[19] + p[4]*v_tmp[21] - 1.0*v_tmp[16]*(p[4] + p[5]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[15];
Jv_tmp[17] = 0.55*p[2]*v_tmp[23] + 1.1*p[2]*v_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*v_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*v_tmp[18];
Jv_tmp[18] = 0.55*p[2]*v_tmp[23] - 1.0*v_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[19] = 0.45*p[2]*v_tmp[23] - 1.0*p[4]*v_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[18];
Jv_tmp[20] = 0.55*p[2]*v_tmp[25] + 0.45*p[2]*v_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[17] - 1.0*v_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*v_tmp[18];
Jv_tmp[21] = p[4]*v_tmp[19] + 0.45*p[2]*v_tmp[23] - 2.0*p[4]*v_tmp[21] + 0.9*p[2]*v_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[20];
Jv_tmp[22] = p[4]*v_tmp[25] - 1.0*p[5]*v_tmp[22] + p[2]*v_tmp[33];
Jv_tmp[23] = p[2]*v_tmp[28];
Jv_tmp[24] = p[2]*v_tmp[28] + 2.0*p[2]*v_tmp[32];
Jv_tmp[25] = 0.45*p[2]*v_tmp[24] + p[2]*v_tmp[27] - 1.0*p[4]*v_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[26];
Jv_tmp[26] = p[2]*v_tmp[12] + 0.55*p[2]*v_tmp[24] - 1.0*v_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[27] = p[0]*v_tmp[34] + 0.45*p[2]*v_tmp[32] - 1.0*v_tmp[27]*(p[1] + p[4]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[12];
Jv_tmp[28] = p[0]*v_tmp[9] - 1.0*p[1]*v_tmp[28];
Jv_tmp[29] = p[2]*v_tmp[7] - 0.0002*v_tmp[29];
Jv_tmp[30] = p[2]*v_tmp[11] - 1.0*p[0]*v_tmp[30] + p[1]*v_tmp[32] + 0.0002*v_tmp[29];
Jv_tmp[31] = p[1]*v_tmp[12] + 0.55*p[2]*v_tmp[30] + 0.0002*v_tmp[0] - 1.0*v_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
Jv_tmp[32] = p[2]*v_tmp[10] + p[0]*v_tmp[30] - 1.0*p[1]*v_tmp[32];
Jv_tmp[33] = p[0]*v_tmp[2] + p[4]*v_tmp[27] - 1.0*v_tmp[33]*(p[1] + p[5]);
Jv_tmp[34] = p[1]*v_tmp[27] + 0.45*p[2]*v_tmp[30] - 1.0*v_tmp[34]*(p[0] + p[4]) + 0.0002*v_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*v_tmp[31];

  for (ix=0; ix<35; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(N_Vector vB, N_Vector JvB, realtype t,
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
JvB_tmp[0] = vB_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) - 0.0002*vB_tmp[31] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[3];
JvB_tmp[1] = (p[5] + 0.0002)*vB_tmp[1] - 0.0002*vB_tmp[2];
JvB_tmp[2] = vB_tmp[2]*(p[0] + p[5]) - 1.0*p[0]*vB_tmp[33];
JvB_tmp[3] = (p[4] + 0.0002)*vB_tmp[3] - 1.0*p[4]*vB_tmp[1] - 0.0002*vB_tmp[34];
JvB_tmp[4] = 0.0004*vB_tmp[4] - 0.0002*vB_tmp[5];
JvB_tmp[5] = (p[0] + 0.0002)*vB_tmp[5] - 1.0*p[0]*vB_tmp[7] - 0.0004*vB_tmp[8];
JvB_tmp[6] = 0.0002*vB_tmp[5] - 0.0002*vB_tmp[4] + 0.0002*vB_tmp[6] - 0.0002*vB_tmp[8] - 0.0002*vB_tmp[9];
JvB_tmp[7] = (p[1] + 0.0002)*vB_tmp[7] - 1.0*p[2]*vB_tmp[29] - 1.0*p[1]*vB_tmp[5] - 0.0002*vB_tmp[11];
JvB_tmp[8] = 2.0*p[0]*vB_tmp[8] - 1.0*p[0]*vB_tmp[11];
JvB_tmp[9] = p[0]*vB_tmp[9] - 1.0*p[0]*vB_tmp[8] - 1.0*p[0]*vB_tmp[10] + p[0]*vB_tmp[11] - 1.0*p[0]*vB_tmp[28];
JvB_tmp[10] = 2.0*p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[11] - 1.0*p[2]*vB_tmp[32];
JvB_tmp[11] = vB_tmp[11]*(p[0] + p[1]) - 2.0*p[0]*vB_tmp[10] - 1.0*p[2]*vB_tmp[30] - 2.0*p[1]*vB_tmp[8];
JvB_tmp[12] = vB_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*p[1]*vB_tmp[31] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[27] - 1.0*p[2]*vB_tmp[26];
JvB_tmp[13] = p[5]*vB_tmp[13] - 1.0*p[5]*vB_tmp[14];
JvB_tmp[14] = 2.0*p[5]*vB_tmp[14];
JvB_tmp[15] = vB_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[16];
JvB_tmp[16] = vB_tmp[16]*(p[4] + p[5]) - 2.0*p[4]*vB_tmp[14];
JvB_tmp[17] = (2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*vB_tmp[17] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[20];
JvB_tmp[18] = vB_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[19] - 1.0*(4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[21] - 1.0*(p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*vB_tmp[17] + 2.0*p[2]*(pow(p[3],2))*vB_tmp[20];
JvB_tmp[19] = p[4]*vB_tmp[16] - 1.0*p[4]*vB_tmp[14] - 1.0*p[4]*vB_tmp[13] + p[4]*vB_tmp[19] - 1.0*p[4]*vB_tmp[21];
JvB_tmp[20] = vB_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[21] - 1.0*p[4]*vB_tmp[15];
JvB_tmp[21] = 2.0*p[4]*vB_tmp[21] - 1.0*p[4]*vB_tmp[16];
JvB_tmp[22] = p[5]*vB_tmp[22] - 0.45*p[2]*vB_tmp[16] - 0.55*p[2]*vB_tmp[15];
JvB_tmp[23] = - 0.55*p[2]*vB_tmp[17] - 0.55*p[2]*vB_tmp[18] - 0.45*p[2]*vB_tmp[19] - 0.45*p[2]*vB_tmp[21];
JvB_tmp[24] = - 0.45*p[2]*vB_tmp[25] - 0.55*p[2]*vB_tmp[26];
JvB_tmp[25] = p[4]*vB_tmp[25] - 0.9*p[2]*vB_tmp[21] - 1.0*p[4]*vB_tmp[22] - 0.55*p[2]*vB_tmp[20];
JvB_tmp[26] = vB_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 0.45*p[2]*vB_tmp[20] - 1.1*p[2]*vB_tmp[17] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[25];
JvB_tmp[27] = vB_tmp[27]*(p[1] + p[4]) - 1.0*p[1]*vB_tmp[34] - 1.0*p[4]*vB_tmp[33] - 1.0*p[2]*vB_tmp[25];
JvB_tmp[28] = p[1]*vB_tmp[11] - 1.0*p[1]*vB_tmp[9] - 1.0*p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[8] - 1.0*p[2]*vB_tmp[23] - 1.0*p[2]*vB_tmp[24] + p[1]*vB_tmp[28];
JvB_tmp[29] = 0.0002*vB_tmp[29] - 0.45*p[2]*vB_tmp[3] - 0.55*p[2]*vB_tmp[0] - 0.0002*vB_tmp[30];
JvB_tmp[30] = p[0]*vB_tmp[30] - 1.0*p[0]*vB_tmp[32] - 0.55*p[2]*vB_tmp[31] - 0.45*p[2]*vB_tmp[34];
JvB_tmp[31] = vB_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*vB_tmp[34] - 1.0*p[0]*vB_tmp[12];
JvB_tmp[32] = p[1]*vB_tmp[32] - 2.0*p[2]*vB_tmp[24] - 0.45*p[2]*vB_tmp[27] - 1.0*p[1]*vB_tmp[30] - 0.55*p[2]*vB_tmp[12];
JvB_tmp[33] = vB_tmp[33]*(p[1] + p[5]) - 1.0*p[2]*vB_tmp[22] - 1.0*p[1]*vB_tmp[2];
JvB_tmp[34] = vB_tmp[34]*(p[0] + p[4]) - 1.0*p[0]*vB_tmp[27] - 1.0*p[4]*vB_tmp[2];

  for (ix=0; ix<35; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(long int N, realtype t, N_Vector x,
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
J->data[0] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2)) - 0.0002;
J->data[3] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[31] = 0.0002;
J->data[36] = - 1.0*p[5] - 0.0002;
J->data[37] = 0.0002;
J->data[72] = - 1.0*p[0] - 1.0*p[5];
J->data[103] = p[0];
J->data[106] = p[4];
J->data[108] = - 1.0*p[4] - 0.0002;
J->data[139] = 0.0002;
J->data[144] = -0.0004;
J->data[145] = 0.0002;
J->data[180] = - 1.0*p[0] - 0.0002;
J->data[182] = p[0];
J->data[183] = 0.0004;
J->data[214] = 0.0002;
J->data[215] = -0.0002;
J->data[216] = -0.0002;
J->data[218] = 0.0002;
J->data[219] = 0.0002;
J->data[250] = p[1];
J->data[252] = - 1.0*p[1] - 0.0002;
J->data[256] = 0.0002;
J->data[274] = p[2];
J->data[288] = -2.0*p[0];
J->data[291] = p[0];
J->data[323] = p[0];
J->data[324] = -1.0*p[0];
J->data[325] = p[0];
J->data[326] = -1.0*p[0];
J->data[343] = p[0];
J->data[360] = -2.0*p[1];
J->data[361] = p[1];
J->data[382] = p[2];
J->data[393] = 2.0*p[1];
J->data[395] = 2.0*p[0];
J->data[396] = - 1.0*p[0] - 1.0*p[1];
J->data[415] = p[2];
J->data[432] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[1] - 1.0*p[2]*(pow(p[3],2));
J->data[446] = p[2];
J->data[447] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[451] = p[1];
J->data[468] = -1.0*p[5];
J->data[469] = p[5];
J->data[504] = -2.0*p[5];
J->data[540] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[5] - 1.0*p[2]*(pow(p[3],2));
J->data[541] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[574] = 2.0*p[4];
J->data[576] = - 1.0*p[4] - 1.0*p[5];
J->data[612] = 2.0*(pow((p[3] - 1.0),2))*p[2] - 2.0*p[2]*(pow(p[3],2));
J->data[615] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[647] = p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2];
J->data[648] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[649] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[650] = -2.0*p[2]*(pow(p[3],2));
J->data[651] = 4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[678] = p[4];
J->data[679] = p[4];
J->data[681] = -1.0*p[4];
J->data[684] = -1.0*p[4];
J->data[686] = p[4];
J->data[715] = p[4];
J->data[720] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[4] - 1.0*p[2]*(pow(p[3],2));
J->data[721] = 4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3];
J->data[751] = p[4];
J->data[756] = -2.0*p[4];
J->data[785] = 0.55*p[2];
J->data[786] = 0.45*p[2];
J->data[792] = -1.0*p[5];
J->data[822] = 0.55*p[2];
J->data[823] = 0.55*p[2];
J->data[824] = 0.45*p[2];
J->data[826] = 0.45*p[2];
J->data[865] = 0.45*p[2];
J->data[866] = 0.55*p[2];
J->data[895] = 0.55*p[2];
J->data[896] = 0.9*p[2];
J->data[897] = p[4];
J->data[900] = -1.0*p[4];
J->data[927] = 1.1*p[2];
J->data[930] = 0.45*p[2];
J->data[935] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[936] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[970] = p[2];
J->data[972] = - 1.0*p[1] - 1.0*p[4];
J->data[978] = p[4];
J->data[979] = p[1];
J->data[988] = p[1];
J->data[989] = p[1];
J->data[990] = p[1];
J->data[991] = -1.0*p[1];
J->data[1003] = p[2];
J->data[1004] = p[2];
J->data[1008] = -1.0*p[1];
J->data[1015] = 0.55*p[2];
J->data[1018] = 0.45*p[2];
J->data[1044] = -0.0002;
J->data[1045] = 0.0002;
J->data[1080] = -1.0*p[0];
J->data[1081] = 0.55*p[2];
J->data[1082] = p[0];
J->data[1084] = 0.45*p[2];
J->data[1097] = p[0];
J->data[1116] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[0] - 1.0*p[2]*(pow(p[3],2));
J->data[1119] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[1132] = 0.55*p[2];
J->data[1144] = 2.0*p[2];
J->data[1147] = 0.45*p[2];
J->data[1150] = p[1];
J->data[1152] = -1.0*p[1];
J->data[1157] = p[1];
J->data[1177] = p[2];
J->data[1188] = - 1.0*p[1] - 1.0*p[5];
J->data[1192] = p[4];
J->data[1217] = p[0];
J->data[1224] = - 1.0*p[0] - 1.0*p[4];

  for (iJ=0; iJ<1225; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(realtype t, N_Vector x,
  	N_Vector xdot, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(J);
  J->rowvals[0] = 0;
  J->rowvals[1] = 3;
  J->rowvals[2] = 31;
  J->rowvals[3] = 1;
  J->rowvals[4] = 2;
  J->rowvals[5] = 2;
  J->rowvals[6] = 33;
  J->rowvals[7] = 1;
  J->rowvals[8] = 3;
  J->rowvals[9] = 34;
  J->rowvals[10] = 4;
  J->rowvals[11] = 5;
  J->rowvals[12] = 5;
  J->rowvals[13] = 7;
  J->rowvals[14] = 8;
  J->rowvals[15] = 4;
  J->rowvals[16] = 5;
  J->rowvals[17] = 6;
  J->rowvals[18] = 8;
  J->rowvals[19] = 9;
  J->rowvals[20] = 5;
  J->rowvals[21] = 7;
  J->rowvals[22] = 11;
  J->rowvals[23] = 29;
  J->rowvals[24] = 8;
  J->rowvals[25] = 11;
  J->rowvals[26] = 8;
  J->rowvals[27] = 9;
  J->rowvals[28] = 10;
  J->rowvals[29] = 11;
  J->rowvals[30] = 28;
  J->rowvals[31] = 10;
  J->rowvals[32] = 11;
  J->rowvals[33] = 32;
  J->rowvals[34] = 8;
  J->rowvals[35] = 10;
  J->rowvals[36] = 11;
  J->rowvals[37] = 30;
  J->rowvals[38] = 12;
  J->rowvals[39] = 26;
  J->rowvals[40] = 27;
  J->rowvals[41] = 31;
  J->rowvals[42] = 13;
  J->rowvals[43] = 14;
  J->rowvals[44] = 14;
  J->rowvals[45] = 15;
  J->rowvals[46] = 16;
  J->rowvals[47] = 14;
  J->rowvals[48] = 16;
  J->rowvals[49] = 17;
  J->rowvals[50] = 20;
  J->rowvals[51] = 17;
  J->rowvals[52] = 18;
  J->rowvals[53] = 19;
  J->rowvals[54] = 20;
  J->rowvals[55] = 21;
  J->rowvals[56] = 13;
  J->rowvals[57] = 14;
  J->rowvals[58] = 16;
  J->rowvals[59] = 19;
  J->rowvals[60] = 21;
  J->rowvals[61] = 15;
  J->rowvals[62] = 20;
  J->rowvals[63] = 21;
  J->rowvals[64] = 16;
  J->rowvals[65] = 21;
  J->rowvals[66] = 15;
  J->rowvals[67] = 16;
  J->rowvals[68] = 22;
  J->rowvals[69] = 17;
  J->rowvals[70] = 18;
  J->rowvals[71] = 19;
  J->rowvals[72] = 21;
  J->rowvals[73] = 25;
  J->rowvals[74] = 26;
  J->rowvals[75] = 20;
  J->rowvals[76] = 21;
  J->rowvals[77] = 22;
  J->rowvals[78] = 25;
  J->rowvals[79] = 17;
  J->rowvals[80] = 20;
  J->rowvals[81] = 25;
  J->rowvals[82] = 26;
  J->rowvals[83] = 25;
  J->rowvals[84] = 27;
  J->rowvals[85] = 33;
  J->rowvals[86] = 34;
  J->rowvals[87] = 8;
  J->rowvals[88] = 9;
  J->rowvals[89] = 10;
  J->rowvals[90] = 11;
  J->rowvals[91] = 23;
  J->rowvals[92] = 24;
  J->rowvals[93] = 28;
  J->rowvals[94] = 0;
  J->rowvals[95] = 3;
  J->rowvals[96] = 29;
  J->rowvals[97] = 30;
  J->rowvals[98] = 30;
  J->rowvals[99] = 31;
  J->rowvals[100] = 32;
  J->rowvals[101] = 34;
  J->rowvals[102] = 12;
  J->rowvals[103] = 31;
  J->rowvals[104] = 34;
  J->rowvals[105] = 12;
  J->rowvals[106] = 24;
  J->rowvals[107] = 27;
  J->rowvals[108] = 30;
  J->rowvals[109] = 32;
  J->rowvals[110] = 2;
  J->rowvals[111] = 22;
  J->rowvals[112] = 33;
  J->rowvals[113] = 2;
  J->rowvals[114] = 27;
  J->rowvals[115] = 34;
  J->colptrs[0] = 0;
  J->colptrs[1] = 3;
  J->colptrs[2] = 5;
  J->colptrs[3] = 7;
  J->colptrs[4] = 10;
  J->colptrs[5] = 12;
  J->colptrs[6] = 15;
  J->colptrs[7] = 20;
  J->colptrs[8] = 24;
  J->colptrs[9] = 26;
  J->colptrs[10] = 31;
  J->colptrs[11] = 34;
  J->colptrs[12] = 38;
  J->colptrs[13] = 42;
  J->colptrs[14] = 44;
  J->colptrs[15] = 45;
  J->colptrs[16] = 47;
  J->colptrs[17] = 49;
  J->colptrs[18] = 51;
  J->colptrs[19] = 56;
  J->colptrs[20] = 61;
  J->colptrs[21] = 64;
  J->colptrs[22] = 66;
  J->colptrs[23] = 69;
  J->colptrs[24] = 73;
  J->colptrs[25] = 75;
  J->colptrs[26] = 79;
  J->colptrs[27] = 83;
  J->colptrs[28] = 87;
  J->colptrs[29] = 94;
  J->colptrs[30] = 98;
  J->colptrs[31] = 102;
  J->colptrs[32] = 105;
  J->colptrs[33] = 110;
  J->colptrs[34] = 113;
  J->colptrs[35] = 116;
J->data[0] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2)) - 0.0002;
J->data[1] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[2] = 0.0002;
J->data[3] = - 1.0*p[5] - 0.0002;
J->data[4] = 0.0002;
J->data[5] = - 1.0*p[0] - 1.0*p[5];
J->data[6] = p[0];
J->data[7] = p[4];
J->data[8] = - 1.0*p[4] - 0.0002;
J->data[9] = 0.0002;
J->data[10] = -0.0004;
J->data[11] = 0.0002;
J->data[12] = - 1.0*p[0] - 0.0002;
J->data[13] = p[0];
J->data[14] = 0.0004;
J->data[15] = 0.0002;
J->data[16] = -0.0002;
J->data[17] = -0.0002;
J->data[18] = 0.0002;
J->data[19] = 0.0002;
J->data[20] = p[1];
J->data[21] = - 1.0*p[1] - 0.0002;
J->data[22] = 0.0002;
J->data[23] = p[2];
J->data[24] = -2.0*p[0];
J->data[25] = p[0];
J->data[26] = p[0];
J->data[27] = -1.0*p[0];
J->data[28] = p[0];
J->data[29] = -1.0*p[0];
J->data[30] = p[0];
J->data[31] = -2.0*p[1];
J->data[32] = p[1];
J->data[33] = p[2];
J->data[34] = 2.0*p[1];
J->data[35] = 2.0*p[0];
J->data[36] = - 1.0*p[0] - 1.0*p[1];
J->data[37] = p[2];
J->data[38] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[1] - 1.0*p[2]*(pow(p[3],2));
J->data[39] = p[2];
J->data[40] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[41] = p[1];
J->data[42] = -1.0*p[5];
J->data[43] = p[5];
J->data[44] = -2.0*p[5];
J->data[45] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[5] - 1.0*p[2]*(pow(p[3],2));
J->data[46] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[47] = 2.0*p[4];
J->data[48] = - 1.0*p[4] - 1.0*p[5];
J->data[49] = 2.0*(pow((p[3] - 1.0),2))*p[2] - 2.0*p[2]*(pow(p[3],2));
J->data[50] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[51] = p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2];
J->data[52] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[53] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[54] = -2.0*p[2]*(pow(p[3],2));
J->data[55] = 4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[56] = p[4];
J->data[57] = p[4];
J->data[58] = -1.0*p[4];
J->data[59] = -1.0*p[4];
J->data[60] = p[4];
J->data[61] = p[4];
J->data[62] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[4] - 1.0*p[2]*(pow(p[3],2));
J->data[63] = 4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3];
J->data[64] = p[4];
J->data[65] = -2.0*p[4];
J->data[66] = 0.55*p[2];
J->data[67] = 0.45*p[2];
J->data[68] = -1.0*p[5];
J->data[69] = 0.55*p[2];
J->data[70] = 0.55*p[2];
J->data[71] = 0.45*p[2];
J->data[72] = 0.45*p[2];
J->data[73] = 0.45*p[2];
J->data[74] = 0.55*p[2];
J->data[75] = 0.55*p[2];
J->data[76] = 0.9*p[2];
J->data[77] = p[4];
J->data[78] = -1.0*p[4];
J->data[79] = 1.1*p[2];
J->data[80] = 0.45*p[2];
J->data[81] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[82] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[2]*(pow(p[3],2));
J->data[83] = p[2];
J->data[84] = - 1.0*p[1] - 1.0*p[4];
J->data[85] = p[4];
J->data[86] = p[1];
J->data[87] = p[1];
J->data[88] = p[1];
J->data[89] = p[1];
J->data[90] = -1.0*p[1];
J->data[91] = p[2];
J->data[92] = p[2];
J->data[93] = -1.0*p[1];
J->data[94] = 0.55*p[2];
J->data[95] = 0.45*p[2];
J->data[96] = -0.0002;
J->data[97] = 0.0002;
J->data[98] = -1.0*p[0];
J->data[99] = 0.55*p[2];
J->data[100] = p[0];
J->data[101] = 0.45*p[2];
J->data[102] = p[0];
J->data[103] = (pow((p[3] - 1.0),2))*p[2] - 1.0*p[0] - 1.0*p[2]*(pow(p[3],2));
J->data[104] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3];
J->data[105] = 0.55*p[2];
J->data[106] = 2.0*p[2];
J->data[107] = 0.45*p[2];
J->data[108] = p[1];
J->data[109] = -1.0*p[1];
J->data[110] = p[1];
J->data[111] = p[2];
J->data[112] = - 1.0*p[1] - 1.0*p[5];
J->data[113] = p[4];
J->data[114] = p[0];
J->data[115] = - 1.0*p[0] - 1.0*p[4];
  return(0);
}


 int JBBand__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(long int N, realtype t, N_Vector x,
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
JB->data[0] = p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002;
JB->data[29] = -0.55*p[2];
JB->data[36] = p[5] + 0.0002;
JB->data[38] = -1.0*p[4];
JB->data[71] = -0.0002;
JB->data[72] = p[0] + p[5];
JB->data[103] = -1.0*p[1];
JB->data[104] = -1.0*p[4];
JB->data[105] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[108] = p[4] + 0.0002;
JB->data[134] = -0.45*p[2];
JB->data[144] = 0.0004;
JB->data[146] = -0.0002;
JB->data[179] = -0.0002;
JB->data[180] = p[0] + 0.0002;
JB->data[181] = 0.0002;
JB->data[182] = -1.0*p[1];
JB->data[216] = 0.0002;
JB->data[250] = -1.0*p[0];
JB->data[252] = p[1] + 0.0002;
JB->data[285] = -0.0004;
JB->data[286] = -0.0002;
JB->data[288] = 2.0*p[0];
JB->data[289] = -1.0*p[0];
JB->data[291] = -2.0*p[1];
JB->data[308] = -1.0*p[1];
JB->data[321] = -0.0002;
JB->data[324] = p[0];
JB->data[343] = -1.0*p[1];
JB->data[359] = -1.0*p[0];
JB->data[360] = 2.0*p[1];
JB->data[361] = -2.0*p[0];
JB->data[378] = -1.0*p[1];
JB->data[392] = -0.0002;
JB->data[393] = -1.0*p[0];
JB->data[394] = p[0];
JB->data[395] = -1.0*p[1];
JB->data[396] = p[0] + p[1];
JB->data[413] = p[1];
JB->data[432] = p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[451] = -1.0*p[0];
JB->data[452] = -0.55*p[2];
JB->data[468] = p[5];
JB->data[474] = -1.0*p[4];
JB->data[503] = -1.0*p[5];
JB->data[504] = 2.0*p[5];
JB->data[506] = -2.0*p[4];
JB->data[509] = -1.0*p[4];
JB->data[540] = p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[545] = -1.0*p[4];
JB->data[547] = -0.55*p[2];
JB->data[575] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[576] = p[4] + p[5];
JB->data[579] = p[4];
JB->data[581] = -1.0*p[4];
JB->data[582] = -0.45*p[2];
JB->data[612] = 2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[613] = - 1.0*p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[618] = -0.55*p[2];
JB->data[621] = -1.1*p[2];
JB->data[648] = p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[653] = -0.55*p[2];
JB->data[683] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[684] = p[4];
JB->data[688] = -0.45*p[2];
JB->data[717] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[718] = 2.0*p[2]*(pow(p[3],2));
JB->data[720] = p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[725] = -0.55*p[2];
JB->data[726] = -0.45*p[2];
JB->data[753] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 4.0*p[2]*(pow(p[3],2));
JB->data[754] = -1.0*p[4];
JB->data[755] = 4.0*(p[3] - 1.0)*p[2]*p[3] - 4.0*p[2]*(pow(p[3],2));
JB->data[756] = 2.0*p[4];
JB->data[758] = -0.45*p[2];
JB->data[760] = -0.9*p[2];
JB->data[792] = p[5];
JB->data[795] = -1.0*p[4];
JB->data[803] = -1.0*p[2];
JB->data[833] = -1.0*p[2];
JB->data[868] = -1.0*p[2];
JB->data[872] = -2.0*p[2];
JB->data[899] = -0.45*p[2];
JB->data[900] = p[4];
JB->data[901] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[902] = -1.0*p[2];
JB->data[922] = -1.0*p[2];
JB->data[934] = -0.55*p[2];
JB->data[936] = p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[957] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[972] = p[1] + p[4];
JB->data[977] = -0.45*p[2];
JB->data[979] = -1.0*p[0];
JB->data[989] = -1.0*p[0];
JB->data[1008] = p[1];
JB->data[1022] = -1.0*p[2];
JB->data[1044] = 0.0002;
JB->data[1061] = -1.0*p[2];
JB->data[1079] = -0.0002;
JB->data[1080] = p[0];
JB->data[1082] = -1.0*p[1];
JB->data[1085] = -0.0002;
JB->data[1097] = -1.0*p[1];
JB->data[1115] = -0.55*p[2];
JB->data[1116] = p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2];
JB->data[1130] = -1.0*p[2];
JB->data[1150] = -1.0*p[0];
JB->data[1152] = p[1];
JB->data[1157] = -1.0*p[0];
JB->data[1182] = -1.0*p[4];
JB->data[1188] = p[1] + p[5];
JB->data[1193] = -0.0002;
JB->data[1217] = -1.0*p[1];
JB->data[1220] = -0.45*p[2];
JB->data[1221] = 2.0*(p[3] - 1.0)*p[2]*p[3] - 2.0*p[2]*(pow(p[3],2));
JB->data[1224] = p[0] + p[4];

  for (iJ=0; iJ<1225; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, SlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(JB);
  JB->rowvals[0] = 0;
  JB->rowvals[1] = 29;
  JB->rowvals[2] = 1;
  JB->rowvals[3] = 3;
  JB->rowvals[4] = 1;
  JB->rowvals[5] = 2;
  JB->rowvals[6] = 33;
  JB->rowvals[7] = 34;
  JB->rowvals[8] = 0;
  JB->rowvals[9] = 3;
  JB->rowvals[10] = 29;
  JB->rowvals[11] = 4;
  JB->rowvals[12] = 6;
  JB->rowvals[13] = 4;
  JB->rowvals[14] = 5;
  JB->rowvals[15] = 6;
  JB->rowvals[16] = 7;
  JB->rowvals[17] = 6;
  JB->rowvals[18] = 5;
  JB->rowvals[19] = 7;
  JB->rowvals[20] = 5;
  JB->rowvals[21] = 6;
  JB->rowvals[22] = 8;
  JB->rowvals[23] = 9;
  JB->rowvals[24] = 11;
  JB->rowvals[25] = 28;
  JB->rowvals[26] = 6;
  JB->rowvals[27] = 9;
  JB->rowvals[28] = 28;
  JB->rowvals[29] = 9;
  JB->rowvals[30] = 10;
  JB->rowvals[31] = 11;
  JB->rowvals[32] = 28;
  JB->rowvals[33] = 7;
  JB->rowvals[34] = 8;
  JB->rowvals[35] = 9;
  JB->rowvals[36] = 10;
  JB->rowvals[37] = 11;
  JB->rowvals[38] = 28;
  JB->rowvals[39] = 12;
  JB->rowvals[40] = 31;
  JB->rowvals[41] = 32;
  JB->rowvals[42] = 13;
  JB->rowvals[43] = 19;
  JB->rowvals[44] = 13;
  JB->rowvals[45] = 14;
  JB->rowvals[46] = 16;
  JB->rowvals[47] = 19;
  JB->rowvals[48] = 15;
  JB->rowvals[49] = 20;
  JB->rowvals[50] = 22;
  JB->rowvals[51] = 15;
  JB->rowvals[52] = 16;
  JB->rowvals[53] = 19;
  JB->rowvals[54] = 21;
  JB->rowvals[55] = 22;
  JB->rowvals[56] = 17;
  JB->rowvals[57] = 18;
  JB->rowvals[58] = 23;
  JB->rowvals[59] = 26;
  JB->rowvals[60] = 18;
  JB->rowvals[61] = 23;
  JB->rowvals[62] = 18;
  JB->rowvals[63] = 19;
  JB->rowvals[64] = 23;
  JB->rowvals[65] = 17;
  JB->rowvals[66] = 18;
  JB->rowvals[67] = 20;
  JB->rowvals[68] = 25;
  JB->rowvals[69] = 26;
  JB->rowvals[70] = 18;
  JB->rowvals[71] = 19;
  JB->rowvals[72] = 20;
  JB->rowvals[73] = 21;
  JB->rowvals[74] = 23;
  JB->rowvals[75] = 25;
  JB->rowvals[76] = 22;
  JB->rowvals[77] = 25;
  JB->rowvals[78] = 33;
  JB->rowvals[79] = 28;
  JB->rowvals[80] = 28;
  JB->rowvals[81] = 32;
  JB->rowvals[82] = 24;
  JB->rowvals[83] = 25;
  JB->rowvals[84] = 26;
  JB->rowvals[85] = 27;
  JB->rowvals[86] = 12;
  JB->rowvals[87] = 24;
  JB->rowvals[88] = 26;
  JB->rowvals[89] = 12;
  JB->rowvals[90] = 27;
  JB->rowvals[91] = 32;
  JB->rowvals[92] = 34;
  JB->rowvals[93] = 9;
  JB->rowvals[94] = 28;
  JB->rowvals[95] = 7;
  JB->rowvals[96] = 29;
  JB->rowvals[97] = 11;
  JB->rowvals[98] = 29;
  JB->rowvals[99] = 30;
  JB->rowvals[100] = 32;
  JB->rowvals[101] = 0;
  JB->rowvals[102] = 12;
  JB->rowvals[103] = 30;
  JB->rowvals[104] = 31;
  JB->rowvals[105] = 10;
  JB->rowvals[106] = 30;
  JB->rowvals[107] = 32;
  JB->rowvals[108] = 2;
  JB->rowvals[109] = 27;
  JB->rowvals[110] = 33;
  JB->rowvals[111] = 3;
  JB->rowvals[112] = 27;
  JB->rowvals[113] = 30;
  JB->rowvals[114] = 31;
  JB->rowvals[115] = 34;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 2;
  JB->colptrs[2] = 4;
  JB->colptrs[3] = 8;
  JB->colptrs[4] = 11;
  JB->colptrs[5] = 13;
  JB->colptrs[6] = 17;
  JB->colptrs[7] = 18;
  JB->colptrs[8] = 20;
  JB->colptrs[9] = 26;
  JB->colptrs[10] = 29;
  JB->colptrs[11] = 33;
  JB->colptrs[12] = 39;
  JB->colptrs[13] = 42;
  JB->colptrs[14] = 44;
  JB->colptrs[15] = 48;
  JB->colptrs[16] = 51;
  JB->colptrs[17] = 56;
  JB->colptrs[18] = 60;
  JB->colptrs[19] = 62;
  JB->colptrs[20] = 65;
  JB->colptrs[21] = 70;
  JB->colptrs[22] = 76;
  JB->colptrs[23] = 79;
  JB->colptrs[24] = 80;
  JB->colptrs[25] = 82;
  JB->colptrs[26] = 86;
  JB->colptrs[27] = 89;
  JB->colptrs[28] = 93;
  JB->colptrs[29] = 95;
  JB->colptrs[30] = 97;
  JB->colptrs[31] = 101;
  JB->colptrs[32] = 105;
  JB->colptrs[33] = 108;
  JB->colptrs[34] = 111;
  JB->colptrs[35] = 116;
  return(0);
}


 int sx__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1] - 1.0*x_tmp[2];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6] - 1.0*x_tmp[5];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7] + x_tmp[5];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6] - 2.0*x_tmp[8] + x_tmp[9];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6] - 1.0*x_tmp[9];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28] + x_tmp[9] + 2.0*x_tmp[11];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7] + x_tmp[8] - 1.0*x_tmp[9] - 1.0*x_tmp[11];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] + x_tmp[31] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) + x_tmp[34] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28] + x_tmp[9];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29] - 1.0*x_tmp[30];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] - 1.0*x_tmp[31] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32] + x_tmp[30];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]) + x_tmp[2];
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] - 1.0*x_tmp[34] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31];

  } break;

  case 1: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1] + x_tmp[33];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6] + x_tmp[7];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7] - 1.0*x_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6] + 2.0*x_tmp[11] + x_tmp[28];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6] + x_tmp[28];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28] - 2.0*x_tmp[10] + x_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7] + x_tmp[10] - 1.0*x_tmp[11] - 1.0*x_tmp[28];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] - 1.0*x_tmp[12] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) - 1.0*x_tmp[27] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28] - 1.0*x_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29] + x_tmp[32];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] + x_tmp[12] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32] - 1.0*x_tmp[32];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]) - 1.0*x_tmp[33];
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] + x_tmp[27] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31];

  } break;

  case 2: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] + (pow((p[3] - 1.0),2))*x_tmp[0] - 1.0*(pow(p[3],2))*x_tmp[0] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) + 0.55*x_tmp[29];
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] + 2.0*(pow(p[3],2))*x_tmp[0] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + 0.45*x_tmp[29] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] + (pow((p[3] - 1.0),2))*x_tmp[12] - 1.0*(pow(p[3],2))*x_tmp[12] + 0.55*x_tmp[32] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] + (pow((p[3] - 1.0),2))*x_tmp[15] - 1.0*(pow(p[3],2))*x_tmp[15] + 0.55*x_tmp[22] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) + 2.0*(pow(p[3],2))*x_tmp[15] + 0.45*x_tmp[22] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] + 2.0*(pow((p[3] - 1.0),2))*x_tmp[17] + (pow((p[3] - 1.0),2))*x_tmp[18] - 2.0*(pow(p[3],2))*x_tmp[17] + (pow(p[3],2))*x_tmp[18] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + 0.55*x_tmp[23] + 1.1*x_tmp[26] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] + (pow((p[3] - 1.0),2))*x_tmp[18] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(pow(p[3],2))*x_tmp[18] + 0.55*x_tmp[23];
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + 2.0*(pow(p[3],2))*x_tmp[18] + 0.45*x_tmp[23] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (pow((p[3] - 1.0),2))*x_tmp[20] + 2.0*(pow(p[3],2))*x_tmp[17] - 2.0*(pow(p[3],2))*x_tmp[18] - 1.0*(pow(p[3],2))*x_tmp[20] + 0.55*x_tmp[25] + 0.45*x_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[17];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + 4.0*(pow(p[3],2))*x_tmp[18] + 4.0*(pow(p[3],2))*x_tmp[20] + 0.45*x_tmp[23] + 0.9*x_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[3]*x_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33] + x_tmp[33];
sxdot_tmp[23] = p[2]*sx_tmp[28] + x_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32] + x_tmp[28] + 2.0*x_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + 2.0*(pow(p[3],2))*x_tmp[26] + 0.45*x_tmp[24] + x_tmp[27] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] + (pow((p[3] - 1.0),2))*x_tmp[26] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 1.0*(pow(p[3],2))*x_tmp[26] + x_tmp[12] + 0.55*x_tmp[24];
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) + 2.0*(pow(p[3],2))*x_tmp[12] + 0.45*x_tmp[32] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29] + x_tmp[7];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29] + x_tmp[11];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + (pow((p[3] - 1.0),2))*x_tmp[31] - 1.0*(pow(p[3],2))*x_tmp[31] + 0.0002*sx_tmp[0] + 0.55*x_tmp[30] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32] + x_tmp[10];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]);
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 2.0*(pow(p[3],2))*x_tmp[31] + 0.0002*sx_tmp[3] + 0.45*x_tmp[30] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31] - 2.0*(p[3] - 1.0)*p[3]*x_tmp[31];

  } break;

  case 3: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002) + (2.0*p[3] - 2.0)*p[2]*x_tmp[0] - 2.0*p[2]*p[3]*x_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[0] + 2.0*p[2]*p[3]*x_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[12] - 2.0*p[2]*p[3]*x_tmp[12];
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[15] - 2.0*p[2]*p[3]*x_tmp[15];
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[15] + 2.0*p[2]*p[3]*x_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18] + 2.0*(2.0*p[3] - 2.0)*p[2]*x_tmp[17] + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 4.0*p[2]*p[3]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[18];
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[18] + 2.0*p[2]*p[3]*x_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[20] - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[17] + 2.0*p[2]*p[3]*x_tmp[17] - 4.0*p[2]*p[3]*x_tmp[18] - 2.0*p[2]*p[3]*x_tmp[20];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[18] - 4.0*(p[3] - 1.0)*p[2]*x_tmp[20] + 6.0*p[2]*p[3]*x_tmp[18] + 4.0*p[2]*p[3]*x_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[26] + 2.0*p[2]*p[3]*x_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[26] - 2.0*p[2]*p[3]*x_tmp[26];
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[12] + 2.0*p[2]*p[3]*x_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]) + (2.0*p[3] - 2.0)*p[2]*x_tmp[31] - 2.0*p[2]*p[3]*x_tmp[31];
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]);
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31] - 2.0*(p[3] - 1.0)*p[2]*x_tmp[31] + 2.0*p[2]*p[3]*x_tmp[31];

  } break;

  case 4: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1] + x_tmp[3];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1] + x_tmp[34];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] - 1.0*x_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13] + x_tmp[19];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19] + 2.0*x_tmp[16] + x_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] + x_tmp[20] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) - 1.0*x_tmp[16] - 1.0*x_tmp[19] + x_tmp[21] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] - 1.0*x_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] - 1.0*x_tmp[20] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + x_tmp[19] - 2.0*x_tmp[21] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33] + x_tmp[25];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] - 1.0*x_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) - 1.0*x_tmp[27] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]) + x_tmp[27];
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] - 1.0*x_tmp[34] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31];

  } break;

  case 5: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1] - 1.0*x_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1] - 1.0*x_tmp[2];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13] - 1.0*x_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19] + x_tmp[13] - 2.0*x_tmp[14];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] - 1.0*x_tmp[15] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) - 1.0*x_tmp[16] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33] - 1.0*x_tmp[22];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]) - 1.0*x_tmp[33];
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31];

  } break;

  case 6: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]);
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31];

  } break;

  case 7: {
sxdot_tmp[0] = 0.55*p[2]*sx_tmp[29] - 1.0*sx_tmp[0]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2] + 0.0002);
sxdot_tmp[1] = p[4]*sx_tmp[3] - 1.0*(p[5] + 0.0002)*sx_tmp[1];
sxdot_tmp[2] = p[1]*sx_tmp[33] + p[4]*sx_tmp[34] - 1.0*sx_tmp[2]*(p[0] + p[5]) + 0.0002*sx_tmp[1];
sxdot_tmp[3] = 0.45*p[2]*sx_tmp[29] - 1.0*(p[4] + 0.0002)*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[0];
sxdot_tmp[4] = 0.0002*sx_tmp[6] - 0.0004*sx_tmp[4];
sxdot_tmp[5] = p[1]*sx_tmp[7] - 1.0*(p[0] + 0.0002)*sx_tmp[5] + 0.0002*sx_tmp[4] - 0.0002*sx_tmp[6];
sxdot_tmp[6] = -0.0002*sx_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[5] - 1.0*(p[1] + 0.0002)*sx_tmp[7];
sxdot_tmp[8] = p[0]*sx_tmp[9] - 2.0*p[0]*sx_tmp[8] + 2.0*p[1]*sx_tmp[11] + p[1]*sx_tmp[28] + 0.0004*sx_tmp[5] + 0.0002*sx_tmp[6];
sxdot_tmp[9] = p[1]*sx_tmp[28] - 1.0*p[0]*sx_tmp[9] + 0.0002*sx_tmp[6];
sxdot_tmp[10] = p[0]*sx_tmp[9] + 2.0*p[0]*sx_tmp[11] - 2.0*p[1]*sx_tmp[10] + p[1]*sx_tmp[28];
sxdot_tmp[11] = p[0]*sx_tmp[8] - 1.0*p[0]*sx_tmp[9] + p[1]*sx_tmp[10] - 1.0*p[1]*sx_tmp[28] - 1.0*sx_tmp[11]*(p[0] + p[1]) + 0.0002*sx_tmp[7];
sxdot_tmp[12] = p[0]*sx_tmp[31] + 0.55*p[2]*sx_tmp[32] - 1.0*sx_tmp[12]*(p[2]*(pow(p[3],2)) + p[1] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[13] = p[4]*sx_tmp[19] - 1.0*p[5]*sx_tmp[13];
sxdot_tmp[14] = p[5]*sx_tmp[13] - 2.0*p[5]*sx_tmp[14] + 2.0*p[4]*sx_tmp[16] + p[4]*sx_tmp[19];
sxdot_tmp[15] = 0.55*p[2]*sx_tmp[22] + p[4]*sx_tmp[20] - 1.0*sx_tmp[15]*(p[2]*(pow(p[3],2)) + p[5] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[16] = 0.45*p[2]*sx_tmp[22] - 1.0*p[4]*sx_tmp[19] + p[4]*sx_tmp[21] - 1.0*sx_tmp[16]*(p[4] + p[5]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[15];
sxdot_tmp[17] = 0.55*p[2]*sx_tmp[23] + 1.1*p[2]*sx_tmp[26] - 1.0*(2.0*p[2]*(pow(p[3],2)) - 2.0*(pow((p[3] - 1.0),2))*p[2])*sx_tmp[17] + (p[2]*(pow(p[3],2)) + (pow((p[3] - 1.0),2))*p[2])*sx_tmp[18];
sxdot_tmp[18] = 0.55*p[2]*sx_tmp[23] - 1.0*sx_tmp[18]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[19] = 0.45*p[2]*sx_tmp[23] - 1.0*p[4]*sx_tmp[19] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18];
sxdot_tmp[20] = 0.55*p[2]*sx_tmp[25] + 0.45*p[2]*sx_tmp[26] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[17] - 1.0*sx_tmp[20]*(p[2]*(pow(p[3],2)) + p[4] - 1.0*(pow((p[3] - 1.0),2))*p[2]) - 2.0*p[2]*(pow(p[3],2))*sx_tmp[18];
sxdot_tmp[21] = p[4]*sx_tmp[19] + 0.45*p[2]*sx_tmp[23] - 2.0*p[4]*sx_tmp[21] + 0.9*p[2]*sx_tmp[25] + (4.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[18] + (4.0*p[2]*(pow(p[3],2)) - 4.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[20];
sxdot_tmp[22] = p[4]*sx_tmp[25] - 1.0*p[5]*sx_tmp[22] + p[2]*sx_tmp[33];
sxdot_tmp[23] = p[2]*sx_tmp[28];
sxdot_tmp[24] = p[2]*sx_tmp[28] + 2.0*p[2]*sx_tmp[32];
sxdot_tmp[25] = 0.45*p[2]*sx_tmp[24] + p[2]*sx_tmp[27] - 1.0*p[4]*sx_tmp[25] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[26];
sxdot_tmp[26] = p[2]*sx_tmp[12] + 0.55*p[2]*sx_tmp[24] - 1.0*sx_tmp[26]*(p[2]*(pow(p[3],2)) - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[27] = p[0]*sx_tmp[34] + 0.45*p[2]*sx_tmp[32] - 1.0*sx_tmp[27]*(p[1] + p[4]) + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[12];
sxdot_tmp[28] = p[0]*sx_tmp[9] - 1.0*p[1]*sx_tmp[28];
sxdot_tmp[29] = p[2]*sx_tmp[7] - 0.0002*sx_tmp[29];
sxdot_tmp[30] = p[2]*sx_tmp[11] - 1.0*p[0]*sx_tmp[30] + p[1]*sx_tmp[32] + 0.0002*sx_tmp[29];
sxdot_tmp[31] = p[1]*sx_tmp[12] + 0.55*p[2]*sx_tmp[30] + 0.0002*sx_tmp[0] - 1.0*sx_tmp[31]*(p[2]*(pow(p[3],2)) + p[0] - 1.0*(pow((p[3] - 1.0),2))*p[2]);
sxdot_tmp[32] = p[2]*sx_tmp[10] + p[0]*sx_tmp[30] - 1.0*p[1]*sx_tmp[32];
sxdot_tmp[33] = p[0]*sx_tmp[2] + p[4]*sx_tmp[27] - 1.0*sx_tmp[33]*(p[1] + p[5]);
sxdot_tmp[34] = p[1]*sx_tmp[27] + 0.45*p[2]*sx_tmp[30] - 1.0*sx_tmp[34]*(p[0] + p[4]) + 0.0002*sx_tmp[3] + (2.0*p[2]*(pow(p[3],2)) - 2.0*(p[3] - 1.0)*p[2]*p[3])*sx_tmp[31];

  } break;

  }
 for (ix=0; ix<35; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(int ip, N_Vector sx0, void *user_data)
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
sx0_tmp[4] = (2.0*p[6] - 1.0)*(k[7] - 1.0);
sx0_tmp[5] = (k[8] - 1.0)*p[7];
sx0_tmp[6] = 1.0 - 1.0*k[0];
sx0_tmp[7] = - 1.0*(k[9] - 1.0)*(p[6] + p[7] - 1.0) - 1.0*(k[9] - 1.0)*p[6];
sx0_tmp[10] = (k[20] - 1.0)*(p[6] + p[7] - 1.0) + (k[20] - 1.0)*(p[6] + p[7]);
sx0_tmp[11] = -1.0*(k[15] - 1.0)*p[7];
sx0_tmp[28] = k[2] - 1.0;

  } break;

  case 7: {
sx0_tmp[5] = (k[8] - 1.0)*p[6];
sx0_tmp[7] = -1.0*(k[9] - 1.0)*p[6];
sx0_tmp[8] = (2.0*p[7] - 1.0)*(k[14] - 1.0);
sx0_tmp[9] = 1.0 - 1.0*k[1];
sx0_tmp[10] = (k[20] - 1.0)*(p[6] + p[7] - 1.0) + (k[20] - 1.0)*(p[6] + p[7]);
sx0_tmp[11] = - 1.0*(k[15] - 1.0)*(p[6] + p[7] - 1.0) - 1.0*(k[15] - 1.0)*p[7];
sx0_tmp[28] = k[2] - 1.0;

  } break;

  }

  return;
}


void y__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*23];
y[it+nt*1] = x[it+nt*18] + x[it+nt*19];
y[it+nt*2] = x[it+nt*13];
y[it+nt*3] = x[it+nt*24];
y[it+nt*4] = x[it+nt*25] + x[it+nt*26];
y[it+nt*5] = x[it+nt*22];
y[it+nt*6] = x[it+nt*17] + 2.0*x[it+nt*20] + x[it+nt*21];
y[it+nt*7] = x[it+nt*15] + x[it+nt*16];
y[it+nt*8] = x[it+nt*14];
    
    return;
}


void dydp__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*315);
dydx[119] = 1.0;
dydx[134] = 1.0;
dydx[142] = 1.0;
dydx[151] = 1.0;
dydx[159] = 1.0;
dydx[163] = 1.0;
dydx[172] = 1.0;
dydx[186] = 2.0;
dydx[195] = 1.0;
dydx[203] = 1.0;
dydx[207] = 1.0;
dydx[219] = 1.0;
dydx[229] = 1.0;
dydx[238] = 1.0;
  
  return;
}


void sy__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(23+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(18+np*nx)] + sx[it+nt*(19+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(24+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(25+np*nx)] + sx[it+nt*(26+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(22+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(17+np*nx)] + 2.0*sx[it+nt*(20+np*nx)] + sx[it+nt*(21+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(15+np*nx)] + sx[it+nt*(16+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];

  } break;

  }
  
  return;
}
int root__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
void deltadisc__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int idisc, N_Vector x, void *user_data){
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
void sdeltadisc__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
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


void dxdotdp__sA_sPB_3o_2B_S_atS__T_atS__B_dtS(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(2+ip*nx)] = -1.0*x[it+nt*2];
dxdotdp[(5+ip*nx)] = -1.0*x[it+nt*5];
dxdotdp[(7+ip*nx)] = x[it+nt*5];
dxdotdp[(8+ip*nx)] = x[it+nt*9] - 2.0*x[it+nt*8];
dxdotdp[(9+ip*nx)] = -1.0*x[it+nt*9];
dxdotdp[(10+ip*nx)] = x[it+nt*9] + 2.0*x[it+nt*11];
dxdotdp[(11+ip*nx)] = x[it+nt*8] - 1.0*x[it+nt*9] - 1.0*x[it+nt*11];
dxdotdp[(12+ip*nx)] = x[it+nt*31];
dxdotdp[(27+ip*nx)] = x[it+nt*34];
dxdotdp[(28+ip*nx)] = x[it+nt*9];
dxdotdp[(30+ip*nx)] = -1.0*x[it+nt*30];
dxdotdp[(31+ip*nx)] = -1.0*x[it+nt*31];
dxdotdp[(32+ip*nx)] = x[it+nt*30];
dxdotdp[(33+ip*nx)] = x[it+nt*2];
dxdotdp[(34+ip*nx)] = -1.0*x[it+nt*34];

  } break;

  case 1: {
dxdotdp[(2+ip*nx)] = x[it+nt*33];
dxdotdp[(5+ip*nx)] = x[it+nt*7];
dxdotdp[(7+ip*nx)] = -1.0*x[it+nt*7];
dxdotdp[(8+ip*nx)] = 2.0*x[it+nt*11] + x[it+nt*28];
dxdotdp[(9+ip*nx)] = x[it+nt*28];
dxdotdp[(10+ip*nx)] = x[it+nt*28] - 2.0*x[it+nt*10];
dxdotdp[(11+ip*nx)] = x[it+nt*10] - 1.0*x[it+nt*11] - 1.0*x[it+nt*28];
dxdotdp[(12+ip*nx)] = -1.0*x[it+nt*12];
dxdotdp[(27+ip*nx)] = -1.0*x[it+nt*27];
dxdotdp[(28+ip*nx)] = -1.0*x[it+nt*28];
dxdotdp[(30+ip*nx)] = x[it+nt*32];
dxdotdp[(31+ip*nx)] = x[it+nt*12];
dxdotdp[(32+ip*nx)] = -1.0*x[it+nt*32];
dxdotdp[(33+ip*nx)] = -1.0*x[it+nt*33];
dxdotdp[(34+ip*nx)] = x[it+nt*27];

  } break;

  case 2: {
dxdotdp[(0+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*0] - 1.0*(pow(p[3],2))*x[it+nt*0] + 0.55*x[it+nt*29];
dxdotdp[(3+ip*nx)] = 2.0*(pow(p[3],2))*x[it+nt*0] + 0.45*x[it+nt*29] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*0];
dxdotdp[(12+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*12] - 1.0*(pow(p[3],2))*x[it+nt*12] + 0.55*x[it+nt*32];
dxdotdp[(15+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*15] - 1.0*(pow(p[3],2))*x[it+nt*15] + 0.55*x[it+nt*22];
dxdotdp[(16+ip*nx)] = 2.0*(pow(p[3],2))*x[it+nt*15] + 0.45*x[it+nt*22] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*15];
dxdotdp[(17+ip*nx)] = 2.0*(pow((p[3] - 1.0),2))*x[it+nt*17] + (pow((p[3] - 1.0),2))*x[it+nt*18] - 2.0*(pow(p[3],2))*x[it+nt*17] + (pow(p[3],2))*x[it+nt*18] + 0.55*x[it+nt*23] + 1.1*x[it+nt*26];
dxdotdp[(18+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*18] - 1.0*(pow(p[3],2))*x[it+nt*18] + 0.55*x[it+nt*23];
dxdotdp[(19+ip*nx)] = 2.0*(pow(p[3],2))*x[it+nt*18] + 0.45*x[it+nt*23] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*18];
dxdotdp[(20+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*20] + 2.0*(pow(p[3],2))*x[it+nt*17] - 2.0*(pow(p[3],2))*x[it+nt*18] - 1.0*(pow(p[3],2))*x[it+nt*20] + 0.55*x[it+nt*25] + 0.45*x[it+nt*26] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*17];
dxdotdp[(21+ip*nx)] = 4.0*(pow(p[3],2))*x[it+nt*18] + 4.0*(pow(p[3],2))*x[it+nt*20] + 0.45*x[it+nt*23] + 0.9*x[it+nt*25] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*18] - 4.0*(p[3] - 1.0)*p[3]*x[it+nt*20];
dxdotdp[(22+ip*nx)] = x[it+nt*33];
dxdotdp[(23+ip*nx)] = x[it+nt*28];
dxdotdp[(24+ip*nx)] = x[it+nt*28] + 2.0*x[it+nt*32];
dxdotdp[(25+ip*nx)] = 2.0*(pow(p[3],2))*x[it+nt*26] + 0.45*x[it+nt*24] + x[it+nt*27] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*26];
dxdotdp[(26+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*26] - 1.0*(pow(p[3],2))*x[it+nt*26] + x[it+nt*12] + 0.55*x[it+nt*24];
dxdotdp[(27+ip*nx)] = 2.0*(pow(p[3],2))*x[it+nt*12] + 0.45*x[it+nt*32] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*12];
dxdotdp[(29+ip*nx)] = x[it+nt*7];
dxdotdp[(30+ip*nx)] = x[it+nt*11];
dxdotdp[(31+ip*nx)] = (pow((p[3] - 1.0),2))*x[it+nt*31] - 1.0*(pow(p[3],2))*x[it+nt*31] + 0.55*x[it+nt*30];
dxdotdp[(32+ip*nx)] = x[it+nt*10];
dxdotdp[(34+ip*nx)] = 2.0*(pow(p[3],2))*x[it+nt*31] + 0.45*x[it+nt*30] - 2.0*(p[3] - 1.0)*p[3]*x[it+nt*31];

  } break;

  case 3: {
dxdotdp[(0+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*0] - 2.0*p[2]*p[3]*x[it+nt*0];
dxdotdp[(3+ip*nx)] = 2.0*p[2]*p[3]*x[it+nt*0] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*0];
dxdotdp[(12+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*12] - 2.0*p[2]*p[3]*x[it+nt*12];
dxdotdp[(15+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*15] - 2.0*p[2]*p[3]*x[it+nt*15];
dxdotdp[(16+ip*nx)] = 2.0*p[2]*p[3]*x[it+nt*15] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*15];
dxdotdp[(17+ip*nx)] = 2.0*(2.0*p[3] - 2.0)*p[2]*x[it+nt*17] + (2.0*p[3] - 2.0)*p[2]*x[it+nt*18] - 4.0*p[2]*p[3]*x[it+nt*17] + 2.0*p[2]*p[3]*x[it+nt*18];
dxdotdp[(18+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*18] - 2.0*p[2]*p[3]*x[it+nt*18];
dxdotdp[(19+ip*nx)] = 2.0*p[2]*p[3]*x[it+nt*18] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*18];
dxdotdp[(20+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*20] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*17] + 2.0*p[2]*p[3]*x[it+nt*17] - 4.0*p[2]*p[3]*x[it+nt*18] - 2.0*p[2]*p[3]*x[it+nt*20];
dxdotdp[(21+ip*nx)] = 6.0*p[2]*p[3]*x[it+nt*18] - 4.0*(p[3] - 1.0)*p[2]*x[it+nt*20] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*18] + 4.0*p[2]*p[3]*x[it+nt*20];
dxdotdp[(25+ip*nx)] = 2.0*p[2]*p[3]*x[it+nt*26] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*26];
dxdotdp[(26+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*26] - 2.0*p[2]*p[3]*x[it+nt*26];
dxdotdp[(27+ip*nx)] = 2.0*p[2]*p[3]*x[it+nt*12] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*12];
dxdotdp[(31+ip*nx)] = (2.0*p[3] - 2.0)*p[2]*x[it+nt*31] - 2.0*p[2]*p[3]*x[it+nt*31];
dxdotdp[(34+ip*nx)] = 2.0*p[2]*p[3]*x[it+nt*31] - 2.0*(p[3] - 1.0)*p[2]*x[it+nt*31];

  } break;

  case 4: {
dxdotdp[(1+ip*nx)] = x[it+nt*3];
dxdotdp[(2+ip*nx)] = x[it+nt*34];
dxdotdp[(3+ip*nx)] = -1.0*x[it+nt*3];
dxdotdp[(13+ip*nx)] = x[it+nt*19];
dxdotdp[(14+ip*nx)] = 2.0*x[it+nt*16] + x[it+nt*19];
dxdotdp[(15+ip*nx)] = x[it+nt*20];
dxdotdp[(16+ip*nx)] = x[it+nt*21] - 1.0*x[it+nt*19] - 1.0*x[it+nt*16];
dxdotdp[(19+ip*nx)] = -1.0*x[it+nt*19];
dxdotdp[(20+ip*nx)] = -1.0*x[it+nt*20];
dxdotdp[(21+ip*nx)] = x[it+nt*19] - 2.0*x[it+nt*21];
dxdotdp[(22+ip*nx)] = x[it+nt*25];
dxdotdp[(25+ip*nx)] = -1.0*x[it+nt*25];
dxdotdp[(27+ip*nx)] = -1.0*x[it+nt*27];
dxdotdp[(33+ip*nx)] = x[it+nt*27];
dxdotdp[(34+ip*nx)] = -1.0*x[it+nt*34];

  } break;

  case 5: {
dxdotdp[(1+ip*nx)] = -1.0*x[it+nt*1];
dxdotdp[(2+ip*nx)] = -1.0*x[it+nt*2];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];
dxdotdp[(14+ip*nx)] = x[it+nt*13] - 2.0*x[it+nt*14];
dxdotdp[(15+ip*nx)] = -1.0*x[it+nt*15];
dxdotdp[(16+ip*nx)] = -1.0*x[it+nt*16];
dxdotdp[(22+ip*nx)] = -1.0*x[it+nt*22];
dxdotdp[(33+ip*nx)] = -1.0*x[it+nt*33];

  } break;

  }
  }
  
  return;
}
