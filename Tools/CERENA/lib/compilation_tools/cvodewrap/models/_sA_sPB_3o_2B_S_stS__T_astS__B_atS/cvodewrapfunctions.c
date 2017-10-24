#include "_sA_sPB_3o_2B_S_stS__T_astS__B_atS.c"
                
                int cvodewrap_init(void *cvode_mem, N_Vector x, double t){
                    return CVodeInit(cvode_mem, xdot__sA_sPB_3o_2B_S_stS__T_astS__B_atS, RCONST(t), x);
                }
                int cvodewrap_binit(void *cvode_mem, int which, N_Vector xB, double t){
                    return CVodeInitB(cvode_mem, which, xBdot__sA_sPB_3o_2B_S_stS__T_astS__B_atS, RCONST(t), xB);
                }
                int cvodewrap_qbinit(void *cvode_mem, int which, N_Vector xQB){
                    return CVodeQuadInitB(cvode_mem, which, xQB__sA_sPB_3o_2B_S_stS__T_astS__B_atS, xQB);
                }
                
                void fx0(N_Vector x0, void *user_data){
                    UserData data = (UserData) user_data;
                    x0__sA_sPB_3o_2B_S_stS__T_astS__B_atS(x0, data);
                }
                
                int cvodewrap_SetDenseJacFn(void *cvode_mem){
                    return CVDlsSetDenseJacFn(cvode_mem, J__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetSparseJacFn(void *cvode_mem){
                    return CVSlsSetSparseJacFn(cvode_mem, JSparse__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetBandJacFn(void *cvode_mem){
                    return CVDlsSetBandJacFn(cvode_mem, JBand__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetJacTimesVecFn(void *cvode_mem){
                    return CVSpilsSetJacTimesVecFn(cvode_mem, Jv__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetDenseJacFnB(void *cvode_mem,int which){
                    return CVDlsSetDenseJacFnB(cvode_mem, which, JB__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetSparseJacFnB(void *cvode_mem, int which){
                    return CVSlsSetSparseJacFnB(cvode_mem, which, JSparseB__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetBandJacFnB(void *cvode_mem,int which){
                    return CVDlsSetBandJacFnB(cvode_mem, which, JBBand__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int cvodewrap_SetJacTimesVecFnB(void *cvode_mem,int which){
                    return CVSpilsSetJacTimesVecFnB(cvode_mem, which, JvB__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                int fJ(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
                    return J__sA_sPB_3o_2B_S_stS__T_astS__B_atS(N,t,x,fx,J,user_data,tmp1,tmp2,tmp3);
                }
                int fJB(long int N, realtype t, N_Vector x, N_Vector xB, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
                    return JB__sA_sPB_3o_2B_S_stS__T_astS__B_atS(N,t,x,xB,fx,J,user_data,tmp1,tmp2,tmp3);
                }
                
                void fsx0(N_Vector *sx0, void *user_data){
                    UserData data = (UserData) user_data;
                    int ip;
                    int *plist = data->plist;
                    int np = *data->np;
                    for (ip=0; ip<np; ip++) {
                       sx0__sA_sPB_3o_2B_S_stS__T_astS__B_atS(plist[ip], sx0[ip], data);
                    }
                }
                
                int cvodewrap_SensInit1(void *cvode_mem, int np, int sensi_meth, N_Vector *sx){
                    return CVodeSensInit1(cvode_mem, np, sensi_meth, sx__sA_sPB_3o_2B_S_stS__T_astS__B_atS, sx);
                }
                int cvodewrap_RootInit(void *cvode_mem, int nr){
                    return CVodeRootInit(cvode_mem, nr, root__sA_sPB_3o_2B_S_stS__T_astS__B_atS);
                }
                void froot(double t, N_Vector x, realtype *gout, void *user_data){
                    root__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, x, gout, user_data);
                }
                
                void fy(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
                    y__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, nt, it, y, p, k, u, x);
                }
                
                void fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data){
                    xdot__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t,x,xdot,user_data);
                }
                void fdydp(double t, int nt, int it,double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
                    dydp__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, nt, it, dydp, y, p, k, u, x, plist, np, ny);
                }
                void fdxdotdp(double t, int nt, int it,double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
                    dxdotdp__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, nt, it, dxdotdp, p, k, u, x, plist, np, nx);
                }
                void fdydx(double t,double *dydx, double *y, double *p, double *k, double *x){
                    dydx__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, dydx, y, p, k, x);
                }
                
                
                void fsy(double t, int nt, int it, int *plist, int nx, int ny, int np, double *sy, double *p, double *k, double *x, double *sx){
                    int ip;
                    for (ip=0; ip<np; ip++) {
                        sy__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, nt, it, plist[ip], ip,  nx, ny,  sy, p, k, x, sx);
                    }
                }
                double fsroot(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(sroot__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, ip, ir, x, sx, user_data));
                }
                double fsrootval(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(srootval__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, ip, ir, x, sx, user_data));
                }
                double fs2root(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(s2root__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, ip, jp, ir, x, sx, user_data));
                }
                double fs2rootval(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(s2rootval__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, ip, jp, ir, x, sx, user_data));
                }
                void deltadisc(double t, int idisc, N_Vector x, void *user_data){
                       deltadisc__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, idisc, x, user_data);
                }
                void sdeltadisc(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
                    UserData data = (UserData) user_data;
                    sdeltadisc__sA_sPB_3o_2B_S_stS__T_astS__B_atS(t, idisc, x, sx, data);
                }
                void cvodewrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data){
                    char buffer [250];
                    sprintf(buffer,"CVODES ERROR: in module %s in function %s : %s ",module,function,msg);
                    mexWarnMsgTxt(buffer);
                }

