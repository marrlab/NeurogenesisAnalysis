% simulate__sA_sPB_3o_2B_S_stS__T_dtS__B_astS.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%
% USAGE:
% ======
% [...] = simulate__sA_sPB_3o_2B_S_stS__T_dtS__B_astS(tout,theta)
% [...] = simulate__sA_sPB_3o_2B_S_stS__T_dtS__B_astS(tout,theta,kappa,options)
% sol = simulate__sA_sPB_3o_2B_S_stS__T_dtS__B_astS(...)
% [status,tout,x,y,sx,sy] = simulate__sA_sPB_3o_2B_S_stS__T_dtS__B_astS(...)
%
% INPUTS:
% =======
% tout ... 1 dimensional vector of timepoints at which a solution to the ODE is desired
% theta ... 1 dimensional parameter vector of parameters for which sensitivities are desired.
%           this corresponds to the specification in model.sym.p
% kappa ... 1 dimensional parameter vector of parameters for which sensitivities are not desired.
%           this corresponds to the specification in model.sym.k
%           Arbitrary initial conditions can be provided in kappa (see ACME/doc/ACME_doc.pdf for detailed instructions).
% options ... additional options to pass to the cvodes solver. Refer to the cvodes guide for more documentation.
%    .cvodes_atol ... absolute tolerance for the solver. default is specified in the user-provided syms function.
%    .cvodes_rtol ... relative tolerance for the solver. default is specified in the user-provided syms function.
%    .cvodes_maxsteps    ... maximal number of integration steps. default is specified in the user-provided syms function.
%    .tstart    ... start of integration. for all timepoints before this, values will be set to initial value.
%    .sens_ind ... 1 dimensional vector of indexes for which sensitivities must be computed.
%           default value is 1:length(theta).
%    .sx0 ... user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters].
%        default is sensitivity initialisation based on the derivative of the state initialisation.%    .lmm    ... linear multistep method for forward problem.
%        1: Adams-Bashford
%        2: BDF (DEFAULT)
%    .iter    ... iteration method for linear multistep.
%        1: Functional
%        2: Newton (DEFAULT)
%    .linsol   ... linear solver module.
%        direct solvers:
%        1: Dense (DEFAULT)
%        2: Band (not implented)
%        3: LAPACK Dense (not implented)
%        4: LAPACK Band  (not implented)
%        5: Diag (not implented)
%        implicit krylov solvers:
%        6: SPGMR
%        7: SPBCG
%        8: SPTFQMR
%        sparse solvers:
%        9: KLU
%    .stldet   ... flag for stability limit detection. this should be turned on for stiff problems.
%        0: OFF
%        1: ON (DEFAULT)
%    .qPositiveX   ... vector of 0 or 1 of same dimension as state vector. 1 enforces positivity of states.
%    .sensi_meth   ... method for sensitivity computation.
%        1: Forward Sensitivity Analysis (DEFAULT)
%        2: Adjoint Sensitivity Analysis
%    .ism   ... only available for sensi_meth == 1. Method for computation of forward sensitivities.
%        1: Simultaneous (DEFAULT)
%        2: Staggered
%        3: Staggered1
%    .Nd   ... only available for sensi_meth == 2. Number of Interpolation nodes for forward solution. 
%              Default is 1000. 
%    .interpType   ... only available for sensi_meth == 2. Interpolation method for forward solution.
%        1: Hermite (DEFAULT)
%        2: Polynomial
%    .lmmB   ... only available for sensi_meth == 2. linear multistep method for backward problem.
%        1: Adams-Bashford
%        2: BDF (DEFAULT)
%    .iterB   ... only available for sensi_meth == 2. iteration method for linear multistep.
%        1: Functional
%        2: Newton (DEFAULT)
%
% Outputs:
% ========
% sol.status ... flag for status of integration. generally status<0 for failed integration
% sol.tout ... vector at which the solution was computed
% sol.x ... time-resolved state vector
% sol.y ... time-resolved output vector
% sol.sx ... time-resolved state sensitivity vector
% sol.sy ... time-resolved output sensitivity vector
% sol.xdot time-resolved right-hand side of differential equation
% sol.rootval value of root at end of simulation time
% sol.srootval value of root at end of simulation time
% sol.root time of events
% sol.sroot value of root at end of simulation time
function varargout = simulate__sA_sPB_3o_2B_S_stS__T_dtS__B_astS(varargin)

% DO NOT CHANGE ANYTHING IN THIS FILE UNLESS YOU ARE VERY SURE ABOUT WHAT YOU ARE DOING
% MANUAL CHANGES TO THIS FILE CAN RESULT IN WRONG SOLUTIONS AND CRASHING OF MATLAB
if(nargin<2)
    error('Not enough input arguments.');
else
    tout=varargin{1};
    phi=varargin{2};
end
if(nargin>=3)
    kappa=varargin{3};
   if(length(kappa)==0)
    kappa(1:70)=0;
   end
else
    kappa=zeros(1,70);
end
if(nargout>1)
    if(nargout>4)
        options_cvode.sensi = 1;
    else
        options_cvode.sensi = 0;
    end
else
    options_cvode.sensi = 1;
end
theta = phi;


if(length(theta)<11)
    error('provided parameter vector is too short');
end
if(length(kappa)<70)
    error('provided constant vector is too short');
end

options_cvode.cvodes_atol = 1e-08;
options_cvode.cvodes_rtol = 1e-08;
options_cvode.cvodes_maxsteps = 10000;
options_cvode.sens_ind = 1:11;
options_cvode.nx = 35; % MUST NOT CHANGE THIS VALUE
options_cvode.ny = 9; % MUST NOT CHANGE THIS VALUE
options_cvode.nr = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.ndisc = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.nnz = 123; % MUST NOT CHANGE THIS VALUE
options_cvode.tstart = 0;
options_cvode.lmm = 2;
options_cvode.iter = 2;
options_cvode.linsol = 9;
options_cvode.stldet = 1;
options_cvode.Nd = 1000;
options_cvode.interpType = 1;
options_cvode.lmmB = 2;
options_cvode.iterB = 2;
options_cvode.ism = 1;
options_cvode.sensi_meth = 1;

options_cvode.nmaxroot = 100;

options_cvode.ubw = 34;

options_cvode.lbw = 31;

options_cvode.ss = 0;
options_cvode.qPositiveX = zeros(length(tout),35);

sol.status = 0;
sol.t = tout;
sol.x = zeros(length(tout),35);
sol.y = zeros(length(tout),9);
sol.xdot = zeros(length(tout),35);
sol.J = zeros(35,35);
sol.dydx = zeros(9,35);
sol.dydp = zeros(9,11);
sol.root = NaN(options_cvode.nmaxroot,0);
sol.rootval = NaN(options_cvode.nmaxroot,0);
sol.numsteps = zeros(length(tout),1);
sol.numrhsevals = zeros(length(tout),1);
sol.numlinsolvsetups = zeros(length(tout),1);
sol.numerrtestfails = zeros(length(tout),1);
sol.order = zeros(length(tout),1);
sol.numnonlinsolviters = zeros(length(tout),1);
sol.numjacevals = zeros(length(tout),1);
sol.numliniters = zeros(length(tout),1);
sol.numconvfails = zeros(length(tout),1);
sol.numprecevals = zeros(length(tout),1);
sol.numprecsolves = zeros(length(tout),1);

sol.numjtimesevals = zeros(length(tout),1);
sol.numstepsS = zeros(length(tout),1);
sol.numrhsevalsS = zeros(length(tout),1);
sol.numlinsolvsetupsS = zeros(length(tout),1);
sol.numerrtestfailsS = zeros(length(tout),1);
sol.orderS = zeros(length(tout),1);
sol.numnonlinsolvitersS = zeros(length(tout),1);
sol.numjacevalsS = zeros(length(tout),1);
sol.numlinitersS = zeros(length(tout),1);
sol.numconvfailsS = zeros(length(tout),1);
sol.numprecevalsS = zeros(length(tout),1);
sol.numprecsolvesS = zeros(length(tout),1);
sol.numjtimesevalsS = zeros(length(tout),1);

pbar = ones(size(theta));
pbar(pbar==0) = 1;
xscale = [];
if(nargin>=4)
    options_cvode = cw_setdefault(varargin{4},options_cvode);
else
end
if(options_cvode.ss>0)
    options_cvode.sensi = 0;
end
options_cvode.np = length(options_cvode.sens_ind); % MUST NOT CHANGE THIS VALUE
sol.dxdotdp = zeros(35,options_cvode.np);
plist = options_cvode.sens_ind-1;
if(options_cvode.sensi>0)
    sol.xS = zeros(length(tout),35,length(options_cvode.sens_ind));
    sol.yS = zeros(length(tout),9,length(options_cvode.sens_ind));
    sol.rootS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
    sol.rootvalS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
end
if(max(options_cvode.sens_ind)>11)
    error('Sensitivity index exceeds parameter dimension!')
end
rt = [7  10  28  25  19  20  14   5   6   8  35   1   4   2   9  12  30  31  32   3  11  29  13  33  34  26  27  24  23  18  21  16  22  17  15];
cw__sA_sPB_3o_2B_S_stS__T_dtS__B_astS(sol,tout,theta(1:11),kappa(1:70),options_cvode,plist,pbar,xscale);
sol.x = sol.x(:,rt);
sol.xdot = sol.xdot(:,rt);
if(options_cvode.sensi>0)
    sol.sx = sol.xS(:,rt,:);
    sol.sy = sol.yS;
    sol.sroot = sol.rootS;
    sol.srootval = sol.rootvalS;
elseif(options_cvode.ss>0)
    sol.sx = -sol.J\sol.dxdotdp;
    sol.sy = sol.dydx*sol.sx + sol.dydp;
    sol.sx = sol.sx(rt,:);
end
if(nargout>1)
    varargout{1} = sol.status;
    varargout{2} = sol.t;
    varargout{3} = sol.x;
    varargout{4} = sol.y; % Moments of species
    if(nargout>4)
        varargout{5} = sol.sx;
        varargout{6} = sol.sy;
    end
else
    sol.theta = theta;
    sol.kappa = kappa;
    varargout{1} = sol;
end
end
