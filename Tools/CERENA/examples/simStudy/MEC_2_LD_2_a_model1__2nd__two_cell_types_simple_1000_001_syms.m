% function [model] = MEC_2_LD_2_a_model1__2nd__two_cell_types_simple_1000_001_syms(f0_user)
function [model] = MEC_2_LD_2_a_model1__2nd__two_cell_types_simple_1000_001_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms mu_1 mu_2 C_1_1 C_1_2 C_2_2

x = [
mu_1, mu_2, C_1_1, C_1_2, C_2_2 ...
];

% PARAMETERS

syms r_divA pA_AA pA_BB 

% KAPPA (constant parameters)

syms indmu1 indmu2 indC1 indC2 indC3 indC4 kmu01 kmu02 kC01 kC02 kC03 kC04 

syms t

p = [r_divA,pA_AA,pA_BB];

k = [indmu1,indmu2,indC1,indC2,indC3,indC4,kmu01,kmu02,kC01,kC02,kC03,kC04];

if nargin > 0
   f0_user = varargin{1};
   if ~isnumeric(f0_user)
      p_user = setdiff(symvar(f0_user),p);
      % ADDITIONAL PARAMETERS IN INITIAL CONDITIONS
      p = [p,p_user];
   end
	fmu01 = f0_user(1); 
	fmu02 = f0_user(2); 
	fC01 = f0_user(3); 
	fC02 = f0_user(4); 
	fC03 = f0_user(5); 
	fC04 = f0_user(6); 
else
	fmu01 = 1; 
	fmu02 = 0; 
	fC01 = 999/1000; 
	fC02 = 0; 
	fC03 = 0; 
	fC04 = 0; 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = mu_1*pA_AA*r_divA - mu_1*pA_BB*r_divA;
xdot(2) = 2*mu_1*pA_BB*r_divA - mu_1*r_divA*(pA_AA + pA_BB - 1);
xdot(3) = 2*C_1_1*pA_AA*r_divA - 2*C_1_1*pA_BB*r_divA + mu_1*pA_AA*r_divA + mu_1*pA_BB*r_divA;
xdot(4) = C_1_2*pA_AA*r_divA - C_1_1*r_divA*(pA_AA + pA_BB - 1) + 2*C_1_1*pA_BB*r_divA - C_1_2*pA_BB*r_divA - 2*mu_1*pA_BB*r_divA;
xdot(5) = 4*C_1_2*pA_BB*r_divA - mu_1*r_divA*(pA_AA + pA_BB - 1) - 2*C_1_2*r_divA*(pA_AA + pA_BB - 1) + 4*mu_1*pA_BB*r_divA;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = indmu1*kmu01 - fmu01*(indmu1 - 1);
x0(2) = indmu2*kmu02 - fmu02*(indmu2 - 1);
x0(3) = indC1*kC01 - fC01*(indC1 - 1);
x0(4) = indC2*kC02 - fC02*(indC2 - 1);
x0(5) = indC3*kC03 - fC03*(indC3 - 1);
x0(6) = indC4*kC04 - fC04*(indC4 - 1);

% OBSERVABLES

y = sym(zeros(5,1));

y(1) = mu_1;
y(2) = mu_2;
y(3) = C_1_1;
y(4) = C_1_2;
y(5) = C_2_2;

% SYSTEM STRUCT

model.sym.nmx = 0;
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 0;
end