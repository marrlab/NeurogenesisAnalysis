function [model] = neuro_syms()

%%
% CVODES OPTIONS

% set the default absolute tolerance
model.atol = 1e-8; 
% set the default relative tolerance
model.rtol = 1e-8; 
% set the default maximum number of integration steps
model.maxsteps = 1e4; 
% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
% model.param = 'log';

%%
% STATES


% create state syms
syms x1 x2 x3 

% create state vector
x = [
x1 x2 x3 
];

%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms p1 p2 p3 p4 p7 p8 p9
% create parameter vector 
p = [p1,p2,p3,p4,p7,p8,p9];

%%
% CONSTANTS

% syms k1 k2 k3
% 
% k=[k1,k2,k3];
k=[];

%%
% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

xdot = sym(zeros(size(x)));

xdot(1) = -p1*x1;
xdot(2) = p1*x1+(p2+p3)*x3-p2*x2;
xdot(3) = p2*x2-(p2+p3)*x3-p4*x3;
 

%%
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));
x0(1) = p7;
x0(2) = p8;
x0(3) = p9;

%%
% OBSERVABLES

y = x;


%%
% SYSTEM STRUCT

model.sym.x = x;
model.sym.k = k;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
% model.event = event;
end
