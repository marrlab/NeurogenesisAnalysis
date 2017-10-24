function [opt] = setPaths()

%directory of the RUN_N_modelselection file
opt.RUN_N_dir=cd;
%directory of toolboxes
cd('../../');
dir1 = cd;
cd(opt.RUN_N_dir);
addpath(genpath([dir1,'/Tools']));

%directory of Cerena toolbox
opt.CERENApath = [dir1,'/Tools/CERENA/examples/neurogenesis/'];

end