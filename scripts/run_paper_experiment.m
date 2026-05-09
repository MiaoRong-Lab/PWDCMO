%RUN_PAPER_EXPERIMENT Standalone PWDCMO experiment entry.
%
% Adjust instance, random seed, and output options as needed before running
% the full experiment. The default time points match the original scripts:
% 0:0.25:20 with 20 independent runs.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Standalone'));

instance = 1;              % DCMOP1..DCMOP8
N = 200;
M = 2;
V = 10;
pps = 100;
G = 1000;
min_range = zeros(1,V);
max_range = [1,2*ones(1,V-1)];

options = struct();
options.seed = [];
options.max_batch = 20;
options.time_points = 0:0.25:20;
options.plot = false;
options.save_figures = false;
options.figure_dir = fullfile(root,'results','figures');
options.verbose = true;

[OP_solution,OP_solution_pf,PF_true,time1,time2,IGD1,IGD2,HV,SP] = ...
    main(instance,N,M,V,pps,G,min_range,max_range,options);

outDir = fullfile(root,'results');
if ~exist(outDir,'dir')
    mkdir(outDir);
end
save(fullfile(outDir,sprintf('PWDCMO_DCMOP%d.mat',instance)), ...
    'OP_solution','OP_solution_pf','PF_true','time1','time2','IGD1','IGD2','HV','SP','options');
