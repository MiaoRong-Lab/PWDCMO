%RUN_SMOKE_TEST Quick reproducibility check for the standalone code.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Standalone'));
rng(1);

options = struct();
options.seed = 1;
options.max_batch = 1;
options.time_points = 0;
options.plot = false;
options.verbose = false;

N = 20;
M = 2;
V = 10;
pps = 10;
min_range = zeros(1,V);
max_range = [1,2*ones(1,V-1)];

[OP_solution,OP_solution_pf,PF_true,time1,time2,IGD1,IGD2,HV,SP] = ...
    main(1,N,M,V,pps,1,min_range,max_range,options);

outDir = fullfile(root,'results');
if ~exist(outDir,'dir')
    mkdir(outDir);
end
save(fullfile(outDir,'smoke_test.mat'), ...
    'OP_solution','OP_solution_pf','PF_true','time1','time2','IGD1','IGD2','HV','SP');
fprintf('Smoke test completed. IGD = %.6f\n',IGD2{1}{1});
