%RUN_PLATEMO_EXAMPLE Example of running PWDCMO through PlatEMO.
%
% Set platemoRoot to the directory that contains the official PlatEMO files.

root = fileparts(fileparts(mfilename('fullpath')));
platemoRoot = fullfile(getenv('TEMP'),'PlatEMO_ref','PlatEMO');
if ~exist(fullfile(platemoRoot,'platemo.m'),'file') && ~exist(fullfile(platemoRoot,'Algorithms'),'dir')
    error('Please set platemoRoot to your local PlatEMO/PlatEMO directory.');
end

addpath(genpath(platemoRoot));
addpath(genpath(fullfile(root,'PlatEMO')));

Problem = PWDCMO_DCMOP('N',100,'D',10,'maxFE',20000,'parameter',{1,20,4});
Algorithm = PWDCMO('save',1,'metName',{'IGD','HV'});
Algorithm.Solve(Problem);
