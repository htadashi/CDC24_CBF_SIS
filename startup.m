%% Load dependencies

% Util functions
addpath(genpath(".\util"));

% MOSEK
MOSEK_path = "C:\Program Files\Mosek\10.1\toolbox\r2017aom";
addpath(genpath(MOSEK_path));

% YALMIP 
YALMIP_path = ".\dependencies\YALMIP-R20230622";
addpath(genpath(YALMIP_path));