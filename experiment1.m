clear all; clc;
yalmip('clear');

subsystem_params = load_sys_parameters(3, 0.001, 1);
simulate_SIR_local(subsystem_params, true);