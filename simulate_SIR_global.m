function [B] = simulate_SIR_global(subsystem_params, debug)

% System parameters
n = subsystem_params.n;
N = subsystem_params.N;
P = subsystem_params.P;
m = subsystem_params.m;
mu = subsystem_params.mu;
beta = subsystem_params.beta;
r = subsystem_params.r;
p = subsystem_params.p;
a = subsystem_params.a;
b = subsystem_params.b;
c = subsystem_params.c;

% State variables
S = sdpvar(n, 1, 'full');
I = sdpvar(n, 1, 'full');
R = sdpvar(n, 1, 'full');
x = [];
for i = 1:n
    x = [x; S(i); I(i); R(i)];
end

% Input variables
u = sdpvar(N, 1,'full');  

% f{j} is the dynamics of the j-th switching mode
f = cell(1, P);

for j = 1:P
    for i = 1:n
        f_ij   = [m(i) - beta(i, j)*S(i)*I(i) - mu(i)*S(i)        + sum(a(i,:) * S); ...
                         beta(i, j)*S(i)*I(i) - (mu(i)+r(i))*I(i) + sum(b(i,:) * I); ...
                         r(i)*I(i)            - mu(i)*R(i)        + sum(c(i,:) * R)];        
        f{j}   = [f{j}; f_ij];
    end
end

% g describes the jump map of the entire state space
g = [];
for i = 1:n
    g = [g;
         (1-p)*S(i); % Next S
               I(i); % Next I
         R(i) + p*S(i)]; % Next R
end


% Define the state space, and the initial, safe and unsafe sets 
infected_threshold = 5;

X           = [100-(sum(S.^2)+sum(I.^2)+sum(R.^2))];
unsafe_set  = [X];
for i = 1:n
    unsafe_set = [unsafe_set, I(i)-infected_threshold];
end
% Note: The safe set must be the complement of the unsafe set
%       and must be a basic semi-algebraic set.
%       If this is not the case, we should replace Xsafe by X or R^n 
%       in the SOS conditions 
%
% i.e. 
% ```
% for i = 1:n
%     safe_set = [safe_set, -I(i)+infected_threshold];
% end
% ```
% is not valid because it is not the correct complement.

% Note: here we do not use initial_set = [safe_set, ...
%  1-(sum(S.^2)+sum(I.^2)+sum(R.^2))]; because the intersection is 
% already given by the array below:
initial_set = [1-(sum(S.^2)+sum(I.^2)+sum(R.^2))];

if(debug)
    disp('X:');
    sdisplay(X);
    disp('initial_set:');
    for i = 1:length(initial_set)
        sdisplay(initial_set(i));
    end
    %disp('safe_set:');
    %for i = 1:length(safe_set)
    %    sdisplay(safe_set(i));
    %end
    disp('unsafe_set:');
    for i = 1:length(unsafe_set)
        sdisplay(unsafe_set(i));
    end
end

sets = struct('space_set',   X, ...
              'initial_set', initial_set, ...
              'safe_set',    X, ...
              'unsafe_set',  unsafe_set);

% Define the optimization parameters 
optimization_parameters = struct( ...
    'deg_B',  2, ...
    'deg_nu',  0, ...
    'deg_tau', 2, ...
    'lambda', -1e-1, ...
    'theta',  1.0, ...
    'set_3d', 'Rn' ...
);

% Define solver options
% http://www.penopt.com/doc/penbmi2_1.pdf
solver_parameters = struct( ...
    'num_eps',       0.001, ...   % to change strict inequalities to non-strict
    'solver',        'mosek', ...
    'debug',         debug, ...
    'warm_start',    0 ...
);


[B, nu, ~] = find_global_barrier(f, g, x, u, sets, ...
                              optimization_parameters, solver_parameters, 'find_only_barrier');
