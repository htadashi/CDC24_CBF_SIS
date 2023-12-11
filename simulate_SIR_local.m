function [B_i, eigvecs] = simulate_SIR_local(subsystem_params, debug)

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

x = cell(n, 1);
for i = 1:n
    x{i} = [S(i);I(i);R(i)];
end

% Input variables
u = sdpvar(n, 1,'full');  

% f{i,j} is the dynamics of the i-th subsystem with j-th switching mode
f = cell(n, P);

for i = 1:n
    for j = 1:P
        f{i,j}   = [m(i) - beta(i, j)*S(i)*I(i) - mu(i)*S(i)        + sum(a(i,:) * S); ...
                           beta(i, j)*S(i)*I(i) - (mu(i)+r(i))*I(i) + sum(b(i,:) * I); ...
                           r(i)*I(i)            - mu(i)*R(i)        + sum(c(i,:) * R)];        
    end
end

% g is a polynomial map describing the jump map of the i-th subsystem
g = cell(n, 1);
for i = 1:n
   g{i} = [(1-p)*S(i); 
                 I(i); 
                 R(i) + p*S(i)];
end

% Define the vector of internal inputs
omega = cell(n, 1);
for i = 1:n
    for j = 1:n
        if (j ~= i) && (a(i,j) ~= 0 || b(i,j) ~= 0 || c(i,j) ~= 0)
            omega{i} = [omega{i}; x{j}];
        end
    end
end

% Define the state space, and the initial, safe and unsafe sets
infected_threshold = 5;

X           = cell(n,1);
safe_set    = cell(n,1);
unsafe_set  = cell(n,1);
initial_set = cell(n,1);

% Note: 
% The relation [initial_set{i} ⊆ safe_set{i}] must be true
for i = 1:n
    X{i}           = [100-(S(i)^2 + I(i)^2 + R(i)^2)];
    safe_set{i}    = [X{i}, -I(i)+infected_threshold];
    unsafe_set{i}  = [X{i},  I(i)-infected_threshold];
    initial_set{i} = [1-(S(i)^2 + I(i)^2 + R(i)^2)];
end

if(debug)
    for i = 1:n
        disp(['X{' num2str(i) '}:']);
        for j = 1:length(X{i})
            sdisplay(X{i}(j));
        end
        disp(['safe_set{' num2str(i) '}:']);
        for j = 1:length(safe_set{i})
            sdisplay(safe_set{i}(j));
        end
        disp(['unsafe_set{' num2str(i) '}:']);
        for j = 1:length(unsafe_set{i})
            sdisplay(unsafe_set{i}(j));
        end
        disp(['initial_set{' num2str(i) '}:']);
        for j = 1:length(initial_set{i})
            sdisplay(initial_set{i}(j));
        end        
    end
end

% Define the space of internal inputs
W = cell(n,1);

for i = 1:n
    for j = 1:n
        if (j ~= i) && (a(i,j) ~= 0 || b(i,j) ~= 0 || c(i,j) ~= 0)
            W{i} = [W{i}; X{j}];
        end
    end
end

sets = cell(n,1);
for i = 1:n
    sets{i} = struct('space_set',   X{i}, ...
                     'W_set',       W{i}, ...
                     'initial_set', initial_set{i}, ...
                     'safe_set',    safe_set{i},    ...
                     'unsafe_set',  unsafe_set{i});
end

% Define the optimization parameters 
% Note: theta must be equal to 1 in the compositional setting
optimization_parameters = cell(n,1);
L = 1e4;
for i = 1:n
    optimization_parameters{i} = struct( ...
        'deg_B',     2, ...
        'deg_nu',    0, ...
        'deg_tau',   2, ...
        'epsilon_1', 0, ...
        'epsilon_2', 0, ...
        'lambda',    -1e-1, ...
        'alpha',     1e-1, ...
        'L',         L, ...
        'c_eps',     1e-8, ...
        'set_9d',    'Rn' ... 
    );
end

% Define solver options
solver_parameters = struct( ...
    'num_eps',       0.001, ...   % to change strict inequalities to non-strict
    'solver',        'mosek', ...
    'debug',         debug, ...
    'warm_start',    0 ...
);

% Search for local barrier certificates
B_i   = cell(n,1);
nu_i  = cell(n,1);
c_i   = cell(n,1);
f_i   = cell(1,P);
solvertime = zeros(n,1);

for i = 1:n
    %Note: f{i,:} unpacks the cell instead of selecting all columns
    [B_i{i}, nu_i{i}, c_i{i}, solvertime(i)] = find_local_barrier(f(i,:), g{i}, x{i}, u(i), omega{i}, sets{i}, ...
                                  optimization_parameters{i}, solver_parameters, 'find_only_barrier');
end

%% Reconstruct global barrier certificate
Lambda_matrix = zeros(n, n);
Gamma_matrix  = zeros(n, n);
for i = 1:n
    Lambda_matrix(i,i) = optimization_parameters{i}.lambda;
    for j = 1:n
        Gamma_matrix(i,j) = c_i{i}*optimization_parameters{j}.alpha;
    end
end
eta_I = min(diag(Lambda_matrix)) * eye(n,n);

% Return largest magnitude eigenvalue and corresponding (left) eigenvector
[eigvecs, max_eigval] = eigs(Lambda_matrix + Gamma_matrix - eta_I, 1);
assert(max_eigval > 0, "The maximal eigenvalue of matrix Λ + Γ - ηI should be positive")

% Display global B (using max formulation)
if(debug)
    disp('Global B:');
    for i=1:n
        disp(['B_' num2str(i) '/v_' num2str(i) ':']);
        sdisplay(B_i{i}/eigvecs(i));
    end
    disp('Right eigenvector:')
    disp(eigvecs)
end

end