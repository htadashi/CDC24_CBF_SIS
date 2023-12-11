function subsystem_params = load_sys_parameters(subsystem_order, random_matrix_std, experiment_n)

    % System parameters
    n  = subsystem_order;   % Number of patches (subsystems)
    N  = 3*n;               % Number of states
    P  = 2;                 % Number of modes
    
    % Infection and vaccination parameters
    % Recruitment rate
    m    = 0.02   * ones(n,1);
    % Death rate
    mu   = 0.0001 * ones(n,1);
    % beta(i,j) is the transmission rate in the i-th city at the j-th mode
    beta = [0.001*ones(n,1), 0.002*ones(n,1)];
    % Recovery rate
    r = 0.05 * ones(n,1);
    % Vaccine success probability
    p = 0.8; 

    % Immigration and emigration rates
    if(experiment_n == 1)
        a = [-0.0007,  0.0014,  0.0005;
              0.0006, -0.0021,  0.0019;
              0.0001,  0.0007, -0.0024];

        b = [-0.0014,  0.0002,  0.0002;
              0.0011, -0.0005,  0.0001;
              0.0003,  0.0003, -0.0003];

        c = [-0.0002,  0.0003,  0.0004;
              0.0001, -0.0006,  0.0013;
              0.0001,  0.0003, -0.0017];
    elseif(experiment_n == 2)
        a = random_ring_matrix(n, random_matrix_std);
        b = random_ring_matrix(n, random_matrix_std);
        c = random_ring_matrix(n, random_matrix_std);
    end

    subsystem_params = struct(...
        'n', n, ...               
        'N', N, ...              
        'P', P, ...     
        'm', m, ...
        'mu', mu, ...           
        'beta', beta, ...          
        'r', r, ...                
        'p', p, ...                
        'a', a, ...       
        'b', b, ...           
        'c', c  ...
    );

end