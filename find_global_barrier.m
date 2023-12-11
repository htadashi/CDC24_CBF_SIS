function [B_f, nu_f, solver_time] = find_global_barrier(f, g, x, u, sets, ...
    optimization_parameters, solver_parameters, problem_type)
 
    %% Extract parameters
    deg_B     = optimization_parameters.deg_B;
    deg_nu    = optimization_parameters.deg_nu;
    deg_tau   = optimization_parameters.deg_tau;
    lambda    = optimization_parameters.lambda;
    theta     = optimization_parameters.theta;
    
    n_eps      = solver_parameters.num_eps;

    X_set       = sets.space_set;
    initial_set = sets.initial_set;
    safe_set    = sets.safe_set;
    unsafe_set  = sets.unsafe_set;

    Constraints = [];
    Variables   = [];

    %% Variables 
    % Barrier function
    [B, B_i] = polynomial(x, deg_B);
    dBdx = jacobian(B, x);
    if(solver_parameters.debug)
        disp('B:')
        sdisplay(B)
        disp('dB/dx:')
        for j = 1:length(x)
            sdisplay(dBdx(j))
        end
    end    
    Variables   = [Variables; B_i(:)];
    % Feedback controller
    if(strcmp(problem_type,'find_only_barrier'))
        % Only search for barrier
        nu = 0;
    else
        nu = double2sdpvar(zeros(length(u),1));
        for i=1:length(u)
            [nu(i), nu_i{i}] = polynomial(x, deg_nu);
            Variables   = [Variables; nu_i{i}(:)];
        end
    end

    % Slack variables (Putinar's Positivstellensatz)
    [cone_X,       vars_X,       constraints_X]       = generatecone(X_set, x, deg_tau);
    [cone_initial, vars_initial, constraints_initial] = generatecone(initial_set, x, deg_tau);
    if(strcmp(optimization_parameters.set_3d, 'Xsafe'))
       [cone_safe,    vars_safe,    constraints_safe]    = generatecone(safe_set, x, deg_tau); 
    end
    [cone_unsafe,  vars_unsafe,  constraints_unsafe]  = generatecone(unsafe_set, x, deg_tau);

    Constraints = [Constraints, constraints_X];  
    Variables   = [Variables; vars_X];
    Constraints = [Constraints, constraints_initial]; 
    Variables   = [Variables; vars_initial];
    if(strcmp(optimization_parameters.set_3d, 'Xsafe'))
        Constraints = [Constraints, constraints_safe]; 
        Variables   = [Variables; vars_safe];
    end
    Constraints = [Constraints, constraints_unsafe]; 
    Variables   = [Variables; vars_unsafe];
    
    %% Constraints    
    % [3a] B(x) ≤ 0 for all x ∈ X₀     
    %      -B(x) - cone_initial is SOS
    condition_init = -B -cone_initial;
    Constraints = [Constraints, sos(condition_init)]; 

    % [3b] B(x) ≥ ϵ for all x ∈ Xunsafe
    %     [B(x) ≥ ϵ ⇒ B(x) > ε₂]
    %      B(x) - ϵ - cone_unsafe is SOS
    condition_unsafe = B -n_eps -cone_unsafe;
    Constraints = [Constraints, sos(condition_unsafe)]; 

    % [3d] B(g(x,u(x),t)) ≤ θ B(x) for all x ∈ Xsafe, t ∈ Ω
    %      θB(x) - B(g(x,u(x),t)) - cone_safe is SOS
    % Notes:
    % 1. if g is the identity map and θ = 1, the constraint must be 
    %    eliminated 
    % 2. if optimization_parameters.set_3d = 'Xsafe' then [3d] will be
    % checked. Otherwise:
    % 
    % if optimization_parameters.set_3d = "X", then the
    %    condition will be checked instead:
    %
    % [3d'] B(g(x,u(x),t)) ≤ θ B(x) for all x ∈ X, t ∈ Ω
    %         θB(x) - B(g(x,u(x),t)) - cone_X is SOS
    %
    % if optimization_parameters.set_3d = "Rn", then the
    %    condition will be checked instead:
    %
    % [3d"] B(g(x,u(x),t)) ≤ θ B(x) for all x ∈ R^(n_x), t ∈ Ω
    %        θB(x) - B(g(x,u(x),t)) is SOS    
    %
    % 3. the following implication diagram is not necessarily commutative: 
    %
    %  [3d]    ←  [3d SOS]
    %    ↑           |
    %  [3d']   ←  [3d' SOS]
    %    ↑           |
    %  [3d"]   ←  [3d" SOS]
    %
    if(~(is_identity_map(g) && theta == 1))
        Bg = replace(B, x, g);
        switch(optimization_parameters.set_3d)
            case 'Rn'
                condition_jump = theta*B - Bg;
            case 'X'
                condition_jump = theta*B - Bg - cone_X;
            case 'Xsafe'
                condition_jump = theta*B - Bg - cone_safe;
        end        
        Constraints = [Constraints, sos(condition_jump)]; 
    end

    % [3c] ∂B/∂x fp(x, u(x)) ≤ λB(x) for all x ∈ X
    condition_dynamics = [];
    for j = 1:length(f)
       % Add feedback u(t) = ν(x(t))
       fv = replace(f{j}, u, nu);
       if(solver_parameters.debug)
            disp('f{j}:')
            sdisplay(f{j})
            disp('fv:')
            sdisplay(fv)
       end
       condition_dynamics = [condition_dynamics, ...
          lambda*B - dBdx*fv - cone_X];
       Constraints = [Constraints, sos(condition_dynamics(j))];
    end    

    %% Debug
    if(solver_parameters.debug)
        disp('Barrier function definition:')
        sdisplay(B)
        disp('Feedback control definition:')
        sdisplay(nu)
        disp('Cones:')
        disp('Cone X:')
        sdisplay(cone_X)
        disp('Cone initial:')
        sdisplay(cone_initial)
        disp('Cone unsafe:')
        sdisplay(cone_unsafe)
        if(strcmp(optimization_parameters.set_3d, 'Xsafe'))
            disp('Cone safe:')
            sdisplay(cone_safe)
        end
        disp('Variables:')
        disp(Variables)
        disp('Constraints:')
        disp(Constraints)

        disp('Initial set condition:')
        sdisplay(condition_init)
        disp('Unsafe set condition:')
        sdisplay(condition_unsafe)
        if exist('condition_jump', 'var')
            disp('Jump set condition:')
            sdisplay(condition_jump)
        end
        disp('Dynamics condition:')
        for j=1:length(f)
            disp([num2str(j) '-th mode:'])
            sdisplay(condition_dynamics(j))
        end
    end

    %% Solve SOS and check residuals
    % - To use default SOS, use sos.model = 1
    % - To use DSOS, use sos.model  = 5  [experimental]
    % - To use SDSOS, use sos.model = 6  [experimental]
    solver_options = sdpsettings( 'solver',              solver_parameters.solver, ...
                                  'sos.model',           1, ...
                                  'verbose',             solver_parameters.debug, ...
                                  'debug',               solver_parameters.debug);
    
    [sol,v,Q,res] = solvesos(Constraints,0,solver_options,Variables);    

    % Analyze SOS constraints
    primal_res = check(Constraints(logical(is(Constraints, 'sos'))));
    min_eigval = [];
    size_Q = [];
    reliable = [];
    for i = 1:length(Q)
        min_eigval = [min_eigval; min(eig(Q{i}))];
        size_Q     = [size_Q; size(Q{i}, 1)];      
        reliable   = [reliable; min(eig(Q{i})) >  size(Q{i},1)*primal_res(i)];
    end
    if(solver_parameters.debug)
        X_labels       = generate_constraint_labels(length(X_set), 'X (g_%d)');
        initial_labels = generate_constraint_labels(length(initial_set), 'X0 (g_%d)');
        if(strcmp(optimization_parameters.set_3d, 'Xsafe'))
            Xsafe_labels = generate_constraint_labels(length(safe_set), 'Xsafe (g_%d)');
        else
            Xsafe_labels = [];
        end
        Xunsafe_labels = generate_constraint_labels(length(unsafe_set), 'Xunsafe (g_%d)');    
        cond_3c_labels = generate_constraint_labels(length(f), '3c (mode %d)');
        if(is_identity_map(g) && theta == 1)
            cond_3d_label = [];
        else
            cond_3d_label = "3d";        
        end
        constraints_labels = [...
            X_labels; ...
            initial_labels; ...
            Xsafe_labels; ...
            Xunsafe_labels; ...        
            "3a"; ...
            "3b"; ...
            cond_3d_label; ... 
            cond_3c_labels ...
        ];
        table(constraints_labels, min_eigval, primal_res, size_Q, reliable)
    end
    assert(all(reliable(i) == 1), 'Barrier computed global conditions is unreliable');

    % Return solutions
    B_f  = replace(B, B_i, double(B_i));
    if(strcmp(problem_type,'find_only_barrier'))
        nu_f = nu;
    else
        nu_f = double2sdpvar(zeros(length(u),1));
        for i=1:length(u)       
            nu_f(i) = replace(nu(i), nu_i{i}, double(nu_i{i}));
        end
    end
    solver_time = sol.solvertime;

    if(solver_parameters.debug)
        disp(sol)
        disp('B:')
        sdisplay(B_f)
        disp('nu:')
        sdisplay(nu_f)
    end
end

function is_id_map = is_identity_map(g)

    numericArray = cellfun(@str2double, sdisplay(jacobian(g)));
    is_id_map = isequal(numericArray, eye(size(numericArray)));    

end

function [multiplicative_cone, extra_variables, SOS_constraints] = ...
    generatecone(polynomial_set, x, degree)
    
    extra_variables = [];
    SOS_constraints = [];

    multiplicative_cone = 0;
    for p = 1:length(polynomial_set)
        [tau{p}, tau_i{p}]  = polynomial(x, degree);
        SOS_constraints     = [SOS_constraints, sos(tau{p})];
        multiplicative_cone = multiplicative_cone + tau{p} * polynomial_set(p);
        extra_variables     = [extra_variables; tau_i{p}(:)];
    end

end

function labels = generate_constraint_labels(n, constraint_str_spec)
    labels = [];
    for i=1:n
        label = sprintf(constraint_str_spec, i);
        labels = [labels; label];
    end
end  