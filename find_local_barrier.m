function [B_f, nu_f, c_f, solver_time] = find_local_barrier(f, g, x, u, omega, sets, ...
    optimization_parameters, solver_parameters, problem_type)
 
    %% Extract parameters
    epsilon_1 = optimization_parameters.epsilon_1;
    epsilon_2 = optimization_parameters.epsilon_2;
    deg_B     = optimization_parameters.deg_B;
    deg_nu    = optimization_parameters.deg_nu;
    L         = optimization_parameters.L;
    c_eps     = optimization_parameters.c_eps;
    alpha     = optimization_parameters.alpha;
    lambda    = optimization_parameters.lambda;

    n_eps      = solver_parameters.num_eps;

    Xi_set      = sets.space_set;
    Wi_set      = sets.W_set;
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

    % Optimize parameter c
    c = sdpvar(1,1);
    Variables = [Variables; c];
    Constraints = [Constraints, c >= c_eps];
     
    % Slack variables (Putinar's Positivstellensatz)
    if(isempty(Xi_set))
        cone_Xi = 0;
    else
        [cone_Xi, vars_Xi, constraints_Xi] = generatecone(Xi_set, x, deg_B);
        Constraints = [Constraints, constraints_Xi];
        Variables   = [Variables; vars_Xi];
    end

    if(isempty(Wi_set))
        cone_Wi = 0;
    else
        [cone_Wi, vars_Wi, constraints_Wi] = generatecone(Wi_set, x, deg_B);
        Constraints = [Constraints, constraints_Wi];
        Variables   = [Variables; vars_Wi];
    end

    [cone_initial, vars_initial, constraints_initial] = generatecone(initial_set, x, deg_B);
    Constraints = [Constraints, constraints_initial];
    Variables   = [Variables; vars_initial];

    [cone_safe,    vars_safe,    constraints_safe]    = generatecone(safe_set, x, deg_B);
    Constraints = [Constraints, constraints_safe];
    Variables   = [Variables; vars_safe];
    
    [cone_unsafe,  vars_unsafe,  constraints_unsafe]  = generatecone(unsafe_set, x, deg_B);
    Constraints = [Constraints, constraints_unsafe];
    Variables   = [Variables; vars_unsafe];
    
    %% Constraints

    % [9a] Bᵢ(xᵢ) ≤ ε₁ᵢ for all xᵢ ∈ Pᵢ(X₀)     
    %       ε₁ᵢ - Bᵢ(xᵢ) - cone_initial is SOS
    condition_init = epsilon_1 -B -cone_initial;
    Constraints = [Constraints, sos(condition_init)];

    % [9b] Bᵢ(xᵢ) - ε₂ᵢ ≥ ϵ for all xᵢ ∈ Pᵢ(Xunsafe)  
    %     [Bᵢ(xᵢ) - ε₂ᵢ ≥ ϵ ⇒ Bᵢ(xᵢ) > ε₂ᵢ]
    %      Bᵢ(xᵢ) - ε₂ᵢ - ϵ - cone_unsafe is SOS
    condition_unsafe = B -epsilon_2 -n_eps -cone_unsafe;
    Constraints = [Constraints, sos(condition_unsafe)];

    % [9d] Bᵢ(gᵢ(xᵢ,ωᵢ,uᵢ(xᵢ),t)) ≤ Bᵢ(xᵢ) for all xᵢ ∈ Pᵢ(Xsafe), ωᵢ ∈ Wᵢ, t ∈ Ωᵢ
    %       Bᵢ(xᵢ) - Bᵢ(gᵢ(xᵢ,ωᵢ,uᵢ(xᵢ),t)) - cone_safe is SOS
    % Notes: 
    % 1. if g is the identity map, the constraint below must be eliminated      
    % 2. if optimization_parameters.set_9d = "Xsafe" then [9d] will be
    % checked. Otherwise:
    % 
    % if optimization_parameters.set_9d = "Xi", then the
    %    condition will be checked instead:
    %
    % [9d'] Bᵢ(gᵢ(xᵢ,ωᵢ,uᵢ(xᵢ),t)) ≤ Bᵢ(xᵢ) for all xᵢ ∈ Xᵢ, ωᵢ ∈ Wᵢ, t ∈ Ωᵢ
    %         Bᵢ(xᵢ) - Bᵢ(gᵢ(xᵢ,ωᵢ,uᵢ(xᵢ),t)) - (cone_Xi+cone_Wi) is SOS
    %
    %    if optimization_parameters.set_9d = "Rn", then the
    %    condition will be checked instead:
    %
    % [9d"] Bᵢ(gᵢ(xᵢ,ωᵢ,uᵢ(xᵢ),t)) ≤ Bᵢ(xᵢ) for all xᵢ ∈ R^(n_x), ωᵢ ∈ R^(#Iᵢ), t ∈ Ωᵢ
    %        Bᵢ(xᵢ) - Bᵢ(gᵢ(xᵢ,ωᵢ,uᵢ(xᵢ),t)) is SOS    
    %
    % 3. the following implication diagram is not necessarily commutative: 
    %
    %  [9d]   ←  [9d SOS]
    %    ↑          |
    %  [9d']  ←  [9d' SOS]
    %    ↑          |
    %  [9d"]  ←  [9d" SOS]
    %
    if(~is_identity_map(g))
        Bg = replace(B, x, g);
        switch(optimization_parameters.set_9d)
            case 'Rn'
                condition_jump = B -Bg;
            case 'Xi'
                condition_jump = B -Bg -(cone_Xi+cone_Wi);
            case 'Xsafe'
                condition_jump = B -Bg -cone_safe;
        end    
        Constraints = [Constraints, sos(condition_jump)];
    end

    % [15a] ‖xᵢ‖² ≤ αᵢBᵢ(xᵢ) + L for all xᵢ ∈ Xᵢ
    %       αᵢBᵢ(xᵢ) + L - ‖xᵢ‖² - cone_Xi is SOS
    condition_norm = alpha*B + L - norm(x)^2 - cone_Xi;
    Constraints = [Constraints, sos(condition_norm)];

    % [15b] (dBᵢ/dxᵢ)fᵢ,ⱼ(xᵢ) ≤ λᵢBᵢ(xᵢ) + cᵢ(‖ωᵢ‖² - L Sᵢ) 
    % for all xᵢ ∈ Xᵢ, ωᵢ ∈ Wᵢ, j = 1, ..., P 
    condition_dynamics = [];
    size_Ii = length(omega)/length(x);
    for j = 1:length(f)
       % Add feedback u(t) = ν(xᵢ(t))
       fv = replace(f{j}, u, nu);
       if(solver_parameters.debug)
            disp('f{j}:')
            sdisplay(f{j})
            disp('fv:')
            sdisplay(fv)
       end
       condition_dynamics = [condition_dynamics, ...
          lambda*B+c*(norm(omega)^2- L*size_Ii)-dBdx*fv-(cone_Xi+cone_Wi)];
       Constraints = [Constraints, sos(condition_dynamics(j))];
    end    

    %% Debug
    if(solver_parameters.debug)
        disp('Barrier function definition:')
        sdisplay(B)
        disp('Feedback control definition:')
        sdisplay(nu)
        disp('Cones:')
        disp('Cone X_i:')
        sdisplay(cone_Xi)
        disp('Cone W_i:')
        sdisplay(cone_Wi)
        disp('Cone initial:')
        sdisplay(cone_initial)
        disp('Cone unsafe:')
        sdisplay(cone_unsafe)
        disp('Cone safe:')
        sdisplay(cone_safe)
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
        disp('Norm condition:')
        sdisplay(condition_norm)
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
    
    [sol,v,Q,res] = solvesos(Constraints,[],solver_options,Variables);

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
        Xi_labels       = generate_constraint_labels(length(Xi_set), 'Xi (g_%d)');
        Wi_labels       = generate_constraint_labels(length(Wi_set), 'Wi (g_%d)');
        Xinitial_labels = generate_constraint_labels(length(initial_set), 'X0i (g_%d)');
        Xsafe_labels    = generate_constraint_labels(length(safe_set), 'Xsafe (g_%d)');
        X_unsafe_labels = generate_constraint_labels(length(unsafe_set), 'Xunsafe (g_%d)');
        cond_15b_labels = generate_constraint_labels(length(f), '15b (mode %d)');
        if(is_identity_map(g))
            cond_9d_label = [];
        else
            cond_9d_label = "9d";
        end
    
        constraints_labels = [...
            Xi_labels; ...
            Wi_labels; ...
            Xinitial_labels; ...
            Xsafe_labels; ...        
            X_unsafe_labels; ...
            "9a"; ...
            "9b"; ...
            cond_9d_label; ...
            "15a"; ...
            cond_15b_labels ...
        ];
        table(constraints_labels, min_eigval, primal_res, size_Q, reliable)
    end
    assert(all(reliable(i) == 1), 'Barrier computed using local barriers is unreliable');

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
    c_f = double(c);
    solver_time = sol.solvertime;

    if(solver_parameters.debug)
        disp(sol)
        disp('B_i:')
        sdisplay(B_f)
        disp('nu_i:')
        sdisplay(nu_f)
        disp('c_i')
        sdisplay(c_f) 
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
    