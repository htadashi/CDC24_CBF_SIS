clc; clearvars;

n_subsystems = 3:6; 
no_experiments_per_order = 20;
no_samples = 5;

% Run local conditions then run global conditions with same parameters
n_sucessful_experiments = [];
local_avg  = [];
global_avg = [];
local_std  = [];
global_std = [];

elapsed_times_local  = {};
elapsed_times_global = {};

for n = n_subsystems
    rng(2);
    elapsed_time_local  = nan(no_experiments_per_order, 1);
    elapsed_time_global = nan(no_experiments_per_order, 1);
    successes_cases = 0;
    for i = 1:no_experiments_per_order  
        if(successes_cases >= no_samples)
            continue
        end
        
        subsystem_params = load_sys_parameters(n, 0.001, 2);
        try 

            yalmip('clear')
            tic
            simulate_SIR_local(subsystem_params, false);
            disp(['Exp ' num2str(i) ': Local condition found reliable solution.']);
            elapsed_time_local(i) = toc;

            yalmip('clear')
            tic 
            simulate_SIR_global(subsystem_params, false);
            disp(['Exp ' num2str(i) ': Global condition found reliable solution.']);
            elapsed_time_global(i) = toc;
            
            successes_cases = successes_cases + 1;
        catch exception
            toc;
            disp(['Exp ' num2str(i) ': ' exception.message]);
        end    
    end

    elapsed_times_local{n} = elapsed_time_local;
    elapsed_times_global{n} = elapsed_time_global;

    elapsed_time_local_wo_NaN  = elapsed_time_local(~isnan(elapsed_time_local));
    elapsed_time_global_wo_NaN = elapsed_time_global(~isnan(elapsed_time_global));
    
    n_sucessful_experiments = [n_sucessful_experiments; length(elapsed_time_local_wo_NaN)];
    local_avg  = [local_avg;  mean(elapsed_time_local_wo_NaN(1:no_samples))];
    global_avg = [global_avg; mean(elapsed_time_global_wo_NaN(1:no_samples))];
    local_std  = [local_std;  std(elapsed_time_local_wo_NaN(1:no_samples))];
    global_std = [global_std; std(elapsed_time_global_wo_NaN(1:no_samples))];

end

% Generate report
table(n_subsystems', n_sucessful_experiments, global_avg, global_std, local_avg, local_std)
