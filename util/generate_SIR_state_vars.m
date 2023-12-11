function [vars] = generate_SIR_state_vars(n)
    prefixes = {'S', 'I', 'R'}; % Variable prefixes
    
    vars = cell(1, n * length(prefixes));
    index = 1;

    for group = 1:n
        for prefix = prefixes
            vars{index} = [prefix{1}, num2str(group)];
            index = index + 1;
        end
    end
end