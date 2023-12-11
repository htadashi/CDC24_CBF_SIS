function [matrix] = random_ring_matrix(n, std)
    matrix = zeros(n, n);
    for i=1:n
        matrix(i, switch_index(mod(i-1, n), n)) = std * abs(randn());
        matrix(i, switch_index(mod(i+1, n), n)) = std * abs(randn());
    end
    column_sum = sum(matrix, 1);
    for i=1:n
        matrix(i,i) = -column_sum(i);
    end
end

function idx = switch_index(index, n)
    if index == 0
        idx = n;
    else
        idx = index;
    end    
end    