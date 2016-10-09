function [sol_positions,sol_exists] = solve_board(B_init)
    global A
    
    n = length(B_init);   % order of initial square board
    nn = n*n;
    
    % Create A^(i,j) and a^(k) matrices/vectors
    A = Create_Aij(n);
    a = Create_a(A);

    % concatenate a^(k) vectors into AA matrix and reshape the initial
    % configuration into a vector
    b_init = Vectorize(B_init);
    AA = zeros(nn,nn);
    for j=1:nn
        AA(:,j)=a{j};
    end

    % Turn off the warning in case of no solution
    nosol_id = 'comm:gflineq:NoSolution';
    warning('off',nosol_id)
    % compute the solution c
    [sol_positions,sol_exists] = gflineq(AA,b_init,2);
    warning('on',nosol_id)

    if sol_exists
        sol_positions = find(sol_positions);
    end
end