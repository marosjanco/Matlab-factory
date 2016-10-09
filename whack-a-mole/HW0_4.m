close all
clear all

n = 10;                           % order of the square matrix (board)
disp_sol_visually = 1;
disp_sol_positions = false;
random_hitting_order = true;     % hit solution positions randomly or 
                                  % in ascending order


% Initiate a random board configuration
B_init = Create_initial(n);

% Solve the initial board by getting all positions to hit (in c)
[sol_pos,solution_exist] = solve_board(B_init);

if solution_exist
    if disp_sol_positions && disp_sol_visually
        fprintf('The sequence of solution positions to be hit:\n')
        disp(sol_pos')
        pause(1)
        disp_whole_sol(B_init,sol_pos,random_hitting_order)
    elseif disp_sol_positions && ~disp_sol_visually
        fprintf('The sequence of solution positions to be hit:\n')
        disp(sol_pos')
    elseif ~disp_sol_positions && disp_sol_visually
        disp_whole_sol(B_init,sol_pos,random_hitting_order)
    else
        fprintf('You opted to show nothing!\n')
    end
else
    fprintf('No solution for the initial setting!\n')
end