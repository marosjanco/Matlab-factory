function [row,col] = Get_coordinates(pos,matrix_order)
    row = matrix_order - floor((pos-1)/matrix_order);
    col = mod((pos-1),matrix_order) + 1;
end