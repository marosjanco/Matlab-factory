function pos = Get_position(row,col,matrix_order)
    pos = matrix_order*(matrix_order-row)+col;
end