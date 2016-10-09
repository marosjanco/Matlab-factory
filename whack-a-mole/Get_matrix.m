function matrix = Get_matrix(vector)
    matrix_order = sqrt(length(vector));
    matrix = flipud(reshape(vector,matrix_order,matrix_order)');
end
