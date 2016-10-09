function a = Create_a(A)
    matrix_order = length(A{1,1});
    Vn = 1:matrix_order;
    a = cell(1,matrix_order^2);
    for i=Vn
        for j=Vn
            a{Get_position(i,j,matrix_order)} = Vectorize(A{i,j});
        end
    end
end