function A = Create_Aij(matrix_order)
    A = cell(matrix_order,matrix_order);
    Vn = 1:matrix_order;

    for i=Vn
        for j=Vn
            A{i,j}=zeros(matrix_order,matrix_order);
            A{i,j}(i,j)=1;
            Indices = [i,j-1;i,j+1;i-1,j;i+1,j];
            for k = 1:4
                    row = Indices(k,1);
                    col = Indices(k,2);
                    if ismember(row,Vn) && ismember(col,Vn)
                        A{i,j}(row,col)=1;
                    end
            end
        end
    end
end