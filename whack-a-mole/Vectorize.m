function vector = Vectorize(matrix)
    vector = reshape((flipud(matrix))',numel(matrix),1);
end