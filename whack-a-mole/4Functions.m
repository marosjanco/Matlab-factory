function Binit = Create_initial(matrix_order)
       Binit = randi(2,matrix_order)-1;
end

function [row,col] = Get_coordinates(pos,matrix_order)
    row = matrix_order - floor((pos-1)/matrix_order);
    col = mod((pos-1),matrix_order) + 1;
end

function pos = Get_position(row,col,matrix_order)
    pos = matrix_order*(matrix_order-row)+col;
end

function vector = Vectorize(matrix)
    vector = reshape((flipud(matrix))',1,numel(matrix));
end

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

function a = Create_a(A)
    matrix_order = length(A{1,1});
    A = cell(matrix_order,matrix_order);
    Vn = 1:matrix_order;
    a = cell(1,matrix_order^2);
    for i=Vn
        for j=Vn
            a{Get_position(i,j,matrix_order)} = Vectorize(A{i,j});
        end
    end
end

function matrix = Get_matrix(vector)
    matrix_order = sqrt(length(vector));
    matrix = flipud(reshape(vector,matrix_order,matrix_order)');
end

function disp_bin_matrix(bin_matrix)
    matrix_order = length(bin_matrix);
    matrix_orderplus = matrix_order+1;
    range = 1:matrix_order;
    range_edges = range + .5;
    
    close all
    imagesc(range_edges,range_edges,bin_matrix);
    hold on
    set(gca,'XTick',range_edges,'YTick',range_edges,'XTicklabel',range,'YTicklabel',range,...
        'XLim',[1 matrix_orderplus],'YLim',[1 matrix_orderplus])
    for i = 2:matrix_order
       plot([.5,matrix_orderplus],[i,i],'g-','linewidth',.5);
       plot([i,i],[.5,matrix_orderplus],'g-','linewidth',.5);
    end
    axis square
    if all(bin_matrix(:)==0)        
        colormap([1 1 1])
    elseif all(bin_matrix(:)==1) 
        colormap([0 0 0])
    else
        colormap([1 1 1;0 0 0]);
    end
    xlabel('column')
    ylabel('row')
    title('Whack-a-mole')
    hold off
end
