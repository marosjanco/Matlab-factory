function disp_bin_matrix(bin_matrix,text)
    matrix_order = length(bin_matrix);
    matrix_orderplus = matrix_order+1;
    range = 1:matrix_order;
    range_edges = range + .5;
    
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
    title(text)
    hold off
end
