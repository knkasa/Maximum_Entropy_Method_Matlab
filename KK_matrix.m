function matrix = KK_matrix( Kdif2mat, params )
    
    matrix=Kdif2mat'*Kdif2mat./(params.error.^2);
    
end

