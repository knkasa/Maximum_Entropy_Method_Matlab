function vector = KG_vector( Kdif2mat, L2vec, params )
    
    vector=Kdif2mat'*L2vec./(params.error.^2);
    
end

