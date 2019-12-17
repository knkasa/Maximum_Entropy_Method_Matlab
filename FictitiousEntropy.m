function scalar = FictitiousEntropy( RATvec,params )
    
    rhovec=rho_vector(RATvec,params);
    scalar=sum(rhovec-params.mvec-rhovec.*RATvec,1);
    
end

