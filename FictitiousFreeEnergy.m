function scalar = FictitiousFreeEnergy( RATvec,Kdif2mat,L2vec,params )
    
    rhovec=rho_vector(RATvec,params);
    entropy=FictitiousEntropy(RATvec,params);
    GGvec=Kdif2mat*rhovec;
    Lfic=sum((L2vec-GGvec).^2,1)./(params.error.^2)./2;
    scalar=Lfic-params.alpha.*entropy;
    
end

