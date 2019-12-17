function vector = errorfun( RATvec,KKmat,KGvec,params )
    
    rhovec=rho_vector(RATvec,params);
    vector=KKmat*rhovec-KGvec+params.alpha.*RATvec;
    
end

