function scalar = P_alpha_Gm( Efic,Svec,params )
    
%     lnZs=params.dim./2.*log(2.*pi./params.alpha);
%     lnZl=params.dim./2.*log(2.*pi)+params.dim.*log(params.error);
    scalar=-Efic+0.5.*sum(log(params.alpha./(ones(size(Svec)).*params.alpha+Svec)));
    
end

