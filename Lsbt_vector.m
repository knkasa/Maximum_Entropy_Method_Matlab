function vector = Lsbt_vector( expdata, params )
    
    % This uses the data with elastic component already subtracted
    % Interpolate experiment data points
    Lvec=zeros(size(params.xvec));
    Lvec(params.xvec>=min(expdata.x).*params.charge&params.xvec<=max(expdata.x).*params.charge)=...
        interp1(expdata.x,expdata.y,params.xvec(params.xvec>=min(expdata.x).*params.charge&params.xvec<=max(expdata.x).*params.charge)./params.charge);
    vector=Lvec;
    
end

