function matrix = Kdif2_matrix( params )
    
%Calculate x value for each matrix element
  xmat=repmat(params.xvec,1,numel(params.xvec)); 
    
   %Calculate phonon energy for each matrix element
   Epmat=repmat(params.xvec',numel(params.xvec),1);
    
    %Numerically obtain second derivative ( function f3 is defined below )
    matrix = params.Nbands.*4.*pi.*params.mstar.*params.m0.*(params.charge./   ...
              (params.hbar.*2.*pi)).^3./(params.delta).^2.*(f3(Epmat,xmat+ones(size(xmat)).* ...
              params.delta)+f3(Epmat,xmat-ones(size(xmat)).*params.delta)-2.*f3(Epmat,xmat));
    
    % define function f3 used in above equation
    function res=f3(Epmat,xmat)     %Epmat: phonon, xmat: bias v
        EFmat=ones(size(xmat)).*params.EF;
        I1=(Epmat<=xmat-EFmat);
        I2=(Epmat>xmat-EFmat & Epmat<=xmat);
        I3=(Epmat>xmat);
        res=zeros(size(xmat));
        res(I1)=EFmat(I1).^2;
        res(I2)=-(Epmat(I2)-xmat(I2)).*(Epmat(I2)-xmat(I2)+2.*EFmat(I2));
        res(I3)=0;
    end
    
end

