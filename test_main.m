function test_main( inputfilepath ) 
% specify CSV file at inputfilepath, which includes various parameters for calculation 
% run code as follows:
% *******  test_main('./result/test_input.csv')   ******   ;

%------------obtaining the parameters from text files------------------
%The format in text file is
%ParameterName,Parameter Value
%These parameters are saved in variable "inputparam"
fid=fopen(inputfilepath,'r');
while(1)
    str=fgetl(fid);
    if(str<0)
        break;
    end
    arr=strsplit(str,',');
    if(numel(arr{1})>0)
        inputparam.(arr{1})=arr{2};
    end
end
fclose(fid);

%correcting the file path format if necessary
if(numel(inputparam.dirpath)>0 && strcmp(inputparam.dirpath(end),'/')~=1)
    inputparam.dirpath=[inputparam.dirpath '/'];
end

display(inputparam)



%Exp_Data includes the data points from experiment using the CSV format, ->
% ->  which is loaded by the custom funtion "loadcsvdata"
Exp_Data=loadcsvdata(inputparam.datasbtfilepath);   %Exp data file

%Initial_Guess includes the first guess for the inverted function
Initial_Guess=loadcsvdata(inputparam.initialguesspath);

%Process input parameters
%Processed parameters are saved in the variable ****"params"*****
%Parameters will be named   "params.xxx"

% Whether or not to display graphs
if(isfield(inputparam,'display') && eval(inputparam.display))
    params.display=1;
else
    params.display=0;
end

% Whether or not to optimize alpha (weighting parameter for entropy component)
% see input file .csv    ( 1=YES   0=No )
if(isfield(inputparam,'alphaoptim') && eval(inputparam.alphaoptim))
    params.alphaoptim=1;
else
    params.alphaoptim=0;
end

%Switch for different way of calculating the target function. Specific to
%my application.   (If model function m(x) /= const, hwnorm should be 1 )
%  Choose hwnorm=0  for general case, unless you know what m(x) should be.
if(isfield(inputparam,'hwnorm') && eval(inputparam.hwnorm))
    params.hwnorm=1;
else
    params.hwnorm=0;
end

%Later the experimental data points are devided by factor to satisfy the
%normalization condition.   see input file  .csv   factor=0.5651(default) ??? 
% in hwnorm=0, doesnt matter. 
if(isfield(inputparam,'factor') && eval(inputparam.factor)>0)
    factor=eval(inputparam.factor);
else
    factor=-1;
end
if(isfield(inputparam,'Nbands') && eval(inputparam.Nbands)>0)
    params.Nbands=eval(inputparam.Nbands);
else
    params.Nbands=3;
end
if(isfield(inputparam,'mstartun') && eval(inputparam.mstartun)>0)
    params.mstartun=eval(inputparam.mstartun);
else
    params.mstartun=eval(inoutparam.mstar);
end

%defining parameters using params
params.hbar=1.05457173E-34;
params.kB=1.3806488E-23;
params.charge=1.60217657E-19;
params.phi0=2.067833758E-15;
params.eps0=8.854187817E-12;
params.m0=9.10938291E-31;

params.mstar=eval(inputparam.mstar);

%Fermi energy of semiconductor in Jule.  This is defined in "input.xxx" 
params.EF=eval(inputparam.EF).*params.charge;     
params.area=eval(inputparam.contactarea);       %Contact area size (m^-2)

params.thickness=eval(inputparam.thickness);      %barrier thickness
params.beta=1;                 %correction from uniform barrier
params.A=(4.*pi.*params.beta.*params.thickness./(params.hbar.*2.*pi)).*sqrt(2.*params.mstar.*params.m0);  %WKB cofficient
params.barrier=params.charge.*eval(inputparam.barrierheight);                %Jule, average barrier height from EF1

%parames.error is the error values for each data points for Exp_Data. This
%is set to constant for current case.  (standard error for exp data)
params.error=eval(inputparam.error);

%Initial guess for alpha (weighting parameter for entropy component). 
%Choosing correct guess is important or it does not converge 
params.alpha=eval(inputparam.alpha);

%Minimun value of rhovec (inverted function * Delta x) below which program thinks is
%negligible  (care needs to be taken
params.minlevel=eval(inputparam.minlevel);

%Interval of x points  ( x value of the inverted function)   ***********
params.delta=params.charge.*eval(inputparam.delta);

% ( x value of the inverted function in units [eV])  
%Minimum value of x used in the calculation (the first component of xvec)
% charge 10^-19 is included
params.xoffset=params.charge.*eval(inputparam.xoffset);

%Defines the number of points and dimension of matrix/vector used in the
%calculation  (# of x value also for the inverted function)
params.dim=eval(inputparam.dimension);

%********************************************
%Vector for x-axis ( note that units is in Joul, not eV)
params.xvec=linspace(params.xoffset,params.xoffset+(params.dim-1).*params.delta,params.dim)';
xscaled=params.xvec./max(params.xvec);
%*******************************************


%"Model function" used in MEM. If integrate over x, it should give 1.
params.mvec=1;
if(params.hwnorm==0)
    params.mvec=params.mvec./(params.xvec./max(params.xvec));
end
params.mvec=params.mvec./sum(params.mvec);

%Path to a temporary file which records realtime calculation status/history
tmpfilepath=[inputparam.dirpath num2str(params.EF/params.charge) '_' num2str(params.barrier/params.charge) '_' num2str(params.thickness/0.382E-9) '_dump.txt'];
fid=fopen(tmpfilepath,'a');
fprintf(fid,'---------------Start New Session----------------\n');
fclose(fid);

% Initial guess of alpha2f? which produces rho vector (initial guess of inverted function) 
a2Fini=ones(size(params.xvec)).*1E+50;
%Interpolate the initial guess
a2Fini(params.xvec>=min(Initial_Guess.x)&params.xvec<=max(Initial_Guess.x))=interp1(Initial_Guess.x,Initial_Guess.y,params.xvec(params.xvec>=min(Initial_Guess.x)&params.xvec<=max(Initial_Guess.x)));
a2Fini(params.xvec<min(Initial_Guess.x)|params.xvec>max(Initial_Guess.x))=min(a2Fini);

%Obtain rho vector (see MEM notes) from alpha2F  
% **** rho is the solution you're looking for *******    
% Note: exp_data = sum( k(n,m)*rho(m) )
if(params.hwnorm==1)
    rhovecini=a2Fini.*params.delta.*params.xvec;
    if(factor<0)
        factor=sum(rhovecini);
    end
    rhovecini=rhovecini./factor;
else
    if(factor<0)   % factor=1 if m(x)=constant  
        factor=1;
    end
    
    %******************* initial guess of inverted function  *********
    rhovecini=a2Fini.*params.delta./factor; % initial guess is constant
    
    %**************************************************************
    
end

%These are modified later during the iteration
errorini=params.error;
newfactor=factor;

if(params.display)
    f1=figure();
    f2=figure();
    f3=figure();
end




relative_err1=0;    %Difference of alpha from the optimum value
relative_err2=0;
alpha1=0;
alpha2=0;
multiplier=0.1;     %Used to calculate next alpha in Dogleg method (see Mem notes pg.44)
j=1;
k=1;
%-----------------------------iteration for alpha--------------------------
flag1=0;
while(1)
    
    %-------------------------iteration for normalization------------------
    flag2=0;
    while(1)
        
        if(k>1)
            %Adjusting rho vector so that log(rhovec) does not diverge
            params.minlevel=1/params.dim*1E-10;
            rhovecini(rhovecini<params.minlevel)=params.minlevel;
        end
        
        %Obtaining RAT vector used in MEM (see MEM notes P42)
        RATvecini=log(rhovecini./params.mvec);
        
        % Uses experiment data with elastic component already subtracted
        %-------------------------------------------------------------------------------------------
        % The Exp_Data has to be normalized in units of A/m^2!!
        L2vec=Lsbt_vector(Exp_Data,params)./factor;
        params.error=errorini./factor;
        %-------------------------------------------------------------------------------------------
        if(params.display)
            figure(f1);
            plot(params.xvec,rhovecini);
            hold on
            %plot(params.xvec,params.mvec);
            %hold off
            figure(f2);
            plot(params.xvec./params.charge,L2vec);
            drawnow
        end
        
        %Obtain kernel of the problem of interest
        if(params.hwnorm==1)
            Kdif2mat=Kdif2_matrix_hwnorm(params);
        else
            Kdif2mat=Kdif2_matrix(params);
        end
        
        %Calculate KK matrix and KG vector used in MEM (see MEM notes P42) using functions
        KKmat=KK_matrix(Kdif2mat,params);
        KGvec=KG_vector(Kdif2mat,L2vec,params);
        err=errorfun(RATvecini,KKmat,KGvec,params); % derivative of free energy
        
        %Solve dim simultaneous equations using MATLAB function
        options=optimoptions('fsolve','display','iter');
        options=optimoptions(options,'OutputFcn',@(x,optimValues,state,varargin) optimmyplot(x,optimValues,state,varargin,params,tmpfilepath));
        try
            
        %  The derivative of Entropy is set to zero to find rho(x)
            [results.x,results.fval]=fsolve(@(x) errorfun(x,KKmat,KGvec,params),RATvecini,options);
       
        catch ME
            % In case of error, save the current data
            save('dump_RATvec.mat','RATvecini');
            save('dump_rhovec.mat','rhovecini');
            save('dump_params.mat','params');
            rethrow(ME);
        end
        
        %Saving calculated results in "results" struct
        results.rhovec=rho_vector(results.x,params);
        results.Jcalc=Kdif2mat*results.rhovec;
        results.Efic=FictitiousFreeEnergy(results.x,Kdif2mat,L2vec,params);
        results.L2vec=L2vec;
        results.params=params;
        results.err=errorfun(results.x,KKmat,KGvec,params);

        %Calculate parameters to optimize the  alpha value (see MEM notes)
        % lambda is the eigenvalue of the matrix
        LAMBDAmat=LAMBDA_matrix(results.rhovec,KKmat,params);
        % svd = singular value decomposition
        Svec=svd(LAMBDAmat);
        % Ng = sum lam/(a+lam)
        Ng=sum(Svec./(ones(size(Svec)).*results.params.alpha+Svec));
        Sfic=FictitiousEntropy(results.x,params); % Sfic = entropy
        Ent=-2.*results.params.alpha.*Sfic;
        
        results.P_a_Gm = P_alpha_Gm(results.Efic,Svec,results.params);
        
        %Processing iteration for the normalization
        %Estimate factor to devide experiment data points (L2vec) to
        %satisfy the normalization condition (integral(alpha2F dx)=1)
        newfactor=sum(results.rhovec);
        factor=factor*newfactor;
        results.factor=factor;
        if(abs(newfactor-1)<0.01)
            flag2=1;
        end
        fprintf('%12s %12s %12.6g %12.6g %12.6g %12.6g\n','','',relative_err1,alpha1,newfactor,factor);
        fid=fopen(tmpfilepath,'a');
        fprintf(fid,'%12s %12s %12.6g %12.6g %12.6g %12.6g\n','','',relative_err1,alpha1,newfactor,factor);
        fclose(fid);
        
        k=k+1;
        
        if(flag2)
            break;
        end
        
        % Use previous results for initial guess of rho vector
        rhovecini=results.rhovec;
        
    % end for inner while    
    end %---------------end-iteration for normalization-end----------------
    
    
    if(params.alphaoptim<1)
        break;
    end
    
    % Processing interation for the alpha value using Dogleg method
    relative_err2=relative_err1;
    relative_err1=(Ng-Ent);
    alpha2=alpha1;
    alpha1=params.alpha;
    params.alpha=alpha1.*(1-multiplier);
    relative_err_hist(j)=relative_err1;
    alpha_hist(j)=alpha1;
    if(params.display)
        figure(f3)
        plot(alpha_hist,relative_err_hist,'o')
        drawnow
    end
    j=j+1;
    if(relative_err_hist(1)*relative_err1<0)
        if(multiplier<0.01)
            flag1=1;
        else
            multiplier=multiplier/10;
            params.alpha=alpha2.*(1-multiplier);
        end
    end
        
    fprintf('%12s %12s %12.6g %12.6g %12.6g %12.6g\n','','',relative_err1,alpha1,newfactor,factor);
    fid=fopen(tmpfilepath,'a');
    fprintf(fid,'%12s %12s %12.6g %12.6g %12.6g %12.6g\n','','',relative_err1,alpha1,newfactor,factor);
    fclose(fid);
    
    if(flag1)
        break;
    end
    
    % If the previous result seems reasonable, use it for the initial guess
    if(nnz(results.rhovec./params.charge<params.minlevel)<numel(results.rhovec)/10)
        rhovecini=rho_vector(results.x,params);
    end
    
end %-------------------end-iteration for alpha-end------------------------

fprintf('%12.6g %12.6g %12.6g\n',results.params.alpha,factor,results.P_a_Gm);
fid=fopen(tmpfilepath,'a');
fprintf(fid,'%12.6g %12.6g %12.6g\n',results.params.alpha,factor,results.P_a_Gm);
fclose(fid);

save([inputparam.dirpath num2str(params.EF/params.charge) '_' num2str(params.barrier/params.charge) '_' num2str(params.thickness/0.382E-9) '_' num2str(params.alpha) '_results.mat'],'results');

if(params.display)
    figure()
    plot(params.xvec./params.charge',results.err')
    figure(f2)
    hold on
    plot(params.xvec./params.charge',results.L2vec')
    plot(params.xvec./params.charge',results.Jcalc')
    hold off
    figure(f1)
    plot(params.xvec./params.charge',results.rhovec./params.delta')
end

end

