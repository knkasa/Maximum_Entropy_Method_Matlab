function stop = optimmyplot(~,optimValues,state,varargin,params,filepath)

stop = false;

switch state
    case 'init'
        fid=fopen(filepath,'a');
        fprintf(fid,'%12s %12s\n','IterNum','Fval');
        fclose(fid);
    case 'iter'
        if(mod(optimValues.iteration,10)==0)
            fid=fopen(filepath,'a');
            fprintf(fid,'%12g %12.6g\n',optimValues.iteration,optimValues.fval'*optimValues.fval);
            fclose(fid);
        end
end

end
