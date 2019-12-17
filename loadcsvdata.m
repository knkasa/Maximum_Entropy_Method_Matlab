function data = loadcsvdata( path )
    
    if(exist(path,'file')==0)
        error(['File ' path ' does not exist.']);
    end
    fid=fopen(path,'r');
    str=fgetl(fid);
    fclose(fid);
    if(isletter(str(1))||isletter(str(2)))
        tmp=csvread(path,1);
    else
        tmp=csvread(path);
    end
    
    data.x=tmp(:,1);
    data.y=tmp(:,2);
    
end

