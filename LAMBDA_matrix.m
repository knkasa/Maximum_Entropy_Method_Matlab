function matrix = LAMBDA_matrix( rhovec,KKmat,params )
    
   % calculating lambda which is simply KKmat(i,j)

    asqvec=sqrt(rhovec);
    matrix=zeros(numel(rhovec));
    for i=1:size(rhovec)
        for j=1:size(rhovec)
            matrix(i,j)=asqvec(i)*KKmat(i,j)*asqvec(j);
        end
    end
    
end

