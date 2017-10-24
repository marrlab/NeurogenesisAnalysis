function Z = cov2corr(Y,opt)
    Z=zeros(size(Y));
    idx_start=sum(opt.outVec~=0)+1;
    Z(:,1:idx_start) = Y(:,1:idx_start); %copy first order moments (means)
    k=sum(opt.outVec~=0);
    idx_var=[];
    %find indices for variances
    while k>0
        idx_var = [idx_var k+1];
        k=k-1;
    end
    idx_var = cumsum(idx_var);
    i=1;
    for idx = idx_start:size(Y,2)
        for j=1:length(idx_var)
            if j>i
                Z(:,idx)=Y(:,idx)./(Y(:,idx_var(i)).*Y(:,idx_var(j)));
            end
        end
        i=1+1;
    end
    
end