function par0 = LHS_constraint(par_guess,par_min,par_max,n_starts,constr,scale,scaleVec)

par0_exp=[];
par0=[];
%LH sampling
par_number = length(par_min);
stopp=1;

while stopp==1;
    theta_0 = [par_guess,...
              bsxfun(@plus,par_min,bsxfun(@times,par_max - par_min,...
                           lhsdesign(n_starts - size(par_guess,2),par_number,'smooth','off')'))];
    theta_0_exp= theta_0;
    if (strcmp(scale,'log') || strcmp(scale,'partly_log'))
        %re-transform: theta_0_exp is linear scale
           theta_0_exp(scaleVec,:)=exp(theta_0(scaleVec,:));
    end

    %check if constraint is fullfilled for starting values:
    if (~isempty(constr.A) && ~isempty(constr.b))
        A=constr.A;
        b=constr.b;
        for i=1:n_starts
            if sum(A*theta_0_exp(:,i)<=b)==length(b) %% all constraints fullfilled
                par0_exp=[par0_exp theta_0_exp(:,i)];
                if (size(par0_exp,2) ==n_starts)
                    stopp=0;
                    break;
                end
            end
        end
    else
        par0_exp=theta_0_exp;
        stopp=0;
    end  
    
end

par0= par0_exp;
if (strcmp(scale,'log') || strcmp(scale,'partly_log'))
       par0(scaleVec,:)=log(par0_exp(scaleVec,:));
end

end