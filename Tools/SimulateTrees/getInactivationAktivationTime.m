function t = getInactivationAktivationTime(SC_cyclingTime,mu1,mu2)
%#codegen
coder.inline('never')

% rn=unifrnd(0,1);
% t_inact = (-1/lambda2)*log(rn);
t_inact = exprnd(mu2,1,1);
    if (t_inact<=SC_cyclingTime)
%         rn2=random('unif',0,1,1,1);
%         rn2 = rand(1);
%         t = t_inact+(-1/lambda1).*log(rn2);
        t=t_inact+exprnd(mu1,1,1);
    else
        t=0;
    end
end