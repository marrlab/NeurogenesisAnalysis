function [t] = drawTime(lambda,optDistr,time_str)

    switch optDistr
        case 'exp'
            t = exprnd(lambda,1);
        case 'erlang'
            switch time_str
                case {'t_act','t_inact'}
                    k=1;
                case 't_div'
                    k=7;
                case 't_B_N'
                    k=4;
                case 't_death'
                    k=10;
            end
            lambda_n = k/lambda; 
            rn_unif = unifrnd(0,1,k,1);
            t = (-1/lambda_n)*log(prod(rn_unif,1));
    end

end