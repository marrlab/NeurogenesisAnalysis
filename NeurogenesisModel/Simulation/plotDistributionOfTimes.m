function [T] = plotDistributionOfTimes(dataSet,rate_str,theta,optDistr,OPT)
n=100;

 for d=1:length(dataSet)
    %% plot distribution of times according to resulting rates:
    times_str = rate_str;
    idx = strwcmp(rate_str,'r_*'); 
    times_str=times_str(idx==1);
    times_str = strrep(times_str,'r','t');
    time_means=1./theta{d}(idx==1);
    t=zeros(length(time_means),n);
    a=length(time_means);
    for i=1:length(time_means)
        lambda=time_means(i);
        switch optDistr
            case 'exp'
                k=1;
                t(i,:) = exprnd(lambda,1,n);
                title_str=[optDistr,' distribution for ',times_str{i},'=\lambda=',num2str(lambda)];
                lambda_n = k/lambda;
                F=5;
                Y=0.01;
            case 'erlang'
                switch i
                    case 1
                        k=1;
                        F=5;
                        Y=0.01;
                    case 2
                        k=1;
                        F=5;
                        Y=0.01;
                    case 3
                        k=10;
                        F=3;
                        Y=0.06;
                    case 4
%                         k=3;
                        k=10;
                        F=3;
                        Y=0.006;
                    case 5
                        k=10;
                        F=3;
                        Y=0.004;
                end
                lambda_n = k/lambda; 
                rn_unif = unifrnd(0,1,k,n);
                t(i,:) = (-1/lambda_n)*log(prod(rn_unif,1));
                title_str=[optDistr,'(',num2str(k),',',times_str{i},') PDF'];
        end
        if k==1
            K = 'exponential';
        else
            K= 'kernel';
        end
        figure(2)
        subplot(1,a,i)
        x=0:0.1:F*lambda;
        switch d
            case 1
                plot(x,(lambda_n^k.*x.^(k-1).*exp(-lambda_n.*x))./factorial(k-1),'k');
                hold on;
                plot([lambda lambda],[0 Y],'k')
            case 2
                plot(x,(lambda_n^k.*x.^(k-1).*exp(-lambda_n.*x))./factorial(k-1),'k--');
                hold on;
                plot([lambda lambda],[0 Y],'k--')
        end
        title(title_str);
        if d==1
            axis([0,F*lambda,0,Y])
        end
    end
    T{d}=t;
end
saveFigs(OPT,['timeDistribution_',optDistr])

end