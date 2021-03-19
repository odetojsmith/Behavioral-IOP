function [min_val_opt,Phi_val_opt,gamma_opt] = golden_search(dim,opt,len)
    a=0.1 * opt.phiuynorm;                            % start of interval
    %b=min(0.95 * 1 / opt.ep,10);                            % end of interval
    b=min(0.95 * 1 / opt.ep,opt.alpha);
    epsilon=0.0000000001;               % accuracy value
    iter= 10;                       % maximum number of iterations
    tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
    k=0;                            % number of iterations

    x1=a+(1-tau)*(b-a);             % computing x values
    x2=a+tau*(b-a);

    [f_x1,~] = minimum_value_evaluation(dim,opt,len,x1);  
    f_x1 = (1 - opt.ep * x1)^(-1) * f_x1;
    [f_x2,~] = minimum_value_evaluation(dim,opt,len,x2); 
    f_x2 = (1 - opt.ep * x2)^(-1) * f_x2;

%     plot(x1,f_x1,'rx')              % plotting x
%     plot(x2,f_x2,'rx')
    while ((abs(b-a)>epsilon) && (k<iter))
        k=k+1;
        if(f_x1<f_x2)
            b=x2;
            x2=x1;
            x1=a+(1-tau)*(b-a);
            
            [f_x1,~] = minimum_value_evaluation(dim,opt,len,x1);  
            f_x1 = (1 - opt.ep * x1)^(-1) * f_x1;
            [f_x2,~] = minimum_value_evaluation(dim,opt,len,x2); 
            f_x2 = (1 - opt.ep * x2)^(-1) * f_x2;           
        else
            a=x1;
            x1=x2;
            x2=a+tau*(b-a);

            [f_x1,~] = minimum_value_evaluation(dim,opt,len,x1);  
            f_x1 = (1 - opt.ep * x1)^(-1) * f_x1;
            [f_x2,~] = minimum_value_evaluation(dim,opt,len,x2); 
            f_x2 = (1 - opt.ep * x2)^(-1) * f_x2;
        end
        k=k+1;
    end
    % chooses minimum point
    if(f_x1<f_x2)
        gamma_opt = x1;
        [min_val_opt,Phi_val_opt] = minimum_value_evaluation(dim,opt,len,x1);
    else
        gamma_opt = x2;
        [min_val_opt,Phi_val_opt] = minimum_value_evaluation(dim,opt,len,x2);
    end
end