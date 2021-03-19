function impulse_mat_tot = SMM_Predict_impulse(u_input,y_output,dim,len,sigma)
    impulse_mat = [];
    u_recent_p = zeros(dim.u * (len.t_p-1),1);
    y_recent_p = zeros(dim.y * (len.t_p-1),1);
% Simulation with input output noises
    [U_d_noisy,Y_d_noisy] = hankel_generation(len,u_input,y_output);
    [U_p,U_f,Y_p,Y_f] = slipt_the_matrix(U_d_noisy,Y_d_noisy,dim,len);
    H1 = [U_p;U_f;Y_p];
    H2 = Y_f;
    Yp = Y_p;
    yvar = sigma^2;
    ypvar = yvar;
    M = size(Y_f,2);
    LL = len.n_horizon;
    L0 = len.t_p;
    U = U_d_noisy;
    Y = Y_d_noisy;
    u0 = u_recent_p;
    y0 = y_recent_p;
    
    niter 	= 5;    
    for i_iter = 1:dim.u
        u_recent_f = zeros(dim.u * (len.n_horizon+1),1);
        u_recent_f(i_iter,1) = 1;
        util = [u_recent_p;u_recent_f;y_recent_p];
        uref = u_recent_f;
        
        ghat1 	= zeros(size(Yp,2),niter+1);
        ghat1(:,1) = H1'/(H1*H1')*util;         % Initialize at the pseudo inverse soln
        for i = 1:niter
            lambda = LL*(yvar*norm(ghat1(:,i))^2+ypvar)/norm(ghat1(:,i))^2+L0*yvar;
            F = lambda*eye(M)+Yp'*Yp;
            % Closed-form
            ghat1(:,i+1) = (inv(F)-F\U'/(U/F*U')*U/F)*Yp'*y0+F\U'/(U/F*U')*[u0;uref];
        end
        yhatres_1 = H2*ghat1(:,end);
        y_impulse = yhatres_1;
        impulse_mat = [impulse_mat y_impulse];
        u_recent_f(i_iter,1) = 0;
    end
    T_ori = impulse_mat;
    Z = z_gen_y(dim,len.n_horizon);
    Z_exp = Z;
    T_mat = T_ori;
    for i = 1 : len.n_horizon-1
        T_mat = [T_mat Z_exp*T_ori];
        Z_exp = Z_exp*Z;
    end
    impulse_mat_tot = T_mat;
%%%%%%%%%%%%%%%%%%%
                            % No iteration is needed for impulse response est
    
end