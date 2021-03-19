function y_free = SMM_Predict_y_free(u_input,y_output,dim,len,sigma)
    u_recent_p = reshape(transpose(u_input(len.t_opt_start - len.t_p+1 : len.t_opt_start-1,:)),[(len.t_p-1)*dim.u 1]);
    u_recent_f = zeros(dim.u * (len.n_horizon+1),1);
    y_recent_p = reshape(transpose(y_output(len.t_opt_start - len.t_p+1 : len.t_opt_start-1,:)),[(len.t_p-1)*dim.y 1]);
% Simulation with input output noises
    [U_d_noisy,Y_d_noisy] = hankel_generation(len,u_input,y_output);
    [U_p,U_f,Y_p,Y_f] = slipt_the_matrix(U_d_noisy,Y_d_noisy,dim,len);
    H1 = [U_p;U_f;Y_p];
    H2 = Y_f;
    Yp = Y_p;
    yvar = sigma^2;
    ypvar = yvar;
    util = [u_recent_p;u_recent_f;y_recent_p];
    M = size(Y_f,2);
    LL = len.n_horizon;
    L0 = len.t_p;
    U = U_d_noisy;
    Y = Y_d_noisy;
    u0 = u_recent_p;
    y0 = y_recent_p;
    uref = u_recent_f;

%%%%%%%%%%%%%%%%%%%
    niter 	= 5;                            % No iteration is needed for impulse response est
    ghat1 	= zeros(size(Yp,2),niter+1);
    ghat1(:,1) = H1'/(H1*H1')*util;         % Initialize at the pseudo inverse soln
    for i = 1:niter
        lambda = LL*(yvar*norm(ghat1(:,i))^2+ypvar)/norm(ghat1(:,i))^2+L0*yvar;
        F = lambda*eye(M)+Yp'*Yp;
        % Closed-form
        ghat1(:,i+1) = (inv(F)-F\U'/(U/F*U')*U/F)*Yp'*y0+F\U'/(U/F*U')*[u0;uref];
    end
    yhatres_1 = H2*ghat1(:,end);
    y_free = yhatres_1;
end