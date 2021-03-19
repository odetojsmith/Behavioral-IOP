function data_mat = data_process(len,dim,U_d,Y_d,u_input,y_output_data)
    [U_p,U_f,Y_p,Y_f] = slipt_the_matrix(U_d,Y_d,dim,len);
    H = [U_p;Y_p;U_f];
    u_recent_p = reshape(transpose(u_input(len.t_opt_start - len.t_p+1 : len.t_opt_start-1,:)),[(len.t_p-1)*dim.u 1]);
    u_recent_f = zeros(dim.u * (len.n_horizon+1),1);
    y_recent_p = reshape(transpose(y_output_data(len.t_opt_start - len.t_p+1 : len.t_opt_start-1,:)),[(len.t_p-1)*dim.y 1]);
    % G_O = pinv(H) * [H_0]
    % Compute a set of generators for the free responses of the system
    data_mat.y_free = Y_f * pinv(H) * [u_recent_p;y_recent_p;u_recent_f];
    % Full rank decomposition by means of SVD

    G_T = pinv([U_p; Y_p; U_f])*[zeros(dim.u*(len.t_p-1), dim.u); zeros(dim.y*(len.t_p-1), dim.u); eye(dim.u); zeros(dim.u*(len.n_horizon), dim.u)];
    T_ori = Y_f * G_T;
    
    Z = z_gen_y(dim,len.n_horizon);
    Z_exp = Z;
    T_mat = T_ori;
    for i = 1 : len.n_horizon-1
        T_mat = [T_mat Z_exp*T_ori];
        Z_exp = Z_exp*Z;
    end
    data_mat.T_mat = T_mat;
    data_mat.bfG = data_mat.T_mat;
end