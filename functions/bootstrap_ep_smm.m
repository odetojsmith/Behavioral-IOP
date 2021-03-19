function ep_max = bootstrap_ep_smm(N,u_input,y_output_data,dim,len,sigma,data_mat)
if sigma == 0
    ep_max = 0;
    return
end
ep_bostr_a = [];
ep_bostr_b = [];
    for i_boot = 1 : N
        noise_u = normrnd(0,sigma,size(u_input));
        u_input_data_with_noise = u_input - noise_u;
        noise_y = normrnd(0,sigma,size(y_output_data));
        y_output_data_with_noise = y_output_data + noise_y;
        
        data_mat_noisy.y_free = SMM_Predict_y_free(u_input_data_with_noise,y_output_data_with_noise,dim,len,sigma);
        data_mat_noisy.bfG = SMM_Predict_impulse(u_input_data_with_noise,y_output_data_with_noise,dim,len,sigma);
    %%%%%%%%%%%%%%%%%%%

        ep_bostr_a = [ep_bostr_a;norm(data_mat.bfG-data_mat_noisy.bfG)];
        ep_bostr_b = [ep_bostr_b;norm(data_mat.y_free-data_mat_noisy.y_free)];
    end
ep_a = prctile(ep_bostr_a,90);
ep_b = prctile(ep_bostr_b,90);
ep_max = max(ep_a, ep_b);
end