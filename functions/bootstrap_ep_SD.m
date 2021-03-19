function ep_max = bootstrap_ep_SD(N,u_input,y_output_data,dim,len,sigma,data_mat_noisy)

ep_bostr_a = [];
ep_bostr_b = [];
y_output_bfG = sim_from_impulse(data_mat_noisy,u_input,dim);
    for i = 1 : N
        noise_u = normrnd(0,sigma,size(u_input));
        u_input_data_with_noise = u_input - noise_u;
        noise_y = normrnd(0,sigma,size(y_output_data));
        y_output_data_with_noise = y_output_bfG + noise_y;
        
        [U_d_noisy,Y_d_noisy] = hankel_generation(len,u_input_data_with_noise,y_output_data_with_noise);
        data_mat_noisy_noisy = data_process(len,dim,U_d_noisy,Y_d_noisy,u_input_data_with_noise,y_output_data_with_noise);
        
        
        ep_bostr_a = [ep_bostr_a;norm(data_mat_noisy_noisy.bfG-data_mat_noisy.bfG)];
        ep_bostr_b = [ep_bostr_b;norm(data_mat_noisy_noisy.y_free-data_mat_noisy.y_free)];
    end
ep_a = prctile(ep_bostr_a,90);
ep_b = prctile(ep_bostr_b,90);
ep_max = max(ep_a, ep_b);

end