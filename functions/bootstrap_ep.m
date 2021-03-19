function ep_max = bootstrap_ep(N,u_input,y_output_data,dim,len,sigma,data_mat)

ep_bostr_a = [];
ep_bostr_b = [];
    for i = 1 : N
        noise_u = normrnd(0,sigma,size(u_input));
        u_input_data_with_noise = u_input - noise_u;
        noise_y = normrnd(0,sigma,size(y_output_data));
        y_output_data_with_noise = y_output_data + noise_y;
    % Simulation with input output noises

        [U_d_noisy,Y_d_noisy] = hankel_generation(len,u_input_data_with_noise,y_output_data_with_noise);
        data_mat_noisy = data_process(len,dim,U_d_noisy,Y_d_noisy,u_input_data_with_noise,y_output_data_with_noise);

        ep_bostr_a = [ep_bostr_a;norm(data_mat.T_mat-data_mat_noisy.T_mat)];
        ep_bostr_b = [ep_bostr_b;norm(data_mat.y_free-data_mat_noisy.y_free)];
    end
ep_a = prctile(ep_bostr_a,90);
ep_b = prctile(ep_bostr_b,90);
ep_max = max(ep_a, ep_b);
end