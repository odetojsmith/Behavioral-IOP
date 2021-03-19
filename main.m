clear all
close all
addpath('./functions')
addpath('./scripts')

% Method control
smm_icon = 0; % 0 with least square BM prediction; 1 with SMM BM prediction

% Load the simulation system
sig = 0.0001;
len = set_lengths();
[dim,sys_state] = sys_gen(len,0.9);
len.control_len = dim.x+3;

% Simulation conditions (for comparison)
x_ini_fix = [1;-1];
rho_list = [0.4 0.5 0.6 0.7 0.8 0.9 0.99];
howmany_systems = length(rho_list);
sigma_list = 10*10^(-4)*(0:100/9:100);
howmany_sigmas = length(sigma_list);

% Create space for lists 
J_hat_list = zeros(howmany_sigmas,howmany_systems);
ep_max_list = zeros(howmany_sigmas,howmany_systems);
J_delta_ratio = zeros(howmany_sigmas,howmany_systems);
J_star_list = zeros(howmany_sigmas,howmany_systems);

% Load the fixed input
if isfile('u_input.mat')
    load u_input.mat
else
    u_input = normrnd(0,1,len.n_length_sim,dim.u);
    save('u_input.mat','u_input');
end

for plot_number = 1:howmany_systems
    % Generate system
    rho = rho_list(plot_number);
    [dim,sys_state] = sys_gen(len,rho);
    % Check controllability and observability
    sanity_check(sys_state);
   
    % Generate appropriate control input to enforce a desired x_ini (for fair comparisons purposes)
    u_input = uinput_design(sys_state,dim,len,x_ini_fix,u_input);

    
    % Simulate the systems to get the states and outputs
    sys_state_onlystate = onlystate(sys_state,dim);
    x_state = lsim(sys_state_onlystate,u_input,len.t,zeros(dim.x,1));
    x_ini = x_state(len.t_opt_start, :);
    x_ini = x_ini';
    y_output_vec = lsim(sys_state,u_input,len.t,zeros(dim.x,1));
    y_output_data = y_output_vec;
    
    %% Calculate the non-robust data-driven SLS problem and compare with LQR
    % Genera the Hankel Matrix and then the Impulse Matrice/y_free
    [U_d,Y_d] = hankel_generation(len,u_input,y_output_data);                  % This is using all data, not only historical?
    data_mat = data_process(len,dim,U_d,Y_d,u_input,y_output_data);
    
    % Define the properties of the optimization problem and the noises
    [opt,noise] = opt_and_noise_less_u(dim,data_mat,len,sig);
    % Derive the corresponding matrices in the optimization problem
    Phi_val = opt_non_robust(dim,len,data_mat,opt);
    
    % system simulation and controller implementation
    x_real = zeros(dim.x,len.n_horizon+1);
    x_real(:,1) = x_ini;
    y_ini = sys_state.C * x_ini + transpose(noise.v(1,:));
    y_real = y_ini;
    
    K = Phi_val.Phi_uy_val * inv(Phi_val.Phi_yy_val);
    y_real_data = transpose(y_ini);
    y_vec = y_ini;
    y_vec = [y_vec;zeros(len.n_horizon*dim.y,1)];
    u_real_data = transpose(K(1:dim.u,1:dim.y) * y_ini);
    
    for i = 2 : len.n_horizon+1
        x_real(:,i) = sys_state.A * x_real(:,i-1) + sys_state.B * (transpose(u_real_data(i-1,:)) + transpose(noise.w(i-1,:)));
        y_real_data(i,:) = transpose(sys_state.C * x_real(:,i)) + noise.v(i,:);
        y_vec((i-1)*dim.y+1 : i*dim.y,1) = transpose(y_real_data(i,:));
        if i < len.n_horizon+1
            u_i_vec = K((i-1)*dim.u+1 : i*dim.u , :) * y_vec;
            u_real_data = [u_real_data;transpose(u_i_vec)];
        end
    end
    J_star = compute_cost(opt.LR,opt.Sigma,K,data_mat.bfG);

    % LQG implementation and plot of comparison
    y_real_data_lqg = LQG_imp(len,dim,sys_state,opt,noise,x_ini);
%     if plot_number == howmany_systems
%         plots_1;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Robust formulation
    N = 2000; % Bootstrap numbers
    J_hat_non_robust_list = []; 
    sigma_ok_list = [];
    norm_bound = norm(Phi_val.Phi_uy_val); % This norm is used to derive alpha
   
    for i = 1 : howmany_sigmas
        fprintf('Iteration %d with sigma = %f\n',i,sigma_list(i))
        sigma = sigma_list(i);
        
        % If one wants strict comparison with the bound, she needs to apply
        % this bound check. See Theorem 3
        %ep_max = bootstrap_ep(N,u_input,y_output_data,dim,len,sigma,data_mat);
        %     if ep_max > 1/5/norm_bound
        %         continue
        %     end

        % apply the noises on the past data (historical and past)
        noise_u = normrnd(0,sigma,size(u_input));
        u_input_data_with_noise = u_input - noise_u;
        noise_y = normrnd(0,sigma,size(y_output_data));
        y_output_data_with_noise = y_output_data + noise_y;
        
        % Generate the matrices for the optimization problem
        [U_d_noisy,Y_d_noisy] = hankel_generation(len,u_input_data_with_noise,y_output_data_with_noise);
        if smm_icon == 0 
            data_mat_noisy = data_process(len,dim,U_d_noisy,Y_d_noisy,u_input_data_with_noise,y_output_data_with_noise);
        else
            if sigma == 0
                data_mat_noisy.y_free = data_mat.y_free;
                data_mat_noisy.bfG = data_mat.bfG;
            else
                data_mat_noisy.y_free = SMM_Predict_y_free(u_input_data_with_noise,y_output_data_with_noise,dim,len,sigma);
                data_mat_noisy.bfG = SMM_Predict_impulse(u_input_data_with_noise,y_output_data_with_noise,dim,len,sigma);
            end
        end
%         error_bfG = norm(data_mat_noisy.bfG-data_mat.bfG);
%         error_y_free = norm(data_mat_noisy.y_free-data_mat.y_free);
%         ep_max_list(i,plot_number) = max(error_bfG,error_y_free);
       
                % run bootstarp
        sigma_ok_list = [sigma_ok_list sigma];
        if smm_icon == 1 
            ep_max_list(i,plot_number) = bootstrap_ep_smm(N,u_input,y_output_data,dim,len,sigma,data_mat);
        else
            ep_max_list(i,plot_number) = bootstrap_ep(N,u_input,y_output_data,dim,len,sigma,data_mat);
        end
        
        
        opt.phiuynorm = norm(Phi_val.Phi_uy_val);
        opt.ep = ep_max_list(i,plot_number);
        opt.alpha = 3*norm_bound;
        opt.G = data_mat_noisy.bfG;
        opt.y_free = data_mat_noisy.y_free;
        opt.scaler_1 = error_scale(opt.ep,opt.alpha,opt.G);
        opt.scaler_2 = error_scale(opt.ep,opt.alpha,opt.y_free);
        % Golden Search for the optimal solution
        [min_val_opt,Phi_val_noisy,gamma_opt] = golden_search(dim,opt,len);

        %%%Compute Khatstar
        K_hatstar = Phi_val_noisy.Phi_uy_val*inv(Phi_val_noisy.Phi_yy_val);
        %%%Compute costs
        J_hat_list(i,plot_number) = compute_cost(opt.LR,opt.Sigma,K_hatstar,data_mat.bfG)
        J_star_list(i,plot_number) = J_star
        J_delta_ratio(i,plot_number) = (J_hat_list(i,plot_number)^2- J_star^2)/J_star^2
        ep_max_list
    end
 
    plots_2;
    %%
    figure(2)
    subplot(111)
    ylabel('$\epsilon$','Interpreter','latex','FontSize',20);
    xlabel('$\sigma$','Interpreter','latex','FontSize',20);
    title('Noise to estimation error','FontSize',15) 
    legend('$\rho$ = 0.4','$\rho$ = 0.5','$\rho$ = 0.6','$\rho$ = 0.7','$\rho$ = 0.8','$\rho$ = 0.9','$\rho$ = 0.99','Interpreter','latex')
    subplot(122)
    xlabel('$\epsilon$','Interpreter','latex','FontSize',20);
    title('Suboptimality','FontSize',15) 
    legend('$\rho$ = 0.4','$\rho$ = 0.5','$\rho$ = 0.6','$\rho$ = 0.7','$\rho$ = 0.8','$\rho$ = 0.9','$\rho$ = 0.99','Interpreter','latex')

    

    
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



