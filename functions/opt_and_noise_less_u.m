function [opt,noise] = opt_and_noise_less_u(dim,data_mat,len,ep)

% Opt: configuration about the optimization system
% Noise parameters and realizations
noise.sigma_w_nu = ep;
noise.sigma_v_nu = ep;
noise.sigma_w = 1;
noise.sigma_v = 1;
L_weight = 1;
R_weight = 1;


noise.Sigma_w = diag_mat_state((noise.sigma_w)^2 * eye(dim.u),len.n_horizon);
noise.Sigma_v = diag_mat_state((noise.sigma_v)^2 * eye(dim.y),len.n_horizon+1);


noise.w = normrnd(0,noise.sigma_w,len.n_horizon,dim.u);
noise.v = normrnd(0,noise.sigma_v,len.n_horizon+1,dim.y); 
noise.w_nu = normrnd(0,noise.sigma_w_nu,len.n_horizon,dim.u);
noise.v_nu = normrnd(0,noise.sigma_w_nu,len.n_horizon+1,dim.y); 

opt.L = kron(eye(len.n_horizon),eye(dim.y));
opt.L = blkdiag(opt.L,L_weight*eye(dim.y));

opt.R = R_weight * diag_mat_state(eye(dim.u),len.n_horizon);
LR_cell = {sqrtm(opt.L),sqrtm(opt.R)};
opt.LR = blkdiag(LR_cell{:});
opt.I_tot_u = diag_mat_state(eye(dim.u),len.n_horizon);
opt.I_tot_y = diag_mat_state(eye(dim.y),len.n_horizon+1);
opt.zero_yu = zeros(dim.y * (len.n_horizon+1), dim.u * len.n_horizon);

% opt.L_uni = diag_mat_state(eye(dim.y),len.n_horizon+1);
% opt.R_uni = diag_mat_state(eye(dim.u),len.n_horizon);
% LR_cell_uni = {sqrtm(opt.L_uni),sqrtm(opt.R_uni)};  
% opt.LR_uni = blkdiag(LR_cell_uni{:}); 


Sigma_cell = {sqrtm(noise.Sigma_v),sqrtm(noise.Sigma_w)};
Sigma = blkdiag(Sigma_cell{:});
y_free_vec = zeros(size(Sigma,1),1);
y_free_vec(1 : dim.y*(len.n_horizon+1)) = data_mat.y_free;
opt.Sigma = [Sigma y_free_vec];

end
