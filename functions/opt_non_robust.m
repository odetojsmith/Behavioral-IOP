function [Phi_val] = opt_non_robust(dim,len,data_mat,opt)


% Define the optimization variables

Phi_yy = sdpvar(dim.y_horizon,dim.y_horizon,'full');
Phi_yu = sdpvar(dim.y_horizon,dim.u_horizon,'full');
Phi_uy = sdpvar(dim.u_horizon,dim.y_horizon,'full');
Phi_uu = sdpvar(dim.u_horizon,dim.u_horizon, 'full');
Phi = [Phi_yy Phi_yu;Phi_uy Phi_uu];

% Define the optimization problem
J_mat = opt.LR * Phi *  opt.Sigma;
objective  = norm(J_mat,'fro');
constraint = [[opt.I_tot_y, -data_mat.bfG] * Phi == [opt.I_tot_y opt.zero_yu]];
constraint = [constraint, Phi * [-data_mat.bfG;opt.I_tot_u] == [opt.zero_yu;opt.I_tot_u]];
for i = 1 : len.n_horizon-1
    constraint = [constraint, Phi_yu((i-1) * dim.y+1 : i * dim.y , i * dim.u+1 : end) == zeros(dim.y,(len.n_horizon-i)*dim.u),...
                              Phi_yy((i-1) * dim.y+1 : i * dim.y , i *dim.y+1 : end) == zeros(dim.y,(len.n_horizon+1-i)*dim.y),...
                              Phi_uu((i-1) * dim.u+1 : i * dim.u , i *dim.u+1 : end) == zeros(dim.u,(len.n_horizon-i)*dim.u),...
                              Phi_uy((i-1) * dim.u+1 : i * dim.u , i *dim.y+1 : end) == zeros(dim.u,(len.n_horizon+1-i)*dim.y)];
end
constraint = [constraint, Phi_yy((len.n_horizon-1) * dim.y+1 : len.n_horizon * dim.y , len.n_horizon *dim.y+1 : end) == zeros(dim.y,dim.y),...
                          Phi_uy((len.n_horizon-1) * dim.y+1 : len.n_horizon * dim.y , len.n_horizon *dim.y+1 : end) == zeros(dim.u,dim.y),];
                      


% Calculate the optimizer
options = sdpsettings('solver','mosek');
sol = optimize(constraint,objective,options);
Phi_val.Phi_yu_val = value(Phi_yu);
Phi_val.Phi_yy_val = value(Phi_yy);
Phi_val.Phi_uu_val = value(Phi_uu);
Phi_val.Phi_uy_val = value(Phi_uy);
end