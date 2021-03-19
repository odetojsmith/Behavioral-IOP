function [min_val,Phi_val] = minimum_value_evaluation(dim,opt,len,gamma)
Phi_yy = sdpvar(dim.y_horizon,dim.y_horizon,'full');
Phi_yu = sdpvar(dim.y_horizon,dim.u_horizon,'full');
Phi_uy = sdpvar(dim.u_horizon,dim.y_horizon,'full');
Phi_uu = sdpvar(dim.u_horizon,dim.u_horizon, 'full');
Phi = [Phi_yy Phi_yu;Phi_uy Phi_uu];

Mat = [sqrt(1+opt.scaler_1+opt.scaler_2)*Phi_yy, Phi_yu, Phi_yy * opt.y_free;
       sqrt(1+opt.scaler_2)*Phi_uy, Phi_uu, Phi_uy * opt.y_free];
objective  = norm(Mat,'fro');
constraint = [[opt.I_tot_y, -opt.G] * Phi == [opt.I_tot_y opt.zero_yu]];
constraint = [constraint, Phi * [-opt.G;opt.I_tot_u] == [opt.zero_yu;opt.I_tot_u]];
for i = 1 : len.n_horizon-1
    constraint = [constraint, Phi_yu((i-1) * dim.y+1 : i * dim.y , i * dim.u+1 : end) == zeros(dim.y,(len.n_horizon-i)*dim.u),...
                              Phi_yy((i-1) * dim.y+1 : i * dim.y , i *dim.y+1 : end) == zeros(dim.y,(len.n_horizon+1-i)*dim.y),...
                              Phi_uu((i-1) * dim.u+1 : i * dim.u , i *dim.u+1 : end) == zeros(dim.u,(len.n_horizon-i)*dim.u),...
                              Phi_uy((i-1) * dim.u+1 : i * dim.u , i *dim.y+1 : end) == zeros(dim.u,(len.n_horizon+1-i)*dim.y)];
end
constraint = [constraint, Phi_yy((len.n_horizon-1) * dim.y+1 : len.n_horizon * dim.y , len.n_horizon *dim.y+1 : end) == zeros(dim.y,dim.y),...
                          Phi_uy((len.n_horizon-1) * dim.y+1 : len.n_horizon * dim.y , len.n_horizon *dim.y+1 : end) == zeros(dim.u,dim.y),];
constraint = [constraint, norm(Phi_uy) <= gamma, norm(Phi_uy) <= opt.alpha];
options = sdpsettings('solver','mosek','verbose',0);
sol = optimize(constraint,objective,options);
Phi_val.Phi_yu_val = value(Phi_yu);
Phi_val.Phi_yy_val = value(Phi_yy);
Phi_val.Phi_uu_val = value(Phi_uu);
Phi_val.Phi_uy_val = value(Phi_uy);
min_val = value(objective);
end