function  y_real_data_lqg = LQG_imp(len,dim,sys_state,opt,noise,x_ini)

LQR_K = cell(len.n_horizon,1);
LQR_P = cell(len.n_horizon+1,1);
y_location_f = len.n_horizon * dim.y + 1 : (len.n_horizon+1) * dim.y;
y_location = 1 : dim.y;
u_location = 1 : dim.u;
LQR_P{len.n_horizon+1,1} = transpose(sys_state.C) * opt.L(y_location_f,y_location_f) * sys_state.C;
Q = opt.L(y_location,y_location);
Q = transpose(sys_state.C) *  Q * sys_state.C;
R = opt.R(u_location,u_location);
for j = 1:len.n_horizon
    i = len.n_horizon+1-j;
    LQR_P{i,1} = P_back_prop(Q,R,LQR_P{i+1,1},sys_state);
    LQR_K{i,1} = K_feedback(Q,R,LQR_P{i+1,1},sys_state);
end

x_real_lqg = zeros(dim.x,len.n_horizon+1);
x_real_lqg(:,1) = x_ini;
y_ini = sys_state.C * x_ini + transpose(noise.v(1,:));
y_real_data_lqg = transpose(y_ini);
u_lqg_data = zeros(len.n_horizon,dim.u);
u_lqg_data(1,:) = transpose(LQR_K{1,1} * x_ini);
x_est = x_ini;
x_real_lqg(:,2) = sys_state.A * x_ini + sys_state.B * (transpose(u_lqg_data(1,:)) + transpose(noise.w(1,:)));
y_real_data_lqg = [y_real_data_lqg;transpose(sys_state.C * x_real_lqg(:,2)) + noise.v(2,:)];

W = sys_state.B * (noise.sigma_w)^2 * transpose(sys_state.B);
V = (noise.sigma_v)^2 * eye(dim.y);
Sigma_KF_forward = zeros(dim.x,dim.x);
for i = 2 : len.n_horizon
    [Sigma_KF_forward,L] = KF_Sigma_update(Sigma_KF_forward,sys_state,W,V);
    e = transpose(y_real_data_lqg(i,:)) - sys_state.C * (sys_state.A * x_est + sys_state.B * transpose(u_lqg_data(i-1,:)));
    x_est = sys_state.A * x_est + sys_state.B * transpose(u_lqg_data(i-1,:)) + L * e;
    u_lqg_data(i,:) = transpose(LQR_K{i,1} * x_est);
    x_real_lqg(:,i+1) = sys_state.A * x_real_lqg(:,i)  + sys_state.B * (transpose(u_lqg_data(i,:)) + transpose(noise.w(i,:)));
    y_real_data_lqg = [y_real_data_lqg;transpose(sys_state.C * x_real_lqg(:,i+1)) + noise.v(i,:)];
end