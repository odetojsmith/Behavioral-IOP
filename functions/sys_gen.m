function [dim,sys_state] = sys_gen(len,rho)
    A = rho*[0.8 0.4;0.8 -0.6];
%     B = [1;2;1.5]; 
%     C = [1 1 1];
     B = [1 0.2;2 0.3]; 
     C = [1 1;0.7 0.2];
    D = zeros(size(C,1),size(B,2));
    sys_state = ss(A,B,C,D,-1);
    dim.x = size(sys_state.A,1);
    dim.u = size(sys_state.B,2);
    dim.y = size(sys_state.C,1);
    dim.y_horizon = dim.y * (len.n_horizon+1);
    dim.u_horizon = dim.u * (len.n_horizon);
end
