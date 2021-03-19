function sys_state_onlystate = onlystate(sys_state,dim)
    C = eye(dim.x);
    D = zeros(dim.x,dim.u);
    sys_state_onlystate = ss(sys_state.A,sys_state.B,C,D,-1);
end