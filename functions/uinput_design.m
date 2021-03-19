function u_input = uinput_design(sys_state,dim,len,x_ini_fix,u_input)

    sys_state_onlystate = onlystate(sys_state,dim);
    x_state = lsim(sys_state_onlystate,u_input,len.t,zeros(dim.x,1));
    x_before_control = transpose(x_state(len.t_opt_start-len.control_len, :));
    sys_A = sys_state.A;
    sys_B = sys_state.B;
    control_mat = sys_B;
    for i = 2 : len.control_len
        sys_B = sys_A * sys_B;
        control_mat = [sys_B control_mat];
    end
    u_vec = pinv(control_mat) * (x_ini_fix-(sys_A)^len.control_len * x_before_control);
    u_mat = transpose(reshape(u_vec,dim.u,length(u_vec)/dim.u));
    u_input(len.t_opt_start - len.control_len:len.t_opt_start-1,:) = u_mat;
