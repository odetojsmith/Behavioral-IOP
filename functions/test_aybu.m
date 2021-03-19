function res = test_aybu(A_y,B_u,y_output,u_input,dim,len)

    u_vec = u_input(len.t_opt_start-dim.x+1:len.t_opt_start+len.n_horizon);
    y_vec = y_output(len.t_opt_start-dim.x+1:len.t_opt_start+len.n_horizon);
    
    res = A_y * y_vec + B_u * u_vec;
end
