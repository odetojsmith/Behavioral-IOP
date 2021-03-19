function sys_arx = space_to_arx(sys_state)
    n_x = size(sys_state.A,1);
    sys_arx = cell(2,n_x);   
    sys_arx_intm = tf(sys_state);   
    mat_u = cell2mat(sys_arx_intm.Numerator);
    mat_y = -cell2mat(sys_arx_intm.Denominator);
    for i = 1 : n_x
        sys_arx{1,i} = mat_y(:,i+1 : n_x+1 : end);
        sys_arx{2,i} = mat_u(:,i+1 : n_x+1 : end);
    end
end
