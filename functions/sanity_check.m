
function sanity_check(sys_state)
    C_matrix = [sys_state.B sys_state.A * sys_state.B (sys_state.A)^2 * sys_state.B];
    O_matrix = [sys_state.C; sys_state.C * sys_state.A; sys_state.C * (sys_state.A)^2];
    c_eig = abs(svds(C_matrix));
    o_eig = abs(svds(O_matrix));
    if min(c_eig)<0.01 || min(o_eig)<0.01
        error("not observable or controllable")
    end
end