
function [U_d,Y_d] = hankel_generation(len,u_data,y_data)
    U_d = [];
    Y_d = [];
    n_row_U_d = len.t_p+len.n_horizon;
    for i = 1 : len.t_d-n_row_U_d+1
        U_vec = [];
        Y_vec = [];
        for j = 1 : n_row_U_d
            U_vec = [U_vec; transpose(u_data(j+i-1,:))];
            Y_vec = [Y_vec; transpose(y_data(j+i-1,:))];
        end
        U_d = [U_d U_vec];
        Y_d = [Y_d Y_vec];
    end
end