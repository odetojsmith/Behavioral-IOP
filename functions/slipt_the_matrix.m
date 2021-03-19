
function [U_p,U_f,Y_p,Y_f] = slipt_the_matrix(U_d,Y_d,dim,len)
    U_p = U_d(1: dim.u * (len.t_p-1) , :);
    U_f = U_d(dim.u * (len.t_p-1)+1 : end, :);
    Y_p = Y_d(1: dim.y * (len.t_p-1) , :);
    Y_f = Y_d(dim.y * (len.t_p-1)+1 : end, :);
end