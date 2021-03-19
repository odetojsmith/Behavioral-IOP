
function Z = z_gen_y(dim,n_horizon)
    n = 1+n_horizon;
    Z_1 = eye((n-1) * dim.y);
    Z_1 = [zeros(dim.y,(n-1) * dim.y);Z_1];
    Z = [Z_1 zeros(n* dim.y,dim.y)];
end
