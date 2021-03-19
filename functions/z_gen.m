function Z = z_gen(dim,n_horizon)
    n = 1+n_horizon;
    Z_1 = eye((n-1) * dim.x);
    Z_1 = [zeros(dim.x,(n-1) * dim.x);Z_1];
    Z = [Z_1 zeros(n* dim.x,dim.x)];
end