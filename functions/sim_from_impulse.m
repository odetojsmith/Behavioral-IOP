function y_output_bfG_data = sim_from_impulse(data_mat_noisy,u_input,dim)
    n = size(u_input,1);
    imp_mat = data_mat_noisy.bfG(:,1:dim.u);
    m = size(imp_mat,1);
    y_output_bfG_data = zeros(n+m,dim.y);
    for i = 1:n
        y_u_1 = imp_mat(:,1) * u_input(i,1);
        y_u_2 = imp_mat(:,2) * u_input(i,2);
        y_u = y_u_1+y_u_2;
        y_u = transpose(reshape(y_u,dim.y,m/dim.y));
        y_output_bfG_data(i:i+m/dim.y-1,:) = y_output_bfG_data(i:i+m/dim.y-1,:)+y_u;
    end
    y_output_bfG_data = y_output_bfG_data(1:n,:);
end