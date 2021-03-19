function len = set_lengths()
    len.n_horizon = 10;
    len.n_length_sim = 300;
    len.t_opt_start = 250;
    len.t_d = 200;
    len.t_p = 30;
    t = 1:len.n_length_sim;
    len.t = t'-1;
end