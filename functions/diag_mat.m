

function A = diag_mat(sys_arx,k,i,rep_times)
    A_1 = sys_arx{i,k};
    ACell = repmat({A_1}, 1, rep_times);
    A = blkdiag(ACell{:});
end