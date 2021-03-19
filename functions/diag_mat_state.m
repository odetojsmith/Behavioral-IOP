
function A = diag_mat_state(A,rep_times)
    ACell = repmat({A}, 1, rep_times);
    A = blkdiag(ACell{:});
end