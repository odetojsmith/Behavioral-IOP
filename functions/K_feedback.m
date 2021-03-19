function K = K_feedback(Q,R,P,sys_state)
    A = sys_state.A;
    B = sys_state.B;
    trans = B' * P * A;
    K = - inv(R + B' * P * B) * trans;
end