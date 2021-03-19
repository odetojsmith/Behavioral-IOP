function P_plus = P_back_prop(Q,R,P,sys_state)
    A = sys_state.A;
    B = sys_state.B;
    trans = B' * P * A;
    P_plus = Q + A' * P * A - trans' * inv(R + B' * P * B) * trans;
end