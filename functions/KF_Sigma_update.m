
function [Sigma_KF_forward,L] = KF_Sigma_update(Sigma_KF_forward,sys,W,V)
    C = sys.C;
    A = sys.A;
    mid = inv(C * Sigma_KF_forward * C'+V);
    trans = C * Sigma_KF_forward;
    Sigma_KF = Sigma_KF_forward - trans' * mid * trans;
    Sigma_KF_forward = A * Sigma_KF * A' + W;
    mid = inv(C * Sigma_KF_forward * C'+V);
    trans = C * Sigma_KF_forward;
    L =  trans' * mid; 
end