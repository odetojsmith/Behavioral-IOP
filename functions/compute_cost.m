function cost = compute_cost(LR,Sigma,K,G)
    Closed_Loops = [inv(eye(size(G,1))-G*K) inv(eye(size(G,1))-G*K)*G;
                    K*inv(eye(size(G,1))-G*K) inv(eye(size(K,1))-K*G)];
    Matrix_cost = LR * Closed_Loops *  Sigma;
    cost = norm(Matrix_cost,'fro');
end
