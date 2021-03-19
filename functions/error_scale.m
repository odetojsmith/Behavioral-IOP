function h = error_scale(ep,alpha,Y)
Y_norm = norm(Y);
h = ep^2*(2+alpha*Y_norm)^2 + 2* ep * norm(Y,'fro') * (2 + alpha * Y_norm);
end
