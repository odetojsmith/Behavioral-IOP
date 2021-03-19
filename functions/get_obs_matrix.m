function [O, x_ini, F] = get_obs_matrix(Ud, Yd, Tini, Tf, n)
    
    % Complexity parameters
    m = size(Ud, 1); p = size(Yd, 1);
    % Length of the collected historical data
    T = min(size(Ud, 2), size(Yd, 2));
    % Number of block rows and columns in the Hankel matrices
    L = Tini + Tf; M = T - L + 1;
    
    % Arrange the historical data into the corresponding Hankel matrices
    Hu = zeros(m*L, M); Hy = zeros(p*L, M); % Pre-allocate memory
    for i = 1:M % Fill the Hankel matrices one column at a time
        tmp = Ud(:, i:(i + L - 1));
        Hu(:, i) = tmp(:); % Vectorize
        tmp = Yd(:, i:(i + L - 1));
        Hy(:, i) = tmp(:);
    end
    % Split the Hankel matrices into the blocks Up, Yp, Uf, Yp
    Up = Hu(1:m*Tini, :); Yp = Hy(1:p*Tini, :); 
    Uf = Hu((m*Tini + 1):m*L, :); Yf = Hy((p*Tini + 1):p*L, :); 
%     % Necessary and sufficient for using an image representation
%     if (rank([Hu; Hy]) < m*L + n)
%         disp('Image representation is not valid');
%     end
    % Compute the least-squares solution for G
    G = pinv([Up; Yp; Uf])*[Up; Yp; zeros(m*Tf, M)];
    % Compute a set of generators for the free responses of the system
    F = Yf*G;
    % Full rank decomposition by means of SVD
    [U, S, V] = svd(F); % By construction, it holds that U*S*V' = F
    O = U(:, 1:n);
    x_ini = S(1:n, 1:n)*V(:, 1:n)';
end   