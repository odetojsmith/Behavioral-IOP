function [T] = get_imp_resp_matrix(Ud, Yd, Tini, Tf, n)

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
    G = pinv([Up; Yp; Uf])*[zeros(m*Tini, m); zeros(p*Tini, m); eye(m); zeros(m*(Tf - 1), m)];
    % Compute a set of generators for the free responses of the system
    T = Yf*G;
end   