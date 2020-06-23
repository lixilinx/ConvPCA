function [y, e] = FreqPCA(x, N, F)
% PCA in the frequency domain
% x: observation 
% N: number of principle components
% F: frame size or FFTSize, comparable to filter length
% y: principle components 
% e: minor components
% 
if mod(F, 2)==1
    F = F + 1;
    warning('Change F to F+1 to facilitate the overlap-and-add operation.');
end
win = coswin(F); % I use the cos window

[M, T] = size(x);
x = [zeros(M,F), x, zeros(M,F)]; % pad zeros at both ends 
T = T + 2*F;

% first, we estimate the power spectral density matrix
Px = zeros(M, M, F);
i = 1;
while i+F-1 <= T 
    X = fft(win.*x(:, i:i+F-1), [], 2);
    for j = 1 : F
        Px(:,:,j) = Px(:,:,j) + X(:,j)*X(:,j)';
    end
    i = i + F/2; % half overlap between successive frames 
end

% then, do PCA at each frequency bin
Uf = zeros(N, M, F); 
Vf = zeros(M-N, M, F);
for i = 1 : F
    [V, D] = eig(Px(:,:,i) + Px(:,:,i)');
    [~, I] = sort(diag(D), 'descend');
    Uf(:,:, i) = V(:, I(1:N))'; % eigenvectors for the principal components
    Vf(:,:, i) = V(:, I(N+1:end))'; % eigenvectors for the minor components
end

% last, do filtering in the frequency domain
y = zeros(N, T);
e = zeros(M-N, T);
i = 1; Y = zeros(N, F); E = zeros(M-N, F);
while i+F-1 <= T
    X = fft(win.*x(:, i:i+F-1), [], 2);
    
    for j = 1 : F
        Y(:,j) = Uf(:,:, j)*X(:,j);
        E(:,j) = Vf(:,:, j)*X(:,j);
    end
    
    y(:, i:i+F-1) = y(:, i:i+F-1) + ifft(Y, [], 2);
    e(:, i:i+F-1) = e(:, i:i+F-1) + ifft(E, [], 2);
    
    i = i + 0.5*F;
end