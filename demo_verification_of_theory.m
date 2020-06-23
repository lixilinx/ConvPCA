clear all; close all; clc
% Let's verify 1) the conservation of variances and 2) W^H*W = eye

x = randn(7, 1000) + sqrt(-1)*randn(7, 1000);
for i = 1 : 5
    x = (randn(7,7) + sqrt(-1)*randn(7,7))*x(:, 1:end-1) + x(:, 2:end);
end

M = 7; N = 2; L = 3;
W = randn(M, N, L+1) + sqrt(-1)*randn(M, N, L+1); W = sqrt(N)*W/norm(W(:));
[W, err] = ConvPCA(x, W, 0.5, 1000, 1);

y = PolyMatFilter(PolyMatH(W), x);
power_x = sum(sum(x.*conj(x)))/size(x,2); power_y = sum(sum(y.*conj(y)))/size(x,2);
(power_y + err)/power_x % this should be 1, the conservation of variances

PolyMatMult(PolyMatH(W),W) % this should be close to eye(N)*z^(-L-1); not exact due to finite filter length