clear all; close all; clc

T = 2100
tau0 = 3
s = filter(1, [1, -0.9], randn(1, T));
x = [s(tau0+1:end); s(1:end-tau0)];
x = x + 0.1*randn(size(x));

[M, T] = size(x); N = 1; L = 10;
W = zeros(M, N, L+1); temp = dctmtx(M); W(:,:,L/2+1) = temp(1:N, :)';
y = zeros(N, T);
err = zeros(1, T);
m1 = zeros(size(W)); %1st momentum: smoothed gradient
m2 = 1; % 2nd momentum: (gradient norm)^2
step_size = sqrt(1e-2/(L+1)/M/N)
for t = L+1 : T-L
    if mod(t-L-1,500)==0
        w12 = squeeze(W)';
        subplot(2,5,6+(t-L-1)/500); 
        stem([0:L], w12(:,1), 'ko')
        hold on; stem([0:L], w12(:,2), 'kx')
        xlim([0, L]); ylim([-0.1, 0.8])
        if t==L+1
            ylabel({'${\bf W}(z)=[w_1(z); w_2(z)]$', '$w_1(z): -$o; $w_2(z): -$x'},'Interpreter','latex')
        end
    end
    
    for i = 0 : L
        y(:,t) = y(:,t) + W(:,:,i+1)'*x(:, t+i);
    end
    e = x(:,t);
    for i = 0 : L
        e = e - W(:,:,i+1)*y(:, t-i);
    end
    err(t) = e'*e;
    grad = zeros(size(W));
    for i = 0 : L
        grad(:,:,i+1) = -e*y(:, t-i)';
        for j = 0 : L
            grad(:,:,i+1) = grad(:,:,i+1) - x(:, t+i-j)*e'*W(:,:, j+1);
        end
    end
    m1 = 0.9*m1 + 0.1*grad;
    m2 = 0.9*m2 + 0.1*sum(grad(:).^2);
    W = W - step_size*(m1./sqrt(m2 + eps) + 1e-1*sign(W));
end
subplot(2,1,1); semilogy(err, 'k-'); xlim([1, 2000]); xlabel('$t$','Interpreter','latex'); ylabel('$\|{\bf e}(t)\|^2$','Interpreter','latex')