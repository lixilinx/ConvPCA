function [W, err] = SparseConvPCA(x, W, step_size, eta, max_num_iter, echo)
%
% The only difference with ConvPCA is the added parameter eta
% It is for penalty eta*sum(abs(W)) to encourage sparse filter coefficients
%
[M, T] = size(x);
[M, N, L] = size(W); L = L - 1;
Rx = zeros(M, M, 4*L+1);
Rx(:,:,2*L+1) = (x*x')/T;
for i = 1 : 2*L
    Rx(:,:, 2*L+1+i) = x(:, i+1:T)*x(:, 1:T-i)'/T;
    Rx(:,:, 2*L+1-i) = Rx(:,:, 2*L+1+i)';
end

scale = sqrt(M/trace(Rx(:,:,2*L+1))); % to scale x to the same average power per component
vareps = eps;
err = 0;
momentum = zeros(size(W));
for iter = 1 : max_num_iter
    % we will use these common terms more than one time
    % I do not early trucate RxW and WhRxW to keep code easier to read
    RxW = PolyMatMult(Rx, W); % exp of z: 2*L --> -3*L
    WhRxW = PolyMatMult(PolyMatH(W), RxW); % exp of z: 3*L --> -3*L
    WhW = PolyMatMult(PolyMatH(W), W); % exp of z: L --> -L
    
    % calculate the fitting error, i.e., the power of minor components
    err = trace(Rx(:,:, 2*L+1)) - 2*trace(WhRxW(:,:, 3*L+1));
    err = err + sum(sum(sum( conj(WhW).*WhRxW(:,:, 2*L+1:4*L+1) ))); % tricks: tr(A*B)=sum(A.*B^T); WhRxW is Hermitian 

    % calculate the trucated gradient
    % only take the parts with exp of z from 0 down to -L
    Grad = -2*RxW(:,:, 2*L+1 : 3*L+1);
    term1 = PolyMatMult(W, WhRxW(:,:, 2*L+1 : 4*L+1)); % exp of z: L --> -2*L
    Grad = Grad + term1(:,:, L+1 : 2*L+1);
    term2 = PolyMatMult(RxW(:,:, L+1:4*L+1), WhW); % exp of z: 2*L --> -3*L
    Grad = Grad + term2(:,:, 2*L+1 : 3*L+1);
   
    % gradient descent with momentum for acceleration
    momentum = 0.9*momentum + 0.1*(scale*scale*Grad + 0.5*eta*W./sqrt(W.*conj(W) + vareps));
    W = W - step_size*momentum;
    
    % display progress
    if echo
        fprintf('iter: %u; total loss: %g. ', iter, scale*scale*err + eta*sum(sum(sum(sqrt(W.*conj(W) + vareps)))));
        if mod(iter, 100)==0
            fprintf('\n');
        end
    end
end