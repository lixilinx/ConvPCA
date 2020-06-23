function [W, err] = ConvPCA(x, W, step_size, max_num_iter, echo)
% z-domain implementation of convolutional PCA  
% INPUTS:
%   x: observations with size (M, T)
%   W: filter initial guess with size (M, N, L+1)
%   step_size: step size for gradient descent
%   max_num_iter: maximum number of iterations
%   echo: 1 (0) for (not) displaying intermediate results 
% OUTPUTS:
%   W: filter
%   err: fitting MSE
[M, T] = size(x);
[M, N, L] = size(W); L = L - 1;
Rx = zeros(M, M, 4*L+1);
Rx(:,:,2*L+1) = (x*x')/T;
for i = 1 : 2*L
    Rx(:,:, 2*L+1+i) = x(:, i+1:T)*x(:, 1:T-i)'/T;
    Rx(:,:, 2*L+1-i) = Rx(:,:, 2*L+1+i)';
end

scale = sqrt(M/trace(Rx(:,:,2*L+1))); % to scale x to the same average power per component
tol = 1e-6;
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
    momentum = 0.9*momentum + 0.1*Grad;
    W = W - (step_size*scale*scale)*momentum;
    
    % display progress
    if echo
        fprintf('iter: %u; err: %g. ', iter, err);
        if mod(iter, 100)==0
            fprintf('\n');
        end
    end
    if scale*scale*max(abs(Grad(:))) < tol
        break;
    end
end
% lastly, rotate W such that Ry(0) is diagonal
[V, ~] = eig(WhRxW(:,:, 3*L+1)+WhRxW(:,:, 3*L+1)');%Ry(0) must be Hermitian
for i = 0 : L
    W(:,:,i+1) = W(:,:,i+1)*V;
end