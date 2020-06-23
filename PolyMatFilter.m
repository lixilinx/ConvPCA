function y = PolyMatFilter(A, x)
% Apply polynomial matrix filter on x
%
% to verify the implementation:
%   a = randn(1,1,3);
%   b = randn(1,5);
%   PolyMatFilter(a,b) - conv(squeeze(a), b) % should be 0
%
[I,J,M] = size(A); [J,N] = size(x);
y = zeros(I, M+N-1);
for i = 0 : M+N-2
    for j = max(0, i+1-N) : min(M-1, i)
        y(:,i+1) = y(:,i+1) + A(:,:,j+1)*x(:,i-j+1);
    end
end