clear all; close all; clc

load('FOETAL_ECG.dat');
x = FOETAL_ECG(:, 2:end)';
for n = 1 : 8
    x(n,:) = x(n,:)/std(x(n,:)); %  standardizing 
end
power_x = sum(x(:).^2)/size(x,2);

% ordinary PCA
[V, D] = eig(x*x');

y = V(:, end:-1:end-2)'*x; % 3 PCs
e = V(:, 1:end-3)'*x;
subplot(6,3,1); plot(y(1,:), 'k'); axis('tight'); axis('off'); title({'Ordinary PCA', '{\itN}=3, 97.7%'})
subplot(6,3,4); plot(y(2,:), 'k'); axis('tight'); axis('off')
subplot(6,3,7); plot(y(3,:), 'k'); axis('tight'); axis('off')
power_y = sum(y(:).^2)/size(x,2); power_e = sum(e(:).^2)/size(x,2);
fprintf('\nPCA. PCs power %g; error power %g; Unexplained power %g\n', power_y/power_x, power_e/power_x, abs(power_y + power_e - power_x)/power_x)

y = V(:, end:-1:end-1)'*x; % 2 PCs
e = V(:, 1:end-2)'*x;
subplot(6,3,10); plot(y(1,:), 'k'); axis('tight'); axis('off'); title('{\itN}=2, 93.0%')
subplot(6,3,13); plot(y(2,:), 'k'); axis('tight'); axis('off')
power_y = sum(y(:).^2)/size(x,2); power_e = sum(e(:).^2)/size(x,2);
fprintf('\nPCA. PCs power %g; error power %g; Unexplained power %g\n', power_y/power_x, power_e/power_x, abs(power_y + power_e - power_x)/power_x)

y = V(:, end)'*x;   % 1 PC
e = V(:, 1:end-1)'*x;
subplot(6,3,16); plot(y(1,:), 'k'); axis('tight'); axis('off'); title('{\itN}=1, 74.8%')
power_y = sum(y(:).^2)/size(x,2); power_e = sum(e(:).^2)/size(x,2);
fprintf('\nPCA. PCs power %g; error power %g; Unexplained power %g\n', power_y/power_x, power_e/power_x, abs(power_y + power_e - power_x)/power_x)

% PCA in the frequency domain
[y, e] = FreqPCA(x, 3, 20); % 3 PCs
subplot(6,3,2); plot(y(1,:), 'k'); axis('tight'); axis('off'); title({'Freq. PCA', '{\itN}=3, 96.8%'})
subplot(6,3,5); plot(y(2,:), 'k'); axis('tight'); axis('off')
subplot(6,3,8); plot(y(3,:), 'k'); axis('tight'); axis('off')
power_y = sum(y(:).^2)/size(x,2); power_e = sum(e(:).^2)/size(x,2);
fprintf('\nFreq PCA. PCs power %g; error power %g; Unexplained power %g\n', power_y/power_x, power_e/power_x, abs(power_y + power_e - power_x)/power_x)

[y, e] = FreqPCA(x, 2, 20); % 2 PCs
subplot(6,3,11); plot(y(1,:), 'k'); axis('tight'); axis('off'); title('{\itN}=2, 94.1%')
subplot(6,3,14); plot(y(2,:), 'k'); axis('tight'); axis('off')
power_y = sum(y(:).^2)/size(x,2); power_e = sum(e(:).^2)/size(x,2);
fprintf('\nFreq PCA. PCs power %g; error power %g; Unexplained power %g\n', power_y/power_x, power_e/power_x, abs(power_y + power_e - power_x)/power_x)

[y, e] = FreqPCA(x, 1, 20); % 1 PC
subplot(6,3,17); plot(y(1,:), 'k'); axis('tight'); axis('off'); title('{\itN}=1, 87.4%')
power_y = sum(y(:).^2)/size(x,2); power_e = sum(e(:).^2)/size(x,2);
fprintf('\nFreq PCA. PCs power %g; error power %g; Unexplained power %g\n\n', power_y/power_x, power_e/power_x, abs(power_y + power_e - power_x)/power_x)

% Convolutional PCA
M = 8; N = 3; L = 10; % 3 PCs
W = zeros(M, N, L+1); temp = dctmtx(M); W(:,:,L/2+1) = temp(1:N, :)';
[W, err] = ConvPCA(x, W, 0.2, 200, 1);
y = PolyMatFilter(PolyMatH(W), x);
subplot(6,3,3); plot(y(1,:), 'k'); axis('tight'); axis('off'); title({'Conv. PCA', '{\itN}=3, 98.8%'})
subplot(6,3,6); plot(y(2,:), 'k'); axis('tight'); axis('off')
subplot(6,3,9); plot(y(3,:), 'k'); axis('tight'); axis('off')
power_y = sum(y(:).^2)/size(x,2); 
fprintf('\nConv. PCA. PCs power %f; error power %f; Unexplained power %f\n\n', power_y/power_x, err/power_x, abs(power_y + err - power_x)/power_x)

M = 8; N = 2; L = 10; % 2 PCs
W = zeros(M, N, L+1); temp = dctmtx(M); W(:,:,L/2+1) = temp(1:N, :)';
[W, err] = ConvPCA(x, W, 0.2, 200, 1);
y = PolyMatFilter(PolyMatH(W), x);
subplot(6,3,12); plot(y(1,:), 'k'); axis('tight'); axis('off'); title('{\itN}=2, 97.0%')
subplot(6,3,15); plot(y(2,:), 'k'); axis('tight'); axis('off')
power_y = sum(y(:).^2)/size(x,2); 
fprintf('\nConv. PCA. PCs power %f; error power %f; Unexplained power %f\n\n', power_y/power_x, err/power_x, abs(power_y + err - power_x)/power_x)

M = 8; N = 1; L = 10; % 1 PC
W = zeros(M, N, L+1); temp = dctmtx(M); W(:,:,L/2+1) = temp(1:N, :)';
[W, err] = ConvPCA(x, W, 0.2, 200, 1);
y = PolyMatFilter(PolyMatH(W), x);
subplot(6,3,18); plot(y(1,:), 'k'); axis('tight'); axis('off'); title('{\itN}=1, 89.7%')
power_y = sum(y(:).^2)/size(x,2); 
fprintf('\nConv. PCA. PCs power %f; error power %f; Unexplained power %f\n\n', power_y/power_x, err/power_x, abs(power_y + err - power_x)/power_x)