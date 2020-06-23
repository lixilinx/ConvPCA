clear all; close all; clc

recordings = audioread('array_2mics.wav');
T = 2048
for trial = 1 : 4
    x = recordings(2000:2000+T-1, :)';
    
    % standard cross-correlation for delay estimation
    r = zeros(1, 33);
    for i = -16:16
        r(i+17) = sum(x(1, 17:end-16).*x(2, 17-i:end-i-16));
    end
    subplot(4,3,3*trial-2);
    plot([-16:16], r, '.k-');
    if trial==1
        title({'Cross-correlation', 'coefficients'},'Interpreter','latex')
    end
    if trial==4
        xlabel('Delay')
    end
    ylabel({'Sample size', ['{\itT}=',num2str(T)]})
    set(gca,'YTick',[])
    axis('tight'); %axis('off');
    
    % generalized cross-correlation phase transform (GCC-PHAT)
    X = fft(x, [], 2);
    r = X(1,:).*conj(X(2,:));
    r = r./abs(r);
    r = ifftshift(ifft(r));
    r = r(T/2+1-16:T/2+1+16);
    subplot(4,3,3*trial-1);
    plot([-16:16], r, '.k-'); axis('tight'); set(gca,'YTick',[]); %axis('off');
    if trial==1
        title({'GCC-PHAT', 'coefficients'},'Interpreter','latex')
    end
    if trial==4
        xlabel('Delay')
    end
    
    % conv PCA for delay estimation
    M = 2; N = 1; L = 32;
    W = zeros(M, N, L+1); temp = dctmtx(M); W(:,:,L/2+1) = temp(1:N, :)';
    [W, ~] = SparseConvPCA(x, W, 0.2, 0.2, 200, 1);
    [W, ~] = SparseConvPCA(x, W, 0.02, 0.2, 100, 1); % a smaller step size for refinement
    w12 = squeeze(W)';
    subplot(4,3,3*trial);
    plot(w12(:,1), '.k-');
    hold on; plot(w12(:,2), 'k-x'); axis('tight'); set(gca,'YTick',[]); %axis('off');
    if trial==1
        title({'Conv. PCA', '$w_1(z): \cdot  \quad w_2(z): \times $'},'Interpreter','latex')
    end
    if trial==4
        xlabel('Tap index')
    end

    T = 2*T;
end