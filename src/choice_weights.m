clear all
close all
clc

%% Generate dataset

load('mri'); clear siz map
HRII = double(reshape(D,[size(D,1),size(D,2),size(D,4)])); clear D

d1 = 4; D1 = eye(size(HRII,1)); D1 = D1(1:d1:end,:);
d2 = 4; D2 = eye(size(HRII,2)); D2 = D2(1:d2:end,:);
d3 = 3; D3 = eye(size(HRII,3)); D3 = D3(1:d3:end,:);

Y1 = tmprod(HRII,D1,1); Y2 = tmprod(HRII,D2,2); Y3 = tmprod(HRII,D3,3);
Y1 = awgn(Y1,25,'measured'); Y2 = awgn(Y2,25,'measured'); Y3 = awgn(Y3,25,'measured');

%% Run RICOTTA and evaluate the RMSE as a function of mu

R = [32 32 22]; alpha = 0.01;

l1=0.1:0.1:1;
l3=0.1:0.1:1;

for n = 1:length(l1)
    for m = 1:length(l3)
        [l1(n) l3(m)]
        [U,V,W,G] = ricotta(Y1,Y2,Y3,D1,D2,D3,R,l1(n),l1(n),l3(m),alpha);
        HRII_hat1 = lmlragen({U,V,W},G);
        snr1(n,m) = r_snr(HRII,HRII_hat1);
        cc1(n,m) = cc(HRII,HRII_hat1);
        rmse1(n,m) = rmse(HRII,HRII_hat1);
    end
end


figure
subplot(1,3,1); imagesc(snr1); title('R-SNR'); xlabel('l_1'); ylabel('l_3');
subplot(1,3,2); imagesc(cc1); title('CC'); xlabel('l_1'); ylabel('l_3');
subplot(1,3,3); imagesc(rmse1); title('RMSE'); xlabel('l_1'); ylabel('l_3');



