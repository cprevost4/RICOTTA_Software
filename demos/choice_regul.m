%        Tensor-based multi-frame super-resolution MRI                    %
%-------------------------------------------------------------------------%

% Copyright (c) 2022 Clemence Prevost, Freddy Odille
% https://github.com/cprevost4/RICOTTA_Software
% Contact: clemence.prevost@univ-lille.fr

% This software reproduces the results from the paper called:
% "MULTI-FRAME SUPER-RESOLUTION MRI USING COUPLED LOW-RANK TUCKER
% APPROXIMATION" - C.Prévost, F. Odille
%
% In order to run the demo, you will need to add to your MATLAB path:
% - Tensorlab 3.0: https://www.tensorlab.net
% - Sharpness Index: http://www.mi.parisdescartes.fr/~moisan/sharpness/
% - Beltrami primal-dual solver: https://math.montana.edu/dzosso/code/
%
%-------------------------------------------------------------------------%
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

l1 = 1; l2 = 1; l3 = 1;

R = [32 32 22];

alpha = [1;2;5].*logspace(-5,0,6);
alpha = alpha(:);

for r=1:length(alpha)
    alpha(r)
    [U,V,W,G] = ricotta(Y1,Y2,Y3,D1,D2,D3,R,l1,l2,l3,alpha(r));
    HRII_hat1 = lmlragen({U,V,W},G);
    snr1(r) = r_snr(HRII,HRII_hat1);
    cc1(r) = cc(HRII,HRII_hat1);
    rmse1(r) = rmse(HRII,HRII_hat1);
end


figure
subplot(1,3,1); semilogx(alpha,snr1); title('R-SNR'); xlabel('R_1'); ylabel('R_3');
subplot(1,3,2); semilogx(alpha,cc1); title('CC'); xlabel('R_1'); ylabel('R_3');
subplot(1,3,3); semilogx(alpha,rmse1); title('RMSE'); xlabel('R_1'); ylabel('R_3');



