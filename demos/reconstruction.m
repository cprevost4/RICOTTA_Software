%        Tensor-based multi-frame super-resolution MRI                    %
%-------------------------------------------------------------------------%

% Copyright (c) 2022 Clemence Prevost, Freddy Odille
% https://github.com/cprevost4/RICOTTA_Software
% Contact: clemence.prevost@univ-lille.fr

% This software reproduces the results from the paper called:
% "MULTI-FRAME SUPER-RESOLUTION MRI USING COUPLED LOW-RANK TUCKER
% APPROXIMATION" - C.Prévost, F. Odille (2022)
%
% In order to run the demo, you will need to add to your MATLAB path:
% - Tensorlab 3.0: https://www.tensorlab.net
% - Sharpness Index: http://www.mi.parisdescartes.fr/~moisan/sharpness/
% - Beltrami primal-dual solver: https://math.montana.edu/dzosso/code/
%
%-------------------------------------------------------------------------%
%                              CONTENT
% - /baseline_algorithms: contains codes and adaptors for other methods
% - /data : contains data for synthetic examples (Section VI.D)
% - /demos : contains demo files that produce tables and figures
% - /metrics : contains the metrics use for comparison in the paper
% - /src : contains helpful files to run the demos
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

%% Run RICOTTA

l1 = 1; l2 = 1; l3 = 1; alpha = 0.01;
R = [32 32 22];

tic;
[U,V,W,G] = ricotta(Y1,Y2,Y3,D1,D2,D3,R,l1,l2,l3,alpha);
t1 = toc;

HRII_hat1 = lmlragen({U,V,W},G);

snr1 = r_snr(HRII,HRII_hat1);
cc1 = cc(HRII,HRII_hat1);

for k=1:size(HRII_hat1,3)
    si_depth(k) = sharpness_index(HRII_hat1(:,:,k));
end
mean_si_depth1 = mean(si_depth);

for k=1:size(HRII_hat1,2)
    si_dim1(k) = sharpness_index(reshape(HRII_hat1(k,:,:),size(HRII_hat1,2),size(HRII_hat1,3)));
    si_dim2(k) = sharpness_index(reshape(HRII_hat1(:,k,:),size(HRII_hat1,1),size(HRII_hat1,3)));
end
mean_si_plane1 = 0.5*mean(si_dim1) + 0.5*mean(si_dim2);

%% Run Block-RICOTTA

l1 = 1; l2 = 1; l3 = 1; alpha = 0.01;
opts.Nblocks = [2 2];
R = [16 16 10];

tic;
[HRII_hat2] = block_ricotta(Y1,Y2,Y3,D1,D2,D3,R,l1,l2,l3,alpha,opts);
t2 = toc;

snr2 = r_snr(HRII,HRII_hat2);
cc2 = cc(HRII,HRII_hat2);

for k=1:size(HRII_hat1,3)
    si_depth(k) = sharpness_index(HRII_hat2(:,:,k));
end
mean_si_depth2 = mean(si_depth);

for k=1:size(HRII_hat1,2)
    si_dim1(k) = sharpness_index(reshape(HRII_hat2(k,:,:),size(HRII_hat1,2),size(HRII_hat1,3)));
    si_dim2(k) = sharpness_index(reshape(HRII_hat2(:,k,:),size(HRII_hat1,1),size(HRII_hat1,3)));
end
mean_si_plane2 = 0.5*mean(si_dim1) + 0.5*mean(si_dim2);

%% RICOTTA w/o regul + Beltrami 3D

l1 = 1; l2 = 1; l3 = 1; alpha = 0.01;
R = [32 32 22];

tic;
[U,V,W,G] = ricotta_noregul(Y1,Y2,Y3,D1,D2,D3,R,l1,l2,l3);
t11 = toc;

HRII_hat11 = lmlragen({U,V,W},G);

s = 5;
beta = 1; tol = 1e-3; lambda = 1/s;
r1 = 0.1; r2 = 0.1/lambda;
Maxit = 100;

tic;
[HRII_hat3,k] = beltrami3D(HRII_hat11,Maxit,beta,lambda,tol,r1,r2);
t22 = toc;
t3 = t11 + t22;

snr3 = r_snr(HRII,HRII_hat3);
cc3 = cc(HRII,HRII_hat3);

for k=1:size(HRII_hat3,3)
    si_depth(k) = sharpness_index(HRII_hat3(:,:,k));
end
mean_si_depth3 = mean(si_depth);

for k=1:size(HRII_hat3,2)
    si_dim1(k) = sharpness_index(reshape(HRII_hat3(k,:,:),size(HRII_hat3,2),size(HRII_hat3,3)));
    si_dim2(k) = sharpness_index(reshape(HRII_hat3(:,k,:),size(HRII_hat3,1),size(HRII_hat3,3)));
end
mean_si_plane3 = 0.5*mean(si_dim1) + 0.5*mean(si_dim2);

%% Figures

clims = [ min(HRII(:)) max(HRII(:)) ];

figure
colormap('gray')
subaxis(4,3,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII(:,:,5), clims ); colorbar; ylabel('Reference');axis off;
subaxis(4,3,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII(:,60,:), clims ); colorbar, axis off;
subaxis(4,3,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII(60,:,:), clims ); colorbar, axis off;
axis off
subaxis(4,3,4, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat1(:,:,5), clims ); colorbar; ylabel('RICOTTA');axis off;
subaxis(4,3,5, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat1(:,60,:), clims ); colorbar, axis off;
subaxis(4,3,6, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat1(60,:,:), clims ); colorbar, axis off;
axis off
subaxis(4,3,7, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat2(:,:,5), clims ); colorbar; ylabel('Block-RICOTTA');axis off;
subaxis(4,3,8, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat2(:,60,:), clims ); colorbar, axis off;
subaxis(4,3,9, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat2(60,:,:), clims ); colorbar, axis off;
axis off
subaxis(4,3,10, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat3(:,:,5), clims ); colorbar; ylabel('RICOTTA+Beltrami');axis off;
subaxis(4,3,11, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat3(:,60,:), clims ); colorbar, axis off;
subaxis(4,3,12, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat3(60,:,:), clims ); colorbar, axis off;
axis off

%% Tables

table = ["Algorithm" "R-SNR" "CC" "SI_{1,2}" "SI_3" "Time (sec)";
"Best" "Infty" "1" "Infty" "Infty" "0";
"RICOTTA" snr1 cc1 mean_si_plane1 mean_si_depth1 t1;
"Block-RICOTTA" snr2 cc2 mean_si_plane2 mean_si_depth2 t2; 
"RICOTTA + Beltrami3D" snr3 cc3 mean_si_plane3 mean_si_depth3 t3]
