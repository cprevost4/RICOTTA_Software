%        Tensor-based multi-frame super-resolution MRI                    %
%-------------------------------------------------------------------------%

% Copyright (c) 2023 Clemence Prevost, Freddy Odille
% https://github.com/cprevost4/RICOTTA_Software
% Contact: clemence.prevost@univ-lille.fr

% This software reproduces the results from the paper called:
% "MULTI-FRAME SUPER-RESOLUTION MRI USING COUPLED LOW-RANK TUCKER
% APPROXIMATION" - C.Prévost, F. Odille (2022)
% and "Non-local tensor sparse coding for multi-image
% super-resolution in magnetic resonance imaging" - C.Prévost, F. Odille (2022)
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

d1 = 4; %D1 = eye(size(HRII,1)); D1 = D1(1:d1:end,:);
% d2 = 4; D2 = eye(size(HRII,2)); D2 = D2(1:d2:end,:);
d3 = 3; %D3 = eye(size(HRII,3)); D3 = D3(1:d3:end,:);

D1 = zeros(size(HRII,1)/d1,size(HRII,1));
for n=1:size(D1,1)
    D1(n,(n-1)*d1+1:n*d1) = 1/d1;
end
D2 = D1;

D3 = zeros(size(HRII,3)/d3,size(HRII,3));
for n=1:size(D3,1)
    D3(n,(n-1)*d3+1:n*d3) = 1/d3;
end

Y1 = tmprod(HRII,D1,1); Y2 = tmprod(HRII,D2,2); Y3 = tmprod(HRII,D3,3);
Y1 = awgn(Y1,25,'measured'); Y2 = awgn(Y2,25,'measured'); Y3 = awgn(Y3,25,'measured');

%% Sparse non-local approach (3 images)

patchsize = 2; overlap = 1; nb_cluster = 10; par.lambda = 1e-5; par.mu = 1e-5;

R = [patchsize^2,size(Y1,3),2000];

EZ = kmeans_misr(Y1,Y2,Y3,D1,D2,D3,patchsize,overlap,nb_cluster,par,R,3);

for k=1:size(EZ,3)
    si_depth(k) = sharpness_index(EZ(:,:,k));
end
mean_si_depth0 = mean(si_depth);

for k=1:size(EZ,2)
    si_dim1(k) = sharpness_index(reshape(EZ(k,:,:),size(EZ,2),size(EZ,3)));
    si_dim2(k) = sharpness_index(reshape(EZ(:,k,:),size(EZ,1),size(EZ,3)));
end
mean_si_plane0 = 0.5*mean(si_dim1) + 0.5*mean(si_dim2);


%% Sparse non-local approach (2 images)

patchsize = 2; overlap = 1; nb_cluster = 10; par.lambda = 1e-5; par.mu = 1e-5;

R = [patchsize^2,size(Y1,3),2000];

EZ2 = kmeans_misr(Y1,Y2,Y3,D1,D2,D3,patchsize,overlap,nb_cluster,par,R,2);

for k=1:size(EZ2,3)
    si_depth(k) = sharpness_index(EZ2(:,:,k));
end
mean_si_depth01 = mean(si_depth);

for k=1:size(EZ,2)
    si_dim1(k) = sharpness_index(reshape(EZ2(k,:,:),size(EZ2,2),size(EZ2,3)));
    si_dim2(k) = sharpness_index(reshape(EZ2(:,k,:),size(EZ2,1),size(EZ2,3)));
end
mean_si_plane01 = 0.5*mean(si_dim1) + 0.5*mean(si_dim2);



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


%% Figures

clims = [min(HRII(:)) max(HRII(:))];


figure
colormap('gray')
subaxis(6,3,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII(:,:,4),clims); colorbar; ylabel('Reference');axis off;
subaxis(6,3,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII(:,60,:),clims); colorbar, axis off;
subaxis(6,3,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII(60,:,:),clims); colorbar, axis off;
axis off
subaxis(6,3,4, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat1(:,:,4),clims); colorbar; ylabel('RICOTTA');axis off;
subaxis(6,3,5, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat1(:,60,:),clims); colorbar, axis off;
subaxis(6,3,6, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat1(60,:,:),clims); colorbar, axis off;
axis off
subaxis(6,3,7, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat2(:,:,4),clims); colorbar; ylabel('Block-RICOTTA');axis off;
subaxis(6,3,8, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat2(:,60,:),clims); colorbar, axis off;
subaxis(6,3,9, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( HRII_hat2(60,:,:),clims); colorbar, axis off;
axis off
subaxis(6,3,13, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( EZ(:,:,5),clims); colorbar; ylabel('Sparse (3im)');axis off;
subaxis(6,3,14, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( EZ(:,60,:),clims); colorbar, axis off;
subaxis(6,3,15, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( EZ(60,:,:),clims); colorbar, axis off;
axis off
subaxis(6,3,16, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( EZ2(:,:,5),clims); colorbar; ylabel('Sparse (2im)');axis off;
subaxis(6,3,17, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( EZ2(:,60,:),clims); colorbar, axis off;
subaxis(6,3,18, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05); imagesc_abs_squeeze( EZ2(60,:,:),clims); colorbar, axis off;
axis off

%% Tables

table = ["Algorithm" "SI_{1,2}" "SI_3";
"Best" "Infty" "Infty";
"Proposed (3 im.)" mean_si_plane0 mean_si_depth0 ;
"Proposed (2 im.)" mean_si_plane01 mean_si_depth01 ;
"RICOTTA" mean_si_plane1 mean_si_depth1;
"Block-RICOTTA" mean_si_plane2 mean_si_depth2]
