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

%% Examine the condition number of the XtX matrix used for solving G

% Find singular vectors and their projections
[U, S1, ~] = svd([tens2mat(Y2,1,[]) tens2mat(Y3,1,[])],0);
[V, S2, ~] = svd([tens2mat(Y1,2,[]) tens2mat(Y3,2,[])],0);
[W, S3, ~] = svd([tens2mat(Y1,3,[]) tens2mat(Y2,3,[])], 0);


PU = D1 * U; PV = D2 * V; PW = D3 * W;
I1 = size(PU,1); J2 = size(PV,1); K3 = size(PW,1);

sigmaU = zeros(I1,2); sigmaV = zeros(J2,2); s1W = zeros(K3,1);
for i=1:size(sigmaU,1), 
  sU = svd(PU(:,1:i)); 
  sigmaU(i,:) = [sU(1) sU(end)]; 
end
for i=1:size(sigmaV,1), 
  sV = svd(PV(:,1:i)); 
  sigmaV(i,:) = [sV(1) sV(end)]; 
end
for i=1:size(s1W,1), s1W(i) = norm(PW(:,1:i),2); end

minIJ_H = min(I1,J2);
condXtX = (sigmaU(1:minIJ_H,2).^2 + sigmaV(1:minIJ_H,2).^2 + s1W(K3)^2) ./ ...
          (sigmaU(1:minIJ_H,2).^2 + sigmaV(1:minIJ_H,2).^2);
      
figure 
subplot(1,2,1)
semilogy(condXtX)
xlabel('$R_1 = R_2$', 'Interpreter', 'latex'); ylabel('$\mathrm{cond}({\bf X}^{\rm T} {\bf X}) $', 'Interpreter', 'latex')
set(gcf, 'Position',  [100, 100, 200, 160])
subplot(1,2,2)
plot(s1W);
xlabel('$R_3$', 'Interpreter', 'latex'); ylabel('$\sigma_1({\bf P}_3 \widehat{\bf W})$', 'Interpreter', 'latex')
set(gcf, 'Position',  [100, 100, 200, 160])
%saveas(gcf,'figures/fig17.fig')

%% Run RICOTTA and evaluate the RMSE as a function of the ranks

l1 = 1; l2 = 1; l3 = 1; alpha = 0.001;


for r1 = 1:128
    for r3 = 1:27
        R = [r1 r1 r3]
        h = size(Y1,1); l = size(Y2,2); p = size(Y3,3);
        if (R(1)>h && R(2)>l && R(3)>p)
            snr1(r1,r3) = NaN; cc1(r1,r3) = NaN; rmse1(r1,r3) = NaN;
        else
            [U,V,W,G] = ricotta(Y1,Y2,Y3,D1,D2,D3,R,l1,l2,l3,alpha);
            HRII_hat1 = lmlragen({U,V,W},G);
            snr1(r1,r3) = r_snr(HRII,HRII_hat1);
            cc1(r1,r3) = cc(HRII,HRII_hat1);
            rmse1(r1,r3) = rmse(HRII,HRII_hat1);
        end
    end
end

figure
subplot(1,3,1); imagesc(snr1); title('R-SNR'); xlabel('R_1'); ylabel('R_3');
subplot(1,3,2); imagesc(cc1); title('CC'); xlabel('R_1'); ylabel('R_3');
subplot(1,3,3); imagesc(rmse1); title('RMSE'); xlabel('R_1'); ylabel('R_3');



