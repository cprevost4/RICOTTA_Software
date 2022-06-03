function [I,k] = beltrami3D(I0,Maxit,beta,lambda,tol,r1,r2)

%% Initialization

I = I0;
beta2 = beta.^2;
phix = zeros(size(I));
phiy = zeros(size(I));
phiz = zeros(size(I));

primal = zeros(Maxit+1,1);

[I_x,I_y,I_z] = Forward(I);

N = sqrt(1 + beta2*( I_x.^2 + I_y.^2 + I_z.^2));

primal(1) = sum(N(:) + lambda/2*(I(:)-I0(:)).^2);

%% Algo

Maxit = 2000;

% maximum Maxit iterations
for k = 1:Maxit
    
    % phi step
    [I_x,I_y,I_z] = Forward(I);
    
    Den = real(beta*sqrt(beta2 - phix.^2 - phiy.^2 - phiz.^2));
    Den(Den<0) = 0;
    phix = phix + r1*(I_x.*Den - phix);
    phiy = phiy + r1*(I_y.*Den - phiy);
    phiz = phiz + r1*(I_z.*Den - phiz);
    
    proj = sqrt( phix.^2 + phiy.^2 + phiz.^2);
    proj( proj < 1 ) = 1;
    phix = phix ./ proj;
    phiy = phiy ./ proj;
    phiz = phiz ./ proj;
    
    % I step
    I = I + r2*lambda*( I0 - I + 1/lambda*(BackwardX(phix) + BackwardY(phiy) + BackwardZ(phiz)));
    
    % convergence?
    N = sqrt(1 + beta2*( I_x.^2 + I_y.^2 + I_z.^2));
    primal(k+1) = sum(N(:) + lambda/2*(I(:)-I0(:)).^2);
    
    if abs(primal(k+1)-primal(k)) < primal(1)*tol
        break;
    end
end


end

