function [Y_hat] = block_ricotta(Y1,Y2,Y3,P1,P2,P3,R,l1,l2,l3,alpha,opts)

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'Nblocks') || isempty(opts.Nblocks)
    opts.Nblocks = [1,1];
end

range1 = [size(Y1,1),size(Y1,2)]; 
range2 = [size(Y2,1),size(Y2,2)]; 
range3 = [size(Y3,1),size(Y3,2)]; 

step1 = ceil(range1 ./ opts.Nblocks); 
step2 = ceil(range2 ./ opts.Nblocks); 
step3 = ceil(range3 ./ opts.Nblocks); 

Y_hat = zeros(size(Y2,1), size(Y3,2), size(Y1,3));


for i1=1:opts.Nblocks(1)
  for i2=1:opts.Nblocks(2)
      
   ind1_min = [i1-1,i2-1].*step1 + 1; ind1_max = min([i1,i2].*step1, range1);
   ind2_min = [i1-1,i2-1].*step2 + 1; ind2_max = min([i1,i2].*step2, range2);
   ind3_min = [i1-1,i2-1].*step3 + 1; ind3_max = min([i1,i2].*step3, range3);
  
   
   [U, ~, ~] = svds([tens2mat(Y2(ind2_min(1):ind2_max(1),ind2_min(2):ind2_max(2),:),1,[]) ...
                    tens2mat(Y3(ind3_min(1):ind3_max(1),ind3_min(2):ind3_max(2),:),1,[])],R(1));
   [V, ~, ~] = svds([tens2mat(Y1(ind1_min(1):ind1_max(1),ind1_min(2):ind1_max(2),:),2,[]) ...
                     tens2mat(Y3(ind3_min(1):ind3_max(1),ind3_min(2):ind3_max(2),:),2,[])],R(2));
   [W, ~, ~] = svds([tens2mat(Y1(ind1_min(1):ind1_max(1),ind1_min(2):ind1_max(2),:),3,[]) ...
                    tens2mat(Y2(ind2_min(1):ind2_max(1),ind2_min(2):ind2_max(2),:),3,[])], R(3));

    U_tilde = P1(ind1_min(1):ind1_max(1),ind2_min(1):ind2_max(1))*U;
    V_tilde = P2(ind2_min(2):ind2_max(2),ind3_min(2):ind3_max(2))*V;
    W_tilde = P3*W;
   
    if R(1)>R(3) || R(2)||R(3)
    A = l1*(U_tilde'*U_tilde);
    A = A + alpha*eye(size(A));
    B = eye(R(2)*R(3));
    C = eye(R(1));
    D = kron(l2*eye(R(3)),V_tilde'*V_tilde)+ kron(W_tilde'*W_tilde,l3*eye(R(2)));
    D = D + alpha*eye(size(D));
    E = l1*tmprod(Y1(ind1_min(1):ind1_max(1),ind1_min(2):ind1_max(2),:),{U_tilde', V', W'},[1,2,3]) ...
        + l2*tmprod(Y2(ind2_min(1):ind2_max(1),ind2_min(2):ind2_max(2),:),{U', V_tilde', W'},[1,2,3]) ...
        + l3*tmprod(Y3(ind3_min(1):ind3_max(1),ind3_min(2):ind3_max(2),:),{U', V', W_tilde'},[1,2,3]);
    E = reshape(E,R(1),R(2)*R(3));
    G = reshape(bartelsStewart(A,[],[],D,E), R);
else
    A = l1*kron(eye(R(2),U_tilde'*U_tilde)) + l2 * kron(V_tilde'*V_tilde,eye(R(1)));
    A = A + alpha*eye(size(A));
    B = eye(R(3));
    C = eye(R(1)*R(2));
    D = l3*(W_tilde'*W_tilde);
    D = D + alpha*eye(size(D));
    E = l1*tmprod(Y1(ind1_min(1):ind1_max(1),ind1_min(2):ind1_max(2),:),{U_tilde', V', W'},[1,2,3]) ...
        + l2*tmprod(Y2(ind2_min(1):ind2_max(1),ind2_min(2):ind2_max(2),:),{U', V_tilde', W'},[1,2,3]) ...
        + l3*tmprod(Y3(ind3_min(1):ind3_max(1),ind3_min(2):ind3_max(2),:),{U', V', W_tilde'},[1,2,3]);
    E = reshape(E,R(3),R(1)*R(2));
    G = reshape(bartelsStewart(A,[],[],D,E), R);
    
end

   
   
   
%    I0 = lmlragen({U,V,W},G); 
%    
%    s = 10;
%     beta = 1; tol = 1e-3; lambda = 1/s;
%     r1 = 0.1; r2 = 0.1/lambda;
%     Maxit = 100;
% 
%     [I2,k] = beltrami3D(I0,Maxit,beta,lambda,tol,r1,r2);
%    
%     Y_hat(ind3_min(1):ind3_max(1),ind3_min(2):ind3_max(2), :) = I2;
    
    Y_hat(ind3_min(1):ind3_max(1),ind3_min(2):ind3_max(2), :) = lmlragen({U,V,W},G);
    
  end
end


end


