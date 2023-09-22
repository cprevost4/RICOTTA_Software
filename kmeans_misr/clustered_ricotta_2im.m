function EZ = clustered_ricotta_2im(tensor_Y1,tensor_Y3,D3,R,patch_nb,par,nb_cluster,bparams)

K = size(tensor_Y1{1},2);
patchsize = sqrt(size(tensor_Y3{1},1));

%storage for the reconstructed patches
Z = zeros(patchsize, patchsize, K, bparams.block_num(1)*bparams.block_num(2));

for k = 1:nb_cluster

         Y31 = tens2mat(tensor_Y3{k},1);
         Y31 = unique(Y31','rows'); Y31 = Y31';
         par.K=min(R(1),min(size(Y31,1), size(Y31,2)));
         U3 = Nonnegative_DL(Y31,par);

         Y12 = tens2mat(tensor_Y1{k},2);
         Y12 = unique(Y12','rows'); Y12 = Y12';
         par.K=min(R(2),min(size([Y12],1), size([Y12],2)));
         V1 = Nonnegative_DL([Y12],par);
% 
         Y33 = tens2mat(tensor_Y3{k},3);
         Y33 = unique(Y33','rows'); Y33 = Y33';
         %par.K = R(3);
         %par.K=min(R(3),min(size(Y33,1), size(Y33,2)));
         par.K=min(size(Y33,1), size(Y33,2));
         W3=Nonnegative_DL(Y33, par );

         %W3 = vca(Y33,par.K);


    % opt. c %
    G = sparse_tucker({U3,D3*V1,W3},tensor_Y3{k},par.mu);

%     s = 10;
%     beta = 1; tol = 1e-3; lambda = 1/s;
%     r1 = 0.1; r2 = 0.1/lambda;
%     Maxit = 2000;
    
     %patch2 = beltrami3D(lmlragen({U3,V1,W3},G),Maxit,beta,lambda,tol,r1,r2);
     %patch2 = BeltramiPD(lmlragen({U3,V1,W3},G),Maxit,beta,lambda,tol,r1,r2);

     patch2 = lmlragen({U3,V1,W3},G);


    Z(:,:,:,patch_nb{k}) = reshape(patch2,patchsize,patchsize,K,[]);





end

EZ = JointBlocks(Z, bparams);

end