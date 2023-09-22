function [tensor_Y1,tensor_Y2,tensor_Y3,patch_nb] = coupled_clustering(Y1,Y2,Y3,patchsize,overlap,nb_cluster)

%parameters for splitting the dataset
bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz = [overlap overlap];
I = size(Y2,1); J = size(Y1,2); K = size(Y1,3); K3 = size(Y3,3); 
d1 = I/size(Y1,1); d2 = J/size(Y2,2);

%number of patches
num1 = floor((I-patchsize)/(patchsize-overlap)+1);
num2 = floor((J-patchsize)/(patchsize-overlap)+1);
bparams.block_num = [num1 num2]; 

%storage for the reconstructed patches
Z = zeros(patchsize, patchsize, K, num1*num2);

%split MSI
patches_ref = ExtractBlocks(Y3, bparams); 
patches_vec = tens2mat(patches_ref,4,[]);

 %clustering
 para.K = nb_cluster;
 cluster_ind = fkmeans(patches_vec,para.K,'careful');

 for k = 1:para.K

     % form blocks of MSI
     patch_nb{k} = find(cluster_ind == k);
     blocks_Y3{k} = patches_ref(:,:,:,patch_nb{k});

     d_HSI = d1;
     c1 = 1; c2 = 1; %counter

    % FORM BLOCKS OF Y1 AND Y2
    for i = 1:length(patch_nb{k})

         [x,y] = ind2sub([num1,num2],patch_nb{k}(i)); %il s'agit du patch numéro (x,y)
         xx = 1 + (x - 1)*(patchsize - overlap);
         yy = 1 + (y - 1)*(patchsize - overlap); % le patch (x,y) correspond au bloc (xx,yy)
       
        if (mod(x,d_HSI)~=0) %&& (mod(n,d_HSI)~=0)
            kk= Y1(ceil(x/d_HSI), yy:yy+patchsize-1,:);
            blocks_Y1{k}(:,:,:,c1)=kk;
            c1=c1+1;
        else
           kk = Y1(ceil(x/d_HSI), yy:yy+patchsize-1,:);
           blocks_Y1{k}(:,:,:,c1)=kk;
           c1=c1+1;
    
           kk= Y1(ceil(x/d_HSI)+1, yy:yy+patchsize-1,:);
           blocks_Y1{k}(:,:,:,c1)=kk;
           c1=c1+1;
        end

        if (mod(y,d_HSI)~=0)
            kk= Y2(xx:xx+patchsize-1,ceil(y/d_HSI),:);
            blocks_Y2{k}(:,:,:,c2)=kk;
            c2=c2+1;
        else
           kk= Y2(xx:xx+patchsize-1,ceil(y/d_HSI),:);
           blocks_Y2{k}(:,:,:,c2)=kk;
           c2=c2+1;
    
           kk= Y2(xx:xx+patchsize-1,ceil(y/d_HSI)+1,:);
           blocks_Y2{k}(:,:,:,c2)=kk;
           c2=c2+1;
        end
    end

    % ------- VECTORIZING THE TWO SPATIAL DIMENSIONS --------------------%
    tensor_Y3{k} = reshape(blocks_Y3{k},[patchsize^2,K3,length(patch_nb{k})]);
    tensor_Y1{k} = reshape(blocks_Y1{k},patchsize,K,[]);
    tensor_Y2{k} = reshape(blocks_Y2{k},patchsize,K,[]);

end