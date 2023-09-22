function EZ = kmeans_misr(Y1,Y2,Y3,D1,D2,D3,patchsize,overlap,nb_cluster,par,R,nIm)


[tensor_Y1,tensor_Y2,tensor_Y3,patch_nb] = coupled_clustering(Y1,Y2,Y3,patchsize,overlap,nb_cluster);

%parameters for splitting the dataset
bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz = [overlap overlap];
I = size(Y2,1); J = size(Y1,2); K = size(Y1,3); K3 = size(Y3,3); 
d1 = I/size(Y1,1); d2 = J/size(Y2,2);

%number of patches
num1 = floor((I-patchsize)/(patchsize-overlap)+1);
num2 = floor((J-patchsize)/(patchsize-overlap)+1);
bparams.block_num = [num1 num2]; 

if nIm == 3
    EZ = clustered_ricotta(tensor_Y1,tensor_Y2,tensor_Y3,D3,R,patch_nb,par,nb_cluster,bparams);
else if nIm == 2
    EZ = clustered_ricotta_2im(tensor_Y1,tensor_Y3,D3,R,patch_nb,par,nb_cluster,bparams);
end


end
