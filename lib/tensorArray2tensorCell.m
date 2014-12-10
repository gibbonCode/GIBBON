function [H]=tensorArray2tensorCell(H_mat,siz)


nDims=sqrt(size(H_mat,2));

H_mat=reshape(H_mat',nDims,nDims,size(H_mat,1));
        
H=reshape(mat2cell(H_mat,nDims,nDims,ones(size(H_mat,3),1)),siz);