function IND=box_indices(siz)

IND_all=reshape(1:1:prod(siz),siz(1),siz(2),siz(3));

r1=IND_all(1,:,:); 
r2=IND_all(end,:,:); 
c1=IND_all(:,1,:); 
c2=IND_all(:,end,:); 
s1=IND_all(:,:,1); 
s2=IND_all(:,:,end); 

IND=unique([r1(:);r2(:);c1(:);c2(:);s1(:);s2(:)]);