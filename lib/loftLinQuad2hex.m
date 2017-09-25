function [varargout]=loftLinQuad2hex(Fq,Vq,Vq2,numSteps)

%Get coordinates
X=linspacen(Vq(:,1),Vq2(:,1),numSteps+1);
Y=linspacen(Vq(:,2),Vq2(:,2),numSteps+1);
Z=linspacen(Vq(:,3),Vq2(:,3),numSteps+1);

%Collect node set
V=[X(:) Y(:) Z(:)];

%Create element matrix
E=repmat(Fq,[numSteps,2]);
E_add=0:size(Vq,1):size(Vq,1)*(numSteps-1); 
E_add=E_add(ones(size(Fq,1),1),:);
E_add=E_add(:);
E_add=E_add(:,ones(4,1));
E_add=[E_add E_add+size(Vq,1)];
E=E+E_add;

%Create boundary face set
Fq1=E(1:size(Fq,1),1:4);
Fq2=E(1+(end-size(Fq,1)):end,5:end);

%Collect output
varargout{1}=E;
varargout{2}=V;
varargout{3}=Fq1;
varargout{4}=Fq2;