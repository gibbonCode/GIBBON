function [varargout]=patchEdgeCrossProduct(F,V)

% function [C,CX,CY,CZ]=patchEdgeCrossProduct(F,V)
% ------------------------------------------------------------------------
% 
%
% ------------------------------------------------------------------------

%%
% Deal with 2D patch data
if size(V,2)==2
    V(:,3)=0;
end

%% Get coordinate components
X=V(:,1); Y=V(:,2); Z=V(:,3);

%Get coordinate components at nodal indices for faces
XF=X(F); YF=Y(F); ZF=Z(F);

%Fix issue for single face input
if size(F,1)==1
    XF=XF'; YF=YF'; ZF=ZF';
end

%% Compute cross-products on patch edges. Notation follows form A X B

%Set up indices
indA=1:size(F,2); %Indices for "A vectors" (e.g. first)
indB=[2:size(F,2) 1]; %Indices for "B vectors" (e.g. second)

%Compute A and B components arrays for all edge points
A1=XF(:,indA); B1=XF(:,indB); %X components
A2=YF(:,indA); B2=YF(:,indB); %Y components
A3=ZF(:,indA); B3=ZF(:,indB); %Z components

%Compute cross product components
CX=(A2.*B3-A3.*B2);
CY=(A3.*B1-A1.*B3);
CZ=(A1.*B2-A2.*B1);

%Sum up cross products to form normal direction
C=[sum(CX,2) sum(CY,2) sum(CZ,2)]/2; 

%% Collect ouput 

varargout{1}=C;
if nargout>1
   varargout{2}=CX/2;
   varargout{3}=CY/2;
   varargout{4}=CZ/2;
end

