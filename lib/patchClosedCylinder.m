function [FT,VT,CT]=patchClosedCylinder(inputStruct)



% 2017/04/18

%%
if isfield(inputStruct,'numDigitKeep')
    numDigitKeep=inputStruct.numDigitKeep;
else
    numDigitKeep=5;     
end

%%
cylRadius=inputStruct.cylRadius;
pointSpacing=inputStruct.pointSpacing;
numRadial=round((2*pi*cylRadius)/pointSpacing);
cylHeight=inputStruct.cylHeight;
numHeight=round(inputStruct.cylHeight/((pointSpacing/2)*sqrt(3)));
numHeight=numHeight+iseven(numHeight); %Force uneven for 'tri' method

%%

t=linspace(0,2*pi,numRadial+1);
t=t(1:end-1);
x=cylRadius*cos(t);
y=cylRadius*sin(t);
Vc=[x(:) y(:)];
Vc(:,3)=0; 

cPar.numSteps=numHeight;
cPar.depth=cylHeight; 
cPar.patchType='tri'; 
cPar.dir=0;
cPar.closeLoopOpt=1; 

[F,V]=polyExtrude(Vc,cPar);
C=1*ones(size(F,1),1); %Color for side faces

%%
% Indices for the top and bottom points can be obtained as follows
indTop=numHeight:numHeight:size(V,1);
indBottom=1:numHeight:size(V,1);

%%
% The top and bottom can be meshed using |regionTriMesh2D|

[Ft,Vt]=regionTriMesh2D({V(indTop,[1 2])},[],0);
Vt(:,3)=mean(V(indTop,3));
Ct=2*ones(size(Ft,1),1); %Color for top faces

[Fb,Vb]=regionTriMesh2D({V(indBottom,[1 2])},[],0);
Vb(:,3)=mean(V(indBottom,3));
Fb=fliplr(Fb);
Cb=3*ones(size(Fb,1),1); %Color for bottom faces 

%%

%Join sets 
[FT,VT,CT]=joinElementSets({F,Ft,Fb},{V,Vt,Vb},{C,Ct,Cb});

%Merges nodes
[~,indUnique,indFaceIndexFix]=unique(pround(VT,numDigitKeep),'rows');
VT=VT(indUnique,:);
FT=indFaceIndexFix(FT);

