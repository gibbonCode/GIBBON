function [Vg]=sweepCurveSmooth(varargin)

switch nargin
    case 5
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numSteps=varargin{5};
        p=[];
        f=0.1;        
    case 6    
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numSteps=varargin{5};
        p=varargin{6};
        f=0.1;
    case 7
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numSteps=varargin{5};
        p=varargin{6};
        f=varargin{7};
end

%%

d=sqrt(sum((V1-V2).^2))*f; %Distance between points

V=[V1-d*n1; V1; V1+d*n1; V2-d*n2; V2; V2+d*n2];

%Compute distance metric used for parametric representation
D=pathLength(V);

%Redefine distance metric for evenly spaced points
Dg=linspace(D(2),D(end-1),numSteps)';

p=csapsPar(V,p); %Scale invariant version of smoothening parameter
W=[1 1e9 1 1 1e9 1]; %Weights vector
Vg=zeros(numSteps,size(V,2)); %Initialize curve
for q=1:1:size(V,2) %Loop across dimensions
    v=V(:,q);
    pp = csaps(D,v,p,[],W)'; %Smoothened ppform
    Vg(:,q)=ppval(pp,Dg); %Evaluate ppform
end

