function [e] = ellipseFit3(varargin)

% function [e] = ellipseFit3(V,optMethod,numSample)


%% Parse input

switch nargin
    case 1
        V=varargin{1};
        optMethod=2;
        numSample=[];
    case 2
        V=varargin{1};
        optMethod=varargin{2};
        numSample=[];
    case 3
        V=varargin{1};
        optMethod=varargin{2};
        numSample=varargin{3};
end

%%

mean_V=mean(V,1);
Vk=V-mean_V(ones(size(V,1),1),:);

[Q,~,~]=svd(Vk',0);
Vf=(Q'*Vk')';

%Fit ellipse
[A] = ellipseFit(Vf,optMethod,numSample);

[R,~]=euler2DCM([0 0 -A(5)]);

R1=eye(4,4);
R1(1:3,1:3)=R';

T1=eye(4,4);
T1(1:3,4)=[A(1) A(2) 0];

R2=eye(4,4);
R2(1:3,1:3)=Q;

T2=eye(4,4);
T2(1:3,4)=mean_V;

M=T2*R2*T1*R1';

%%

e.centre=M(1:3,4)'; 
e.radii=A(3:4);
e.axes=M(1:3,1:3);
e.tform=M;
