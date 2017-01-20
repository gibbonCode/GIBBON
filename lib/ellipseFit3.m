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

%Centre on mean
V_offset=mean(V,1);
Vk=V-V_offset(ones(size(V,1),1),:);

%Rotate to 2D problem
[Q,~,~]=svd(Vk',0);
Vf=Vk*Q; %Rotate to XY plane

%Fit ellipse
[A] = ellipseFit(Vf,optMethod,numSample);
[R,~]=euler2DCM([0 0 -A(5)]);
T1=eye(4,4);
T1(1:2,4)=A(1:2);
R1=eye(4,4);
R1(1:3,1:3)=R;

T2=eye(4,4);
T2(1:3,4)=V_offset(:);
R2=eye(4,4);
R2(1:3,1:3)=inv(Q);

M=T1*R1*T2*R2;
e_centre=M(1:3,4)';
Q=M(1:3,1:3);

%%

e.centre=e_centre; 
e.radii=A(3:4);
e.axes=Q;



