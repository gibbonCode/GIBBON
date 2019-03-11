function [V_coil_rep]=gaborCoil(varargin)

% function [V_coil_rep]=gaborCoil(V,E,optStruct)
%-------------------------------------------------------------------------
% This function generates a so-called Gabor coil, i.e. a spiral which is
% modulated similar to a Gabor filter. 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2017/01/01
% 2018/03/27 
%------------------------------------------------------------------------
%%

defaultOpt.numSteps=100;
defaultOpt.numTwist=4;
defaultOpt.coilAmplitude=[];
defaultOpt.coilAmplitudeFactor=1/10;
defaultOpt.f=4;
defaultOpt.funcMethod=2;

switch nargin
    case 2
        V=varargin{1};
        E=varargin{2};
        optStruct=defaultOpt;
    case 3
        V=varargin{1};
        E=varargin{2};
        optStruct=varargin{3};
end

%%

[optStruct]=structComplete(optStruct,defaultOpt,0);

numSteps=optStruct.numSteps;
numTwist=optStruct.numTwist;
coilAmplitude=optStruct.coilAmplitude;
coilAmplitudeFactor=optStruct.coilAmplitudeFactor;
f=optStruct.f;
funcMethod=optStruct.funcMethod;

%%

V_edgeOrigins=V(E(:,1),:);

A=V(E(:,2),:)-V(E(:,1),:); %Edge Vectors

b=[1 0 0]; %Coil standard vector
B=b(ones(size(A,1),1),:);
c=[0 0 1];

A_mag=sqrt(sum(A.^2,2));

thetaRot=acos(dot(A,B,2)./A_mag);

vecRot=cross(A,B);
vecRot_mag=sqrt(sum(vecRot.^2,2));
logicReplace=vecRot_mag<eps(pi);
vecRot(logicReplace,:)=c(ones(nnz(logicReplace),1),:);
vecRot=vecnormalize(vecRot);

if isempty(coilAmplitude) %If coil amplitudes are empty
    coilAmplitude=A_mag*coilAmplitudeFactor; %Set as local edge length multiplied by amplitude factor
end

%% Create standard coil

%Create periodic coil component
t=linspace(0,2*pi*numTwist,numSteps);
x=linspace(0,1,numSteps);
y=cos(t);
z=sin(t);

%Create Gaussian envelope 
xg = linspace(-f,f,numSteps);
hg=exp(-(xg.^2)./2);

%Modulate x and y coordinates to create Gabor coil
y=y.*hg;
z=z.*hg;

%The coil coordinate set
V_coil=[x(:) y(:) z(:)];

%Resample evenly
[V_coil] = evenlySampleCurve(V_coil,numSteps,'pchip',0); 

%Force ends to match edge vertices
V_coil(1,[2 3])=0;
V_coil(end,[2 3])=0;

%% 
        
numEdges=size(E,1);
V_coil_rep=repmat(V_coil,[1,1,numEdges]); %Coils copied in 3rd dim

if numel(coilAmplitude)==1
    coilAmplitude=coilAmplitude*ones(numEdges,1);
end

switch funcMethod
    case 1 %For loop     
        for q=1:1:numEdges
            R=vecAngle2Rot(thetaRot(q),vecRot(q,:));
            V_coil_now=V_coil;
            V_coil_now(:,1)=V_coil_now(:,1)*A_mag(q); %Scaling lenght
            V_coil_now(:,[2 3])=V_coil_now(:,[2 3])*coilAmplitude(q); %Scaling width
            V_coil_now=V_coil_now*R;
            V_coil_rep(:,:,q)=V_coil_now+V_edgeOrigins(q*ones(numSteps,1),:);
        end
    case 2 %Vectorized
        
        % Scaling coil lengths        
        coilScale=repmat(permute(A_mag,[3 2 1]),[numSteps,1,1]); %Coil scaling factors
        V_coil_rep(:,1,:)=V_coil_rep(:,1,:).*coilScale; %Scaling
        
        % Scaling coil amplitudes       
        coilScale=repmat(permute(coilAmplitude,[3 2 1]),[numSteps,2,1]); %Coil scaling factors
        V_coil_rep(:,[2 3],:)=V_coil_rep(:,[2 3],:).*coilScale; %Scaling
        
        % Rotating coils same form as R=vecAngle2Rot(theta,w); V_rot=V*R       
        x=V_coil_rep(:,1,:);
        y=V_coil_rep(:,2,:);
        z=V_coil_rep(:,3,:);
        theta=repmat(permute(thetaRot,[3 2 1]),[numSteps,1,1]); %Coil angles
        w1=repmat(permute(vecRot(:,1),[3 2 1]),[numSteps,1,1]); %Coil rotation vector x
        w2=repmat(permute(vecRot(:,2),[3 2 1]),[numSteps,1,1]); %Coil rotation vector y
        w3=repmat(permute(vecRot(:,3),[3 2 1]),[numSteps,1,1]); %Coil rotation vector z
        
        X=y.*(w3.*sin(theta) - w1.*w2.*(cos(theta) - 1)) - z.*(w2.*sin(theta) + w1.*w3.*(cos(theta) - 1)) + x.*((w2.^2 + w3.^2).*(cos(theta) - 1) + 1);
        Y=z.*(w1.*sin(theta) - w2.*w3.*(cos(theta) - 1)) - x.*(w3.*sin(theta) + w1.*w2.*(cos(theta) - 1)) + y.*((w1.^2 + w3.^2).*(cos(theta) - 1) + 1);
        Z=x.*(w2.*sin(theta) - w1.*w3.*(cos(theta) - 1)) - y.*(w1.*sin(theta) + w2.*w3.*(cos(theta) - 1)) + z.*((w1.^2 + w2.^2).*(cos(theta) - 1) + 1);
        
        V_coil_rep(:,1,:)=X;
        V_coil_rep(:,2,:)=Y;
        V_coil_rep(:,3,:)=Z;
        
        % Translating coils
        
        V_edgeOrigins_rep=repmat(permute(V_edgeOrigins,[3 2 1]),[numSteps,1,1]);
        V_coil_rep=V_coil_rep+V_edgeOrigins_rep;
end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
