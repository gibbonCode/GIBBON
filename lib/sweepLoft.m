function [varargout]=sweepLoft(varargin)

% function [F,V,C,S]=sweepLoft(V1,V2,n1,n2,Vg,numSteps,numTwist,plotOn)
%------------------------------------------------------------------------
%
% The |sweepLoft| function creates a swept loft (as in CAD terminology for
% a shape formed by merging a set of sketches towards each other allong a
% given path or guide curve). The inputs to the function are the start and
% end sketchs (V1 and V2), the start and end normal directions (n1 and n2),
% and the guid curve (Vg). Optional additional inputs are: the number of
% steps for the loft feature (numSteps, same as number of points in guide
% curve if not provided), the number of twists (numTwists, default is zero)
% the shape undergoes around the guide curve, and finally plotOn (default
% is 0, i.e. off) which is a logic to turn on or off plotting within the
% function. The function outputs are patch data i.e. faces (F), the
% vertices (V) and face colors (C, denoting step in lofting process).
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2016/08/01 Created
% 2017/06/27 Fixed bug in relation to compensatory/corrective rotation
% 2018/11/07 Updated so users have the option of using an input structure instead
% 2018/11/07 Added close section option, enabling open or closed sections
% sweeps
%------------------------------------------------------------------------

%% Parse input

defaultOptionStruct.numSteps=[];
defaultOptionStruct.numTwist=0;
defaultOptionStruct.plotOn=0;
defaultOptionStruct.closeSection=1; 

switch nargin
    case 1 %Assume an input structure is used
        optionStruct=varargin{1};
        
        %Check optionStruct against default
        [optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty
        
        V1=optionStruct.V1;
        V2=optionStruct.V2;
        n1=optionStruct.n1;
        n2=optionStruct.n2;
        Vg=optionStruct.Vg;
        numSteps=optionStruct.numSteps;
        numTwist=optionStruct.numTwist;
        plotOn=optionStruct.plotOn;
        closeSection=optionStruct.closeSection;
        
    otherwise %Assume inputs are provided seperately
        if nargin<5
            error('Insufficient input parameters');
        end        
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        Vg=varargin{5};
        numSteps=defaultOptionStruct.numSteps;
        numTwist=defaultOptionStruct.numTwist;
        plotOn=defaultOptionStruct.plotOn;     
        closeSection=defaultOptionStruct.closeSection;
        switch nargin
            case 6
                Vg=varargin{5};
                numSteps=varargin{6};
            case 7
                Vg=varargin{5};
                numSteps=varargin{6};
                numTwist=varargin{7};
            case 8
                numSteps=varargin{6};
                numTwist=varargin{7};
                plotOn=varargin{8};
        end
end

if isempty(numSteps)
    numSteps=size(Vg,1);
end

p1=Vg(1,:);
p2=Vg(end,:);

%%
%Resample guide curve evenly using numSteps
[Vg] = evenlySampleCurve(Vg,numSteps,'linear',0);

stepSize=sqrt(sum((Vg(1,:)-Vg(2,:)).^2,2));
%%

if plotOn
    lineWidth1=6;
    lineWidth2=2;    
    lineWidth3=1; 
    markerSize=25;
    triadSize=stepSize*7; 
    fontSize=20;
    
    hf=cFigure; hold on;
    ht=gtitle('Mapping coordinate systems along guide curve',fontSize);    
    hp(1)=plotV(V1,'r-','LineWidth',lineWidth1);
    hp(2)=plotV(V2,'g-','LineWidth',lineWidth1);
    hp(3)=plotV(Vg,'k--','LineWidth',lineWidth2);
    hp(4)=quiverVec(p1,n1,triadSize,'r');        
    hp(5)=quiverVec(p2,n2,triadSize,'g');
    
    axisGeom(gca,fontSize); camlight headlight; 
    axis manual; 
    colormap gjet;
    ha=gca;
    ha.GridAlpha=0.25;
    ha.LineWidth=1;
    drawnow;
end

%% Define allong curve coordinate systems

Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
U=Uf;
U(:,:,2)=Ub;
U=gnanmean(U,3);
U=vecnormalize(U);

mean_V1=mean(V1,1);
V1m=vecnormalize(V1-mean_V1(ones(size(V1,1),1),:));
e3=n1;
e1=vecnormalize(V1m(1,:));
e2=vecnormalize(cross(e3,e1));
e1=vecnormalize(cross(e2,e3));
R1=[e1;e2;e3];

mean_V2=mean(V2,1);
V2m=vecnormalize(V2-mean_V2(ones(size(V2,1),1),:));
e3=n2;
e1=vecnormalize(V2m(1,:));
e2=vecnormalize(cross(e3,e1));
e1=vecnormalize(cross(e2,e3));
R2=[e1;e2;e3];

R_curve=repmat(eye(3,3),[1,1,size(Vg,1)]);
R_curve(:,:,1)=R1;

R=R1;

if plotOn==1
    hTriad=gobjects(1,1);
end

for q=2:size(Vg,1)
    a=U(q-1,:);
    b=U(q,:);
    
    theta=real(acos(dot(a,b))); %Complex if dot product is out of range [-1.1] due to precission issues
    w=vecnormalize(cross(b,a));
    
    if norm(w)>0.5
        [Rn]=vecAngle2Rot(theta,w);
        R=R*Rn;
        if plotOn
            figure(hf);
            delete(hTriad);
            hTriad=quiverTriad(Vg(q,:),R',triadSize);
            drawnow;
        end
    end
    R_curve(:,:,q)=R;
end

%% Create coordinate matrices

V1s=V1-p1(ones(size(V1,1),1),:);
V1s=V1s*R1';

V2s=V2-p2(ones(size(V2,1),1),:);
V2s=V2s*R2';

% Create coordinate matrices
X=linspacen(V1s(:,1),V2s(:,1),numSteps)';
Y=linspacen(V1s(:,2),V2s(:,2),numSteps)';
Z=linspacen(V1s(:,3),V2s(:,3),numSteps)';

Vs=[X(:) Y(:) Z(:)];

X=reshape(Vs(:,1),size(X));
Y=reshape(Vs(:,2),size(Y));
Z=reshape(Vs(:,3),size(Z));

%%

if plotOn==1
    h1=gobjects(1,numSteps);
end

for q=1:1:numSteps
    
    V2p=[X(q,:)' Y(q,:)' Z(q,:)'];
    V2p=(R_curve(:,:,q)'*V2p')';
    V2p=V2p+Vg(q*ones(size(V2p,1),1),:);
    
    if plotOn==1
        figure(hf);
        ht.String='Initial morphing and sweeping of sections along guide curve';
        h1(q)=plotV(V2p([1:end 1],:),'k-','LineWidth',lineWidth1);        
        if q>1
            h1(q-1).LineWidth=lineWidth3;
        end
        drawnow;
    end
    X(q,:)=V2p(:,1);
    Y(q,:)=V2p(:,2);
    Z(q,:)=V2p(:,3);
end

%%

V2p=[X(end,:)' Y(end,:)' Z(end,:)'];

if size(V2,1)==3
    [V2_temp] = evenlySampleCurve(V2,size(V2,1)*3,'linear',0);
    [V2p_temp] = evenlySampleCurve(V2p,size(V2,1)*3,'linear',0);
    [~,Rc]=rigidTransformationMatrixDirect(V2_temp,V2p_temp);
else
    [~,Rc]=rigidTransformationMatrixDirect(V2,V2p);
end

[theta,w]=rot2VecAngle(Rc);
% theta=0; w=n2';

%%
W=repmat(w(:)',[numSteps,1]);
wt=w;
if plotOn==1
    hVec=gobjects(1,1);
    delete(hTriad);
end
for q=numSteps-1:-1:1
    mean_V_now=mean([X(q,:)' Y(q,:)' Z(q,:)'],1);
    wt=(wt'*R_curve(:,:,q+1)'*R_curve(:,:,q))';
    if plotOn==1
        figure(hf);        
        
        delete(h1(q));
        ht.String='Back tracking and fixing initial orientation mismatch';
        delete(hVec);
        hVec=quiverVec(mean_V_now,wt',triadSize,'k'); drawnow;
    end
    W(q,:)=wt(:)';
end

%%

V2p_c=mean(V2p,1);
V2p=V2p-V2p_c(ones(size(V2p,1),1),:);
V2p=(Rc'*V2p')';
V2p=V2p+V2p_c(ones(size(V2p,1),1),:);

if plotOn==1
    figure(hf);
    delete(hVec);
    delete(h1(end));
    plotV(V2p,'k--','LineWidth',lineWidth1);    
    drawnow;
end

%%

theta_w=linspace(0,1,size(Vg,1));
theta_step=theta*theta_w;

if numTwist>0
    theta_step_twist=linspace(0,2*pi*numTwist,size(Vg,1));
end

if plotOn==1
    h2=gobjects(1,size(Vg,1));   
    hTriad2=gobjects(1,1);
end

for q=1:1:size(Vg,1)
    Vn=[X(q,:)' Y(q,:)' Z(q,:)'];
    Vn_mean=mean(Vn,1);
    Vn=Vn-Vn_mean(ones(size(Vn,1),1),:);
    
    [Rc]=vecAngle2Rot(theta_step(q),W(q,:));
    Vn=(Rc'*Vn')';
    
    if numTwist>0
        [Rc]=vecAngle2Rot(theta_step_twist(q),U(q,:));
        Vn=(Rc'*Vn')';
    end
    
    Vn=Vn+Vn_mean(ones(size(Vn,1),1),:);
    
    if plotOn==1        
        figure(hf);
        ht.String='Correct orientations and add potential twist';
        h2(q)=plotV(Vn([1:end 1],:),'k-','LineWidth',lineWidth1);        
        if q>1
            h2(q-1).LineWidth=lineWidth3;
        end        
        delete(hTriad2);
        hTriad2=quiverTriad(Vg(q,:),(R_curve(:,:,q)*Rc)',triadSize);
        drawnow;
    end
    
    X(q,:)=Vn(:,1);
    Y(q,:)=Vn(:,2);
    Z(q,:)=Vn(:,3);
end

if plotOn==1
    delete(hTriad2);
end

%% Override start and end with input curves

X(1,:)=V1(:,1);
Y(1,:)=V1(:,2);
Z(1,:)=V1(:,3);

X(end,:)=V2(:,1);
Y(end,:)=V2(:,2);
Z(end,:)=V2(:,3);

%% Create patch data

%Color data
c=(1:1:size(Z,1))';
C=c(:,ones(1,size(Z,2)));

%Create quad patch data
[F,V,C] = surf2patch(X,Y,Z,C);

%Close section if required
if closeSection==1
    I=[(1:size(Z,1)-1)' (1:size(Z,1)-1)' (2:size(Z,1))' (2:size(Z,1))' ];
    J=[size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1)];
    F_sub=sub2ind(size(Z),I,J);
    F=[F;F_sub];    
    C(end-size(F_sub,1):end,:)=C(end-size(F_sub,1):end,:)+0.5;   
end
[C]=vertexToFaceMeasure(F,C); %Convert vertex colors to face colors
C=round(C)-1;

if plotOn==1    
    figure(hf);
    delete(h2);
    ht.String='Loft surface';
    gpatch(F,V,C,'k',1);        
    drawnow;
end

%%

varargout{1}=F;
varargout{2}=V;
switch nargout
    case 3
        varargout{3}=C;
    case 4
        varargout{3}=C;
        
        S.X=X;
        S.Y=Y;
        S.Z=Z;
        varargout{4}=S;
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
