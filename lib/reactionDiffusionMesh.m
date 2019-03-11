function [A_out,B_out]=reactionDiffusionMesh(varargin)

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        controlPar=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        controlPar=varargin{3};
end

% Initial values
L=V(:,3)>mean(V(:,3));
defaultControlPar.A=double(~L); 
defaultControlPar.B=double(L); 

% Kill rates
defaultControlPar.f=0.055;
defaultControlPar.k=0.062;

% Diffusion rates
defaultControlPar.da=1;
defaultControlPar.db=0.5;

% Time stepping parameters
defaultControlPar.timeTotal = 1000; %Final time
defaultControlPar.dt=1; %Time step size

%
defaultControlPar.waitbar=1;
defaultControlPar.numSaveSteps=25;

defaultControlPar.maxAngleDeviation=45*(pi/180);
defaultControlPar.selectionMethod='best'; %'Random'
defaultControlPar.triangleConvert=1;
defaultControlPar.fourConnectConvert=1;
[controlPar]=structComplete(controlPar,defaultControlPar,1); %Complement provided with default if missing or empty

%%

A=controlPar.A;
B=controlPar.B;
da=controlPar.da;
db=controlPar.db;
timeTotal=controlPar.timeTotal;
dt=controlPar.dt;
waitbarOption=controlPar.waitbar;
f=controlPar.f;
k=controlPar.k;
numSaveSteps=controlPar.numSaveSteps;

%%

[meshConnectivity]=patchConnectivity(F,V);
vertex_vertex_connectivity=meshConnectivity.vertex.vertex;
logicInvalid=vertex_vertex_connectivity==0;
vertex_vertex_connectivity(logicInvalid)=size(V,1)+1; %Additional point to cope with invalid
N=sum(~logicInvalid,2);

%%

nSteps=round(timeTotal/dt);
if numSaveSteps>1
    saveSteps=round(linspace(1,nSteps,numSaveSteps));
else
    saveSteps=nSteps;
end

if waitbarOption
    hw=waitbar(0,'Reaction diffusion computation... ');
end

if numSaveSteps>1
    A_out=zeros(size(A,1),numel(saveSteps));
    B_out=zeros(size(A,1),numel(saveSteps));
    A_out(:,1)=A;
    B_out(:,1)=B;    
    c=2;
else    
    c=1;
end

for q=2:1:nSteps        
    
    An=[A; 0]; %Add zero to cope with invalid points
    Bn=[B; 0]; %Add zero to cope with invalid points
    
    LA=(sum(An(vertex_vertex_connectivity),2)./N)-A;    
    LB=(sum(Bn(vertex_vertex_connectivity),2)./N)-B;
    
    A = A + (da*LA - A.*B.^2 + f*(1-A))*dt;
    B = B + (db*LB + A.*B.^2 - (k+f)*B)*dt; 
    
    if waitbarOption
        waitbar(q/nSteps,hw,['Reaction diffusion computation... ',num2str(round((q/nSteps)*100)),'%']);
    end      
    
    if any(saveSteps==q)
       A_out(:,c)=A;
       B_out(:,c)=B;
       c=c+1;
    end
    
end

if waitbarOption
    close(hw);
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
