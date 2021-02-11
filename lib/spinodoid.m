function [varargout]=spinodoid(inputStruct)

% -----------------------------------------------------------------------
% This function generates Spinodoid microstructures which are non-periodic
% and stochastic bi-continous microstructures that approximate the
% topology observed in spinodal decomposition. The microstructures are
% generated using a Gaussian random field (GRF) and can be made anisotropic
% by tuning the anisotropy of the GRF.
%
% Based on / how to cite: 
% Kumar, S., Tan, S., Zheng, L., Kochmann, D.M. Inverse-designed spinodoid 
% metamaterials. npj Comput Mater 6, 73 (2020). 
% https://doi.org/10.1038/s41524-020-0341-6
%
% Input structure and default values:
% inputStruct.isocap=true; % option to cap the isosurface
% inputStruct.domainSize=1; % domain size
% inputStruct.resolution=100; % resolution for sampling GRF
% inputStruct.waveNumber=15*pi; % GRF wave number
% inputStruct.numWaves=1000; % number of waves in GRF
% inputStruct.relativeDensity=0.5; % relative density: between [0.3,1]
% inputStruct.thetas=[15 15 0]; % conical half angles (in degrees) along 
%                                xyz axes for controlling the anisotropy. 
%                                Note: each entry must be either 0 or
%                                between [15,90] degrees.
% inputStruct.R = eye(3); % Rotate the GRF, R must be SO(3)
%
%
% Original author: Siddhant Kumar, September 2020
% (contact: siddhantk41@gmail.com)
%
% Change log: 
% 2020/09/19 Siddhant Kumar: Created original code
% -----------------------------------------------------------------------

%% Parse input

%Create default structure
defaultInputStruct.isocap=true; % option to cap the isosurface
defaultInputStruct.domainSize=1; % domain size
defaultInputStruct.resolution=100; % resolution for sampling GRF
defaultInputStruct.waveNumber=15*pi; % GRF wave number
defaultInputStruct.numWaves=1000; % number of waves in GRF
defaultInputStruct.relativeDensity=0.5; % relative density: between [0.3,1]
defaultInputStruct.thetas=[15 15 0]; % conical half angles (in degrees) 
%                                     along xyz axes for controlling the 
%                                     anisotropy. Note: each entry must be 
%                                     either 0 or between [15,90] degrees.
defaultInputStruct.R = eye(3); % Rotate the GRF, R must be SO(3)
defaultInputStruct.ignoreChecks = false; %Ignore checks on parameter values

%Complete input with default if incomplete
[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

%Get parameters from input structure
isocap = inputStruct.isocap; % option to cap the isosurface
domainSize = inputStruct.domainSize; % domain size
resolution = inputStruct.resolution; % resolution for sampling GRF
waveNumber = inputStruct.waveNumber; % GRF wave number
numWaves = inputStruct.numWaves; % number of waves in GRF
relativeDensity = inputStruct.relativeDensity; % relative density: [0.3,1]
thetas= inputStruct.thetas; % conical half angles (in degrees)
R = inputStruct.R; % Rotate the GRF, R must be SO(3)
ignoreChecks = inputStruct.ignoreChecks; %Ignore checks on parameter values

%% Input checks
if(ignoreChecks)
    warning(['Ignoring all checks on parameter values. ',...
        'Unreasonable parameters may give wrong topologies.']);
else
    if((relativeDensity<0.3) || (relativeDensity>1.0))
        error('relativeDensity must be between [0.3,1]')
    end
    if((any(thetas<0)) || (any(thetas>90)))
        error('thetas must be either 0 or between [15,90] degrees')
    end
    for i=1:3
        if((thetas(i)>0) && (thetas(i)<15))
            error('thetas must be either 0 or between [15,90] degrees')
        end
    end

    if(size(R,1)~=size(R,2))
        error('Rotation matrix is not square')
    end
    if(size(R,1)~=3)
        error('Rotation matrix must be 3x3 in size')
    end
    if(abs(det(R)-1)>1e-8)
        error('Rotation matrix: det(R)~=1')
    end
    if(norm(R'*R-eye(3))>1e-8)
        error('Rotation matrix is not orthogonal')
    end
end

%% Define rotated axes
axes1 = [1,0,0];
axes2 = [0,1,0];
axes3 = [0,0,1];

axes1 = (R*axes1')';
axes2 = (R*axes2')';
axes3 = (R*axes3')';

%% Generate wave directions for GRF
%array of all wave directions

waveDirections = zeros(numWaves,3);
for i=1:numWaves
    flag = true; %keep trying until candidate wave vector is found
    while(flag)
        % generate isotropically sampled candidate wave
        candidate = randn(1,3);
        candidate = candidate/norm(candidate);
        % check for allowed wave vector directions
        % angle along first axis
        angle1 = min(...
            acosd(dot(candidate,axes1)),...
            acosd(dot(candidate,-axes1)));
        % angle along second axis
        angle2 = min(...
            acosd(dot(candidate,axes2)),...
            acosd(dot(candidate,-axes2)));
        % angle along third axis
        angle3 = min(...
            acosd(dot(candidate,axes3)),...
            acosd(dot(candidate,-axes3)));
        % check
        if(any([angle1,angle2,angle3]<thetas))
            flag = false;
            break;
        end
    end
    waveDirections(i,:)=candidate;
end

%% Generate wave phase angles for GRF
wavePhases = rand_angle([numWaves,1]); %2*pi*rand(numWaves,1);

%% Discretize the domain
discretization = linspace(0,domainSize,resolution);
[X,Y,Z] = meshgrid(discretization,discretization,discretization);

%% Evaluate GRF on sampling points

GRF = zeros(size(X));
for i=1:numWaves
    dotProduct = waveDirections(i,1)*X ...
               + waveDirections(i,2)*Y ...
               + waveDirections(i,3)*Z;
             
    GRF = GRF+sqrt(2/numWaves)*cos(dotProduct*waveNumber+wavePhases(i));
end

%% Apply levelset
levelset = sqrt(2)*erfinv(2*relativeDensity-1);

%% Compute isosurface and isocaps

%Compute isosurface
[f,v] = isosurface(X,Y,Z,GRF,-levelset);
c=zeros(size(f,1),1);

%Isocaps
if isocap==1
    %Compute isocaps
    [fc,vc] = isocaps(X,Y,Z,GRF,-levelset);

    nc=patchNormal(fc,vc);
    cc=zeros(size(fc,1),1);
    cc(nc(:,1)<-0.5)=1;
    cc(nc(:,1)>0.5)=2;
    cc(nc(:,2)<-0.5)=3;
    cc(nc(:,2)>0.5)=4;
    cc(nc(:,3)<-0.5)=5;
    cc(nc(:,3)>0.5)=6;    
    
    %Join sets
    [f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc});
    
end

%% Merge nodes and clean-up mesh 

%Merge nodes
[f,v]=mergeVertices(f,v); 

%Check for unique faces
[~,indUni,~]=unique(sort(f,2),'rows');
f=f(indUni,:); %Keep unique faces
c=c(indUni);

%Remove collapsed faces
[f,logicKeep]=patchRemoveCollapsed(f); 
c=c(logicKeep);

%Remove unused points
[f,v]=patchCleanUnused(f,v); 

%Invert faces
f=fliplr(f); 

%% Collect output
varargout{1}=f;
varargout{2}=v;
varargout{3}=c;
varargout{4}=GRF;

end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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


