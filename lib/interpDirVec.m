function [varargout]=interpDirVec(interpStruct,Vf,P_I)

% function [Vf_I,interpFun]=interpDirVec(interpStruct,Vf,P_I)
% ------------------------------------------------------------------------
% This function interpolates a direction vector field. Direction vectors
% are here defined as vectors indicating for instance a fibre direction
% such that for a vector 'v' actually defines the same fibre direction as
% the vector '-v'. Hence for interpolation special care must be taken to
% treat this property. 
% Firstly structure tensors are composed of the form V=v'*v. This tensor
% field is then interpolated (in 6 steps due to symmetry). The 1st
% principal eigenvectors of the interpolated tensors then forms the
% interpolated direction vector set. 
% The input is the a variable here called interpStruct. If interpStruct is
% a structure array it and may contain: 
% interpStruct.Points 
% interpStruct.Method (default is 'linear' if not provided)
% interpStruct.ExtrapolationMethod (default is 'none' if not provided)
% If however interpStruct is of the scatteredInterpolant class it will in
% addition contain:
% interpStruct.Values
% If interpStruct is not of the scatteredInterpolant class such a
% representation will be constructed.  
% The inputs are the position vectors interpStruct.Points (coordinates of
% vector origins), the direction or fibre direction vectors Vf, the
% position vectors P_I defining the interpolation coordinates,
% interpStruct.Method defining an interpolation method (default is
% 'linear') for the scatteredInterpolant function
% ('linear','natural','nearest'), and similarly
% interpStruct.ExtrapolationMethod defining the extrapolation method
% (default is 'none'). 
% 
% See also: scatteredInterpolant
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/12/10
%------------------------------------------------------------------------

%% Parse input

isInterpFun= isa(interpStruct,'scatteredInterpolant');

if ~isInterpFun
    if ~isfield(interpStruct,'Method')
        interpStruct.Method='linear';    
    end
    if ~isfield(interpStruct,'ExtrapolationMethod')
        interpStruct.ExtrapolationMethod='none';
    end
end

%%

%Create structure tensors through diadic product in nx6 form (A= v'*v)
A=[Vf.^2 Vf(:,2).*Vf(:,3) Vf(:,1).*Vf(:,3) Vf(:,1).*Vf(:,2)];

%Initialise fitted values
A_fit=zeros(size(P_I,1),6); 

if isInterpFun
    %Use interpStruct as interpFun
    interpFun=interpStruct;
    interpFun.Values=A(:,1); %Update values in interpolation structure
else
    %Initialize interpolation function on first scalar entry
    interpFun=scatteredInterpolant(interpStruct.Points,A(:,1),interpStruct.Method,interpStruct.ExtrapolationMethod);
end

%Interpolate for each of the 6 entries
for q=1:1:6
    if q>1
        interpFun.Values=A(:,q); %Update values in interpolation structure
    end
    A_fit(:,q)=interpFun(P_I); %Interpolate
end

%Create nx9 tensor representation
A_fit_tensor_array=[A_fit(:,1) A_fit(:,6) A_fit(:,5) A_fit(:,6) A_fit(:,2) A_fit(:,4) A_fit(:,5) A_fit(:,4) A_fit(:,3)];

%Convert to mxn cell array containing the 3x3 tensors
[A_fit_tensor_cell]=tensorArray2tensorCell(A_fit_tensor_array,[size(A_fit,1) 1]);

%Eigen decomposition on each cell entry
[Av,Ad]=cellEig(A_fit_tensor_cell);
% [Av,Ad]=cellfun(@eig,A_fit_tensor_cell,'UniformOutput',0);

%Find index of maximum eigenvalue
[~,indMax]=cellfun(@(x) max(diag(x)),Ad,'UniformOutput',0);
indMax=cell2mat(indMax);

%Get matching eigenvectors
Av_max=zeros(size(Ad,1),3);
for q=1:1:size(Ad,1)
    Av_max(q,:)=Av{q}(:,indMax(q));
end

%% Check result

% OLD CODE 
% %Compute trace (should sum to 1)
% Ad_sum_cell=cellfun(@trace,Ad,'UniformOutput',0);
% Ad_sum=cell2mat(Ad_sum_cell);
% logicReplace=(Ad_sum<(1-1e-3));
% Av_max(logicReplace,:)=nan(nnz(logicReplace),3);

Vf_I=Av_max;
[Vf_I]=vecnormalize(Vf_I);

% Invalid vectors have zero magnitude after normalisation (tested as <0.5 here)
H=sqrt(sum(Vf_I.^2,2)); 
logicReplace=(H<0.5);
Vf_I(logicReplace,:)=nan(nnz(logicReplace),3); %Replace invalid vectors witn NaNs

%% Create output

switch nargout
    case 1 
        varargout{1}=Vf_I; 
    case 2
        varargout{1}=Vf_I;
        interpFun.Values=zeros(size(interpFun.Values)); %Set to zeros
        varargout{2}=interpFun; 
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
