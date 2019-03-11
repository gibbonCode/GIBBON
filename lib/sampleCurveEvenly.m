function [Vg] = sampleCurveEvenly(V,cPar)

% [Vg] = sampleCurveEvenly(V,cPar)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 24/10/2013
% 2017/04/28 Links to evenlySampleCurve instead, added warning to replace
% function in the future
%------------------------------------------------------------------------

%% Parse input

warning('This function will be replaced by evenlySampleCurve');

%Check number of input arguments
if nargin ~= 2
    error('myApp:argChk', 'Wrong number of input arguments')
end

%Check control parameter structure for entries
if isempty(cPar) %Assume resampling of open curve
    cPar.nd=size(V,1);
    cPar.interpMethod='pchip';
    cPar.typeOpt='num';
    cPar.closeLoopOpt=0;
else %structure exists but may not contain all fields    
    %Check resampling type
    if ~isfield(cPar,'typeOpt')
        cPar.typeOpt='num'; %Default is the same number of points
    end
    
    %Check nd which defines number of points or curve length increments
    if ~isfield(cPar,'nd')
        cPar.nd=size(V,1); %Default for cPar.typeOpt='num' is the same number of points
    end
    
    %Check interpolation method
    if ~isfield(cPar,'interpMethod')
        cPar.interpMethod='pchip'; %Default is pchip interpolation
    end
    
    %Check closed loop[ option
    if ~isfield(cPar,'closeLoopOpt')
        cPar.closeLoopOpt=0; %Default is pchip interpolation
    end
end

%%

%Redefine distance metric for evenly spaced points
switch cPar.typeOpt
    case 'dist'
        %Compute distance metric used for parametric representation
        D=pathLength(V);
        D_total=max(D(:));
        cPar.nd=round(D_total./cPar.nd); %Calculate n based on length increment
end

[Vg] = evenlySampleCurve(V,cPar.nd,cPar.interpMethod,cPar.closeLoopOpt);

%% 
% 
% sizV=size(V);
% if cPar.closeLoopOpt==1
%     
%     V=[V;V(1,:)]; %Close loop by adding start to end
%     
%     if ~ischar(cPar.interpMethod) %CSAPS SMOOTHEN, interpret cPar.interpMethod as smoothening parameter
%                 
%         midInd=round(size(V,1)/2);
%         %Close loop and aid periodicity by adding start and end points
%         V_addBefore=V(midInd:end-1,:);
%         V_addAfter=V(2:midInd,:);
%         V=[V_addBefore; V; V_addAfter]; 
%              
%         startInd=1+size(V_addBefore,1);
%         endInd=sizV(1)+size(V_addBefore,1)+1;
%     else
%         startInd=1;
%         endInd=size(V,1);
%     end
%     
% else
%     startInd=1;
%     endInd=size(V,1);
% end
% 
% %Compute distance metric used for parametric representation
% D=pathLength(V);
% 
% %Redefine distance metric for evenly spaced points
% switch cPar.typeOpt
%     case 'num'       
%         n=cPar.nd; %Use input as number of points
%         if cPar.closeLoopOpt==1
%             n=n+1;
%         end
%     case 'dist'
%         D_total=max(D(:));
%         n=round(D_total./cPar.nd); %Calculate n based on length increment
% end
% Dg=linspace(D(startInd),D(endInd),n)';
% 
% %Interpolate required x values for even spacing allong curve patch
% Vg=zeros(n,size(V,2));
% 
% %Interpolate using parametric representation
% for q=1:size(V,2)    
%     if ~ischar(cPar.interpMethod) %CSAPS SMOOTHEN, interpret cPar.interpMethod as smoothening parameter p
%         Vg(:,q)=csaps(D,V(:,q),cPar.interpMethod,Dg); %smoothen                
%     else %NORMAL METHODS
%         Vg(:,q)=interp1(D,V(:,q),Dg,cPar.interpMethod);       
%     end
% end
% 
% if cPar.closeLoopOpt==1
%     Vg=Vg(1:end-1,:);
% end
% 
% %Resampling again if smoothening was used since this influences point
% %spacing
% if ~ischar(cPar.interpMethod)
%     cPar.interpMethod='linear';
%     [Vg] = sampleCurveEvenly(Vg,cPar);
% end




 
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
