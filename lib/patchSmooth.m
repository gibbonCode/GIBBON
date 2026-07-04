function [P]=patchSmooth(F,V,IND_V,optionStruct)

% function [P]=patchSmooth(F,V,IND_V,optionStruct)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% 2014/06/02
% 2023/04/26 Switch to use structComplete for input handling, which is
% clearer and more concise
%------------------------------------------------------------------------

%%

%Set default structure
defaultOptionStruct.n=1; %Number of smoothing iterations
defaultOptionStruct.Method='LAP'; %Smoothing method
defaultOptionStruct.RigidConstraints=[]; %Indicices for nodes to hold on to

%Complement input with default if missing
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

if isempty(IND_V)
    [~,IND_V]=patchIND(F,V,2);
%     if isa(F,'cell')        
%         IND_V=[];
%         for q=1:1:numel(F)
%             C=patchConnectivity(F{q},V,'vv');
%             IND_V=[IND_V C.vertex.vertex];            
%         end
%     else
%         C=patchConnectivity(F,V,'vv');
%         IND_V=C.vertex.vertex;
%     end
end

%%

smoothMethod=optionStruct.Method;

%Smooth
switch smoothMethod
    case 'LAP' %Laplacian
        [P]=tesSmooth_LAP(F,V,IND_V,optionStruct);
    case 'HC' %Humphreys Classes
        [P]=tesSmooth_HC(F,V,IND_V,optionStruct);
    case 'tLAP' %Tangent Laplacian       
        
        %Invert face orientation if required
        [logicFlip]=isGlobalSurfDirOutward(F,V);
        if ~logicFlip
            F=fliplr(F); %Flip faces
        end
        
        %Set control parameters for Laplacian smoothening iterations
        optionStruct_t=optionStruct;
        optionStruct_t.n=1;
        optionStruct_t.Tolerance=[];
        P=V;
        for q=1:1:optionStruct.n
            [Ps]=tesSmooth_LAP(F,P,IND_V,optionStruct_t);  %The Laplacian smoothened coordinate set
            D=Ps-P; %smoothening intended displacement vectors
            [Dt]=patchVectorTangent(F,P,D,[]); %Tangential component of displacement
            P=P+Dt; %Displace mesh
        end
        
    case 'tHC' %Tangent HC
                
        %Invert face orientation if required
        [logicFlip]=isGlobalSurfDirOutward(F,V);
        if ~logicFlip
            F=fliplr(F); %Flip faces
        end
        
        %Set control parameters for Laplacian smoothening iterations
        optionStruct_t=optionStruct;
        optionStruct_t.n=1;
        optionStruct_t.Tolerance=[];
        P=V;
        for q=1:1:optionStruct.n
            [Ps]=tesSmooth_HC(F,P,IND_V,optionStruct_t);  %The smoothened coordinate set
            D=Ps-P; %smoothening intended displacement vectors
            [Dt]=patchVectorTangent(F,P,D,[]); %Tangential component of displacement
            P=P+Dt; %Displace mesh
        end
        
    otherwise
        error('Invalid smooth method specified');        
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
