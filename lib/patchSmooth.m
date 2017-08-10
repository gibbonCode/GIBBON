function [P]=patchSmooth(F,V,IND_V,cPar)

% function [P]=patchSmooth(F,V,IND_V,cPar)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/06/02
%------------------------------------------------------------------------

%%

%Get/set method
if isfield(cPar,'Method')
    smoothMethod=cPar.Method;
else
    smoothMethod='LAP'; %DEFAULT
end

if ~isfield(cPar,'n')
    cPar.n=1; %DEFAULT
end

if isempty(IND_V)
    [~,IND_V]=patchIND(F,V,2);
end

%Smooth
switch smoothMethod
    case 'LAP' %Laplacian
        [P]=tesSmooth_LAP(F,V,IND_V,cPar);
    case 'HC' %Humphreys Classes
        [P]=tesSmooth_HC(F,V,IND_V,cPar);
    case 'tLAP' %Tangent Laplacian       
        
        %Invert face orientation if required
        [logicFlip]=isGlobalSurfDirOutward(F,V);
        if ~logicFlip
            F=fliplr(F); %Flip faces
        end
        
        %Set control parameters for Laplacian smoothening iterations
        cPar_t=cPar;
        cPar_t.n=1;
        cPar_t.Tolerance=[];
        P=V;
        for q=1:1:cPar.n
            [Ps]=tesSmooth_LAP(F,P,IND_V,cPar_t);  %The Laplacian smoothened coordinate set
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
        cPar_t=cPar;
        cPar_t.n=1;
        cPar_t.Tolerance=[];
        P=V;
        for q=1:1:cPar.n
            [Ps]=tesSmooth_HC(F,P,IND_V,cPar_t);  %The smoothened coordinate set
            D=Ps-P; %smoothening intended displacement vectors
            [Dt]=patchVectorTangent(F,P,D,[]); %Tangential component of displacement
            P=P+Dt; %Displace mesh
        end
        
    otherwise
        error('Invalid smooth method specified');        
end



 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
