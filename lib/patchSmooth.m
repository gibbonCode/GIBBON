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



