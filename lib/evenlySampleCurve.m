function [Vg] = evenlySampleCurve(V,n,interpMethod,closeLoopOpt)

% [Vg] = evenlySampleCurve(V,n,interpMethod)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 11/07/2013
%------------------------------------------------------------------------

%% 
n_in=n;

sizV=size(V);
if closeLoopOpt==1
    n=n+1;
    
    V=[V;V(1,:)]; %Close loop by adding start to end
    
    if ~ischar(interpMethod); %CSAPS SMOOTHEN, interpret interpMethod as smoothening parameter
                
        midInd=round(size(V,1)/2);
        %Close loop and aid periodicity by adding start and end points
        V_addBefore=V(midInd:end-1,:);
        V_addAfter=V(2:midInd,:);
        V=[V_addBefore; V; V_addAfter]; 
             
        startInd=1+size(V_addBefore,1);
        endInd=sizV(1)+size(V_addBefore,1)+1;
    else
        startInd=1;
        endInd=size(V,1);
    end
    
else
    startInd=1;
    endInd=size(V,1);
end

%Compute distance metric used for parametric representation
D=pathLength(V);

%Redefine distance metric for evenly spaced points
Dg=linspace(D(startInd),D(endInd),n)';

%Interpolate required x values for even spacing allong curve patch
Vg=zeros(n,size(V,2));

%Interpolate using parametric representation
for q=1:size(V,2);    
    if ~ischar(interpMethod); %CSAPS SMOOTHEN, interpret interpMethod as smoothening parameter
        Vg(:,q)=csaps(D,V(:,q),interpMethod,Dg);        
    else %NORMAL METHODS
        Vg(:,q)=interp1(D,V(:,q),Dg,interpMethod);       
    end
end

if closeLoopOpt==1
    Vg=Vg(1:end-1,:);
end

if ~ischar(interpMethod); %CSAPS
    %resample since smoothening alters spacing of points
    [Vg] = evenlySampleCurve(Vg,n_in,'pchip',closeLoopOpt);
end


