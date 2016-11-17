function [Vg] = evenlySampleCurve(varargin)

% function [Vg] = evenlySampleCurve(V,n,interpPar,closeLoopOpt)
% ------------------------------------------------------------------------
% 
% The evenlySampleCurve function samples a curve evenly in n points. The
% curve is parameterized using curve distance and can be closed loop if
% closeLoopOpt==1 (default=0). The resampling is performed using
% interpolation based on the method specified by interpPar. 
% Available methods are those associated with interp1 i.e.: 'linear',
% 'nearest', 'next', 'previous', 'spline', 'pchip' (default), 'cubic'.
% Alternatively interpPar my be set as a scalar in the range 0-1 to use the
% csaps method for cubic spline based smoothening.
% 
% See also: interp1, csaps
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 11/07/2013: Created
% 2013-2015: No log
% 05/05/2015: Improved csaps based smoothening
% 05/05/2015: Updated help and documentation
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        V=varargin{1};
        n=size(V,1);
        interpPar='pchip';
        closeLoopOpt=0;
    case 2
        V=varargin{1};
        n=varargin{2};
        interpPar='pchip';
        closeLoopOpt=0;
    case 3
        V=varargin{1};
        n=varargin{2};
        interpPar=varargin{3};
        closeLoopOpt=0;
    case 4
        V=varargin{1};
        n=varargin{2};
        interpPar=varargin{3};
        closeLoopOpt=varargin{4};
    otherwise
        error('Wrong number of input arguments');
end

if isempty(interpPar)
     interpPar='pchip';
end

if isempty(n)
     n=size(V,1);
end

if isempty(closeLoopOpt)
     closeLoopOpt=0;
end

%% 
n_in=n;

sizV=size(V);
if closeLoopOpt==1
    n=n+1;
    
    V=[V;V(1,:)]; %Close loop by adding start to end
    
    if ~ischar(interpPar) %CSAPS SMOOTHEN, interpret interpPar as smoothening parameter
                
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
% for q=1:size(V,2);    
%     if ~ischar(interpPar); %CSAPS SMOOTHEN, interpret interpPar as smoothening parameter
%         Vg(:,q)=csaps(D,V(:,q),interpPar,Dg);        
%     else %NORMAL METHODS
%         Vg(:,q)=interp1(D,V(:,q),Dg,interpPar);       
%     end
% end
if ~ischar(interpPar); %CSAPS SMOOTHEN, interpret interpPar as smoothening parameter
    Vg = csaps(D,V',interpPar,Dg)'; %Smoothened ppform
else %NORMAL METHODS
    for q=1:size(V,2);        
        Vg(:,q)=interp1(D,V(:,q),Dg,interpPar);
    end
end

if closeLoopOpt==1
    Vg=Vg(1:end-1,:);
end

if ~ischar(interpPar); %CSAPS
    %resample since smoothening alters spacing of points
    [Vg] = evenlySampleCurve(Vg,n_in,'pchip',closeLoopOpt);
end


