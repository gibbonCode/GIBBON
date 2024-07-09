function [VC]=filletCurve(VP,rMax,np,closedLoopOpt)

% [V]=filletCurve(V,r,np,closedLoopOpt)
% ------------------------------------------------------------------------
% This function fillets a curve based on the input radius r using np points
% per fillet arc. If closedLoopOpt==1 then closed end conditions are used
% such that the end and start regions are also filleted.
%
% %% EXAMPLE:
% Vt=[0 0 0; 10 0 0; 5 10 0; 10 0 10; 0 10 10; ];
% r=2; %Fillet radius
% np=25; %Number of points used to construct each fillet edge
% closedLoopOption=0; %Use 1 if curve represents a closed loop but containes unique points
% [VN]=filletCurve(Vt,r,np,closedLoopOption);
%
% figure; hold on;
% plotV(Vt,'k.-.','lineWidth',2,'MarkerSize',25);
% plotV(VN,'r.-','lineWidth',3);
% axis equal; view(3); axis tight;  grid on;
% drawnow;
% %%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/03/19
% 2017/06/03 Fixed bug in relation to angles being 0 or 180 degrees
% 2023/05/08 Extended to allow for multiple radii to be specified
% 2023/05/08 Extended to allow for zero radius (skips fillet)
%------------------------------------------------------------------------

%% Parse input

eps_level=1e-6;

%Cope with 2D input
if size(VP,2)==2
    VP(:,3)=0;
    nDim=2;
else
    nDim=3;
end

if isscalar(rMax) % true for length 1, false for >1 and empty.    
    rOpt = 0;
elseif isempty(rMax)
    rOpt = -1;
elseif length(rMax) == size(VP,1)
    rMax = mcol(rMax);
    rOpt = 1;
else
    error('rMax should either be empty, match VP in length, or be a single number');
end

%%
numPoints=size(VP,1);
VC = [];
if numPoints>2

    if closedLoopOpt==1
        VP = [VP(end,:); VP; VP(1,:)]; % Add start and end
        if rOpt == 1
            rMax = [rMax(end); rMax; rMax(1,:)]; % Add start and end
        end
    end

    if closedLoopOpt~=1
        VC(end+1,:) = VP(1,:);
    end

    % lp = 0;
    i_last = size(VP,1)-1;
    L = sqrt(sum(diff(VP).^2,2));
    for i = 2:1:i_last
        if rOpt == 1
            r = rMax(i);
        else
            r = rMax;
        end
        if i == 2
            e1 = VP(i,:)-VP(i-1,:);
        else
            e1 = e2;
        end
        e2 = VP(i+1,:)-VP(i,:);

        n1 = vecnormalize(e1);
        n2 = vecnormalize(e2);
        n3 = vecnormalize(cross(n1,n2));
        a = pi - acos(dot(n1,n2));
        if abs(a-pi) > eps_level
            b = pi/2 - a/2;

            n1_e = vecnormalize(cross(n3,n1));
            n2_e = vecnormalize(cross(n3,n2));
            m1 = vecnormalize(n1_e + n2_e);
            l1 = L(i-1);
            l2 = L(i);

            l = min(l1,l2)/2;

            rFit = l/tan(b);          
            if rOpt ==-1
                r = rFit;
            end

            if rOpt~=-1 && r<=rFit
                rNow = r;
                lNow = rNow * tan(b);
                if lNow == l
                    fullRound = 1;
                else
                    fullRound = 0;
                end
            else
                fullRound = 1;
                rNow = rFit;
                lNow = l;
            end

            if rNow > 0
                d = abs(rNow/cos(b));
                Q = [-m1; vecnormalize(cross(n3,-m1)); n3]';
                m1p = -m1*Q;
                a = atan(m1p(2)/m1p(1));
                vc = VP(i,:)+ m1.*d;

                t = linspace(a-b,a+b,np)';
                
                Vc = vc + ([rNow.*cos(t) rNow.*sin(t) zeros(size(t))]*Q');

                if closedLoopOpt ~= 1 && i==2 && fullRound==1 && norm(Vc(1,:)-VC(1,:))<eps_level
                    Vc = Vc(2:end,:);
                end

                if i>2 && fullRound==1 && norm(Vc(1,:)-VC(end,:))<eps_level
                    Vc = Vc(2:end,:);
                end

                if closedLoopOpt==1 && i == i_last && fullRound==1 && norm(Vc(end,:)-VC(1,:))<eps_level
                    Vc = Vc(1:end-1,:);
                end
                VC=[VC;Vc];
            else
                VC(end+1,:) = VP(i,:);
            end
        else
            VC(end+1,:) = VP(i,:); % No rounding so add point
            lNow = 0.0;
        end
    end
    if closedLoopOpt~=1
        if norm(VC(end,:)-VP(end,:))<eps_level
            VC(end,:) = VP(end,:); % overwrite end
        else
            VC(end+1,:)=VP(end,:); % add end
        end
    end
else
    VC = VP;
end

%Remove 3rd dimension when 2D input was provided
if nDim==2
    VC=VC(:,1:2);
end

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
