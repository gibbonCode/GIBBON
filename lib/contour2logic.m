function [varargout]=contour2logic(M,v,Vcs)

% function [varargout]=contour2logic(M,v,Vcs)
%-------------------------------------------------------------------------
% This function converts the set of contours contained in the input cell
% Vcs to a logic array. The logic array matches the size of the input image
% M with the voxel size v.
%
% Change log:
% 2018/05/23 Changed to have logicIn and logicOn
% 2018/05/23 Added varargout output for handling both logicIn and logicOn
% 2018/05/23 Added alternative method
%-------------------------------------------------------------------------

%%

methodOption=2;
switch methodOption
    case 1 %Using inpolygon command
        [logicIn,logicOn,M]=inContourPolygon(M,v,Vcs);
    case 2 %Using custom function
        [logicIn,logicOn,M]=inContour(M,v,Vcs);
end

varargout{1}=logicIn;
varargout{2}=logicOn;
varargout{3}=M;

end

%%

function [logicIn,logicOn,M]=inContour(M,v,Vcs)

%Get size of the image
if isvector(M) %Assume M is really size(M)
    siz=M;
    M=zeros(siz); %Overwrite as zeros
else %Compute size(M) from M
    siz=size(M); 
end

hw=waitbar(0,'Computing logic... '); %Create waitbar

try    
    %Initialize logic arrays
    logicIn=false(siz);
    logicOn=false(siz);
    
    %Check empty slices
    try
        logicEmpty=cellfun(@(x) isempty(x{1}),Vcs);
    catch
        logicEmpty=cellfun(@(x) isempty(x),Vcs);
    end
    
    sliceRange=find(~logicEmpty);
    c=1;
    numSteps=numel(sliceRange);
    for qSlice=sliceRange
        numSubContours=numel(Vcs{qSlice}); %Number of sub-contours for the current slice
        
        logicOn_Now=false(siz(1:2));
        
        for qSub=1:1:numSubContours
            Vc=Vcs{qSlice}{qSub}; %Current contour
            if ~isempty(Vc)
                n=size(Vc,1);
                while 1
                    logicVertices_previous=logicOn_Now;
                    Vcc = evenlySampleCurve(Vc,n,'linear',1);
                    [I,J,~]=cart2im(Vcc(:,1),Vcc(:,2),Vcc(:,3),v);
                    I=round(I);
                    
                    if any(I<1) || any(J<1) || any(I>size(logicOn_Now,1)) || any(J>size(logicOn_Now,2))
                        %warning('Contour contains points outside of image space');
                        I(I<1)=1; I(I>size(logicOn_Now,1))=size(logicOn_Now,1);
                        J(J<1)=1; J(J>size(logicOn_Now,2))=size(logicOn_Now,2);
                    end
                    
                    J=round(J);
                    IND=sub2indn(size(logicOn_Now),[I(:) J(:)]);
                    logicOn_Now(IND)=1;
                    if all(logicVertices_previous(:)==logicOn_Now(:)) %If no change
                        break %Terminate while loop
                    else
                        n=n*2; %Increase curve sampling
                    end
                end
            end
        end
        
        %%
        
        %Create boundary, interior and exterior image
        labeledImage = bwlabel(~logicOn_Now,4); %Get labels for ~L image which will segment interior, exterior and boundary
        uniqueLabels=unique(labeledImage(:));
        
        labelsBoundary=labeledImage(logicOn_Now); %The label numbers for the boundary
        
        indExteriorVoxel=1; %First is outside since image is at least a voxel too big on all sides
        labelExterior=labeledImage(indExteriorVoxel); %The label number for the exterior
        
        labelsInterior=uniqueLabels(~ismember(uniqueLabels,[labelsBoundary(:); labelExterior])); %Labels for the interior (possibly multiple)
        
        logicIn_Now=ismember(labeledImage,labelsInterior);
        logicOn_Now=logicOn_Now;
        
        m=zeros(siz(1:2)); %The exterior is set to 0
        m(logicOn_Now)=1; %The boundary is set to 1
        m(logicIn_Now)=2; %Interior is set to 2
        
        M(:,:,qSlice)=m;
        logicIn(:,:,qSlice)=logicIn_Now;
        logicOn(:,:,qSlice)=logicOn_Now;
        
        waitbar(c/numSteps,hw,['Computing logic. ',num2str(round((c/numSteps)*100)),'%']);
        c=c+1;
    end
    
catch ME
    close(hw);
    rethrow(ME);
end

close(hw);

end


%%
function [logicIn,logicOn,M]=inContourPolygon(M,v,Vcs)

%Get size of the image
if isvector(M) %Assume M is really size(M)
    siz=M;
else %Compute size(M) from M
    siz=size(M); 
end

l=true(siz(1:2));
[E,V,~]=im2patch(l,l,'h',v);

%% Loop over contours

hw=waitbar(0,'Computing logic... ');
logicIn=false(siz);
logicOn=false(siz);
logicEmpty=cellfun(@(x) isempty(x{1}),Vcs);

sliceRange=find(~logicEmpty);
c=1;
numSteps=numel(sliceRange);
for qSlice=sliceRange
    numSubContours=numel(Vcs{qSlice}); %Number of sub-contours for the current slice
    for qSub=1:1:numSubContours
        Vd=Vcs{qSlice}{qSub}; %Current contour
        if ~isempty(Vd) %If it is not empty a levelset is computed
            logicIn_Now = inpolygon(V(:,1),V(:,2),Vd(:,1),Vd(:,2)); %Get logic for voxels inside contour
            logicOn_Now=logicIn_Now;
            logicIn_Now=reshape(all(logicIn_Now(E),2),siz(1:2));
            logicOn_Now=reshape(any(logicOn_Now(E),2),siz(1:2)) & ~logicIn_Now;
            logicIn(:,:,qSlice)=logicIn(:,:,qSlice) | logicIn_Now; %Add to current slice
            logicOn(:,:,qSlice)=logicOn(:,:,qSlice) | logicOn_Now; %Add to current slice
        end
    end
    waitbar(c/numSteps,hw,['Computing logic. ',num2str(round((c/numSteps)*100)),'%']);
    c=c+1;
end
close(hw);

M=zeros(siz); %The exterior is set to 0
M(logicOn)=1; %The boundary is set to 1
M(logicIn)=2; %Interior is set to 2

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
