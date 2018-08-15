function [F,V,logicValid,indFix2]=patchEdgeCollapse(F,V,E,logicKeep,meanOption)

% function [F,V,logicValid,indFix2]=patchEdgeCollapse(F,V,E,logicKeep,meanOption)
%-------------------------------------------------------------------------
%
%
%-------------------------------------------------------------------------

%% Parse input

if isempty(logicKeep)
   logicKeep=false(size(E));
   logicKeep(:,1)=1;
end
%% 
%

numEdges=size(E,1);
for q=1:1:numEdges
    edgeCollapseNow=E(q,:); %Current edge to collaps
    logicKeepNow=logicKeep(q,:);
    if all(logicKeepNow==0) || all(logicKeepNow==1)
        error('Invalid logic for keeping edge points provide, cannot remove or keep both points');
    end
    if nnz(ismember(edgeCollapseNow,F))==2 && (edgeCollapseNow(1)~=edgeCollapseNow(2))
        
        F(F==edgeCollapseNow(logicKeepNow==0))=edgeCollapseNow(logicKeepNow==1); %Replace index in face array
        E(E==edgeCollapseNow(logicKeepNow==0))=edgeCollapseNow(logicKeepNow==1); %Replace index in edge array
        if meanOption==1
            V(edgeCollapseNow(1),:)=(V(edgeCollapseNow(1),:)+V(edgeCollapseNow(2),:))/2; %Replace coordinate by mean of the two edge points
        end
        
    end
end

logicValid=~any(diff(sort(F,2),[],2)==0,2);
F=F(logicValid,:); %Keep only valid faces

%Remote unused points and update faces matrix
[F,V,indFix2]=patchCleanUnused(F,V);

%% Collect output

