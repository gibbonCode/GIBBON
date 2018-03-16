function [varargout]=patchNormalFix(varargin)

% function [F,L]=patchNormalFix(F,logicGroup,waitbarOn)
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        F=varargin{1};
        logicGroup=[];
        waitbarOn=0;
    case 2
        F=varargin{1};
        logicGroup=varargin{2};
        waitbarOn=0;
    case 3
        F=varargin{1};
        logicGroup=varargin{2};
        waitbarOn=varargin{3};
end

%%
if isempty(logicGroup)
    logicGroup=false(size(F,1),1);
    logicGroup(1)=1;
elseif nnz(logicGroup)==0    
    logicGroup(1)=1;
end

if ~islogical(logicGroup) 
    indGroup=logicGroup;
    logicGroup=false(size(F,1),1);
    logicGroup(indGroup(:))=1;
end

%%

[E_ind,E1,E2]=getEdgeSets(F);

if waitbarOn
    c=1;
    numSteps=size(F,1);
    hw=waitbar(c/numSteps,['Reorienting faces coherently...',num2str(round(100.*c/numSteps)),'%']);
end

L=false(size(F,1),1);
while 1

    logicFriends=any(ismember(E_ind,E_ind(logicGroup,:)),2) & ~logicGroup;    
    indCheck1=find(logicFriends,1);    
    indCheckEdge1=find(ismember(E_ind(indCheck1,:),E_ind(logicGroup,:)),1);
    edgeCheck1= [E1(indCheck1,indCheckEdge1) E2(indCheck1,indCheckEdge1)];
    
    if isempty(indCheck1)
        break
    end
    
    logicCheckFriends=any(ismember(E_ind,E_ind(indCheck1,indCheckEdge1)),2) & logicGroup;    
    indCheck2=find(logicCheckFriends,1);
    indCheckEdge2=find(E_ind(indCheck2,:)==E_ind(indCheck1,indCheckEdge1),1);
    edgeCheck2= [E1(indCheck2,indCheckEdge2) E2(indCheck2,indCheckEdge2)];

    if all(edgeCheck1==edgeCheck2)        
        F(indCheck1,:)=fliplr(F(indCheck1,:));                
        L(indCheck1)=1;
        [E_ind,E1,E2]=getEdgeSets(F);        
    end    
    logicGroup(indCheck1)=1;        
    
    if waitbarOn
        c=nnz(logicGroup);
        waitbar(c/numSteps,hw,['Reorienting faces coherently...',num2str(round(100.*c/numSteps)),'%']);
    end
end
if waitbarOn
    close(hw);
end

varargout{1}=F;
varargout{2}=L;

end

function [E_ind,E1,E2]=getEdgeSets(F)

E=patchEdges(F,0);
E_sort=sort(E,2);
E_ind=reshape(sub2ind([max(E_sort(:)) max(E_sort(:))],E_sort(:,1),E_sort(:,2)),[size(F,2),size(F,1)])';
E1=reshape(E(:,1),[size(F,2),size(F,1)])';
E2=reshape(E(:,2),[size(F,2),size(F,1)])';

end