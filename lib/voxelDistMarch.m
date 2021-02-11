function [varargout]=voxelDistMarch(varargin)

%% Parse input

%Default option structure
defaultOptionStruct.toleranceLevel=0;
defaultOptionStruct.waitBarOn=false(1,1);
defaultOptionStruct.unitEdgeOn=false(1,1);
defaultOptionStruct.W=[];
defaultOptionStruct.Ws=[];

switch nargin
    case 1
        M=varargin{1};
        indStart=1;
        voxelSize=[1 1 1];
        optionStruct=defaultOptionStruct;
    case 2
        M=varargin{1};
        indStart=varargin{2};
        voxelSize=[1 1 1];
        optionStruct=defaultOptionStruct;
    case 3
        M=varargin{1};
        indStart=varargin{2};
        voxelSize=varargin{3};
        optionStruct=defaultOptionStruct;
    case 4
        M=varargin{1};
        indStart=varargin{2};
        voxelSize=varargin{3};
        optionStruct=varargin{4};
end

if numel(voxelSize)==1
    voxelSize=voxelSize*ones(1,3);
end
defaultOptionStruct.numSeeds=numel(indStart);

[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

if nargout==2
    computeSeedIndex=true(1,1);
else
    computeSeedIndex=false(1,1);
end

% Get optionional inputs
toleranceLevel=optionStruct.toleranceLevel;
numSeeds=optionStruct.numSeeds;
waitBarOn=optionStruct.waitBarOn;
unitEdgeOn=optionStruct.unitEdgeOn;

numVoxels=numel(M); %Number of vertices

W=optionStruct.W;
if isempty(W)
    W=ones(numVoxels,1);
end

Ws=optionStruct.Ws;
if isempty(Ws)
    Ws=ones(numel(indStart),1);
end

numStart=numel(indStart);
if numSeeds<numStart
    warning('Number of seeds is smaller than number of start points, assuming these should be additional points')
    numSeeds=numSeeds+numStart;
    optionStruct.numSeeds=numSeeds; %override input in case function is called recursively
end
numSteps=(numSeeds-numStart)+1;


%% Variables and computations outside of iterative loop

if computeSeedIndex
    indAll=(1:1:numVoxels)'; %Row indices for all points
end

% k=sphereIm(1);
% ind=find(k); 
% [ii,jj,kk]=ind2sub(size(k),ind);
% ii=ii([1 2 3 5 6 7]);
% jj=jj([1 2 3 5 6 7]);
% kk=kk([1 2 3 5 6 7]);
% ii=ii-2;
% jj=jj-2;
% kk=kk-2;

k=ones(3,3,3);
ind=find(k); 
[ii,jj,kk]=ind2sub(size(k),ind);
ii=ii([1:13 15:end]);
jj=jj([1:13 15:end]);
kk=kk([1:13 15:end]);
ii=ii-2;
jj=jj-2;
kk=kk-2;

dd=sqrt((ii.*voxelSize(1)).^2+(jj.*voxelSize(2)).^2+(kk.*voxelSize(3)).^2);

C=maskfind(M,1:numel(M),ii,jj,kk);
siz=size(C); %Size of C
L=C>0; %Logic for valid indices in C

%Calculate edge lengths
if unitEdgeOn
    DE=nan(siz);
    DE(L)=1; %Unit edge lenghts
else
    DE=repmat(dd,[1 numel(M)])';
    DE(~L)=0;
end

%Initialize distance vector
d=nan(numVoxels,1);
d(indStart)=0;
DD=nan(siz);
IND_d=C(L);

%Initialize seed index vector
if computeSeedIndex
    i=nan(numVoxels,1);
    i(indStart)=indStart;
    II=nan(siz);
end

d_previous=d;
i_previous=i;

%% Propagate and iteratively update distances
indSeed=nan(1,numSeeds);
indSeed(1:numStart)=indStart;

if numSteps>1      
    numLoopSteps=numSteps+1;
else
    numLoopSteps=1;
end

if waitBarOn
    hw=waitbar(1/numLoopSteps,['meshDistMarch...',num2str(round(100.*1/numLoopSteps)),'%']);
end

for q=1:1:numLoopSteps
    indStart=indSeed(~isnan(indSeed));
    while 1

        %Update distance data
        DD(L)=d(IND_d); %The distance data currently at neighbours
        if isempty(Ws)
            D_check=DD+DE; %Distance data plus edge lenghts
        else
            WS=ones(numVoxels,1);
            WS(indStart)=Ws;
            WS(~isnan(i))=WS(i(~isnan(i)));
            D_check=DD+DE./WS(:,ones(size(DE,2),1)); %Distance data plus edge lenghts
        end
        
        if computeSeedIndex
            [d,indMin]=min(D_check,[],2,'omitnan'); %Assign minimum distance
            
            %Update index data
            II(L)=i(IND_d); %The seed indices currently at neighbours
            ind=sub2ind(siz,indAll,indMin); %Linear indices for minimum used
            i=II(ind); %Get seed index at point chosen
            i(indStart)=indStart; %Override start seed indices
        else
            d=min(D_check,[],2,'omitnan'); %Assign minimum distance
        end
        d(indStart)=0; %Override starts to be zero

        logicFix=d_previous<d;
        d(logicFix)=d_previous(logicFix);
        if computeSeedIndex
            i(logicFix)=i_previous(logicFix);
        end
        
%         sum((d_previous-d).^2,'omitnan')
%         [~,indMax]=nanmax((d_previous-d).^2)
        
        %Check convergence
        if nnz(isnan(d))==nnz(isnan(d_previous)) %Check if number of NaN's changed
            if sum((d_previous-d).^2,'omitnan')<=toleranceLevel %Check if converged to within tolerance
                break %Stop while loop
            end
        end
        
        d_previous=d; %Store previous
        if computeSeedIndex
            i_previous=i;
        end
        
    end
    if q>1 && q<=numSteps
        [~,indSeed(numStart+(q-1))]=max(d.*W,[],'omitnan');
    end
    
    if waitBarOn
        waitbar(q/numLoopSteps,hw,['meshDistMarch...',num2str(round(100.*q/numLoopSteps)),'%']);
    end    
end

if waitBarOn
    close(hw);
end

%% Store output
varargout{1}=reshape(d,size(M));
if nargout==2
    varargout{2}=reshape(i,size(M));
end

end