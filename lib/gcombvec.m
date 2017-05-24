function [B]=gcombvec(varargin)

numElements=cellfun(@numel,varargin);
numVecs=numel(varargin);
siz=prod(numElements);
B=zeros(numVecs,siz);

for q=1:1:numVecs
     vecNow=varargin{q};
     if q==1
         nRep1=1;
     else
         nRep1=prod(numElements(1:q-1));
     end     
     if q==numVecs
         nRep2=1;
     else
         nRep2=prod(numElements(q+1:end));         
     end
     w=repmat(vecNow,[nRep1 1]);    
     w=repmat(w(:),[nRep2 1]);          
     B(q,:)=w;
end