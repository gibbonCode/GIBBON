function [Fs,Vs,Cs]=scalePatch(F,V,C,scaleFactor)

% --------------------------------------------------------------------
% function [Fs,Vs,Cs]=scalePatch(F,V,C,scaleFactor)
%
% CHANGE LOG: 
% 19/12/2013 Fixed error related to single face entry, see if statement
% related to size(F,1)
%
% Kevin Mattheus Moerman (kevinmoerman@hotmail.com)
% --------------------------------------------------------------------
%%

%If vertices are shared they need to be disconnected by adding vertices
% if numel(F)~=numel(V)
Vn=zeros(numel(F),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q); %The coordinate set
    XF=X(F); %The coordinates for each face vertex
    XFt=XF'; %transpose to prepare for column array
    Vn(:,q)=XFt(:); %Add column as part of new vertex matrix
end

%Fix indices in face matrix
indVn=1:numel(F);
% Ft=F';
% indFt=Ft(:);
Fs=reshape(indVn(:),size(F,2),size(F,1))';
Cs=C;

% else
%     Vn=V;
%     Fn=F;
% end
%Derive face means to shift vertices around mean
Vs=zeros(size(Vn));
for q=1:1:size(V,2)
    X=Vn(:,q);
    XF=X(Fs); 
    if size(F,1)==1
        XF=XF';
    end   
    meanXF=mean(XF,2)*ones(1,size(F,2));     
    meanXFt=meanXF'; 
    meanV(:,q)=meanXFt(:);   
end
Vn_shift=Vn-meanV; %Shift vertices
Vn_shift_scale=Vn_shift*scaleFactor; %Scale shifted vertices with respect to face center 
Vs=Vn_shift_scale+meanV; %Shift scaled vertices back
