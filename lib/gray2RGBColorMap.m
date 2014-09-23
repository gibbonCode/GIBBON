function [C]=gray2RGBColorMap(G,cMap,cLim)

if isempty(cLim)
   cLim=[min(G(:)) max(G(:))]; 
end

G=G(:); %Column

%Snap to color limits
G(G<cLim(1))=cLim(1);
G(G>cLim(2))=cLim(2);

%Normalized to range [0-1]
G=G-min(G(:)); 
max_G=max(G(:));
if max_G>0
    G=G./max(G(:));
end
        
numColorLevels=size(cMap,1); %Number of color levels
linearIndex_G=1+round(G.*(numColorLevels-1));
C=cMap(linearIndex_G,:);

