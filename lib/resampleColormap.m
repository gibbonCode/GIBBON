function [cmap_i]=resampleColormap(cmap,n)

% function [cmap_i]=resampleColormap(cmap,n)
% ------------------------------------------------------------------------
%
% This function resamples the input colormap cmap using n steps. Resampling
% is based on linear interpolation. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

ind=(1:1:size(cmap,1))';
ind_i=linspace(1,size(cmap,1),n)';
cmap_i=zeros(n,size(cmap,2));

%Interpolate color data
for q=1:1:size(cmap,2)
    cmap_i(:,q)=interp1(ind,cmap(:,q),ind_i,'linear');
end
