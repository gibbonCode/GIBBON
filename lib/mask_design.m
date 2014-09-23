function [W]=mask_design(Z)

% function [Im Jm Km Wm]=mask_design(I, J, K, Z)
% ------------------------------------------------------------------------
% This function designs the weights for a mask. It assumes the vector 'Z'
% contains the normalised intensities found in a region of interest. The
% mask weight vector 'W' is calculated such that the value sum(W.*Z) is
% minimum when the mask is centered on the region of interest and maximum
% when it is not.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 13/08/2008
% ------------------------------------------------------------------------


%%

W=Z;
W(Z==1)=-(10./(numel(Z(Z==1))*(Z(Z==1))));
W((Z>0 & Z<1))=(10./(numel(Z(Z>0 & Z<1))*(Z((Z>0 & Z<1)))));
W(Z==0)=999;

%% END