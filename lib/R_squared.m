function [R_sq]=R_squared(yi,fi)  

% function [R_squared]=R_squared()  
% ------------------------------------------------------------------------
% This function calculates R^2 using the data in vector yi and the modelled
% data in the vector fi. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 05/10/2008
% ------------------------------------------------------------------------

y_bar=nanmean(yi);
SS_err=nansum((yi-fi).^2);
SS_tot=nansum((yi-y_bar).^2);
R_sq=1-(SS_err./SS_tot);