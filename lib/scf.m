function scf(h)

%function scf(h)
%SCF Get handle to current figure.
%   SCF(h) Makes the figure specified by the handle h the current figure
%   (without making it appear). 
%
%   See also FIGURE, CLOSE, GCF, CLF, GCA, GCBO, GCO, GCBF.


set(0,'CurrentFigure',h);