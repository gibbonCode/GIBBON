function mfv(varargin)

% function mfv(H)
% %Makes figure or figures visible depending on the input handles. If none
% are provide all are searched and made visible. 


%% Parse input

switch nargin
    case 0 %No handles provided so assume all
        H=findall(0,'type','figure');
    case 1
        H=varargin{1};
    otherwise
        error('Wrong number of input arguments');
end

%% Make figure(s) visible
v=version; 
isOld=~isempty(strfind(v,'R2011')) || ~isempty(strfind(v,'R2012')) || ~isempty(strfind(v,'R2013')) || ~isempty(strfind(v,'R2014a')); 

for q=1:1:numel(H)
    h=H(q);
    
    if isOld
        set(h,'Visible','On');
    else %ASSUMED NEWER VERSION
        h.Visible='On';
    end

    drawnow;

end
