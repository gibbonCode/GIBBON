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
v_end=v(end-2:end-1);
for q=1:1:numel(H)
    h=H(q);
%     try %new        
%         h.Visible='On';
%     catch %old        
%         set(h,'Visible','On');
%     end
switch v_end    
    case {'3a','3b'}
        set(h,'Visible','On');
%     case {'4b','5a'}
%         h.Visible='On';
    otherwise
        h.Visible='On';        
end
drawnow; 

end
