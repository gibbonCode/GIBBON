function axLim=axisLim(V_DEF)

% function axLim=axisLim(V_DEF)
% ------------------------------------------------------------------------
% Computes tight axis limits for the input coordinates V. The coordinates
% may be a 3-dimensional array where the 3rd direction could reflect
% coordinates as a function of time for instance. 
% 
% ------------------------------------------------------------------------

try
    axLim=[min(V_DEF,[],[1 3]); max(V_DEF,[],[1 3])];
catch 
    axLim=[min(min(V_DEF,[],1),[],3); max(max(V_DEF,[],1),[],3)];
end
axLim=axLim(:)';
