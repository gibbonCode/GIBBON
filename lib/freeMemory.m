function [numFreeBytes]=freeMemory

if ispc
    mem_struct=memory;
    numFreeBytes=mem_struct.MaxPossibleArrayBytes;
elseif isunix
    [~,w] = unix('free | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    numFreeBytes = (stats(3)+stats(end));
else
    error('Your platform does not seem to be supported. Code your own solution or contact support.');
end