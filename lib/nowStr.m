function S=nowStr(f)

switch nargin
    case 0
        S = datestr(now,'yyyy_mm_dd');
    case 1
        S = datestr(now,f);
end