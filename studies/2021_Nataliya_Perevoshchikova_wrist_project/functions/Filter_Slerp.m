function [qout] = Filter_Slerp(trajNoisy,hlpf)
%This function will use SLERP  to lowpass filter a noisy trajectory which
%is represented as quaternions.

% The new filter state should be moved toward the current input  by a step
% size that is proportional to the distance between the current input and
% the current filter state. To do this, both the dist and slerp functions
% are used. The dist function returns a measurement in radians of the
% difference in rotation applied by two quaternions. The range of the dist
% function is the half-open interval [0,pi).

%Set parameters - The interpolation parameter to slerp is in the
%closed-interval [0,1] so the output of dist must be re-normalized to this
%range. However, the full range of [0,1] for the interpolation parameter
%gives poor performance so it is limited to a smaller range hrange centered
%at hbias.

% All options have been kept the same as the matlab example...
% Wrtten by Jayishni Maharaj (UQ)
% 15 Feb 2018

hrange = 0.7;%0.4
hbias = 0.7;%0.4
%Limit low and high to the interval [0, 1].
low  = max(min(hbias - (hrange./2), 1), 0);
high = max(min(hbias + (hrange./2), 1), 0);
hrangeLimited = high - low;

%Initialize the filter and preallocate outputs.
y = trajNoisy(1); % initial filter state
qout = zeros(size(y), 'like', y); % preallocate filter output
qout(1) = y;

for ii=2:numel(trajNoisy)
    x = trajNoisy(ii); clear d
    d = dist(y, x);

    % Renormalize dist output to the range [low, high]
    %hlpf = (d./pi).*hrangeLimited + low;
    y = slerp(y,x,hlpf);
    qout(ii,:) = y;
end
end

