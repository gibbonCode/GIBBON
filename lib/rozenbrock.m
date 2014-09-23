function [Z]=rozenbrock(X,Y)

% -----------------------------------------------------------------------
% function [Z]=rozenbrock(X,Y)
%
% "In mathematical optimization, the Rosenbrock function is a non-convex
% function used as a test problem for optimization algorithms. It is also
% known as Rosenbrock's valley or Rosenbrock's banana function. This
% function is often used to test performance of optimization algorithms.
% The global minimum is inside a long, narrow, parabolic shaped flat
% valley. To find the valley is trivial, however to converge to the global
% minimum is difficult.It is defined by Z=(1-X.^2)+100.*((Y-(X.^2)).^2)." 
% 
% From: http://en.wikipedia.org/wiki/Rosenbrock_function
%
% Kevin Moerman
% kevinmoerman@hotmail.com
% -----------------------------------------------------------------------

Z=(1-X.^2)+100.*((Y-(X.^2)).^2);