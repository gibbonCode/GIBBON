function [I,J,K]=immesh(M)

[J,I, K]=meshgrid(0.5:1:(size(M,2)+0.5),0.5:1:(size(M,1)+0.5), 1:1:(size(M,3)));