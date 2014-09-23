function [x]=bias_nodes1d(x,f)

min_x=min(x,[],2)*ones(1,size(x,2)); 
x=x-min_x; 
max_x=max(x,[],2)*ones(1,size(x,2)); 
x=max_x.*((x.^f)./(max_x.^f)); 
x=x+min_x;

end