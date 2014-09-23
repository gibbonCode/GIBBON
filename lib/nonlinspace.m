function y=nonlinspace(range_in,range_out,fform,n)

x=linspace(range_in(1),range_in(2),n); 
eval(['y=',fform,';']);

y=y-y(1);
y=range_out(1)+(abs(range_out(1)-range_out(2))*(y./y(end)));

end