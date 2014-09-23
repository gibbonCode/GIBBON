function P=inv_parbound(Pb,Ps)

w=abs(Ps.ub-Ps.lb); %Interval width
b=Ps.f*w; %Point at which the output interval width equals t*w
a=tan(Ps.t*0.5*pi);
P=Ps.c+((b./a)*tan((pi*(Pb-Ps.lb)./w)-(pi/2)));




