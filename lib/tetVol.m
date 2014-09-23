function VE=tetVol(E,V)

X=V(:,1); Y=V(:,2); Z=V(:,3);
XE=X(E); YE=Y(E); ZE=Z(E);

A=[XE(:,1) YE(:,1) ZE(:,1)];
B=[XE(:,2) YE(:,2) ZE(:,2)];
C=[XE(:,3) YE(:,3) ZE(:,3)];
D=[XE(:,4) YE(:,4) ZE(:,4)];

VE=abs(dot((A-D),cross((B-D),(C-D),2),2))./6;
