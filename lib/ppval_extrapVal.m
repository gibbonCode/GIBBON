function v=ppval_extrapVal(pp,xx,interpLim,extrapVal)

%Evaluate piecewise polynomial
v=ppval(pp,xx);

%Set extrapolation values 
v(xx<interpLim(1))=extrapVal(1);
v(xx>interpLim(2))=extrapVal(2);
