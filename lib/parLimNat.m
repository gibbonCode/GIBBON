function [varargout]=parLimNat(xx_c,xxlim,x)

% [xx,S]=parLimNat(xx_c,[xx_min xx_max],x);

%%

wn=abs(xx_c-xxlim(1));
wp=abs(xxlim(2)-xx_c);

%%

nc=6; 
xtn=linspace(xx_c-2*wn,xx_c,nc);
xtp=linspace(xx_c,xx_c+2*wp,nc);
xt=[xtn xtp(2:end)];
xtlim=[min(xt(:)) max(xt(:))];

xxt=xt;

Ln=xt<=xx_c;
Lp=xt>xx_c;
xxt(Ln)=erf(1/wn*2*(xt(Ln)-xx_c));
xxt(Ln)=xxt(Ln)*wn;
xxt(Lp)=erf(1/wp*2*(xt(Lp)-xx_c));
xxt(Lp)=xxt(Lp)*wp;
xxt=xxt+xx_c;

%% Construct and evaluate Piecewise Cubic Hermite Interpolating Polynomial 

pp = pchip(xt,xxt); 
xx=ppval_extrapVal(pp,x,xtlim,xxlim);

%% Create output
varargout{1}=xx; 
switch nargout
    case 2      
        %Create additional output structure
        S.pp=pp;
        S.interplim=xtlim;
        S.extrapVal=xxlim;        
       
        varargout{2}=S;
end

