function Cmapped=cmaperise(C,cmap,clim)

C=C(:);
p=(C-clim(1))./(clim(2)-clim(1));
p(p<0)=0;
p(p>1)=1;
IND_cmap=round((p*(size(cmap,1)-1))+1);
Cmapped=cmap(IND_cmap,:);

end