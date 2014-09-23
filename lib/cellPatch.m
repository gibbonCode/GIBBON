function HP=cellPatch(Vv,Fv_cell,Cv)

Vv(isinf(Vv))=nan;
L_plot=cellfun(@(X) ~any(isnan(Vv(X(:)))),Fv_cell);

plotRange=find(L_plot);
HP=nan(numel(plotRange),1);
qI=1;
for q=plotRange'
    fv=Fv_cell{q};    
    hp=patch('faces',fv,'vertices',Vv,'faceColor','r');
    HP(qI)=hp; %Store handle
    if nargin>2 %Color specified
        if size(Cv,2)==1 %Colormap driven
            set(hp,'faceColor','flat','CData',Cv(q,:));
        elseif size(Cv,2)==3 %RGB driven
            set(hp,'faceColor',Cv(q,:));
        end
    end
  qI=qI+1;  
end

end