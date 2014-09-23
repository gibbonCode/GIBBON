function [Ri]=interp_polar(TZ,radiusData,TZi,interp_method)

if size(TZ,2)==1 %interp1 function
    thetaData=TZ;
    thetaData=[(thetaData-2.*pi); thetaData; (thetaData+2.*pi)];
    
    thetaData=thetaData.*min(radiusData(:));
    TZi=TZi.*min(radiusData(:));
    
    radiusData=repmat(radiusData,[3,1]);
    if ~ischar(interp_method);
        Ri = csaps(thetaData,radiusData,interp_method,TZi);
    else
        Ri=interp1(thetaData,radiusData,TZi,interp_method);
    end
elseif size(TZ,2)==2 
    thetaData=[(TZ(:,1)-2.*pi); TZ(:,1); (TZ(:,1)+2.*pi)];
    radiusData=repmat(radiusData,[3,1]);
    ZData=repmat(TZ(:,2),[3,1]);
    
    %Scaling angular units to length units using min(R(:))
    thetaData=thetaData.*min(radiusData(:));
    TZi(:,1)=TZi(:,1).*min(radiusData(:));
    
    %Interpolating
    if strcmp(interp_method,'natural') || strcmp(interp_method,'linear') || strcmp(interp_method,'nearest') %TriScatterdInterp function
        F=scatteredInterpolant([thetaData ZData],radiusData,interp_method);
        Ri=F(TZi);
    else %Griddata function
        Ri = griddata(thetaData,ZData,radiusData,TZi(:,1),TZi(:,2),interp_method);
    end
end

end
