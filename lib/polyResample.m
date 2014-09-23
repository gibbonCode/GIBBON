function [Vi]=polyResample(V,d_n,interpOpt,interpMethod)

%Parameterize curve using curve length
D = polyCurveLength(V); %Create distance space

if nnz(D>max(eps(V),[],2))<(size(D,1)-1)
    warning('Resampling not possible due to non-unique points (zero distances detected)');
    Vi=[];
else
    
    %Creating even space steps allong curve
    switch interpOpt
        case 1 %d_n is taken to be number of desired points
            D_geo=linspace(0,D(end),d_n);
        case 2 %d_n is taken to be desired point spacing
            d_n_fix=D(end)/round(D(end)/d_n);
            if d_n_fix~=d_n
                disp(['Point spacing was set to: ', num2str(d_n_fix)])
            end
            D_geo=0:d_n_fix:D(end);
    end
    
    %Interpolating curve points
    Vi=zeros(numel(D_geo),size(V,2)); %Initialising Vi
    for qd=1:1:size(V,2)
        Vi(:,qd) = interp1(D,V(:,qd),D_geo,interpMethod);
    end
end




