function [L]=allchulltouch(DT3_1,DT3_2,DT3_3,V1,V2,V3)

L=0;
Lc123=~isnan(tsearchn(DT3_1.Points,DT3_1.ConnectivityList,[V2; V3]));
Lc1=any(Lc123(1:size(V2,1))) & any(Lc123(size(V2,1)+1:end));
if Lc1
    Lc213=~isnan(tsearchn(DT3_2.Points,DT3_2.ConnectivityList,[V1; V3]));
    Lc2=any(Lc213(1:size(V1,1))) & any(Lc213(size(V1,1)+1:end));
    if Lc2
        Lc312=~isnan(tsearchn(DT3_3.Points,DT3_3.ConnectivityList,[V1; V2]));
        Lc3=any(Lc312(1:size(V1,1))) & any(Lc312(size(V1,1)+1:end));
        if Lc3;
            L=1;
        end
    end
end

% if Lc12 %V2 is touching chull 1
%     Lc13=any(~isnan(tsearchn(DT3_1.Points,DT3_1.ConnectivityList,V3)));
%     if Lc13 %V3 is touching chull 1
%         Lc21=any(~isnan(tsearchn(DT3_2.Points,DT3_2.ConnectivityList,V1)));
%         if Lc21 %V1 is touching chull 2
%             Lc23=any(~isnan(tsearchn(DT3_2.Points,DT3_2.ConnectivityList,V3)));
%             if Lc23 %V3 is touching chull 2
%                 Lc31=any(~isnan(tsearchn(DT3_3.Points,DT3_3.ConnectivityList,V1)));
%                 if Lc31 %V1 is touching chull 3
%                     Lc32=any(~isnan(tsearchn(DT3_3.Points,DT3_3.ConnectivityList,V2)));
%                     if Lc32 %V2 is touching chull 3
%                         L=1;
%                     end
%                 end
%             end
%         end
%     end
% end


% L=0;
% Lc12=any(~isnan(tsearchn(DT3_1.Points,DT3_1.ConnectivityList,V2)));
% if Lc12 %V2 is touching chull 1
%     Lc13=any(~isnan(tsearchn(DT3_1.Points,DT3_1.ConnectivityList,V3)));
%     if Lc13 %V3 is touching chull 1
%         Lc21=any(~isnan(tsearchn(DT3_2.Points,DT3_2.ConnectivityList,V1)));
%         if Lc21 %V1 is touching chull 2
%             Lc23=any(~isnan(tsearchn(DT3_2.Points,DT3_2.ConnectivityList,V3)));
%             if Lc23 %V3 is touching chull 2
%                 Lc31=any(~isnan(tsearchn(DT3_3.Points,DT3_3.ConnectivityList,V1)));
%                 if Lc31 %V1 is touching chull 3
%                     Lc32=any(~isnan(tsearchn(DT3_3.Points,DT3_3.ConnectivityList,V2)));
%                     if Lc32 %V2 is touching chull 3
%                         L=1;
%                     end
%                 end
%             end
%         end
%     end
% end

end