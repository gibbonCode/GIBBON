function [C]=textureCoord2FaceColor(F_uv,ij_M,m,colorMappingOpt)

% function [C]=textureCoord2FaceColor(F_uv,ij_M,m,colorType)
% ------------------------------------------------------------------------
% 
% 2023/06/02 KMM: Created for OBJ texture mapping to patch model data
% ------------------------------------------------------------------------

%%
M=double(m)./double(intmax(class(m)));

C=F_uv;
switch colorMappingOpt
    case 'vertex-averaged' %Vertex averaged
        ij_M=round(ij_M);
        if isa(F_uv,'cell')
            for q=1:1:numel(F_uv)
                ft=F_uv{q};
                RGB_now=zeros(size(ft,1),3);
                for qc=1:1:size(ft,2)
                    ij_M_now=ij_M(ft(:,qc),:);
                    ind_M_now_R=sub2ind(size(M),ij_M_now(:,1),ij_M_now(:,2),1*ones(size(ft(:,qc))));
                    ind_M_now_G=sub2ind(size(M),ij_M_now(:,1),ij_M_now(:,2),2*ones(size(ft(:,qc))));
                    ind_M_now_B=sub2ind(size(M),ij_M_now(:,1),ij_M_now(:,2),3*ones(size(ft(:,qc))));
                    RGB_now=RGB_now+[M(ind_M_now_R) M(ind_M_now_G) M(ind_M_now_B)]/size(ft,2);
                end
                C{q}=RGB_now;
            end
        else
            ft=F_uv;
            RGB_now=zeros(size(ft,1),3);
            for qc=1:1:size(ft,2)
                ij_M_now=ij_M(ft(:,qc),:);
                ind_M_now_R=sub2ind(size(M),ij_M_now(:,1),ij_M_now(:,2),1*ones(size(ft(:,qc))));
                ind_M_now_G=sub2ind(size(M),ij_M_now(:,1),ij_M_now(:,2),2*ones(size(ft(:,qc))));
                ind_M_now_B=sub2ind(size(M),ij_M_now(:,1),ij_M_now(:,2),3*ones(size(ft(:,qc))));
                RGB_now=RGB_now+[M(ind_M_now_R) M(ind_M_now_G) M(ind_M_now_B)]/size(ft,2);
            end
            C=RGB_now;
        end
    case 'face-center' %Face center based color
        if isa(F_uv,'cell')
            for q=1:1:numel(F_uv)
                ft=F_uv{q};
                I=ij_M(:,1);
                J=ij_M(:,2);
                IJ_mean=round([mean(I(ft),2) mean(J(ft),2)]);
                ind_M_now_R=sub2ind(size(M),IJ_mean(:,1),IJ_mean(:,2),1*ones(size(IJ_mean,1),1));
                ind_M_now_G=sub2ind(size(M),IJ_mean(:,1),IJ_mean(:,2),2*ones(size(IJ_mean,1),1));
                ind_M_now_B=sub2ind(size(M),IJ_mean(:,1),IJ_mean(:,2),3*ones(size(IJ_mean,1),1));
                C{q}=[M(ind_M_now_R) M(ind_M_now_G) M(ind_M_now_B)];
            end
        else
            ft=F_uv;
            I=ij_M(:,1);
            J=ij_M(:,2);
            IJ_mean=round([mean(I(ft),2) mean(J(ft),2)]);
            ind_M_now_R=sub2ind(size(M),IJ_mean(:,1),IJ_mean(:,2),1*ones(size(IJ_mean,1),1));
            ind_M_now_G=sub2ind(size(M),IJ_mean(:,1),IJ_mean(:,2),2*ones(size(IJ_mean,1),1));
            ind_M_now_B=sub2ind(size(M),IJ_mean(:,1),IJ_mean(:,2),3*ones(size(IJ_mean,1),1));
            C=[M(ind_M_now_R) M(ind_M_now_G) M(ind_M_now_B)];
        end
    otherwise
        error('Unknown color mapping option');
end

end
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
