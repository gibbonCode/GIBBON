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