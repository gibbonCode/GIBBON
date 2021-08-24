function h = GB_plotboneVRML_New(bone_model,RT,faceColor,lightsOn,parentName,faceAlpha)

bone_conn1 = bone_model.conn;
bone_pts   = bone_model.pts;

if ~exist('lightsOn','var')
    lightsOn = 1;
end

if ~exist('faceColor','var')
    faceColor =  [0.9882 0.9490 0.9020];%[0.2 0.5 0.2];
end
if ~exist('RT','var')
    RT = eye(4,4);
end
if ~exist('faceAlpha','var')
    faceAlpha = 0.4;
end

new_pts = (RT*[bone_pts ones(length(bone_pts),1)]')';

bone_pts = new_pts;

% % % Check Connections
% % if min(min(bone_conn1(1,1:3))) == 1
% %     bone_conn1(:,1:3) = bone_conn1(:,1:3)-1;
% % end


if ~exist('parentName','var')
    h = trisurf(bone_conn1(:,1:3)+1, bone_pts(:,1),bone_pts(:,2),bone_pts(:,3),...
        'EdgeColor','none','LineStyle','none','FaceLighting','phong','FaceColor',faceColor,'FaceAlpha',faceAlpha); hold on
    if lightsOn ~= 0
        GB_giveMeSomeLights(lightsOn);
    end
    
else
    hold(parentName,'on');
    h = trisurf(bone_conn1(:,1:3)+1, bone_pts(:,1),bone_pts(:,2),bone_pts(:,3),...
        'Parent', parentName,'EdgeColor','none','LineStyle','none','FaceLighting','phong','FaceColor',faceColor,'FaceAlpha',faceAlpha); 
    if lightsOn ~= 0
        GB_giveMeSomeLights(lightsOn,parentName);
    end
    
end




