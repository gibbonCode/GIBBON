function [data_XYZ] = GB_fix_gimbal_lock(data_XYZ,plotOn)

if ~exist('plotOn','var')
    plotOn = 1;
end
if size(data_XYZ,2) ~= 3
    error('Cannot use gimbal lock... Give it XYZ angles...');
end

nFrames = size(data_XYZ,1);

dRx = diff(data_XYZ(:,1));
dRy = diff(data_XYZ(:,2));
dRz = diff(data_XYZ(:,3));
if plotOn == 1
    figure(1531); clf; subplot(2,1,1);
    plot((data_XYZ(:,1)),'r','LineWidth',2);hold on
    plot((data_XYZ(:,2)),'g','LineWidth',2)
    plot((data_XYZ(:,3)),'b','LineWidth',2)
    title('Gimbal Lock')
    subplot(2,1,2); hold on;
    plot(dRx); 
    plot(dRy)
    plot(dRz)
    title('Derivative: To detect gimbal lock')

end
dRx_gmbl_fix = abs(dRx)>300; dRx_ind = find(dRx_gmbl_fix == 1);
dRy_gmbl_fix = abs(dRy)>300; dRy_ind = find(dRy_gmbl_fix == 1);
dRz_gmbl_fix = abs(dRz)>300; dRz_ind = find(dRz_gmbl_fix == 1);
fixed.x = 0;
fixed.y = 0;
fixed.z = 0;

if mod(length(dRx_ind),2) == 1
    dRx_ind = [dRx_ind;nFrames];
end
for iX = 1:2:length(dRx_ind)
    first_index = dRx_ind(iX)+1;
    if iX+1 <= length(dRx_ind)
        last_index  = dRx_ind(iX+1);
%         data_XYZ(first_index:last_index,1) = 360 + data_XYZ(first_index:last_index,1);

        %%%%%%%%%% Negative Check
%         disp('IF GIMBAL LOCK IS NOT WORKING, CHECK LINE 39 of CODE. I ADDED THIS FOR NEGATIVE CHECK, BUT DID NOT VALIDATE FOR POSITIVE');
        if (data_XYZ(first_index-1,1) < 0)
            data_XYZ(first_index:last_index,1) = data_XYZ(first_index:last_index,1)-360;

        else
            data_XYZ(first_index:last_index,1) = 360 + data_XYZ(first_index:last_index,1);
        end
        %%%%%%%%%%
    else
        data_XYZ(first_index:end,1) = 360 + data_XYZ(first_index:end,1);
    end
    fixed.x = 1;
end
if mod(length(dRy_ind),2) == 1
    dRy_ind = [dRy_ind;nFrames];
end

for iY = 1:2:length(dRy_ind)
    first_index = dRy_ind(iY)+1;
    if iY+1 <= length(dRy_ind)
        last_index  = dRy_ind(iY+1);
%         data_XYZ(first_index:last_index,2) = 360 + data_XYZ(first_index:last_index,2);
        %%%%%%%%%% Negative Check
%         disp('IF GIMBAL LOCK IS NOT WORKING, CHECK LINE 39 of CODE. I ADDED THIS FOR NEGATIVE CHECK, BUT DID NOT VALIDATE FOR POSITIVE');
        if (data_XYZ(first_index-1,2) < 0)
            data_XYZ(first_index:last_index,2) = data_XYZ(first_index:last_index,2)-360;
        else
            data_XYZ(first_index:last_index,2) = 360 + data_XYZ(first_index:last_index,2);
        end
        %%%%%%%%%%
    else
        data_XYZ(first_index:end,2) = 360 + data_XYZ(first_index:end,2);
    end
    fixed.y = 1;
end

if mod(length(dRz_ind),2) == 1
    dRz_ind = [dRz_ind;nFrames];
end

for iZ = 1:2:length(dRz_ind)
    first_index = dRz_ind(iZ)+1;
    if iZ+1 <= length(dRz_ind)
        last_index  = dRz_ind(iZ+1);
%         data_XYZ(first_index:last_index,3) = 360 + data_XYZ(first_index:last_index,3);
        %%%%%%%%%% Negative Check
%         disp('IF GIMBAL LOCK IS NOT WORKING, CHECK LINE 39 of CODE. I ADDED THIS FOR NEGATIVE CHECK, BUT DID NOT VALIDATE FOR POSITIVE');
        if (data_XYZ(first_index-1,3) < 0)
            data_XYZ(first_index:last_index,3) = data_XYZ(first_index:last_index,3)-360;
        else
            data_XYZ(first_index:last_index,3) = 360 + data_XYZ(first_index:last_index,3);
        end
        %%%%%%%%%%
    else
        data_XYZ(first_index:end,3) = 360 + data_XYZ(first_index:end,3);
    end
    fixed.z = 1;
end

if plotOn == 1
    figure(1531); subplot(2,1,1);
    if fixed.x == 1; plot((data_XYZ(:,1)),'b--','LineWidth',2); end
    if fixed.y == 1; plot((data_XYZ(:,2)),'r--','LineWidth',2); end
    if fixed.z == 1; plot((data_XYZ(:,3)),'g--','LineWidth',2); end
end
