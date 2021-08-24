function [smoothed_16x1,error_all_rmse,not_smoothed] = GB_smooth_aut_kin(input_16x1,plotOn,my_span,plotGimbalFix,myEulSeq)

if ~exist('plotOn','var')
    plotOn = 0; 
end
if ~exist('my_span','var')
    my_span = 5; 
end
if ~exist('plotGimbalFix','var')
    plotGimbalFix = 0; 
end
if ~exist('myEulSeq','var')
    myEulSeq = 'XYZ'; 
end
not_smoothed = 0;

nFrames = length(input_16x1);

%% Pad Data
pad_size = 20;
input_16x1_pad1 = [input_16x1(pad_size:-1:1,:);input_16x1];
input_16x1_pad2 = [input_16x1_pad1;input_16x1_pad1(end:-1:end-pad_size+1,:)];

nFrames_padded = nFrames + 2*pad_size;
%% Transform Data to 6x1 so we can smooth
for iFrame=1:nFrames_padded
    ref_RT_XCS{iFrame,1} = reshape(input_16x1_pad2(iFrame,:),[4 4])';
    
    RT_jt_format_1x6(iFrame,1:3) = GB_rad2deg(tform2eul(ref_RT_XCS{iFrame,1},myEulSeq));
    RT_jt_format_1x6(iFrame,4:6) = ref_RT_XCS{iFrame,1}(1:3,4);
end

%% Fix Gimbal Lock
[RT_jt_format_1x6(:,1:3)] = GB_fix_gimbal_lock(RT_jt_format_1x6(:,1:3),plotGimbalFix);

%% Smoothing
% Smoothed_JT_Vector = RT_jt_format_1x6;
Smoothed_JT_Vector = smooth(RT_jt_format_1x6,my_span,'moving');
Smoothed_JT = reshape(Smoothed_JT_Vector,nFrames_padded,[]);


%% Fix Gimbal Lock
[Smoothed_JT(:,1:3)] = GB_fix_gimbal_lock(Smoothed_JT(:,1:3),0);

%% Put Back in RT_4x4
counter = 1;
smoothed_16x1 = nan(nFrames,16);
for iFrame=(1+pad_size):(nFrames_padded-pad_size)
    new_RT_4x4{counter,1} = eul2tform(GB_deg2rad(Smoothed_JT(iFrame,1:3)),myEulSeq);
    new_RT_4x4{counter,1}(1:3,4) = Smoothed_JT(iFrame,4:6);
    
    smoothed_16x1(counter,:) = reshape(new_RT_4x4{counter,1}',[1 16]);

    counter = counter + 1;
end

%% Calculate RMSE of Smoothing
output_range = (1+pad_size):(nFrames_padded-pad_size);
% error_all = Smoothed_JT(output_range,:) - RT_jt_format_1x6(output_range,:);
[error_all_r2,error_all_rmse] = GB_rsquare(Smoothed_JT(output_range,:),RT_jt_format_1x6(output_range,:));
if error_all_r2 < 0.98 || error_all_rmse>1
    warning('Filtering did not work...');
    not_smoothed = 1;
    smoothed_16x1 = input_16x1;
end
disp(['RMSE of filtering is: ',num2str(error_all_rmse)]);

%% Plot Smoothing
if plotOn == 1
    fig_idx = randi(1000);
    figure(fig_idx);
    for i=1:6
        subplot(2,3,i)
        plot(RT_jt_format_1x6(output_range,i),'b');hold on
        plot(Smoothed_JT(output_range,i),'r');
        grid on
        xlim([1 nFrames])
        xlabel('Frame #')
        switch i
            case 1
                ylabel('X Tran');
                title('X Tran');
            case 2
                ylabel('Y Tran');
                title({['RMSE of filtering is: ',num2str(error_all_rmse)],'Y Tran'});
            case 3
                ylabel('Z Tran');
                title('Z Tran');
            case 4
                ylabel('Z Rot');
                title('Z Rot');
            case 5
                ylabel('X Rot');
                title('X Rot');
            case 6
                ylabel('Y Rot');
                title('Y Rot');
            otherwise
                disp('Check Input!');
        end
        legend('First','Smoothed')
    end
end




