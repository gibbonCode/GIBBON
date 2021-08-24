function eul_6x1_ZYX = GB_RT_to_eul_4x4(input_RT4x4)

if iscell(input_RT4x4)
    nData = size(input_RT4x4,1);
    eul_6x1_ZYX = nan(nData,6);
    
    for iD = 1:nData
        [R,eul_tra] = fX4_to_RT(input_RT4x4{iD});
        
        eul_rot = rad2deg(rotm2eul(R)); % Export in ZYX
        
        eul_6x1_ZYX(iD,:) = [eul_tra eul_rot(3) eul_rot(2) eul_rot(1)];
    end
else
    [R,eul_tra] = fX4_to_RT(input_RT4x4);
    
    eul_rot = rad2deg(rotm2eul(R)); % Export in ZYX
    
    eul_6x1_ZYX = [eul_tra eul_rot(3) eul_rot(2) eul_rot(1)]; % Fix the order: x-y-z translation, x-y-z rotation    
end
