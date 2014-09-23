function dispStartTitleGibbonCode(msgText)

lineSep=repmat('%',1,45);
disp(' '); %Empty line
disp(lineSep); %Stripe
disp(['--- ',msgText,' --- ',datestr(now)]); %Start message
