function [F_uni,IND_F,IND_F_2]=uniqueIntegerRow(F)

Fs=sort(F,2); %Sort so faces with same nodes have the same rows
[~,IND_F,IND_F_2]=unique(Fs,'rows');
F_uni=F(IND_F,:);


