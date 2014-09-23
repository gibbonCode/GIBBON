function [INP,EL]=read_abaqus_inp_3D(fid,n)

% clear all; close all; clc;

% files_to_load=cell(1);
% pathname = 'C:\Users\Kevin\Desktop\MATLAB_files';
% filename = 'Job_Mod_26_imp_v04997_01.inp';
% fid=fopen([pathname,'\', filename],'r');

% n=1;

%[DATA]=read_abaqus_inp_3D(fid,n)
% ------------------------------------------------------------------------
% This functions reads the initial nodal positions for the n-th part (e.g.
% n=2 for the second set of data) from an ABAQUS '.inp' file.
%
% The function assumes the .inp file has the following structure:
% *Heading
% ** Job name: job_Model_5_IMP_exp_match Model name: Model_5_implicit_exp_match
% *Preprint, echo=NO, model=NO, history=NO, contact=NO
% **
% ** PARTS
% **
% *Part, name=Piston
% *Node
%       1,          22.,          10.
%       2,          22.,          -7.
% *Element, type=RAX2
%  1,  1,  6
%  2,  6,  7
% *Node
%      86,  3.55271368e-15,          10.,  3.88450315e-16
% *Nset, nset=Piston-RefPt_, internal
% 86,
% *Elset, elset=_Surf-1_SPOS, internal
%   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16
%  17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32
%  48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63
%  64, 65, 66, 67, 68, 69, 70, 71, 72
% *Surface, type=ELEMENT, name=Surf-1
% _Surf-1_SPOS, SPOS
% *Elset, elset=Piston, generate
%   1,  85,   1
% *End Part
% **
% *Part, name=Sylgard_gel
% *Node
%       1,   22.3250542,   50.9665527
%       2,          32.,    62.739151
% *Element, type=CAX4H
%   1,  517,   59,   60,  591
%   2,  500,  639,  669,  717
%
% The function starts aquiring lines of numerical data from the .inp file
% after '*Part' has appeared at the start of a line for the n-th time. It then keeps
% aquiring data until '*Element' is encountered.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/10/2008
% ------------------------------------------------------------------------

%% FINDING START POINT

break_find=0;
found_start=0;
target='*PART';
while found_start==0;
    l=fgetl(fid);
    if (length(l)>length(target)) && (strcmp(l(1:5),target))
        break_find=break_find+1;
    end
    if break_find==n;
        found_start=1;
    end
end

%Skip one line
l=fgetl(fid);

%% READING INITIAL NODAL POSITIONS

end_point='*ELEMENTS';
found_end=0;
INP=[ ];
while found_end==0;
    l=fgetl(fid);
    if (length(l)>length(end_point)) && strcmp(l(1:8),end_point)==0;
        read_data=sscanf(l, '%f %*s %f %*s %f %*s %f', [1, inf]);
        if (~isempty(read_data))
            INP=[INP;read_data];
        end
    else
        found_end=1;
    end
end

INP=INP(:,2:end);

%% READING ELEMENTS

end_point='*Nset';
found_end=0;
EL=[ ];
while found_end==0;
    l=fgetl(fid);
    if (length(l)>length(end_point)) && strcmp(l(1:5),end_point)==0;
        %         5200,  6215,  6216,  6237,  6236,  5648,  5649,  5670,  5669, 23789, 23852, 23851, 23849, 22136, 22199, 22198,22196, 23788, 23791, 23853, 23850

        read_data_line_part1=sscanf(l, '%f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s', [1, inf]);
        l=fgetl(fid);
        read_data_line_part2=sscanf(l, '%f %*s %f %*s %f %*s %f %*s %f %*s %f %*s', [1, inf]);
        read_data=[read_data_line_part1 read_data_line_part2];

        if (~isempty(read_data))
            EL=[EL;read_data];
        end

    else
        found_end=1;
    end
end

EL=EL(:,2:end);

% figure;fig=gcf; clf(fig); colordef (fig, 'white'); set(fig,'Color',[1 1 1]); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% xlabel('X-J');ylabel('Y-I');zlabel('Z-K'); hold on; grid on;
% hp= patch('Faces',EL,'Vertices',INP);
% set(hp,'FaceColor','b','EdgeColor','k','FaceAlpha',1);


%%
fclose(fid);

%% END