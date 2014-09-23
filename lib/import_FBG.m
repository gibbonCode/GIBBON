function D=import_FBG(loadName)
%-------------------------------------------------------------------------
%
%function D=import_FBG(loadName)
%
%A simple function to import FBG data from text file (specified by
%loadName) and output the data fields (columns) in a structure array D.
%
%
%Kevin Mattheus Moerman
%kevinmoerman@hotmail.com
%2014/01/07 - Created function
%-------------------------------------------------------------------------

%% PARSE TEXT FILE

fid=fopen(loadName);
[cell_out] = textscan(fid,'%s %f %f %f %f %f %f\n', 'delimiter', ',');
fclose(fid);

%% FORMULATE OUTPUT

D.test_date=cell_out{1};
D.cycle_no=cell_out{2};
D.FBG_temp=cell_out{3};
D.FBG_strain=cell_out{4};
D.ind_speed=cell_out{5};
D.ind_depth=cell_out{6};
D.TTL_logic=cell_out{7};

