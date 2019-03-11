function structOut = importAnyBodyOutput(fileName)

% function structOut = importAnyBodyOutput(fileName)
% ------------------------------------------------------------------------
% This function imports an AnyBody output file
%
% DISCLAIMER: This is a beta version of the AnyBody output reading tool! At
% this point, it reads the output for a single segment. Please use at your
% own risk and verify results!
%
% Change log: 
% 2014/09/25 Pavel E. Galibarov (AnyBody), created initial code
% 2014/09/25 Kevin M. Moerman, fixed minor bugs and added to GIBBON library
% with help and documentation
% 2014/10/08 Kevin M. Moerman, renamed function and altered output
% structure lay-out
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

%% Import

fid = fopen(fileName);
cntr = 0;
tline = fgetl(fid);
F=[];
M=[];
while ischar(tline)
   matches = regexp(tline, '^t=');
   if (length(matches)) > 0
      time=tline(matches(1)+2:length(tline)-1);
      cntr=cntr+1;
      poscnt=0; %reset positon index at zero
      structOut(cntr).t=str2num(time);
   end
   
   matches = regexp(tline, '^Pos=');
   if (length(matches)) > 0
      coords=tline(matches(1)+4:length(tline)-1);
      poscnt=poscnt+1; %Increment position index
      structOut(cntr).Pos(poscnt,:)=str2num(coords)';     
      F = zeros(1,3);
      M = zeros(1,3);
   end

   m1 = regexp(tline, '^F\[*[0-9]\]');
   if (length(m1)) > 0
      m1 = regexp(tline, '=');
      F0=str2num(tline(m1(1)+1:length(tline)-1))';
      F=F+F0;
      structOut(cntr).F(poscnt,:)=F;      
   end

   m1 = regexp(tline, '^M\[*[0-9]\]');
   if (length(m1)) > 0
      m1 = regexp(tline, '=');
      M0=str2num(tline(m1(1)+1:length(tline)-1))';
      M=M+M0;
      structOut(cntr).M(poscnt,:)=M;
   end
   tline = fgetl(fid);
end
fclose(fid);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
