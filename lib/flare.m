function [cMap]=flare(varargin)
% function [cMap]=flare(n)
% ------------------------------------------------------------------------
% Inspired by the seaborn flare colormap
%
% ------------------------------------------------------------------------

%%

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[ 236 175 128; ... 
 236 174 127; ... 
 236 173 127; ... 
 236 172 126; ... 
 236 171 126; ... 
 236 170 125; ... 
 236 169 124; ... 
 236 168 124; ... 
 236 167 123; ... 
 236 166 123; ... 
 236 166 122; ... 
 236 165 121; ... 
 235 164 121; ... 
 235 163 120; ... 
 235 162 120; ... 
 235 161 119; ... 
 235 160 118; ... 
 235 159 118; ... 
 235 158 117; ... 
 235 157 117; ... 
 235 156 116; ... 
 235 155 116; ... 
 235 154 115; ... 
 234 153 114; ... 
 234 152 114; ... 
 234 151 113; ... 
 234 150 113; ... 
 234 149 112; ... 
 234 148 111; ... 
 234 147 111; ... 
 234 146 110; ... 
 234 145 110; ... 
 233 144 109; ... 
 233 143 109; ... 
 233 142 108; ... 
 233 141 107; ... 
 233 140 107; ... 
 233 139 106; ... 
 233 138 106; ... 
 233 137 105; ... 
 232 136 105; ... 
 232 135 104; ... 
 232 134 103; ... 
 232 133 103; ... 
 232 132 102; ... 
 232 131 102; ... 
 232 130 101; ... 
 232 129 101; ... 
 231 128 100; ... 
 231 127 100; ... 
 231 126  99; ... 
 231 125  99; ... 
 231 124  98; ... 
 231 123  98; ... 
 230 122  98; ... 
 230 121  97; ... 
 230 120  97; ... 
 230 119  96; ... 
 230 118  96; ... 
 229 118  96; ... 
 229 117  95; ... 
 229 116  95; ... 
 229 115  95; ... 
 229 114  94; ... 
 228 113  94; ... 
 228 112  94; ... 
 228 111  93; ... 
 228 110  93; ... 
 227 109  93; ... 
 227 108  93; ... 
 227 107  92; ... 
 227 106  92; ... 
 226 105  92; ... 
 226 104  92; ... 
 226 103  92; ... 
 225 102  92; ... 
 225 101  91; ... 
 225 100  91; ... 
 224  99  91; ... 
 224  98  91; ... 
 224  97  91; ... 
 223  96  91; ... 
 223  95  91; ... 
 223  94  91; ... 
 222  93  91; ... 
 222  92  91; ... 
 221  91  91; ... 
 221  90  91; ... 
 221  89  91; ... 
 220  89  92; ... 
 220  88  92; ... 
 219  87  92; ... 
 219  86  92; ... 
 218  85  92; ... 
 218  84  92; ... 
 217  83  93; ... 
 216  82  93; ... 
 215  81  93; ... 
 215  80  94; ... 
 214  79  94; ... 
 213  78  94; ... 
 212  77  95; ... 
 211  76  95; ... 
 211  75  96; ... 
 210  75  96; ... 
 209  74  96; ... 
 208  73  97; ... 
 207  72  97; ... 
 206  71  98; ... 
 205  70  98; ... 
 204  70  99; ... 
 203  69  99; ... 
 202  69  99; ... 
 202  68 100; ... 
 201  68 100; ... 
 200  67 101; ... 
 199  67 101; ... 
 198  66 101; ... 
 197  66 102; ... 
 196  65 102; ... 
 195  65 103; ... 
 194  65 103; ... 
 193  64 103; ... 
 192  64 104; ... 
 191  63 104; ... 
 190  63 104; ... 
 189  63 104; ... 
 188  63 105; ... 
 187  62 105; ... 
 186  62 105; ... 
 185  61 106; ... 
 184  61 106; ... 
 183  61 106; ... 
 182  61 106; ... 
 181  60 107; ... 
 180  60 107; ... 
 179  60 107; ... 
 178  60 107; ... 
 177  59 108; ... 
 176  59 108; ... 
 175  59 108; ... 
 174  59 108; ... 
 173  58 108; ... 
 172  58 109; ... 
 171  58 109; ... 
 170  58 109; ... 
 169  57 109; ... 
 168  57 109; ... 
 168  57 110; ... 
 167  57 110; ... 
 166  56 110; ... 
 165  56 110; ... 
 164  56 110; ... 
 163  56 110; ... 
 162  55 110; ... 
 161  55 111; ... 
 160  55 111; ... 
 159  55 111; ... 
 158  54 111; ... 
 157  54 111; ... 
 156  54 111; ... 
 155  54 111; ... 
 154  53 111; ... 
 153  53 112; ... 
 152  53 112; ... 
 151  52 112; ... 
 150  52 112; ... 
 149  52 112; ... 
 148  52 112; ... 
 147  51 112; ... 
 146  51 112; ... 
 145  51 112; ... 
 144  51 112; ... 
 143  50 112; ... 
 142  50 112; ... 
 141  50 112; ... 
 140  49 112; ... 
 139  49 112; ... 
 138  49 112; ... 
 137  49 112; ... 
 136  48 112; ... 
 135  48 112; ... 
 134  48 112; ... 
 133  48 112; ... 
 132  48 112; ... 
 131  47 112; ... 
 130  47 112; ... 
 129  47 112; ... 
 128  47 112; ... 
 127  47 112; ... 
 126  46 112; ... 
 125  46 111; ... 
 124  46 111; ... 
 123  46 111; ... 
 122  45 111; ... 
 121  45 111; ... 
 120  45 111; ... 
 119  45 111; ... 
 118  45 110; ... 
 117  44 110; ... 
 116  44 110; ... 
 115  44 110; ... 
 114  44 110; ... 
 113  44 109; ... 
 112  44 109; ... 
 111  43 109; ... 
 110  43 109; ... 
 109  43 109; ... 
 108  43 108; ... 
 107  43 108; ... 
 106  42 108; ... 
 105  42 107; ... 
 104  42 107; ... 
 103  42 107; ... 
 102  42 107; ... 
 101  41 106; ... 
 100  41 106; ... 
  99  41 106; ... 
  98  41 105; ... 
  97  41 105; ... 
  96  40 105; ... 
  95  40 104; ... 
  94  40 104; ... 
  93  40 104; ... 
  92  39 103; ... 
  91  39 103; ... 
  90  39 103; ... 
  89  39 102; ... 
  88  38 102; ... 
  87  38 102; ... 
  86  38 101; ... 
  85  38 101; ... 
  84  38 101; ... 
  83  37 100; ... 
  82  37 100; ... 
  81  36 100; ... 
  80  36  99; ... 
  79  36  99; ... 
  78  36  99; ... 
  77  35  98; ... 
  76  35  98; ... 
  75  35  98; ... 
  74  34  98
    ];

cMap=cMap./255;

[cMap]=resampleColormap(cMap,n);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
