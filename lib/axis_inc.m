function h=axis_inc(s)

%This function widens the axis limits by h_inc=w*(s-1), where w is the
%length of the diagonal of the axis window and s the scaling factor

%%

%Get axis limits
h=axis;

%Calculate increments
wx=h(2)-h(1); 
h_inc_x=wx*(s-1); %Axis limit increment

wy=h(4)-h(3); 
h_inc_y=wy*(s-1); %Axis limit increment

if numel(h)==6
    wz=h(6)-h(5);
    h_inc_z=wz*(s-1); %Axis limit increments
end

%Adjusting extrema
h(1)=h(1)-h_inc_x;
h(2)=h(2)+h_inc_x;
h(3)=h(3)-h_inc_y;
h(4)=h(4)+h_inc_y;
if numel(h)==6
    h(5)=h(5)-h_inc_z;
    h(6)=h(6)+h_inc_z;
end

%Setting axis limits
axis(h);

end