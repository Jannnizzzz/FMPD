
%
%   direct.m
%
%   Function: wrapper around 'vreckon' that takes radians as an input and
%   outputs in radians
%   
%   Flight Management and Procedure Design
%
%   Copyright (c) 2020 TU Berlin FF
%
function [lat2, lon2, crs21] = direct(lat1,lon1,dist,crs)

%   Convert from radians to degrees
lat1 = lat1*180.0/pi;
lon1 = lon1*180.0/pi;
crs  = crs*180.0/pi;

%   Perform direct algorithm
%[lat2, lon2, crs21] = vreckon(lat1, lon1, dist, crs);
[lat2, lon2, crs21] = geodreckon(lat1, lon1, dist, crs);

%   Convert from degrees to radians
crs21 = crs21*pi/180.0;
lat2 = lat2*pi/180.0;
lon2 = lon2*pi/180.0;

crs21 = crs21 + pi;
if crs21 > 2*pi 
    crs21 = crs21 - 2*pi;
end
