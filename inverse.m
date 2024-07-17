%
%   inverse.m
%
%   Function: wrapper around 'vdist' that takes radians as an input and
%   outputs in radians
%   
%   Flight Management and Procedure Design
%
%   Copyright (c) 2020 TU Berlin FF
%
function [dist12, crs12, crs21] = inverse(lat1,lon1,lat2,lon2)

%   Convert from radians to degrees
lat1 = lat1*180.0/pi;
lon1 = lon1*180.0/pi;
lat2 = lat2*180.0/pi;
lon2 = lon2*180.0/pi;

%   Perform inverse algorithm
%[dist12, crs12, crs21] = vdist(lat1,lon1,lat2,lon2);
[dist12, crs12, crs21] = geoddistance(lat1,lon1,lat2,lon2);

%   Convert from degrees to radians
crs12 = crs12*pi/180.0;
crs21 = crs21*pi/180.0;

crs21 = crs21 + pi;
if crs21 > 2*pi 
    crs21 = crs21 - 2*pi;
end

