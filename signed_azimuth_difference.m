
%
%   signed_azimuth_difference.m
%
%   Function:   Calculate the signed angular difference in azimuth 
%               between two geodesics at the point where they intersect
%   
%   Based on FAA Notice 8260.58A function "signedAzimuthDifference"
%
%   Flight Management and Procedure Design
%
%   Copyright (c) 2020 TU Berlin FF
%
function [da] = signed_azimuth_difference(a1, a2)


da = mod(a1-a2+pi,2*pi) - pi;

return


