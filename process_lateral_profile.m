

%
%   process_lateral_profile.m  - TEMPLATE 
%
%   Flight Management and Procedure Design
%
%   Input - all_wpts: structure with the following members
%
%           line            = line from the navigation database file
%           leg_type        = path termination type (e.g. IF, TF etc.)
%           name            = name of the waypoint;
%           lat             = latitude of the waypoint if available [deg];
%           lon             = longitude of the waypoint if available [deg];
%           alt_constraint  = altitude constraint 
%                             (1: at, 2: above, 3:below, 4: between)
%           alt_top         = top altitude of constraint (top if it is a
%                             between constraint)
%           alt_bottom      = bottom altitude of the altitude window (if between)
%           spd_constraint  = speed constraint (1 = at orn below)
%           spd_value       = value of constraint
%           crs             = course in case of an CF leg
%           center_lat      = latitude of center of radius in case of a RF leg
%           center_lon      = longitude of center of radius in case of a RF leg
%
%   Copyright (c) 2021 TU Berlin FF
%
function [all_wpts_out, lat, lon, lat_wpts, lon_wpts] = process_lateral_profile(all_wpts)
N = 100;
phi = deg2rad(18); % 18 deg bank
FT2M = 0.3048;
KTS2MPS	= (1.852/3.6);

% Compute air density, air pressure, temperature and speed of sound
[rho, p, ~, ~] = Atmos(5000*FT2M);
% Compute the true airspeed based on height and CAS schedule
V = CAStoTAS(220,rho,p)*KTS2MPS;
% Nominal bank angle of 22 degrees
R = V^2/(9.80665*tan(phi));

all_wpts_out = all_wpts;

lat_wpts = cell2mat(all_wpts.lat);
lon_wpts = cell2mat(all_wpts.lon);

%   Compute the distances between the waypoints, 
%   except for ARC --> then NaN
for ind = 1:length(all_wpts.lat)-1
    
    pt1.lat = lat_wpts(ind)*pi/180;
    pt1.lon = lon_wpts(ind)*pi/180;
    
    pt2.lat = lat_wpts(ind+1)*pi/180;
    pt2.lon = lon_wpts(ind+1)*pi/180;
    
    [all_wpts_out.dist_to_nxt{ind}, all_wpts_out.crs_to_nxt{ind}, crs21] = inverse(pt1.lat,pt1.lon,pt2.lat,pt2.lon);
    all_wpts_out.crs_from_prev{ind+1} = crs21 + pi;
    if ind == 1
        all_wpts_out.crs_from_prev{1} = all_wpts_out.crs_to_nxt{ind};
    end
end

% pre-calculate start and end positions as well as center point of fly-by
% maneuvre
for ind = 1:length(all_wpts.lat)
    if all_wpts.leg_type{ind} ~= "RF" & ind > 1 & ind < length(all_wpts.lat) & abs(signed_azimuth_difference(all_wpts_out.crs_to_nxt{ind}, all_wpts_out.crs_from_prev{ind})) > 1e-3 & all_wpts.leg_type{ind+1} ~= "RF" 
        pt1.lat = lat_wpts(ind-1)*pi/180;
        pt1.lon = lon_wpts(ind-1)*pi/180;
        pt3.lat = lat_wpts(ind+1)*pi/180;
        pt3.lon = lon_wpts(ind+1)*pi/180;

        [centerPt, startPt, endPt, ~] = wgs84_tangent_fixed_radius_arc(pt1, all_wpts_out.crs_to_nxt{ind-1}, pt3, all_wpts_out.crs_from_prev{ind+1}, R);
        all_wpts_out.start_fly_by{ind}.lat = startPt.lat * 180/pi;
        all_wpts_out.start_fly_by{ind}.lon = startPt.lon * 180/pi;
        all_wpts_out.end_fly_by{ind}.lat = endPt.lat * 180/pi;
        all_wpts_out.end_fly_by{ind}.lon = endPt.lon * 180/pi;
        all_wpts_out.center_lat{ind} = centerPt.lat * 180/pi;
        all_wpts_out.center_lon{ind} = centerPt.lon * 180/pi;
    else
        all_wpts_out.start_fly_by{ind}.lat = all_wpts.lat{ind};
        all_wpts_out.start_fly_by{ind}.lon = all_wpts.lon{ind};
        all_wpts_out.end_fly_by{ind}.lat = all_wpts.lat{ind};
        all_wpts_out.end_fly_by{ind}.lon = all_wpts.lon{ind};
    end
end

ind = length(all_wpts.lat);
all_wpts_out.crs_to_nxt{ind} = NaN;
all_wpts_out.dist_to_nxt{ind} = NaN;

lat = [];
lon = [];

for ind = 2:length(all_wpts.lat)
    pt1.lat = lat_wpts(ind-1)*pi/180;
    pt1.lon = lon_wpts(ind-1)*pi/180;
    pt2.lat = lat_wpts(ind)*pi/180;
    pt2.lon = lon_wpts(ind)*pi/180;

    if all_wpts.leg_type{ind} == 'TF' | all_wpts.leg_type{ind} == 'CF' | all_wpts.leg_type{ind} == 'IF' | all_wpts.leg_type{ind} == 'DF'
        if isnan(pt1.lat)
            continue
        end

        [dist_start, ~, ~] = inverse(pt1.lat, pt1.lon, all_wpts_out.end_fly_by{ind-1}.lat * pi/180, all_wpts_out.end_fly_by{ind-1}.lon * pi/180);
        [dist_end, ~, ~] = inverse(pt1.lat, pt1.lon, all_wpts_out.start_fly_by{ind}.lat * pi/180, all_wpts_out.start_fly_by{ind}.lon * pi/180);
        for ii = 0:N-1
            [lat(end+1), lon(end+1), ~] = direct(pt1.lat, pt1.lon, dist_start + (dist_end - dist_start) * ii/(N-1), all_wpts_out.crs_to_nxt{ind-1});
        end

        if ~isnan(all_wpts_out.center_lat{ind})
            [~,crs_cs,~] = inverse(all_wpts_out.center_lat{ind}*pi/180, all_wpts_out.center_lon{ind}*pi/180, all_wpts_out.start_fly_by{ind}.lat*pi/180, all_wpts_out.start_fly_by{ind}.lon*pi/180);
            [~,crs_ce,~] = inverse(all_wpts_out.center_lat{ind}*pi/180, all_wpts_out.center_lon{ind}*pi/180, all_wpts_out.end_fly_by{ind}.lat*pi/180, all_wpts_out.end_fly_by{ind}.lon*pi/180);
            d_azimuth = crs_ce - crs_cs;
            if d_azimuth < -pi
                d_azimuth = d_azimuth + 2*pi;
            end
            if d_azimuth > pi
                d_azimuth = d_azimuth - 2*pi;
            end
            for ii = 1:N
                [lat(end+1), lon(end+1), ~] = direct(all_wpts_out.center_lat{ind}*pi/180, all_wpts_out.center_lon{ind}*pi/180, R, crs_cs + d_azimuth * ii/N);
            end
        end
        continue
    end
    if all_wpts.leg_type{ind} == 'RF'
        pt3.lat = all_wpts_out.center_lat{ind}*pi/180;
        pt3.lon = all_wpts_out.center_lon{ind}*pi/180;
        [radius, crs31, ~] = inverse(pt3.lat, pt3.lon, pt1.lat, pt1.lon);
        [~, crs32, ~] = inverse(pt3.lat, pt3.lon, pt2.lat, pt2.lon);
        for ii = 1:N
            [lat(end+1), lon(end+1), crs12] = direct(pt3.lat, pt3.lon, radius, crs31 + (crs32-crs31) * ii/N);
        end
        continue
    end
    disp(all_wpts.leg_type{ind})
end

lat = rad2deg(lat);
lon = rad2deg(lon);
