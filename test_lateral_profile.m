clear all
close all
clc

%%
all_wpts = load("all_wpts.mat").all_wpts;
[all_wpts_out, lat, lon, lat_wpts, lon_wpts] = process_lateral_profile(all_wpts);

figure(1)
hold on
plot(lon, lat, '-')
% geoplot(lat, lon, '-v')
plot(cell2mat(all_wpts_out.lon), cell2mat(all_wpts_out.lat), '*')

for i = 1:length(all_wpts.lon)
    % text(all_wpts_out.start_fly_by{i}.lon, all_wpts_out.start_fly_by{i}.lat, num2str(i))
    % text(all_wpts_out.end_fly_by{i}.lon, all_wpts_out.end_fly_by{i}.lat, num2str(i))
    plot(all_wpts_out.start_fly_by{i}.lon, all_wpts_out.start_fly_by{i}.lat, '+g')
    plot(all_wpts_out.end_fly_by{i}.lon, all_wpts_out.end_fly_by{i}.lat, '+r')
end

for i = 2:length(all_wpts.lon)
    text(all_wpts.lon{i}, all_wpts.lat{i}, [num2str(i), ', ', all_wpts.leg_type{i}])
    plot(all_wpts_out.center_lon{i}, all_wpts_out.center_lat{i}, 'kh')
end

axis equal
xlim([8, 14])

figure(2)
geoplot(lat, lon, '-')
geobasemap topographic