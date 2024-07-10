clear all
close all
clc

%%
all_wpts = load("all_wpts.mat").all_wpts;
[all_wpts_out, lat, lon, lat_wpts, lon_wpts] = process_lateral_profile(all_wpts);

figure(1)
hold on
plot(lon, lat)
plot(cell2mat(all_wpts_out.lon), cell2mat(all_wpts_out.lat), '*')