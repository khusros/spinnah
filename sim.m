% Angle is measured from y-axis
close all;

%%
% Configuration
%
fc = 1e6; % Hz
fs = 10e6; % symbols/sec
symbols_per_frame = 100000; % symbols/frame
frame_interval = symbols_per_frame/fs; % sec/frame

reference_snr = 16; % dB
range_at_reference_snr = 1; % km
Ptx = 1; % W
Gtx = 1; % dB
Grx = 3; % dB

antenna_rotation_rate = 2*pi/12; % rad/sec
antenna_init_angle = 0 % rad
target_xy_coord = [0; 2]; % km
target_speed = 0; % km/h
sensor_init_xy_coord = [0; 0]; % km
sensor_speed = 60; % km/h
antenna_beamwidth = 'very wide'; % options: 'very wide', 'wide', 'less wide', 'narrow'
sensor_trajectory_segments = [[90; 2], [0; 2], [-135; sqrt(8)]]; %[[90; 1],[60; 0.5],[30; 0.5],[0; 0.5],[-50; 0.1],[-90; 2],[-100; 0.5]]; % hdg (deg), distance (km), num samples (filled later)

%%
% Check antenna beamwidth selection
%
antenna_beamwidth_options = {'very wide', 'wide', 'less wide', 'narrow'};
num_elem_options = [2 3 4 8];
valid_antenna_beamwidth_option = false;
for k = 1:length(antenna_beamwidth_options)
  if (strcmp(antenna_beamwidth, antenna_beamwidth_options(k)))
    valid_antenna_beamwidth_option = true;
    num_elem = num_elem_options(k);
  endif
endfor
if (valid_antenna_beamwidth_option)
  num_elem
else
  error('Invalid BW option');
endif

%%
% Calculate array factor (af) and beamwidth
%
lambda = 3e8/fc;
element_separation = lambda/2;
theta = linspace(-pi/2, pi/2, 1000);
N = length(theta);
theta = theta(1:N-1);
af = abs(sin(num_elem * (pi * element_separation / lambda) * sin(theta))./ ...
      (num_elem * sin((pi * element_separation / lambda) * sin(theta))));
antenna_beamwidth = theta(max(find(af > sqrt(0.5)))) - theta(min(find(af > sqrt(0.5))));
antenna_beamwidth_formula = 0.886 * lambda / (num_elem * element_separation);

rotation_per_frame = rad2deg(frame_interval * antenna_rotation_rate);

figure;
polar(theta, abs(af));
title(sprintf("HPBW = %.1f deg (using formula = %.1f deg)\nframe-interval = %.1f deg, ratio=%.1f%%", ...
        rad2deg(antenna_beamwidth), rad2deg(antenna_beamwidth_formula), ...
        rotation_per_frame, rotation_per_frame/rad2deg(antenna_beamwidth)*100));

[~, num_segments] = size(sensor_trajectory_segments);
sensor_trajectory_segments = [sensor_trajectory_segments; zeros(1, num_segments)];
total_sim_ticks = 0;
for k = 1:num_segments
  segment_duration_seconds = (sensor_trajectory_segments(2,k)/sensor_speed)*3600;
  samples = ceil(segment_duration_seconds/frame_interval);
  total_sim_ticks = total_sim_ticks + samples;
  sensor_trajectory_segments(3,k) = samples;
endfor
total_sim_ticks += 1;

sensor_xy_coord = zeros(2, total_sim_ticks);
sensor_xy_coord(:,1) = sensor_init_xy_coord;
p = 2;
for k = 1:num_segments
  for s = 1:sensor_trajectory_segments(3,k);
    sensor_xy_coord(:, p) = sensor_xy_coord(:, p - 1) + ...
          sensor_speed/3600*frame_interval*...
          [sin(deg2rad(sensor_trajectory_segments(1,k))); ...
              cos(deg2rad(sensor_trajectory_segments(1,k)))];
    p = p + 1;
  endfor
endfor

Pn = Ptx * 10^(Gtx/10) * 10^(Grx/10) * lambda^2 / ...
          ((2 * pi * range_at_reference_snr*1000)^2 * 10^(reference_snr/10)); %SNR(@R) = Prx(@R)/Pn

R = zeros(1,total_sim_ticks);
SNR_rx = zeros(1, total_sim_ticks);
SNR_rx_with_ant_rotation = zeros(size(SNR_rx));
for k = 1:total_sim_ticks
  R(k) = sqrt(sum((sensor_xy_coord(:,k) - target_xy_coord).^2));
  Prx = Ptx * 10^(Gtx/10) * 10^(Grx/10) * lambda^2 ./ (2* pi * R(k)*1000).^2;
  SNR_rx(k) = 10*log10(Prx/Pn);
endfor

figure;
plot(sensor_xy_coord(1,:),sensor_xy_coord(2,:),'o', target_xy_coord(1), target_xy_coord(2), 'r+')
title('Sensor xy trajectory and target location');
xlabel('x (km)');
ylabel('y (km)');
xmax = max([max(sensor_xy_coord(1,:)) max(target_xy_coord(1))]);
xmin = min([min(sensor_xy_coord(1,:)) min(target_xy_coord(1))]);
ymax = max([max(sensor_xy_coord(2,:)) max(target_xy_coord(2))]);
ymin = min([min(sensor_xy_coord(2,:)) min(target_xy_coord(2))]);
xlim([xmin-0.1*abs(xmin-xmax) xmax+0.1*abs(xmin-xmax)]);
ylim([ymin-0.1*abs(ymin-ymax) ymax+0.1*abs(ymin-ymax)]);

figure;
%subplot(211);
plot(frame_interval*(0:total_sim_ticks-1), R);
title('Range to target');
xlabel('sec');
ylabel('km');

%subplot(212);
%plot(frame_interval*(0:total_sim_ticks-1), SNR_rx);
%title('SNR');
%xlabel('sec');
%ylabel('dB');

%num_samples = floor(run_time/frame_interval);
%sensor_xy_coord = zeros(2, num_samples);
%sensor_xy_coord(:,1) = sensor_init_xy_coord;
true_angle = zeros(1, total_sim_ticks); % relative true north
true_angle(1) = atan2(-sensor_xy_coord(1,1)+target_xy_coord(1), ...
                -sensor_xy_coord(2,1)+target_xy_coord(2));
antenna_angle = zeros(1,total_sim_ticks); % relative true north
antenna_angle(1) = antenna_init_angle;


boresight_offset_angle = zeros(1, total_sim_ticks);
boresight_offset_angle(1) = antenna_angle(1) - true_angle(1);

antenna_af = zeros(1, total_sim_ticks);
if (boresight_offset_angle(1) > pi/2 && boresight_offset_angle(1) < 3*pi/2)
  antenna_af(1) = 0;
elseif (boresight_offset_angle(1) == pi/2)
  antenna_af(1) = 1;
else
  theta = boresight_offset_angle(1);
  antenna_af(1) = abs(sin(num_elem * (pi * element_separation / lambda) * sin(theta))./(num_elem * sin((pi * element_separation / lambda) * sin(theta))));
endif

rx_data = zeros(1, total_sim_ticks);
for k = 2:total_sim_ticks

  true_angle(k) = atan2(-sensor_xy_coord(1,k)+target_xy_coord(1), ...
                -sensor_xy_coord(2,k)+target_xy_coord(2));
  antenna_angle(k) = antenna_angle(k-1) + antenna_rotation_rate * frame_interval;
  if (antenna_angle(k) > pi)
    antenna_angle(k) = antenna_angle(k) - 2*pi;
  endif

  if (true_angle(k) > 0 && antenna_angle(k) > 0)
    boresight_offset_angle(k) = true_angle(k) - antenna_angle(k);
  elseif (true_angle(k) > 0 && antenna_angle(k) < 0)
    if (true_angle(k) - antenna_angle(k) <= pi)
      boresight_offset_angle(k) = true_angle(k) - antenna_angle(k);
    else
      boresight_offset_angle(k) = true_angle(k) - antenna_angle(k) - 2*pi;
    endif
  elseif (true_angle(k) < 0 && antenna_angle(k) < 0)
    boresight_offset_angle(k) = true_angle(k) - antenna_angle(k);
  elseif (true_angle(k) < 0 && antenna_angle(k) > 0)
    if (true_angle(k) - antenna_angle(k) > -pi)
      boresight_offset_angle(k) = true_angle(k) - antenna_angle(k);
    else
      boresight_offset_angle(k) = true_angle(k) - antenna_angle(k) + 2*pi;
    endif
  endif

  if (boresight_offset_angle(k) > pi/2 || boresight_offset_angle(k) < -pi/2)
    antenna_af(k) = 0;
  elseif (boresight_offset_angle(k) == 0)
    antenna_af(k) = 1;
  else
    theta = boresight_offset_angle(k);
    antenna_af(k) = abs(sin(num_elem * (pi * element_separation / lambda) * sin(theta))./(num_elem * sin((pi * element_separation / lambda) * sin(theta))));
  endif
  antenna_af(k) = antenna_af(k)^2 * 10^(Grx/10);
  if (mod(k,1) == 0)
    k
    Prx = Ptx * 10^(Gtx/10) * antenna_af(k) * lambda^2 ./ (2* pi * R(k)*1000).^2;
    SNR_rx_with_ant_rotation(k) = 10*log10(Prx/Pn);
    %rx_data(k) = mean(sqrt(Prx/Pn) + randn(1,symbols_per_frame));
    rx_data(k) = mean(sqrt(Prx/Pn) + randn(1));
  endif
end

figure;
plot(rx_data(10:10:end));

figure;
subplot(411);
%plot(frame_interval*(0:total_sim_ticks-1), rad2deg(true_angle), 'x-');
plot(rad2deg(true_angle), 'x-');
title('Angle to target');
xlabel('sec');
ylabel('deg');

subplot(412);
%plot(frame_interval*(0:total_sim_ticks-1), rad2deg(antenna_angle), 'x-');
plot(rad2deg(antenna_angle), 'x-');
title('Antenna angle');
xlabel('sec');
ylabel('deg');

subplot(413);
%plot(frame_interval*(0:total_sim_ticks-1), rad2deg(boresight_offset_angle), 'x-');
plot(rad2deg(boresight_offset_angle), 'x-');
title('Boresight offset angle');
xlabel('sec');
ylabel('deg');

subplot(414);
%plot(frame_interval*(0:total_sim_ticks-1), 10*log10(antenna_af), 'x-');
plot(10*log10(antenna_af), 'x-');
title('Antenna factor');
xlabel('sec');
ylabel('dB');



