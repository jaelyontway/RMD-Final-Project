clc; clear; close all;

%% Input Variables
l1 = 20;
l2 = 42.1;
l3 = 11.8;
l4 = 45.6;
l_p = 7.8;
l_vec = [l1, l2, l3, l4];

w_2 = -2; %rad/s CW
alpha_2 = 0; %Constant w_2 

th_1 = 0;
spring_start_angle = 60; %deg
spring_end_angle = 110; %deg 

F_weight = 0.150*9.81 %
%% Begin Analysis
th_2_vec = linspace(180-spring_end_angle, 180-spring_start_angle);

delta = -1; 
delta_ap = 45;

for thIdx = 1:length(th_2_vec)
    % Position Analysis
    th_vec_out(:,thIdx) = fourbarpos(l_vec,th_1,th_2_vec(:,thIdx),delta);

    % Velocity Analysis
    [w_vec_out(:,thIdx),VA(:,thIdx),VAx(:,thIdx),VAy(:,thIdx)] = fourbarvel(l_vec,th_vec_out(:,thIdx),w_2);
    VPAx = -l_p * w_vec_out(3,thIdx) * sind(th_vec_out(3,thIdx) + delta_ap);
    VPAy = l_p * w_vec_out(3,thIdx) * cosd(th_vec_out(3,thIdx) + delta_ap);

    VPx = VPAx + VAx(thIdx);
    VPy = VPAy + VAy(thIdx);

    % Acceleration Analysis
    [a_vec_out(:,thIdx),AA(:,thIdx),Ax(:,thIdx),Ay(:,thIdx)] = fourbaraccel(l_vec,th_vec_out(:,thIdx),w_vec_out(:,thIdx),alpha_2);
    APAx = -w_vec_out(3,thIdx)^2 * l_p * (cosd(th_vec_out(3,thIdx) + delta_ap)) - a_vec_out(3,thIdx) * l_p * sind(th_vec_out(3,thIdx) + delta_ap);
    APAy = -w_vec_out(3,thIdx)^2 * l_p * (sind(th_vec_out(3,thIdx) + delta_ap)) + a_vec_out(3,thIdx) * l_p * cosd(th_vec_out(3,thIdx) + delta_ap);
    %APA = -w_vec_out(3,thIdx)^2 * l_ap *(cosd(th_vec_out(3,thIdx) + delta_ap) + 1i * sind(th_vec_out(3,thIdx) + delta_ap)) + a_vec_out(3,thIdx)* 1i * l_ap * (cosd(th_vec_out(3,thIdx) + delta_ap) + 1i * sind(th_vec_out(3,thIdx) + delta_ap));
    
    APx = APAx + Ax(thIdx);
    APy = APAy + Ay(thIdx);

end

th_3_vec = th_vec_out(3, :); 
th_4_vec = th_vec_out(4, :);

%% Motor Info
motor_stall_torque = 1.3; %kg/cm
gravity = 9.81;  
motor_stall_torque_converted = motor_stall_torque * gravity * 10; 

%% Mechanical Advantage Cal
r_input = 42.1; %Assumed 
r_output = l_vec(4);
MA = (l_vec(4)/l_vec(2))*(r_input/r_output) * (sind(th_4_vec - th_3_vec) ./ sind(th_2_vec - th_3_vec));

%% Convert torque_in to torque_out to Force out 
k = 555; % N*mm/rad, Assumed torsional spring constant 

spring_angle_vec = linspace(spring_start_angle, spring_end_angle, 100);  % 100 steps
spring_angle_change = deg2rad(spring_angle_vec - spring_start_angle); 
tau_input = k*spring_angle_change;

tau_output = MA .* tau_input .* r_output / r_input; 
F_out = tau_output / r_output;

% 2nd way to find F_out
% F_in = tau_input / r_input;
% F_out = F_in * MA;

%% Transfer F_out to F_normal to the ground and F_transverse
global_to_local = 240; %deg CCW 
th_3_global = th_3_vec + global_to_local;
th_3_calculate_normal_force = th_3_global - 180;

F_normal = F_out .* cosd(th_3_calculate_normal_force) + F_weight;
F_transverse = F_out .* sind(th_3_calculate_normal_force);

%% Velocity at the foot 
V_normal = VPy * cosd(th_3_calculate_normal_force) + VPx * cosd(90 - th_3_calculate_normal_force);
V_transverse = VPy * sind(th_3_calculate_normal_force) + VPx * sind(90 - th_3_calculate_normal_force);


%% Acceleration at the foot 
A_normal = APy * cosd(th_3_calculate_normal_force) + APx * cosd(90 - th_3_calculate_normal_force);
A_transverse = APy * sind(th_3_calculate_normal_force) + APx * sind(90 - th_3_calculate_normal_force);


%% Spring Torque (Tau_in) vs. Force output at the foot normal to the ground
figure;
subplot(3,1,1);
plot(tau_input, F_normal);
title('Normal Force (N) vs. Input Torque (N*mm)', 'LineWidth', 5);
xlabel('Input Torque (N*mm)');
ylabel('Normal Force (N)');
grid on;

%% Spring Torque (Tau_in) vs. Velocity at the foot 
subplot(3,1,2);
plot(tau_input, V_normal);
title('Velocity (mm/s) vs. Input Torque (N*mm)', 'LineWidth', 5);
xlabel('Input Torque (N*mm)');
ylabel('Velocity (mm/s)');
grid on;

%% Spring Torque (Tau_in) vs. Acceleration at the foot 
subplot(3,1,3);
plot(tau_input, A_normal);
title('Acceleration (mm/s^2) vs. Input Torque (N*mm)', 'LineWidth', 5);
xlabel('Input Torque (N*mm)');
ylabel('Acceleration (mm/s^2)');
grid on;

%% Plot of vertical velocity at the foot 
% figure;
% plot(th_2_vec, VP_mag);
% title('Velocity Magnitude of the Toe Off Point Relative to Grounded Femur', 'LineWidth', 5);
% xlabel('\theta_2 (degrees)');
% ylabel('velocity (mm/s)');
% hold on

%% Do not really need to show 
% figure;
% plot(th_2_vec, MA);
% title('Mechanical Advantage Relative to Local \theta_2', 'LineWidth', 5);
% xlabel('Local \theta_2 (degrees)');
% ylabel('Mechanical Advantages');
% hold on

%% Plot of normal force of the foot to the ground 
% figure;
% plot(th_2_vec, F_out);
% title('Output Force Relative to Local \theta_2', 'LineWidth', 5);
% xlabel('Local \theta_2 (degrees)');
% ylabel('Output Force (F*mm)');
% hold on