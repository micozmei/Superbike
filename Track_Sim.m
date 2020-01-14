%% Script Information
% Washington Superbike
% Track Simulation 

%% INPUT
clear; close all; clc;

% Motorcycle Input Data
m = 170; % motorcycle mass (kg)
p = 1.3; % wheelbase (m)
h = 0.65; % Vertical cg location (m)
b = 0.65; % back wheel contact point to cg (m)
Cd = 0.55; % Drag coefficent (assume constant)
Cl = -0.19; % Lift coefficent (assume constant)
A = 0.32; % Drag reference area (m^2)
mu = 1.1; % coeifficent of friction betwen tires and road
P_motor_max = 30e3; % Constant Motor Power Input (kW)

% Environment Input Data
g = 9.81; % gravity
rho = 1.225; % air density (kg/m^3)

% Track Data (READ !!!!!)
% You need to have two .mat data files, one containing x cord and the other
% containing the y cord with the following names: 'x_cord_Aragon.mat' and
% 'y_cord_Aragon.mat'. This data is converted into data in the section
% below and plots a map of the track. Be sure to verify that the track map
% looks reasonable. Make sure the two data files are in the same working
% directory (folder) as the script. 
%% TRACK MAP
Track_data_x_raw_1 = load('x_cord_Aragon.mat');
Track_data_y_raw_1 = load('y_cord_Aragon.mat');

Track_data_x_raw = struct2cell(Track_data_x_raw_1);
Track_data_y_raw = struct2cell(Track_data_y_raw_1);

Track_data_x = cell2mat(Track_data_x_raw);
Track_data_y = cell2mat(Track_data_y_raw);

Track_data = zeros(length(Track_data_x),2);
Track_sections = zeros(length(Track_data_x)-2,2);

for j = 1: length(Track_data_x)
Track_data(j,1) = Track_data_x(1,j); % Cordinates of track nodes (x,y) units should be in meters
Track_data(j,2) = Track_data_y(1,j);
end

LIN_1 = linspace(1,length(Track_data_x)-2,length(Track_data_x)-2);
LIN_2 = linspace(2,length(Track_data_x)-1,length(Track_data_x)-2);

Track_sections = [LIN_1; LIN_2]';

% PLOT of TRACK MAP
    figure(11)
    plot(Track_data_x,Track_data_y)
    xlabel('meters')
    ylabel('meters')
    title('Motorland Aragon')
    axis equal
    
%% TRACK SIM
num_nodes = length(Track_data);
num_sections = length(Track_sections);

for j = 1:num_sections
    % Establish Track Node ID's
    n1 = Track_sections(j,1); 
    n2 = Track_sections(j,2);
    
    % Track Node Coordinates
    x1 = Track_data(n1,1);
    y1 = Track_data(n1,2);
    x2 = Track_data(n2,1);
    y2 = Track_data(n2,2);
    
    % Track Section Lengths
    L(j) = sqrt((x2-x1)^2 + (y2-y1)^2);
    
    % Track Section Angles
    section_angle(j) = atand((y2-y1)/(x2-x1));  
end

Rc_a = zeros(1,num_nodes); % Initalize Corner Radius Vector

for j = 1:num_sections - 2
 % Calculating Corner Radius
    x1 = Track_data(j,1);
    y1 = Track_data(j,2);
    x2 = Track_data(j+1,1);
    y2 = Track_data(j+1,2);
    x3 = Track_data(j+2,1);
    y3 = Track_data(j+2,2);
    
    a = sqrt((x3-x1)^2 + (y3-y1)^2);
    b = sqrt((x3-x2)^2 + (y3-y2)^2);
    c = sqrt((x2-x1)^2 + (y2-y1)^2);
    
    Ac = acos((b^2 + c^2 - a^2)/(2*c*b));
    Rc_a(j + 1) = abs(a/(2*sin(pi - Ac)));
end 


%Inital Conditions
v = 1;
a = 0;

for j = 1: num_sections 
  velocity_map(j) = v;
  Fd = 0.5.*rho*Cd*A*v^2;
  Fa = P_motor_max/v - Fd;
  a_max_straight = ((mu*g)*((p-b)/p)/(1-mu*(h/p))) - (Fd)/m;
  a_max_wheelie = g*(b/h) - Fd/m;
  if Rc_a(j) > 10000 % Condition A
    a = Fa/m;
    if a < a_max_straight && a < a_max_wheelie% Condition B
        v = (v^2 + 2*a*L(j))^0.5; 
    elseif a_max_straight < a_max_wheelie
        a = a_max_straight;  
        v = (v^2 + 2*a*L(j))^0.5;
    elseif a_max_wheelie < a_max_straight
        a = a_max_wheelie ;
        v = (v^2 + 2*a*L(j))^0.5;
    end    
  else
  V_max_turn = (((mu*m*g)^2)./((m^2./Rc_a(j).^2)-0.25.*rho.*A^2.*Cd^2)).^0.25;
    if v < V_max_turn && a < a_max_straight && a < a_max_wheelie
        a = Fa/m;
        v = (v^2 + 2*a*L(j))^0.5;
    elseif v > V_max_turn
        v = V_max_turn;
    end
  end
end 

figure(12)
plot(linspace(1,length(Track_sections),length(Track_sections)),velocity_map*2.23)

for j = 1:num_sections - 2
    if velocity_map(num_sections - j) < velocity_map(num_sections - 1 - j)
        ab_flip = ((p-b)/h)*g;
        Fb_max = 0.5*rho*A*velocity_map(num_sections - 1 - j) + m*g*mu;
        a_b_max = Fb_max/m;
        if a_b_max < ab_flip
            velocity_map(num_sections - 1 - j) = ((velocity_map(num_sections - j).^2 + 2*a_b_max.*L(num_sections - j)).^0.5);
        elseif ab_flip < a_b_max
            velocity_map(num_sections - 1 - j) = ((velocity_map(num_sections - j).^2 + 2*ab_flip.*L(num_sections - j)).^0.5);    
        end
    else
    velocity_map(num_sections - 1 - j) = velocity_map(num_sections - 1 - j);
    end
end

hold on
plot(linspace(1,length(Track_sections),length(Track_sections)),velocity_map*2.23)
xlabel('distance (m)')
ylabel('Velocity (mph)')
title('Track Velocity Profile')
grid on
legend('No Braking','Braking')

t = num_sections/mean(velocity_map);
t_plot = linspace(0,t,length(velocity_map));

figure(13)
plot(t_plot,velocity_map*2.23)
xlabel('time (s)')
ylabel('Velocity (mph)')
title('laptime')
grid on