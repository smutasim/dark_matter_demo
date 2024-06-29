clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script calculates the gravitational force and velocity profile for 
% a galaxy with a given mass density distribution. This demonstration uses
% arbitrary values and a numerical approach to demonstrate how our current
% understanding of gravity and mass does not align with our observations.
% In reality, there is a theory that dark matter changes the mass distribution 
% of galaxies to allow for farther objects to travel faster. 
%
% The model generates a galaxy of 'stars' by placing point-masses along 
% many concentric circles of increasing diameter. The galaxy is flat, so
% there is no out-of-plane forces being modeled. Along the x-axis, where 
% y=0, the force on that particular point is towards the center of the
% galaxy. The magnitude will depend on the mass distribution. The force in
% the y-direction is zero, because the mass is symmetric about the x-axis.
% Then, using classical mechanics, the velocity required to keep a star in
% a circular orbit is calculated. This is what we should expect during
% observations.
%
% The user can tweak the parameters in the first code block below to
% simulate different mass density distribution equations or use
% non-arbitrary values to get real numbers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
G = 1; % Gravitational constant (arbitrary units)
radius_g = 5; % Radius of the galaxy
num_circles = 100; % Number of concentric circles
num_points_per_circle = 360;

% Mass density distribution function 
mass_distr = @(x) exp(-x.^2); % Example: Gaussian
%mass_distr = @(x) 0*x + 1; % Example: uniform

%% Generate Galaxy Model
x_list = linspace(0, radius_g, num_circles); % Distances from the center
x_list = x_list(2:end); % remove 0 

% Initialize arrays to store positions and masses
x_positions = 0;
y_positions = 0;
masses = mass_distr(0);

% Generate positions and masses for each circle
for j = 1:numel(x_list)
    radius = x_list(j);
    theta = linspace(0, 2 * pi, num_points_per_circle + 1); % Angle spacing
    theta = theta(1:end-1); % Remove the last element to avoid duplication
    
    % Generate positions for this circle
    circle_x = radius * cos(theta);
    circle_y = radius * sin(theta);
    
    % Calculate mass for each point based on the mass density function
    density = mass_distr(radius);
    total_mass = density * 2 * pi * radius; % Total mass for the circle
    point_masses = total_mass / num_points_per_circle; % Mass per point
    
    % Append to arrays
    x_positions = [x_positions, circle_x];
    y_positions = [y_positions, circle_y];
    masses = [masses, point_masses * ones(1, num_points_per_circle)];
end
clear radius_g circle_x circle_y total_mass point_masses theta radius
clear density num_circles num_points_per_circle
%% Calculate Gravitational Forces
Fx_list = zeros(numel(x_list), 1);
Fy_list = zeros(numel(x_list), 1);

% Compute gravitational forces
for j = 1:numel(x_list)
    x_target = x_list(j);
    y_target = 0;
    Fx = 0;
    Fy = 0;
    
    % Calculate forces from each point mass
    for i = 1:numel(x_positions)
        % Position of point i
        xi = x_positions(i);
        yi = y_positions(i);
        
        % Distance from point i to the target
        dx = xi - x_target;
        dy = yi - y_target;
        distance = sqrt(dx^2 + dy^2);
        
        % Skip if distance is zero (self-force)
        if distance == 0
            continue;
        end
        
        % Gravitational force components (F = G * m1 * m2 / r^2)
        m1 = masses(i);
        m2 = mass_distr(x_target);
        force_magnitude = G * m1 * m2 / distance^2;
        
        % Accumulate force components
        Fx = Fx + force_magnitude * (dx / distance);
        Fy = Fy + force_magnitude * (dy / distance);
    end
    
    % Store total force components for each radius
    Fx_list(j) = Fx;
    Fy_list(j) = Fy;
end
clear distance dx dy force_magnitude Fx Fy m1 m2 xi yi x_target y_target G
%% Calculate Circular Orbit Velocities
v = zeros(numel(x_list), 1);

for j = 1:numel(x_list)
    % Velocity required for circular orbit
    % PE + KE = 0
    % F + mv^2/r = 0
    % F = -mv^2/r
    % v = sqrt(-F*r/m)
    v(j) = sqrt(-Fx_list(j) * x_list(j) / mass_distr(x_list(j)));
end
clear i j Fx_list Fy_list masses x_positions y_positions
%% Plot Results
figure(1)
plot(x_list, mass_distr(x_list))
title('Galaxy Model with Mass Density Distribution')
xlabel('Distance from Galactic Center')
ylabel('Mass Density')
grid on

figure(2)
plot(x_list, v)
title('Expected Velocity at Galactic Location')
xlabel('Distance from Galactic Center')
ylabel('Velocity')
grid on
