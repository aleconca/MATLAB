% Parameters
numIterations = 1000;    % Number of iterations
initialX = 0.1;         % Initial value for x
initialY = 0.1;         % Initial value for y
initialZ = 0.1; 
perturbation = 0.0001;  % Small perturbation to initial conditions

% Initialize arrays to store the results
x = zeros(1, numIterations);
y = zeros(1, numIterations);
z = zeros(1, numIterations);

% Initialize the initial conditions
x(1) = initialX;
y(1) = initialY;
z(1) = initialY;

% Simulate the system:
%The Lorenz system, also known as the Lorenz attractor, is a set of three nonlinear differential equations 
% that describe a chaotic, deterministic dynamical system. It was introduced by the meteorologist Edward Lorenz in 1963 
% as a simplified mathematical model of atmospheric convection. The Lorenz system has since become a classic example 
% of chaotic behavior in a dynamical system and is widely studied in the field of chaos theory.

%The three differential equations that make up the Lorenz system are as follows:

%dx/dt = σ(y - x)
%dy/dt = x(ρ - z) - y
%dz/dt = xy - βz
%Here are the meanings of the variables and parameters:

%-x, y, and z are the state variables that represent the system's state in three-dimensional space.
%-t is time.
%-σ, ρ, and β are the system parameters, which control the behavior of the system. 

% These parameters determine the system's behavior and can be adjusted to create different chaotic or non-chaotic behaviors.

for t = 2:numIterations
    % Update the system equations (Lorenz attractor)
    dx = 10 * (y(t-1) - x(t-1));
    dy = x(t-1) * (28 - z(t-1)) - y(t-1);
    dz = x(t-1) * y(t-1) - 8/3 * z(t-1);
    
    % Apply a small perturbation to initial conditions
    if t == 2
        x(1) = initialX + perturbation;
    end
    
    % Update the variables
    x(t) = x(t-1) + dx * 0.01; 
    y(t) = y(t-1) + dy * 0.01;
    z(t) = z(t-1) + dz * 0.01;
end

% 2D plot
figure;
plot(x, y);
title('Lorenz Attractor - Butterfly Effect');
xlabel('X');
ylabel('Y');









%3D Plot:
%Having defined the Lorenz system equations and parameters as discussed above:

%1. Define the Lorenz system equations and parameters (σ, ρ, β, initial conditions).
%2. Create a time vector to represent the time points at which you want to evaluate the system.
%3. Use a numerical solver like the Runge-Kutta method or MATLAB's built-in ODE solver, ode45, 
% to integrate the Lorenz equations over time to obtain the trajectory of the system.
% we integrate the equations to numerically compute the system's state at different points in time. 
% Integrating the equations essentially means finding a sequence of discrete states that approximate 
% the continuous evolution of the system over time.
%The 2D plot of the Lorenz attractor provided earlier was essentially a cross-section or projection 
% of the 3D chaotic trajectory onto the (x, y) plane. In this case, we sampled the (x, y) points 
% at discrete time intervals and plotted them in 2D. 
% This approach allows you to visualize a subset of the Lorenz attractor without explicitly simulating 
% the entire time evolution.
%4.Extract the x, y, and z coordinates from the solution.
%5.Finally, create a 3D plot of the Lorenz attractor.

% Define Lorenz system parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Define initial conditions
initialX = 1;
initialY = 0;
initialZ = 0;

% Create a time vector
tspan = linspace(0, 50, 5000); % Adjust the time span and resolution as needed

% Define the Lorenz system as a set of differential equations
lorenzEquations = @(t, XYZ) [sigma * (XYZ(2) - XYZ(1)); ...
                             XYZ(1) * (rho - XYZ(3)) - XYZ(2); ...
                             XYZ(1) * XYZ(2) - beta * XYZ(3)];

% Use the ode45 solver to integrate the system
[t, XYZ] = ode45(lorenzEquations, tspan, [initialX, initialY, initialZ]);%Runge cutta

% Extract the x, y, and z coordinates
x = XYZ(:, 1);
y = XYZ(:, 2);
z = XYZ(:, 3);

% Create a 3D plot
figure;
plot3(x, y, z);
title('Lorenz Attractor');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;



