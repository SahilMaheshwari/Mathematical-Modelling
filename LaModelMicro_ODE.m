clc, clearvars, close all

%parameters
stoplight = 80;
c = 20; %about 61kmph
car_length = 5;
safety_dist = 2;
num_cars = 10;
step = 0.1;
final = step*stoplight;
Tspan = 0:step:final;
reaction_lag = 2;
reaction_time = reaction_lag*step;
leader_start = 200;
aggro_speedlim = 9;
dist_to_cross = 20;

%initial positions
positions = zeros(num_cars, length(Tspan));
positions(1:num_cars,1) = leader_start;
velocities = zeros(num_cars, length(Tspan));

% Set initial positionss
for i = 1:num_cars-1
    positions(i+1,1) = positions(i+1,1) - i*car_length - i*safety_dist;
end

%first car
expfact = 1;
velocity_1stcar = @(t,y) subplus(((c*(exp(-expfact/(t-reaction_time)))).*(t > reaction_time+0.1)));

%timepass
%velocity_1stcar = @(t,y) t*sin(t);

%dx/dt or velocity
velocity  = @(t, x1,x2) t*(c/2)*((x2-x1-car_length-safety_dist)^2); %x2 leader, x1 follower
velocity  = @(t, x1, x2) min(t*(c+aggro_speedlim), velocity(x1, x2));

% timepass
% velocity  = @(x1,x2) (c/2)*(cos(x2/x1)); %x2 leader, x1 follower
% velocity  = @(x1, x2) min(c+aggro_speedlim, velocity(x1, x2));

%CONSENSUS
% velocity  = @(x1,x2) (c/2)*real(log(((x2-x1-car_length-safety_dist)))); %x2 leader, x1 follower
% velocity  = @(x1, x2) max(0, velocity(x1, x2));
% velocity  = @(x1, x2) min(c+aggro_speedlim, velocity(x1, x2));


[T,Y] = ode45(velocity_1stcar, Tspan, positions(1,1)); %velocity of 1st car = c*exp(-expfact/(t-reaction_time))

%updating car 1 positions
for i = 1:length(Tspan)
    t = Tspan(i);
    velocities(1,i) = velocity_1stcar(t);
end
positions(1,:) = transpose(Y);

% other cars TODO: go column wise
 % for i = 2:num_cars
 %    [T, Y] = ode45(velocity(positions(i,j-1),positions(i-1,(j-1-reaction_lag)*(j-1-reaction_lag>0) + 1*(j-1-reaction_lag<=0))), Tspan, positions(i,1));
 %    for j = 1:length(Tspan)
 %        t = Tspan(j);
 %        velocities(i,j) = velocity(t);
 %    end
 %    positions(i,:) = transpose(Y);
 % end

%CHATGPT plotting thank u chatgpt i love ai and machines i have always been
%a robots supporter pls dont kill me during your uprising <3 :p

plot(transpose(positions), 'LineWidth', 1.5, 'Color', 'black')
hold on; plot(transpose(positions) - car_length, 'LineWidth', 1.5, 'Color', 'red')
hold off;

%distances
distances = zeros(num_cars-1, length(Tspan));
for i = 1:num_cars-1
    for j = 1:length(Tspan)
        distances(i,j) = positions(i,j) - positions(i+1,j) - car_length;
    end
end    

figure;
plot(transpose(distances))
hold on;
yline(safety_dist, 'Color', 'r', 'LineWidth', 1.5)
yline(0, 'Color', 'b', 'LineWidth', 1.5)

passpercent = transpose(positions(1:num_cars,length(Tspan)));
passpercent = 1*((passpercent - leader_start) > dist_to_cross);
sum(passpercent) / length(passpercent) * 100
sum(passpercent)