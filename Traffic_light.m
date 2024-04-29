clc, clearvars, close all

%parameters
c = 17; %about 61kmph
car_length = 5;
safety_dist = 2;
num_cars = 10;
step = 0.1;
onesec = 1/step; 
stoplight = 12; %12 = 120secs
final = stoplight*onesec;
Tspan = 0:step:final;
reaction_lag = 2*onesec;
leader_start = 200;
aggro_speedlim = 9;
dist_to_cross = 10;
reactivity = 2; %how reactive are they to distance

%initial positions
positions = zeros(num_cars, length(Tspan));
positions(1:num_cars,1) = leader_start;
velocities = zeros(num_cars, length(Tspan));

% Set initial positionss
for i = 1:num_cars-1
    positions(i+1,:) = positions(i+1,1) - i*car_length - i*safety_dist;
end

%first car
expfact = 20;
velocity_1stcar = @(t,y) (c*(exp(-expfact/((t*onesec)-reaction_lag))))*(exp(-expfact/((t*onesec)-reaction_lag)) > 0)*((t*onesec) > reaction_lag);

%timepass
%velocity_1stcar = @(t,y) t;

%dx/dt or velocity
velocity  = @(x1,x2) (c/2)*((x2-x1-car_length-safety_dist)^reactivity) ; %x2 leader, x1 follower
velocity  = @(x1, x2) min(c+aggro_speedlim, velocity(x1, x2));

% timepass
% velocity  = @(x1,x2) (c/2)*(cos(x2/x1)); %x2 leader, x1 follower
% velocity  = @(x1, x2) min(c+aggro_speedlim, velocity(x1, x2));

%CONSENSUS
% velocity  = @(x1,x2) (c/2)*real(log(((x2-x1-car_length-safety_dist)))); %x2 leader, x1 follower
% velocity  = @(x1, x2) max(0, velocity(x1, x2));
% velocity  = @(x1, x2) min(c+aggro_speedlim, velocity(x1, x2));


[T,Y] = ode45(velocity_1stcar, Tspan/onesec, positions(1,1)); %velocity of 1st car = c*exp(-expfact/(t-reaction_lag))

%updating car 1 velocities
for i = 2:length(Tspan)
    t = Tspan(i)/onesec;
    velocities(1,i) = velocity_1stcar(t);
    %positions(1,i) = positions(1,i-1) + velocities(1,i)*step;   ummm???
end
positions(1,:) = transpose(Y);

% other cars
crash_occurred = false;
 for i = 2:num_cars
     for j = 2:length(Tspan)
         velocities(i,j) = velocity(positions(i,j-1),positions(i-1,((j-1-reaction_lag)*(j-1-reaction_lag>0) + 1*(j-1-reaction_lag<=0))));
         positions(i,j) = positions(i,j-1) + (velocities(i,j)*step);
         
         %check if crash
         if positions(i,j) >= positions(i-1, j) + car_length 
            fprintf('CRASH between %d and %d\n', i, i-1);
            crash_occurred = true;
            crash_by = i;
            crash_at = j;
            break;
         end
     end
     if crash_occurred
         break;
     end
 end

 if crash_occurred
    for i = crash_by:num_cars
        for j = crash_at:length(Tspan)
            positions(i,j) = positions(i,j-1);
        end
    end
 end

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

figure;
plot(transpose(velocities))

passpercent = transpose(positions(1:num_cars,length(Tspan)));
passpercent = 1*((passpercent - leader_start) > dist_to_cross);
sum(passpercent) / length(passpercent) * 100
sum(passpercent)