clc, clearvars, close all

global leader_start dist_to_cross c car_length safety_dist num_cars reaction_lag step onesec;

leader_start = 200;
dist_to_cross = 0;
c = 15; %about 61kmph
car_length = 5;
safety_dist = 2;
num_cars = 10;
reaction_lag = 12*onesec;
step = 0.01;
onesec = 1/step; 


[positions, velocities] = initializer();
[positions, velocities] = greenlight(positions, velocities, 10);
[positions, velocities] = yellowlight(positions, velocities, 10, 1.5);

function [positions, velocities] = initializer()
    global leader_start car_length safety_dist num_cars onesec step;
    final = 10*onesec;
    Tspan = 0:step:final;

    %initial positions
    positions = zeros(num_cars, length(Tspan));
    positions(1:num_cars,1) = leader_start;
    velocities = zeros(num_cars, length(Tspan));
    
    % Set initial positionss
    for i = 1:num_cars-1
        positions(i+1,:) = positions(i+1,1) - i*car_length - i*safety_dist;
    end
end

function [positions, velocities] = greenlight(positions, velocities, stoplight)
    global leader_start dist_to_cross c car_length safety_dist num_cars reaction_lag step onesec;
    final = stoplight*onesec;
    Tspan = 0:step:final;

    %first car
    expfact = 20;
    velocity_1stcar = @(t,y) ((c/2)*(exp(-expfact/((t*onesec)-(onesec/6)))))*(exp(-expfact/((t*onesec)-(onesec/6))) > 0)*((t*onesec) > onesec/6);
    
    %dx/dt or velocity
    velocity  = @(x1,x2) ((c/2)*((x2-x1-car_length-safety_dist)^2))*(x2-x1-car_length-safety_dist > 0) ; %x2 leader, x1 follower
    %velocity  = @(x1, x2) min(c+aggro_speedlim, velocity(x1, x2));
    
    [T,Y] = ode45(velocity_1stcar, Tspan/onesec, positions(1,1)); %velocity of 1st car = c*exp(-expfact/(t-reaction_lag))
    
    %updating car 1 velocities
    for i = 2:length(Tspan)
        t = Tspan(i)/onesec;
        velocities(1,i) = velocity_1stcar(t);
    end
    positions(1,:) = transpose(Y);
    
    % other cars
    crash_occurred = false;
     for i = 2:num_cars
         for j = 2:length(Tspan)
             velocities(i,j) = velocity(positions(i,j-1),positions(i-1,((max((floor(j-1-reaction_lag))*(j-1-reaction_lag>0),1) + 1*(j-1-reaction_lag<=0)))));
             positions(i,j) = positions(i,j-1) + (velocities(i,j)*step*step);
             
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
    
    %plotting    
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

end

function [positions, velocities] = yellowlight(positions, velocities, stoplight, carrate)
    global c car_length safety_dist reaction_lag step onesec; 
    final = stoplight*onesec;
    Tspan = 0:step:final;
    aggrospeed = 5.0;

    %removing all cars post 200m mark length
    rows_to_remove = positions(:, end) > 200;
    positions(rows_to_remove, :) = [];
    velocities(rows_to_remove, :) = [];

    %new number of cars
    num_cars = floor(carrate*stoplight/2);
    positions = padarray(positions, [num_cars, 0], 0, 'post');
    velocities = padarray(velocities, [num_cars, 0], 0.5, 'post');
    old_cars = size(positions,1) - num_cars;
    num_cars = size(positions,1);

    %setting cars who havent travelled past stoplight at col=1
    for i = 1:size(positions, 1)
        endval = positions(i,end);
        positions(i,:) = zeros(1, size(positions, 2));
        positions(i,1) = endval;
    end

    %positioning of new cars
    for i = old_cars+1:num_cars
        positions(i,1) = positions(i-1,1) - car_length - 4*safety_dist;
    end

    %same thing as above but for velocity
    for i = 1:size(velocities, 1)
        endval = velocities(i,end);
        velocities(i,:) = zeros(1, size(velocities, 2));
        velocities(i,1) = endval;
    end
    
    function [newtopcar] = topcar (car_num)
        newtopcar = false;
        fastatstart = true;
        carpassed = true;

        %velocity and position top car
        dvdt = @(t,y) [y(2); 0.03+log(200-1*y(1))-2];
        [T,Y] = ode45(dvdt, Tspan/onesec, [positions(car_num,1), velocities(car_num,1)]);
        Y = real(Y);
    
        positions(car_num,:) = transpose(Y(:,1));
        velocities(car_num,:) = transpose(Y(:,2));

        for i = 1:length(positions)
            if positions(car_num,i) >= 200
                if velocities(car_num,i) < aggrospeed
                    fastatstart = false;
                    carpassed = false;
                    break
                end
                break
            end
        end
        if carpassed
            newtopcar = true;
        end
    end

    %making copies
    positions_copy  = positions;
    velocities_copy = velocities;

    first_car_pass = topcar(1);
    car_not_pass = 2;
    for i = 2:num_cars
        if first_car_pass
            first_car_pass = topcar(i);
            car_not_pass = i+1;
        else
            car_not_pass = i;
            break
        end
    end

    if car_not_pass <= num_cars
        car_not_pass

    %for the first car to not pass
    positions(car_not_pass,:) = positions_copy (car_not_pass,:);
    velocities(car_not_pass,:) = velocities_copy(car_not_pass,:);
    
    %to critically dampen
    c = 2;
    b = 2 * sqrt(c);

    dvdt = @(t,y) [y(2); 200-b*y(2)-c*y(1)];
    [T,Y] = ode45(dvdt, Tspan/onesec, [positions(car_not_pass,1), velocities(car_not_pass,1)]);
    positions(car_not_pass,:) = real(transpose(Y(:,1)));
    velocities(car_not_pass,:) = real(transpose(Y(:,2)));
    
        %other cars follow suit
        velocity  = @(x1,x2) ((c/2)*((x2-x1-car_length-safety_dist)^2))*(x2-x1-car_length-safety_dist > 0);
         for i = car_not_pass+1:num_cars
             for j = 2:length(Tspan)
                 velocities(i,j) = velocity(positions(i,j-1),positions(i-1,((max((floor(j-1-reaction_lag))*(j-1-reaction_lag>0),1) + 1*(j-1-reaction_lag<=0)))));
                 positions(i,j) = positions(i,j-1) + (velocities(i,j)*step*step);
                 
             end
         end
    end

    distances = zeros(num_cars-1, length(Tspan));
    for i = 1:num_cars-1
        for j = 1:length(Tspan)
            distances(i,j) = positions(i,j) - positions(i+1,j) - car_length;
        end
    end    
    
    plot(transpose(positions), 'LineWidth', 1.5, 'Color', 'black')
    hold on; plot(transpose(positions) - car_length, 'LineWidth', 1.5, 'Color', 'red')
    hold off;

    figure;
    plot(transpose(distances))
    hold on;
    yline(safety_dist, 'Color', 'r', 'LineWidth', 1.5)
    yline(0, 'Color', 'b', 'LineWidth', 1.5)
    
    figure;
    plot(transpose(velocities))

end
