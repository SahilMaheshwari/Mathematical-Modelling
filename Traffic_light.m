clc, clearvars, close all

global leader_start dist_to_cross c car_length safety_dist num_cars reaction_lag step onesec;

%TOSHHOW
%Animate greenlight
%Anitmate yellow light
%

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
[positions, velocities] = yellowlight(positions, velocities, 10, 2);
moviemaker(positions)

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
    velocity_1stcar = @(t,y) ((c/2)*(exp(-expfact/((t*onesec)-(reaction_lag/6)))))*(exp(-expfact/((t*onesec)-(reaction_lag/6))) > 0)*((t*onesec) > reaction_lag/6);
    %
    %velocity_1stcar = @(t,y) custom func here

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
    % plot(transpose(positions), 'LineWidth', 1.5, 'Color', 'black')
    % hold on; yline(200, 'Color', 'green', LineWidth=1)
    % hold on; plot(transpose(positions) - car_length, 'LineWidth', 1.5, 'Color', 'red')
    % hold off;
    
    %distances
    distances = zeros(num_cars-1, length(Tspan));
    for i = 1:num_cars-1
        for j = 1:length(Tspan)
            distances(i,j) = positions(i,j) - positions(i+1,j) - car_length;
        end
    end    
    
    % figure;
    % plot(transpose(distances))
    % hold on;
    % yline(safety_dist, 'Color', 'r', 'LineWidth', 1.5)
    % yline(0, 'Color', 'b', 'LineWidth', 1.5)
    % 
    % figure;
    % plot(transpose(velocities))
    
    passpercent = transpose(positions(1:num_cars,length(Tspan)));
    passpercent = 1*((passpercent - leader_start) > dist_to_cross);
    sum(passpercent) / length(passpercent) * 100
    sum(passpercent)

end

function [positions, velocities] = yellowlight(positions, velocities, stoplight, carrate)
    global c car_length safety_dist reaction_lag step onesec; 
    reaction_lag = 0;
    final = stoplight*onesec;
    Tspan = 0:step:final;
    aggrospeed = 4.5; %bigger this number, lesser number go through

    %removing all cars post 200m mark length
    rows_to_remove = positions(:, end) > 200;
    positions(rows_to_remove, :) = [];
    velocities(rows_to_remove, :) = [];

    %new number of cars
    num_cars = floor(carrate*stoplight/2);
    positions = padarray(positions, [num_cars, 0], 0, 'post');
    velocities = padarray(velocities, [num_cars, 0], 6.5, 'post');
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
        positions(i,1) = positions(i-1,1) - car_length - 2*safety_dist;
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
        dvdt = @(t,y) [y(2); ((c/2)*sign(y(1)-200)*exp(-1/((y(1)-200)*sign(y(1)-200)))*(1-(180/y(1))))];
        [T,Y] = ode45(dvdt, Tspan/onesec, [positions(car_num,1), velocities(car_num,1)]);
        Y = real(Y);
    
        positions(car_num,:) = transpose(Y(:,1));
        velocities(car_num,:) = transpose(Y(:,2));
        
        %checks if youll be fast enough to cross
        for j = 1:length(positions)
            if positions(car_num,j) >= 200
                if velocities(car_num,j) < aggrospeed
                    fastatstart = false;
                    carpassed = false;
                    break
                end
                break
            end
        end
        
        %checks if you will crash into someone
        if car_num > 1
            for j = 1:length(positions)
                if positions(car_num, j) >= positions(car_num-1, j) - car_length - safety_dist
                    carpassed = false;
                end
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
            car_not_pass = i;
        else
            car_not_pass = i-1;
            break
        end
    end

    if car_not_pass <= num_cars
        car_not_pass

    %for the first car to not pass
    positions(car_not_pass,:) = positions_copy (car_not_pass,:);
    velocities(car_not_pass,:) = velocities_copy(car_not_pass,:);
    
    %to critically dampen
    k = 0.92;
    b = 2 * sqrt(c);

    dvdt = @(t,y) [y(2); 200-b*y(2)-k*y(1)];
    [T,Y] = ode45(dvdt, Tspan/onesec, [positions(car_not_pass,1), velocities(car_not_pass,1)]);
    positions(car_not_pass,:) = real(transpose(Y(:,1)));
    velocities(car_not_pass,:) = real(transpose(Y(:,2)));
    
        %other cars follow suit
        velocity  = @(x1,x2) ((c/2)*((x2-x1-car_length-1*safety_dist)^2))*(x2-x1-car_length-1*safety_dist > 0);
         for i = car_not_pass+1:num_cars
             for j = 2:length(Tspan)
                 velocities(i,j) = (velocity(positions(i,j-1),positions(i-1,((max((floor(j-1-reaction_lag))*(j-1-reaction_lag>0),1) + 1*(j-1-reaction_lag<=0))))));
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
    hold on; yline(200, 'Color', 'green', LineWidth=1)
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

function [] = moviemaker(thematrix)
    
    thematrix = thematrix(:, 1:size(thematrix,2)-1);

    % Define compression factor
    compressionFactor = 100;
    
    % Calculate the number of columns in the compressed matrix
    numCompressedColumns = size(thematrix, 2) / compressionFactor;
    
    % Preallocate memory for the compressed matrix
    compressedMatrix = zeros(size(thematrix, 1), numCompressedColumns);
    
    % Calculate averages of each block
    for i = 1:numCompressedColumns
        % Define start and end indices of the block
        startIndex = (i - 1) * compressionFactor + 1;
        endIndex = i * compressionFactor;
        
        % Extract the block from the original matrix
        block = thematrix(:, startIndex:endIndex);
        
        % Calculate the average of the block and store it in the compressed matrix
        compressedMatrix(:, i) = mean(block, 2);
    end


    thematrix = compressedMatrix;

    myVideo = VideoWriter('movement on traffic light off'); %open video file
    myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    open(myVideo)

    % Plot the initial point
    cars = thematrix(1,1);
    p1 = plot(1, thematrix(1,1), 'o', 'MarkerFaceColor', 'red');
    hold on
    p2 = plot(1, thematrix(2,1), 'o', 'MarkerFaceColor', 'red');
    p3 = plot(1, thematrix(3,1), 'o', 'MarkerFaceColor', 'red');
    p4 = plot(1, thematrix(4,1), 'o', 'MarkerFaceColor', 'red');
    p5 = plot(1, thematrix(5,1), 'o', 'MarkerFaceColor', 'red');
    p6 = plot(1, thematrix(6,1), 'o', 'MarkerFaceColor', 'red');
    p7 = plot(1, thematrix(7,1), 'o', 'MarkerFaceColor', 'red');
    p8 = plot(1, thematrix(8,1), 'o', 'MarkerFaceColor', 'red');
    p9 = plot(1, thematrix(9,1), 'o', 'MarkerFaceColor', 'red');
    p10 = plot(1, thematrix(10,1), 'o', 'MarkerFaceColor', 'red');

    zoom(0.5)
    grid on
    axis manual
    
    % Set initial y-axis limits
    ymin = min(thematrix(1, :));
    ymax = max(thematrix(1, :));
    ylim([ymin-100, ymax+100]);
    xlim([0.5,1.5])
    
    % Iterate through each column of thematrix after the first one
    for k = 2:size(thematrix, 2)
        % Update the YData of the point to animate its movement along the y-axis
        p1.YData = thematrix(1,k);
        p2.YData = thematrix(2,k);
        p3.YData = thematrix(3,k);
        p4.YData = thematrix(4,k);
        p5.YData = thematrix(5,k);
        p6.YData = thematrix(6,k);
        p7.YData = thematrix(7,k);
        p8.YData = thematrix(8,k);
        p9.YData = thematrix(9,k);
        p10.YData = thematrix(10,k);
        drawnow
        pause(0.001) %Pause and grab frame
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
    end
    close(myVideo)
end
