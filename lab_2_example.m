close all
% load the workspace from the prac first.
data = load("matlab_lab_3.mat");

disp(data)
set(0, 'DefaultLineLineWidth', 1.25); % make png plots look better


OP_3 = data.simOut_20250909_145639; % replace with your simout struct for operating point 1



OP_1 = data.simOut_20250909_143502; % replace with your simout struct for operating point 2

OP_2 = data.simOut_20250909_144508; % replace with your simout struct for operating point 3



time_op1 = OP_1.tout;
h1_op1 = OP_1.measurements.Tank_1_Level__m_.Data;
h2_op1 = OP_1.measurements.Tank_2_Level__m_.Data;
step_op1 = OP_1.ref_signal.Data;

time_op2 = OP_2.tout;
h1_op2 = OP_2.measurements.Tank_1_Level__m_.Data;
h2_op2 = OP_2.measurements.Tank_2_Level__m_.Data;
step_op2 = OP_2.ref_signal.Data;




time_op3 = OP_3.tout; % time
h1_op3 = OP_3.measurements.Tank_1_Level__m_.Data; % tank 1 height 
h2_op3 = OP_3.measurements.Tank_2_Level__m_.Data; % tank 2 height
step_op3 = OP_3.ref_signal.Data; % input voltage



Fig_1 = figure
plot(time_op3, smooth(h1_op3)); % may need smoothing - depending on how noisy the signal is
hold on;
plot(time_op3, smooth(h2_op3));
hold on;
plot(time_op3, step_op3);
title("Operating Point 3")
legend("Tank 1", "Tank 2", "Voltage", 'Location', 'southeast')
saveas(gcf, 'ENG331_Lab_2_OP_3.png'); % matlab crashes if i render it. So save it to a png

Fig_2 = figure
plot(time_op1, smooth(h1_op1));
hold on;
plot(time_op1, smooth(h2_op1));
hold on;
plot(time_op1, step_op1);
title("Operating Point 1")
legend("Tank 1", "Tank 2", "Voltage", 'Location','southeast')
saveas(gcf, 'ENG331_Lab_2_OP_1.png');


Fig_3 = figure
plot(time_op2, smooth(h1_op2));
hold on;
plot(time_op2, smooth(h2_op2));
hold on;
plot(time_op2, step_op2);
title("Operating Point 2")
legend("Tank 1", "Tank 2", "Voltage", 'Location','southeast')
saveas(gcf, 'ENG331_Lab_2_OP_2.png');

Fig_4 = figure
plot(time_op3, smooth(h1_op3)); hold on;
plot(time_op3, smooth(h2_op3));
plot(time_op3, step_op3);
title("2nd Order Estimations vs Empirical Model");
legend("h1", "h2", "Step Input",'Location','southeast');

% --- Estimated transfer function ---
num = 0.1344;
den = [1 0.4704 0.0384];
G = tf(num, den);

% Build Δu relative to first input value
du = step_op3 - step_op3(1);

% Use Δu → Δy and add initial operating point of h1
dy_hat = lsim(G, du, time_op3);   % modelled delta output
y_hat = 14.5+dy_hat ;

% Overlay model response
plot(time_op3, y_hat, 'k', 'LineWidth', 1.5, 'DisplayName', '2nd Order Estimation');

% Control Toolbox Estimation
% Build Δu relative to first input value
du = step_op3 - step_op3(1);

% Use Δu → Δy and add initial operating point of h1
num2=0.01126;
den2=[1 0.129 0.005497];
tf1=tf(num2,den2);
dy_hat2 = lsim(tf1, du, time_op3);   % modelled delta output
y_hat2 = 14.5+dy_hat2 ;

% Overlay model response
plot(time_op3, y_hat2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Control Systems Toolbox Estimation');

saveas(gcf, 'ENG331_Lab_2_Comparison.png');

Fig_5 = figure
plot(time_op1, smooth(h1_op1)); 
hold on;
plot(time_op1, smooth(h2_op1));
hold on;
plot(time_op1, step_op1);
hold on;
title("2nd Order Estimations vs Operating Point 3");
legend("h1", "h2", "Step Input",'Location','southeast');
du1 = step_op1 - step_op1(1);
dy_hat3 = lsim(G, du1, time_op1);   % modelled delta output
y_hat3 = 19.5+dy_hat3 ;
plot(time_op1, y_hat3, 'k', 'LineWidth', 1.5, 'DisplayName', 'Theoretical Estimation');
dy_hat4 = lsim(tf1, du1, time_op1);   % modelled delta output
y_hat4 = 19.5+dy_hat4 ;
plot(time_op1, y_hat4, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Control Systems Toolbox Estimation');

saveas(gcf, 'ENG331_Lab_2_Comparison_OP_3.png');

% find final value, gain, settling time, rise time ....


[Final_Value_op1_h1, Gain_op1_h1, Settling_Time_op1_h1] = find_vals(smooth(h1_op1), step_op1, time_op1)
[Final_Value_op2_h1, Gain_op2_h1, Settling_Time_op2_h1] = find_vals(smooth(h1_op2), step_op2, time_op2)
[Final_Value_op3_h1, Gain_op3_h1, Settling_Time_op3_h1] = find_vals(smooth(h1_op3), step_op3, time_op3)

[Final_Value_op1_h2, Gain_op1_h2, Settling_Time_op1_h2] = find_vals(smooth(h2_op1), step_op1, time_op1)
[Final_Value_op2_h2, Gain_op2_h2, Settling_Time_op2_h2] = find_vals(smooth(h2_op2), step_op2, time_op2)
[Final_Value_op3_h2, Gain_op3_h2, Settling_Time_op3_h2] = find_vals(smooth(h2_op3), step_op3, time_op3)


% Non linear simulink model
t = 0:0.1:100;
steps = [9, 8; 9, 10; 5, 4; 5, 6; 7, 6; 7, 8];

Non_linear_data = []; % array of the data that the simulink model will output
open_system('non_linear_h1_h2');

[numRows, numCols] = size(steps);
for idx = 1:numRows
    i = steps(idx, 1)*ones(size(t));
    i(t >= 1) = steps(idx, 2);
    sim_input = timeseries(i, t);
    assignin('base', 'sim_input', sim_input)

    simOut = sim('non_linear_h1_h2');
    results = simOut.get('logsout');
    Non_linear_data = [Non_linear_data, results];
end

fprintf("h1 non-linear:\n");
fprintf("| Voltage | step | Final Value | Gain | Settling Time|\n");
for idx = 1:length(Non_linear_data)
    d = Non_linear_data(idx);
    step = d{1}.Values.Data(1, 1, :);
    h1 = d{2}.Values(:, 1).Data;
    time1 = d{2}.Values(:, 1).Time;
    h2 = d{3}.Values(:, 1).Data;
    time2 = d{3}.Values(:, 1).Time;
    %fprintf("For %i to %i\n", step(1), step(end));
    [FV, G, ST] = find_vals(h1, step, time1);
    fprintf("| %i | $%iV$ | %i | %i | %i|\n", step(1), step(end) - step(1), FV, G, ST);
end

fprintf("h2 non-linear:\n");
fprintf("| Voltage | step | Final Value | Gain | Settling Time|\n");
for idx = 1:length(Non_linear_data)
    d = Non_linear_data(idx);
    step = d{1}.Values.Data(1, 1, :);
    h1 = d{2}.Values(:, 1).Data;
    time1 = d{2}.Values(:, 1).Time;
    h2 = d{3}.Values(:, 1).Data;
    time2 = d{3}.Values(:, 1).Time;
    [FV, G, ST] = find_vals(h2, step, time2);
    %fprintf("For %i to %i\n", step(1), step(end));
    fprintf("| %i | $%iV$ | %i | %i | %i|\n", step(1), step(end) - step(1), FV, G, ST);
end

% linearsied simulink model
t = 0:0.1:100;
steps = [9, 8; 9, 10; 5, 4; 5, 6; 7, 6; 7, 8];

linear_data = []; % array of the data that the simulink model will output
open_system('linear_h1');

[numRows, numCols] = size(steps);
for idx = 1:numRows
    i = steps(idx, 1)*ones(size(t));
    i(t >= 1) = steps(idx, 2);
    sim_input = timeseries(i, t);
    assignin('base', 'sim_input', sim_input)

    simOut = sim('linear_h1');
    results = simOut.get('logsout');
    linear_data = [linear_data, results];
end

fprintf("h1 linear:\n");
fprintf("| Voltage | step | Final Value | Gain | Settling Time|\n");
for idx = 1:length(linear_data)
    d = linear_data(idx);
    step = d{2}.Values.Data(1, 1, :);
    h1 = d{1}.Values(:, 1).Data;
    time1 = d{2}.Values(:, 1).Time;
    %fprintf("For %i to %i\n", step(1), step(end));
    [FV, G, ST] = find_vals(h1, step, time1);

    fprintf("| %i | $%iV$ | %i | %i | %i|\n", step(1), step(end) - step(1), FV, G, ST);
end

%h2 linear
t = 0:0.1:200;
steps = [9, 8; 9, 10; 5, 4; 5, 6; 7, 6; 7, 8];

linear_data = []; % array of the data that the simulink model will output
open_system('linear_h2');

[numRows, numCols] = size(steps);
for idx = 1:numRows
    i = steps(idx, 1)*ones(size(t));
    i(t >= 100) = steps(idx, 2);
    sim_input = timeseries(i, t);
    assignin('base', 'sim_input', sim_input)

    simOut = sim('linear_h2');
    results = simOut.get('logsout');
    linear_data = [linear_data, results];
end


fprintf("h2 linear:\n");
fprintf("| Voltage | step | Final Value | Gain | Settling Time|\n");
for idx = 1:length(linear_data)
    d = linear_data(idx);
    step = d{2}.Values.Data(1, 1, :);
    h1 = d{1}.Values(:, 1).Data;
    time1 = d{2}.Values(:, 1).Time;
    %fprintf("For %i to %i\n", step(1), step(end));
    [OP, FV, G, ST] = find_vals(h1, step, time1);
    a = 0;
    for idx = 1:length(FV)
        for jdx = 1:length(ST)
            for kdx = 1:length(G)
                fprintf("%i | (%i, %i) | $%iV$ | %i | %i | %i|\n", a, OP , step(1), step(end) - step(1), FV(idx), G(kdx), ST(jdx));
                a = a + 1;
            end
        end
    end
end

function [OP, FV_array, G_array, ST_array] = find_vals(h, step, time)
    dif = diff(step);
    OP_array = [];
    FV_array = [];
    G_array = [];
    ST_array = [];
    for i = 1:length(dif)
        if (dif(i) ~= 0)
            start_pos = i;
            end_pos = length(dif);
            for j = (start_pos + 1):length(dif)
                if (dif(j) ~= 0)
                    end_pos = j;
                    break;
                end
            end
            i = end_pos - 1;
            ST_pos = settling_time(h, start_pos, end_pos);
            FV = final_value(h, ST_pos, end_pos);
            G = (FV - h(start_pos))/(step(end_pos) - step(start_pos));
            
            OP = inital_value(h, start_pos);
            OP_array = [OP_array, OP];
            FV_array = [FV_array, FV];
            G_array = [G_array, G];
            ST_array = [ST_array, time(ST_pos) - time(start_pos)];
        end
    end

    
end

function ST = settling_time(h, start_pos, end_pos) % This function isnt great
    percent = 2/100;
    IV = h(start_pos);
    FV = final_value(h, start_pos, end_pos);
    tol = percent*FV;
    LB = FV - tol;
    UB = FV + tol;
    ST = NaN;
    for i = start_pos:end_pos
        if (h(i) >= LB && h(i) <= UB)
            ST = i;
            break;
            %if (all(h(i:end_pos) >= LB & h(i:end_pos) <= LB))
            %    ST = i;
            %    break;
            %end
        end
    end
end

function FV = final_value(h, ST_pos, end_pos)
    FV = mean(h(ST_pos:end_pos));
end

function OP = inital_value(h, OP_pos)
    if OP_pos > 50
        OP = mean(h(((OP_pos - 50):OP_pos)));
    else
        OP = h(OP_pos);
    end
end
