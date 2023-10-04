SOCOCV = importdata("Fine Murata VTC6 SOC OCV Curve.txt"); %Import the cell SOC OCV curve

%Import the Michigan Accel power draw data
Larosa1 = importdata("UT23 Power Draw\Michigan Larosa Accel 1.csv");
Larosa2 = importdata("UT23 Power Draw\Michigan Larosa Accel 2.csv");
Ball1 = importdata("UT23 Power Draw\Michigan Ball Accel 1.csv");
Ball2 = importdata("UT23 Power Draw\Michigan Ball Accel 2.csv");

%Variables
Scount = 115;       %Pack cell series count
Pcount = 5;         %Pack cell parallel count
R_cell = 0.0225;     %Cell internal resistance in Ohm
R_busbars = 0.15;    %Resistance of busbars and other components in the high current path in Ohms
SOC_init = 85;      %Initial SOC of the pack

%Pack parameters
R_pack = R_cell * Scount/Pcount + R_busbars %Total pack internal resistance in Ohms
[value, idx] = min(abs(SOCOCV(:,1)-SOC_init/100)); %Retrieve the index of the closest SOC-OCV point
Cell_OCV = SOCOCV(idx,2); %Retrieve the cell open circuit voltage
Pack_OCV = Scount * Cell_OCV; %Calculate the pack OCV using the cell OCV and Scount

Larosa1_results = zeros(length(Larosa1),3);

for t = 1:length(Larosa1)
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * Larosa1(t,2)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.8
        disp("Cell undervoltage fault at " + string(Larosa1(t)) + "seconds in Larosa 1")
    end
    Larosa1_results(t,1) = Larosa1(t,1);
    Larosa1_results(t,2) = V_cell;
    Larosa1_results(t,3) = I_pack;
end

Larosa2_results = zeros(length(Larosa2),3);

for t = 1:length(Larosa2)
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * Larosa2(t,2)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.8
        disp("Cell undervoltage fault at " + string(Larosa2(t)) + "seconds in Larosa 2")
    end
    Larosa2_results(t,1) = Larosa2(t,1);
    Larosa2_results(t,2) = V_cell;
    Larosa2_results(t,3) = I_pack;
end

Ball1_results = zeros(length(Ball1),3);

for t = 1:length(Ball1)
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * Ball1(t,2)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.8
        disp("Cell undervoltage fault at " + string(Ball1(t)) + "seconds in Ball 1")
    end
    Ball1_results(t,1) = Ball1(t,1);
    Ball1_results(t,2) = V_cell;
    Ball1_results(t,3) = I_pack;
end

Ball2_results = zeros(length(Ball2),3);

for t = 1:length(Ball2)
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * Ball2(t,2)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.8
        disp("Cell undervoltage fault at " + string(Ball2(t)) + "seconds in Ball 2")
    end
    Ball2_results(t,1) = Ball2(t,1);
    Ball2_results(t,2) = V_cell;
    Ball2_results(t,3) = I_pack;
end

%% 

%Electrical Plots

voltage_plot = figure('visible','off');
plot(Larosa1_results(:,1),Larosa1_results(:,2),Larosa2_results(:,1),Larosa2_results(:,2),Ball1_results(:,1),Ball1_results(:,2),Ball2_results(:,1),Ball2_results(:,2));
title("Cell Voltage in Michigan Accel")
xlabel("Time (seconds)")
ylabel("Cell voltage (V)")
saveas(voltage_plot,"Plots/Michigan Acceleration/" + string(SOC_init) + "%SOC " + string(Scount) + "S " + string(R_pack) + "ohm MI Accel Voltage Plot.png")

current_plot = figure('visible','off');
plot(Larosa1_results(:,1),Larosa1_results(:,3),Larosa2_results(:,1),Larosa2_results(:,3),Ball1_results(:,1),Ball1_results(:,3),Ball2_results(:,1),Ball2_results(:,3));
title("Pack Current in Accel")
xlabel("Time (seconds)")
ylabel("Pack Current (A)")
saveas(current_plot,"Plots/Michigan Acceleration/" + string(SOC_init) + "%SOC " + string(Scount) + "S " + string(R_pack) + "ohm MI Accel Current Plot.png")