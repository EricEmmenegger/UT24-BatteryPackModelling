%Import the cell SOC OCV curve and Michigan Endurance power draw data
SOCOCV = importdata("Fine Murata VTC6 SOC OCV Curve.txt");
endurance_data = importdata("UT23 Power Draw\Michigan Endurance.csv");

%Variables
Scount = 115;       %Pack cell series count
Pcount = 5;         %Pack cell parallel count
R_cell = 0.0225;    %Cell internal resistance in Ohm
R_busbars = 0.15;   %Resistance of busbars and other components in the high current path in Ohms
SOC_init = 94;      %Initial SOC of the pack, this was made to match 2023 data
Cp_batt = 960;      %Cell heat capacity in J/(kg*K)
m_batt = 46.6/1000; %Cell mass in kilograms
T_init = 30;        %Initial cell temperature

%Pack parameters
R_pack = R_cell * Scount/Pcount + R_busbars         %Total pack internal resistance in Ohms
Q_batt = SOC_init/100 * 3000;                       %Calculate the initial cell capacity in Ah
SOC = SOC_init;                                     %Initialize SOC variable
T_cell = T_init;                                         %Initialize Temperature variable
[value, idx] = min(abs(SOCOCV(:,1)-SOC/100));       %Retrieve the index of the closest SOC-OCV point
Cell_OCV = SOCOCV(idx,2);                           %Retrieve the cell open circuit voltage
Pack_OCV = Scount*Cell_OCV;                         %Calculate the initial pack open circuit voltage

endurance_results = zeros(length(endurance_data),7);

for t=1:length(endurance_data)
    %Calculating the pack current and cell voltage under load
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * endurance_data(t,3)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.8
        disp("Cell undervoltage fault at " + string(endurance_data(t)) + "seconds")
    end
    %Calculating the cell heat generation and temperature over time
    Qgen_cell = R_cell*(I_pack/Pcount)^2;
    Qaccum_cell = Qgen_cell * 0.05;
    T_cell = T_cell + Qaccum_cell/(Cp_batt*m_batt);
    %This section of the loop records the results in a table
    endurance_results(t,1) = endurance_data(t);
    endurance_results(t,2) = V_cell;
    endurance_results(t,3) = I_pack;
    endurance_results(t,4) = Qgen_cell;
    endurance_results(t,5) = SOC;
    endurance_results(t,6) = Cell_OCV;
    endurance_results(t,7) = T_cell;
    %This section of the loop subtracts the amount of SOC used and
    %determines the new cell OCV
    Q_batt = Q_batt - 0.05*I_pack/(3.6*Pcount);
    SOC = Q_batt/30;
    [value, idx] = min(abs(SOCOCV(:,1)-SOC/100));
    Cell_OCV = SOCOCV(idx,2);
    Pack_OCV = Scount*Cell_OCV;
end

%% 

%Electrical Plots

pack_voltage_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,2)*Scount);
title("Pack Voltage in Endurance")
xlabel("Time (seconds)")
ylabel("Pack voltage (V)")
saveas(pack_voltage_plot,"Electrical Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance Voltage Plot.png")

cell_voltage_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,2));
title("Cell Voltage in Endurance")
xlabel("Time (seconds)")
ylabel("Cell voltage (V)")
saveas(cell_voltage_plot,"Electrical Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance Cell Voltage Plot.png")

ocv_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,6));
title("Cell OCV in Endurance")
xlabel("Time (seconds)")
ylabel("OCV")
saveas(ocv_plot,"Electrical Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance OCV Plot.png")

current_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,3));
title("Pack Current in Endurance")
xlabel("Time (seconds)")
ylabel("Current (A)")
saveas(current_plot,"Electrical Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance Pack Current Plot.png")

soc_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,5));
title("Pack SOC in Endurance")
xlabel("Time (seconds)")
ylabel("SOC (%)")
saveas(soc_plot,"Electrical Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance Pack SOC Plot.png")

%% 

%Heating Plots

qgen_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,4));
title("Heat Generation per Cell in MI Endurance");
xlabel("Time (seconds")
ylabel("Heat Generation (W)")
saveas(qgen_plot,"Heating Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance Cell Heat Generation Plot.png")

T_plot = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
plot(endurance_results(:,1),endurance_results(:,7));
title("Cell Temperature in MI Endurance");
xlabel("Time (seconds")
ylabel("Cell Temperature (Celsius)")
saveas(T_plot,"Heating Plots/Michigan Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm MI Endurance Cell Temperature Plot.png")