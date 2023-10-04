%% Data Calculation

%Import the cell SOC OCV curve and Michigan Endurance power draw data
SOCOCV = importdata("Fine Murata VTC6 SOC OCV Curve.txt");
endurance_data = importdata("UT23 Power Draw\Czech Endurance.csv");

%Variables
Scount = 115;       %Pack cell series count
Pcount = 5;         %Pack cell parallel count
R_cell = 0.0225;    %Cell internal resistance in Ohm
R_busbars = 0.15;   %Resistance of busbars and other components in the high current path in Ohms
SOC_init = 94;      %Initial SOC of the pack, this was made to match 2023 data
Cp_batt = 960;      %Cell heat capacity in J/(kg*K)
m_batt = 46.6/1000; %Cell mass in kilograms
T_init = 30;        %Initial cell temperature
T_amb = 30;         %Ambient temperature
htc = 10;           %Heat transfer coefficient in W/(m^2 * K)
A_ht = pi * 0.018 * 0.0535 / 2; %Area of cell that transfers heat in m^2
%The above calculation assumes that 1/2 of the circumference of the cell
%gets enough airflow to effectively transfer heat, and that only the
%area of the cells not covered by the end caps of the Enepaq brick
%transfers heat.

%Pack parameters
R_pack = R_cell * Scount/Pcount + R_busbars         %Total pack internal resistance in Ohms
Q_batt = SOC_init/100 * 3000;                       %Calculate the initial cell capacity in Ah
SOC = SOC_init;                                     %Initialize SOC variable
T_cell_adiabatic = T_init;                          %Initialize Adiabatic Temperature variable
T_cell_cool = T_init;                               %Initialize cooled temperature variable
[value, idx] = min(abs(SOCOCV(:,1)-SOC/100));       %Retrieve the index of the closest SOC-OCV point
Cell_OCV = SOCOCV(idx,2);                           %Retrieve the cell open circuit voltage
Pack_OCV = Scount*Cell_OCV;                         %Calculate the initial pack open circuit voltage

endurance_results = zeros(length(endurance_data),8);

for t=1:length(endurance_data)
    %Calculating the pack current and cell voltage under load
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * endurance_data(t,2)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.8
        disp("Cell undervoltage fault at " + string(endurance_data(t)) + "seconds")
    end
    
    %Calculating the cell heat generation, heat accumulation and heat loss
    Qgen_cell = R_cell*(I_pack/Pcount)^2;               %Heat generated (W)
    Qcool_cell = htc * A_ht * (T_cell_cool - T_amb);    %Cooling (W)
    Qaccum_cell_adiabatic = Qgen_cell * 0.05;           %Heat accumulated in the cell under adiabatic conditions (J)
    Qaccum_cell_cool = (Qgen_cell - Qcool_cell) * 0.05; %Heat accumulated in the cell with cooling (J)

    %Calculating the cell temperature in degrees Celsius
    T_cell_adiabatic = T_cell_adiabatic + Qaccum_cell_adiabatic/(Cp_batt*m_batt);
    T_cell_cool = T_cell_cool + Qaccum_cell_cool/(Cp_batt*m_batt);

    %This section of the loop records the results in a table
    endurance_results(t,1) = endurance_data(t);
    endurance_results(t,2) = V_cell;
    endurance_results(t,3) = I_pack;
    endurance_results(t,4) = Qgen_cell;
    endurance_results(t,5) = SOC;
    endurance_results(t,6) = Cell_OCV;
    endurance_results(t,7) = T_cell_adiabatic;
    endurance_results(t,8) = T_cell_cool;

    %This section of the loop subtracts the amount of SOC used and
    %determines the new cell OCV
    Q_batt = Q_batt - 0.05*I_pack/(3.6*Pcount);
    SOC = Q_batt/30;
    [value, idx] = min(abs(SOCOCV(:,1)-SOC/100));
    Cell_OCV = SOCOCV(idx,2);
    Pack_OCV = Scount*Cell_OCV;
end

%% Plotting Section

plot_titles =  ["Cell Voltage (V)" "Pack Current (A)" "Cell Heat Generation (W)" "State of Charge (%)" "Cell OCV (V)" "Adiabatic Cell Temperature (deg C)" "Cooled Cell Temperature (deg C)"];
time = endurance_results(:,1);

for i=1:7
    current_figure = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
    data = endurance_results(:,i+1);
    plot(time, data);
    title(plot_titles(i) + " in Czech Endurance")
    xlabel("Time (seconds)")
    ylabel(plot_titles(i))
    saveas(current_figure, "Plots/Czech Endurance/" + string(Scount) + "S " + string(R_pack) + " ohm CZ Endurance " + plot_titles(i) + " Plot.png")
end