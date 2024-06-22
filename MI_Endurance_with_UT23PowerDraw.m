%% Data Calculations

%Import the cell SOC OCV curve and Michigan Endurance power draw data
SOCOCV = importdata("Fine Murata VTC6 SOC OCV Curve.txt");
endurance_data = importdata("UT23 Power Draw\Michigan Endurance.csv");
DCIR_LUT = importdata("DCIR Lookup Table.csv");

%Variables
Scount = 115;       %Pack cell series count
Pcount = 5;         %Pack cell parallel count
R_busbars = 0.15;   %Resistance of busbars and other components in the high current path in Ohms
SOC_init = 94;      %Initial SOC of the pack, this was made to match 2023 data
Cp_batt = 960;      %Cell heat capacity in J/(kg*K)
m_batt = 46.6/1000; %Cell mass in kilograms
T_init = 30;        %Initial cell temperature
T_amb = 30;         %Ambient temperature
htc = 50;           %Heat transfer coefficient in W/(m^2 * K)
A_ht = pi * 0.018 * 0.0535 / 2; %Area of cell that transfers heat in m^2
%The above calculation assumes that 1/2 of the circumference of the cell
%gets enough airflow to effectively transfer heat, and that only the
%area of the cells not covered by the end caps of the Enepaq brick
%transfers heat.
Cp_air = 1005;      %Specific heat capacity of air (J/kg*K)
mdot_air = 0.075;   %Mass Flow rate of air in kg/s


%Cell DC-IR Lookup
[closest_SOC,closest_SOC_id] = min(abs(SOC_init-DCIR_LUT(:,1)));
[closest_temp,closest_temp_id] = min(abs(T_init-DCIR_LUT(1,:)));
R_cell = DCIR_LUT(closest_SOC_id,closest_temp_id) / 1000;

%% Pack Parameters

Q_batt = SOC_init/100 * 3000;                       %Calculate the initial cell capacity in mAh
SOC = SOC_init;                                     %Initialize SOC variable
T_cell_adiabatic = T_init;                          %Initialize Adiabatic Temperature variable
T_cell_cool = T_init;                               %Initialize cooled temperature variable
T_between12 = T_init;                               %Initialize temperature between segment 1 and 2
T_cell_2 = T_init;                                  %Initialize segment 2 temperature
T_between23 = T_init;                               %Initialize temperature between segment 2 and 3
T_cell_3 = T_init;                                  %Initialize segment 3 temperature
T_between34 = T_init;                               %Initialize temperature between segment 3 and 4
T_cell_4 = T_init;                                  %Initialize segment 4 temperature
T_between45 = T_init;                               %Initialize temperature between segment 4 and 5
T_cell_5 = T_init;                                  %Initialize segment 5 temperature
[value, idx] = min(abs(SOCOCV(:,1)-SOC/100));       %Retrieve the index of the closest SOC-OCV point
Cell_OCV = SOCOCV(idx,2);                           %Retrieve the cell open circuit voltage
Pack_OCV = Scount*Cell_OCV;                         %Calculate the initial pack open circuit voltage

endurance_results = zeros(length(endurance_data),12);

%% Endurance for loop

for t=1:length(endurance_data)
    %Calculating pack resistance at that point in time

    R_pack = R_cell * Scount/Pcount + R_busbars;        %Total pack internal resistance in Ohms
    
    %Calculating the pack current and cell voltage under load
    I_pack = (Pack_OCV - sqrt(Pack_OCV^2 - 4000 * R_pack * endurance_data(t,3)))/(2*R_pack);
    V_cell = Cell_OCV - I_pack/Pcount * R_cell;
    if V_cell < 2.5
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

    %Calculate the temperatures for the next 4 segments
    %Segment 2
    T_between12 = T_amb + Qcool_cell/(mdot_air*Cp_air);
    Qcool_cell_2 = htc * A_ht * (T_cell_2 - T_between12);
    Qaccum_cell_2 = (Qgen_cell - Qcool_cell_2) * 0.05;
    T_cell_2 = T_cell_2 + Qaccum_cell_2/(Cp_batt*m_batt);
    %Segment 3
    T_between23 = T_between12 + Qcool_cell_2/(mdot_air*Cp_air);
    Qcool_cell_3 = htc * A_ht * (T_cell_3 - T_between23);
    Qaccum_cell_3 = (Qgen_cell - Qcool_cell_3) * 0.05;
    T_cell_3 = T_cell_3 + Qaccum_cell_3/(Cp_batt*m_batt);
    %Segment 4
    T_between34 = T_between23 + Qcool_cell_3/(mdot_air*Cp_air);
    Qcool_cell_4 = htc * A_ht * (T_cell_4 - T_between34);
    Qaccum_cell_4 = (Qgen_cell - Qcool_cell_4) * 0.05;
    T_cell_4 = T_cell_4 + Qaccum_cell_4/(Cp_batt*m_batt);
    %Segment 5
    T_between45 = T_between34 + Qcool_cell_4/(mdot_air*Cp_air);
    Qcool_cell_5 = htc * A_ht * (T_cell_5 - T_between45);
    Qaccum_cell_5 = (Qgen_cell - Qcool_cell_5) * 0.05;
    T_cell_5 = T_cell_5 + Qaccum_cell_5/(Cp_batt*m_batt);




    %This section of the loop records the results in a table
    endurance_results(t,1) = endurance_data(t);
    endurance_results(t,2) = V_cell;
    endurance_results(t,3) = I_pack;
    endurance_results(t,4) = Qgen_cell;
    endurance_results(t,5) = SOC;
    endurance_results(t,6) = Cell_OCV;
    endurance_results(t,7) = T_cell_adiabatic;
    endurance_results(t,8) = T_cell_cool;
    endurance_results(t,9) = T_cell_2;
    endurance_results(t,10) = T_cell_3;
    endurance_results(t,11) = T_cell_4;
    endurance_results(t,12) = T_cell_5;

    %This section of the loop subtracts the amount of SOC used and
    %determines the new cell OCV and DCIR
    Q_batt = Q_batt - 0.05*I_pack/(3.6*Pcount);
    SOC = Q_batt/30;
    [value, idx] = min(abs(SOCOCV(:,1)-SOC/100));
    Cell_OCV = SOCOCV(idx,2);
    Pack_OCV = Scount*Cell_OCV;
    [closest_SOC,closest_SOC_id] = min(abs(SOC-DCIR_LUT(:,1)));
    [closest_temp,closest_temp_id] = min(abs(T_cell_3-DCIR_LUT(1,:)));
    R_cell = DCIR_LUT(closest_SOC_id,closest_temp_id) / 1000;
end

%% Plotting Section

plot_titles =  ["Cell Voltage (V)" "Pack Current (A)" "Cell Heat Generation (W)" "State of Charge (%)" "Cell OCV (V)" "Adiabatic Cell Temperature (deg C)" "Cooled Cell Temperature (deg C)"];
time = endurance_results(:,1);

for i=1:7
    current_figure = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
    data = endurance_results(:,i+1);
    plot(time, data);
    title(plot_titles(i) + " in Michigan Endurance")
    xlabel("Time (seconds)")
    ylabel(plot_titles(i))
    saveas(current_figure, "Plots/Michigan Endurance/DCIR LUT/" + plot_titles(i) + " Plot.png")
end

current_figure = figure('visible','off','Units','centimeters','Position',[0 0 20 15]);
hold on
plot(endurance_results(:,1),endurance_results(:,8))
plot(endurance_results(:,1),endurance_results(:,9),'y')
plot(endurance_results(:,1),endurance_results(:,10),'b')
plot(endurance_results(:,1),endurance_results(:,11),'g')
plot(endurance_results(:,1),endurance_results(:,12),'r')
xlabel("Time (seconds)")
ylabel(plot_titles(i))
saveas(current_figure,"Plots/Michigan Endurance/DCIR LUT/Segment Temperature Differential.png")