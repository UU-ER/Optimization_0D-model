
%% 0D model for TVSA process
function [Ex_spec_CH4_all, WC_CH4_all, purity_CH4_all, recovery_CH4_all, Ex_spec_CH4_Evac, WC_CH4_Evac, purity_CH4_Evac, recovery_CH4_Evac, Ex_spec_CO2_all, WC_CO2_all, purity_CO2_all, recovery_CO2_all, Ex_spec_CO2_Evac, WC_CO2_Evac, purity_CO2_Evac, recovery_CO2_Evac] = run_0D_model(data)
% PLOT
data.plot = false;
plotAll = false;
process_performance = false;
plot_results = true;

% GENERAL DATA
data.general.gasconstant = 8.314; % universal gas constant, J/mol K = m3*Pa/mol/K
data.general.adiabaticConstant = 1.4;
data.general.pumpEfficiency = 0.72;
data.general.airBlowerEfficiency = 0.72;
data.general.densityAir = 1.23; % kg/m3
data.general.MMCO2 = 44.01; % Molar mass CO2, g/mol
data.general.MMN2 = 28.01; % g/mol
data.general.MMCH4 = 16.04; % g/mol

% PROCESS SPECIFIC DATA
data.process.adsorbentMass = 1; % Mass of the adsorbent, kg
data.process.voidFraction = (data.sorbent(data.currentSorbent).MaterialDensity- data.sorbent(data.currentSorbent).Density)/data.sorbent(data.currentSorbent).MaterialDensity; % column void fraction
data.process.heatResistance = 16.8; % global heat transfer resistance, (J/(m2 K))
data.process.heightFrame = 0.05; % height of sorbent containing frame, (m)
data.process.area = data.process.adsorbentMass/...
    data.sorbent(data.currentSorbent).Density/data.process.heightFrame/2; % column cross section
data.process.noSteps = 100; % number of steps per unit
data.process.Vol = data.process.adsorbentMass/(data.sorbent(data.currentSorbent).Density); % volume column (times void = volume void space)

% Saturation degree and starting concentration BD
full_saturation = true;
for m = 1:size(data.process.Tdes,1)
    % particle density, V_feed, Tdes, pvac
    if full_saturation
        yCO2_sat(m) = data.feed.yCO2;
    else
        load('netSat','netNN_sat');
        Vfeed = data.process.Vfeed;
        Tdes = data.process.Tdes;
        pvac = data.process.pvac ;
        rho = data.sorbent(data.currentSorbent).Density;
        sat(m) = netNN_sat([rho;Vfeed(m);Tdes(m);pvac(m)*10]);
        if sat(m) > 1
            sat(m)=1;
        end
        data.degreeSaturation = sat(m);

        % calculate CO2 concentration at saturation level
        yCO2_sat(m) = calcSat(data,sat(m));
    end
    startingCondBD(m,:) = [yCO2_sat(m),(1-yCO2_sat(m)-data.feed.yCH4),data.feed.yCH4];
end
data.startingCondBD = startingCondBD;

%% PRESSURE AND TEMPEARTURE VECTOR
[pressureVector, timeBD, TimeBDProfilePlot] = pressureProfile(data,data.process.pvac);
pressureVector=linspace(data.process.pamb,data.process.pvac,100);
data.pressureVector = pressureVector; % MPa
data.process.noStepsBD = round(timeBD+2);
data.process.TimeBDProfilePlot = TimeBDProfilePlot; % profile vector reduced to max. amount of steps for plotting
t_BD = data.process.noStepsBD;

[THeatProfileStep, timeHeating, TimeHeatProfilePlot] = temperatureProfile(data,data.process.Tdes);
data.THeatProfileStep = THeatProfileStep; %
data.THeatProfileStep=linspace(data.process.Tads,data.process.Tdes,100);
data.process.noStepsHeating = round(timeHeating+1);
data.process.TimeHeatProfilePlot = TimeHeatProfilePlot; % profile vector reduced to max. amount of steps for plotting
t_Heating = data.process.noStepsHeating;

data.process.coolingTime = 350;
[TCoolProfileStep, TimeCoolProfilePlot, pressureVectorCool] = coolingProfile_1D(data,data.process.Tdes,data.process.pvac); % with profile there were humidity problems during cooling step: coolingProfile_1D(data,data.process.Tdes);
TCoolProfileStep=linspace(data.process.Tdes,data.process.Tads,100);
data.TCoolProfileStepPlot = TCoolProfileStep;
data.TCoolProfileStep = TCoolProfileStep;
pressureVectorCool=linspace(data.process.pvac,data.process.pamb,101);
data.PressureCoolProfileStep = pressureVectorCool;
data.process.TimeCoolProfilePlot = TimeCoolProfilePlot; % profile vector reduced to max. amount of steps for plotting

%% ADSORPTION FEED
data.general.MMFeed = data.feed.yCO2*data.general.MMCO2+data.feed.yN2*data.general.MMN2+data.feed.yCH4*data.general.MMCH4; % g/mol
m_feed_step = 100*ones(size(data.process.Tdes,1),1); % mass feed per step, g
data.process.adsFeedPerSecond = m_feed_step./data.general.MMFeed; % feed per step, mol
data.process.FeedAdsorptionStep = m_feed_step./data.general.MMFeed; % feed per step, mol
data.process.noStepMaxAdsorption = 100000; % maximum steps

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN THE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%------------------------ blowdown and evacuation ------------------------%
% blowdown and evacuation step
% output: pressure, molefraction CO2, solid loadings of CO2 and H2O,
% total moles desorbed, moles of CO2 and H2O desorbed and the energy consumption
% Matrix structure: [Pressure MoleFracCO2 SolidLoadingCO2 SolidLoadingH2O
% MolesDesorbedTotal MolesDesorbedCO2 MolesDesorbedH2O EnergyConsumption]
try
outputBlowEvac = simulateBlowdownEvacuation(data);

%-------------------------------- cooling --------------------------------%
outputCooling = simulateCooling(data,outputBlowEvac);

if outputCooling.yCO2(end)<data.feed.yCO2 && outputCooling.yCH4(end)<data.feed.yCH4
    %------------------------------ adsorption -------------------------------%
    outputAdsorption = simulateAdsorption(data,outputCooling);


    %% results
    for m = 1:length(data.process.Tdes)
        mfeed = sum(outputAdsorption(m).Nfeed_sum)*data.general.MMFeed/1000; % kg
        t_ads(m) = mfeed/data.process.airVelocity(m)/data.process.area/data.general.densityAir; % s

        mfeed1 = sum(outputAdsorption(m).Nfeed_sum(1))*data.general.MMFeed/1000; % kg
        t_ads1(m) = mfeed1/data.process.airVelocity(m)/data.process.area/data.general.densityAir; % s

        if plotAll
            plot_all = 'all_steps';%'all_steps';
        	%     separate_plot_file;
            separate_plot_file_new2;
        end
        Tdes = data.process.Tdes;
        Tamb = data.process.Tamb;
        Tads = data.process.Tads;
        % Heating
        QHeating(m) = outputBlowEvac(m).Q_sum(end); % QHeating, kJ thermal heating
        ExHeating(m) = QHeating*(1-Tamb/Tdes); % Heating Exergy, KJ
        % Cooling
        QCooling_Solid(m) = abs(outputCooling(m).Q_sum(end)); %QCooling, KJ thermal Energy 
        ExCooling_Solid(m) = QCooling_Solid*(Tamb/Tads-1); % Exergy of the Cooling the Solid, KJ
        ExCooling_Air(m) = outputAdsorption(m).E_Refrig+outputCooling(m).E_Refrig; % Exergy of the Cooling the inlet air, KJ
        ExCooling = ExCooling_Air + ExCooling_Solid;
        %Vacuum pump
        E_Vacuum(m) = (outputBlowEvac(m).E_total(end)); % kJ vacuum pump
        %Total Ex: Heating+Cooling+Vacuum
        ExTotal(m) = ExHeating+ExCooling+E_Vacuum(m); % Total Exergy, KJ


        N_CO2_product_Evac(m) = outputBlowEvac(m).Nout_CO2(end)-outputBlowEvac(m).Nout_CO2(data.process.noSteps(m)); % mol only Evacuation, BD not included
        N_CO2_all(m) = outputBlowEvac(m).Nout_CO2(end); % mol BD+Evacution
        N_CH4_product_Evac(m) = outputBlowEvac(m).Nout_CH4(end)-outputBlowEvac(m).Nout_CH4(data.process.noSteps(m)); % mol only Evacuation, BD not included
        N_CH4_all(m) = outputBlowEvac(m).Nout_CH4(end); % mol BD+Evacution
        N_CH4_in_tot=data.feed.yCH4*(sum(outputCooling.Nin_sum)+sum(outputAdsorption.Nfeed_sum));
        N_CO2_in_tot=data.feed.yCO2*(sum(outputCooling.Nin_sum)+sum(outputAdsorption.Nfeed_sum));
        N_N2_in_tot=data.feed.yN2*(sum(outputCooling.Nin_sum)+sum(outputAdsorption.Nfeed_sum));
        %% performance
        % electricity
        E_spec_CO2(m) = E_Vacuum(m)/(N_CO2_product_Evac(m)*data.general.MMCO2/1000)/1000; % MJ/kg
        E_spec_CH4(m) = E_Vacuum(m)/(N_CH4_product_Evac(m)*data.general.MMCH4/1000)/1000; % MJ/kg
        E_spec_eq(m) = E_Vacuum(m)/(N_CO2_product_Evac(m)*data.general.MMCO2/1000+28*N_CH4_product_Evac(m)*data.general.MMCH4/1000)/1000; % MJ/kgCO2eq GHG(CH4)=28 over 100 years
        % heat
        QTotal=QHeating;
        Q_spec_CO2(m) = QTotal(m)/(N_CO2_product_Evac(m)*data.general.MMCO2/1000)/1000; % MJ/kg
        E_spec_CH4(m) = QTotal(m)/(N_CH4_product_Evac(m)*data.general.MMCH4/1000)/1000; % MJ/kg
        Q_spec_eq(m) = QTotal(m)/(N_CO2_product_Evac(m)*data.general.MMCO2/1000+28*N_CH4_product_Evac(m)*data.general.MMCH4/1000)/1000; % MJ/kgCO2eq GHG(CH4)=28 over 100 years
        Q_spec_eq_all(m) = QTotal(m)/(N_CO2_all(m)*data.general.MMCO2/1000+28*N_CH4_all(m)*data.general.MMCH4/1000)/1000; % MJ/kgCO2eq GHG(CH4)=28 over 100 years

        % Exergy
        Ex_spec_CH4_all(m) = ExTotal(m)/(N_CH4_all(m)*data.general.MMCH4/1000)/1000; % MJ/KgCH4 produced during BD+Evac
        Ex_spec_CH4_Evac(m) = ExTotal(m)/(N_CH4_product_Evac(m)*data.general.MMCH4/1000)/1000; % MJ/KgCH4 produced during Evac
        Ex_spec_eq_Evac(m) = ExTotal(m)/(N_CO2_product_Evac(m)*data.general.MMCO2/1000+28*N_CH4_product_Evac(m)*data.general.MMCH4/1000)/1000; % MJ/kgCO2eq GHG(CH4)=28 over 100 years
        Ex_spec_eq_all(m) = ExTotal(m)/(N_CO2_all(m)*data.general.MMCO2/1000+28*N_CH4_all(m)*data.general.MMCH4/1000)/1000; % MJ/kgCO2eq GHG(CH4)=28 over 100 years
        Ex_spec_CO2_all(m) = ExTotal(m)/(N_CO2_all(m)*data.general.MMCO2/1000)/1000; % MJ/KgCO2 produced during BD+Evac
        Ex_spec_CO2_Evac(m) = ExTotal(m)/(N_CO2_product_Evac(m)*data.general.MMCO2/1000)/1000; % MJ/KgCO2 produced during Evac

        % productivity
%         Pr(m) = N_CO2_product(m)*data.general.MMCO2/1000/(data.process.adsorbentMass/data.sorbent(data.currentSorbent).Density)/((data.process.coolingTime+t_ads(m)+data.process.noStepsBD(m)+data.process.noStepsHeating(m))/3600); % kg_CO2/m3_sorbent/h
        %Purity
        purity_CO2_Evac(m) = N_CO2_product_Evac(m)/(outputBlowEvac(m).Nout_sum(end)-outputBlowEvac(m).Nout_sum(data.process.noSteps(m))); % Nout_save_CO2_evac/(Nout_save_total_evac-Nout_save_H2O_evac)
        purity_CO2_all(m) = N_CO2_all(m)/(outputBlowEvac(m).Nout_sum(end));
        purity_CH4_Evac(m) = N_CH4_product_Evac(m)/(outputBlowEvac(m).Nout_sum(end)-outputBlowEvac(m).Nout_sum(data.process.noSteps(m))); % product only from heating step not BD
        purity_CH4_all(m) = N_CH4_all(m)/(outputBlowEvac(m).Nout_sum(end));
        %recovery
        recovery_CO2_Evac(m) = N_CO2_product_Evac(m)/(data.feed.yCO2*(sum(outputAdsorption(m).Nfeed_sum) + sum(outputCooling(m).Nin_sum)));
        recovery_CO2_all(m) = N_CO2_all(m)/(data.feed.yCO2*(sum(outputAdsorption(m).Nfeed_sum) + sum(outputCooling(m).Nin_sum)));
        recovery_CH4_Evac(m) = N_CH4_product_Evac(m)/(data.feed.yCH4*(sum(outputAdsorption(m).Nfeed_sum) + sum(outputCooling(m).Nin_sum)));
        recovery_CH4_all(m) = N_CH4_all(m)/(data.feed.yCH4*(sum(outputAdsorption(m).Nfeed_sum) + sum(outputCooling(m).Nin_sum)));
        %balance:
        Nin_total = sum(outputCooling(m).Nin_sum)+sum(outputAdsorption.Nfeed_sum);
        Nout_total = outputBlowEvac.Nout_sum(end)+sum(outputAdsorption.Nwaste_sum);

        % cyclic working capacity
        WC_CO2_Evac(m) = N_CO2_product_Evac(m)/data.process.adsorbentMass; % mol_CO2_recovered/kg_sorbent in Evacution
        WC_CO2_all(m) = N_CO2_all(m)/data.process.adsorbentMass; % % mol_CO2_recovered/kg_sorbent in Evacution+ BD
        WC_CH4_Evac(m) = N_CH4_product_Evac(m)/data.process.adsorbentMass;
        WC_CH4_all(m) = N_CH4_all(m)/data.process.adsorbentMass;
        sprintf(' Specific Exergy all: %0.2f MJ/kgCH4,\n Working Capacity all: %0.2f KgCH4/kgSorbent,\n purity CH4 all: %0.2f ,\n purity: %0.2f,\n capture rate: %0.2f,\n capture rate_all: %0.2f, \n specific working capacity: %0.2f molCH4/kgSorbent.',Ex_spec_CH4_all, WC_CH4_all, purity_CH4_all, recovery_CH4_all, Ex_spec_CH4_Evac, WC_CH4_Evac, purity_CH4_Evac, recovery_CH4_Evac)

        % compare process performance with 1D
        if process_performance
            plot_process_performance;
        end

        if plot_results
            % 	   figure(2)
            % 	   plot(Pr,(Q_spec),'bo')
        end
    end
else
   Ex_spec_CH4_all=[]; WC_CH4_all=[]; purity_CH4_all=[]; recovery_CH4_all=[]; Ex_spec_CH4_Evac=[]; WC_CH4_Evac=[]; purity_CH4_Evac=[]; recovery_CH4_Evac=[]; Ex_spec_CO2_all=[]; WC_CO2_all=[]; purity_CO2_all=[]; recovery_CO2_all=[]; Ex_spec_CO2_Evac=[]; WC_CO2_Evac=[]; purity_CO2_Evac=[]; recovery_CO2_Evac=[];
end
catch
   Ex_spec_CH4_all=[]; WC_CH4_all=[]; purity_CH4_all=[]; recovery_CH4_all=[]; Ex_spec_CH4_Evac=[]; WC_CH4_Evac=[]; purity_CH4_Evac=[]; recovery_CH4_Evac=[]; Ex_spec_CO2_all=[]; WC_CO2_all=[]; purity_CO2_all=[]; recovery_CO2_all=[]; Ex_spec_CO2_Evac=[]; WC_CO2_Evac=[]; purity_CO2_Evac=[]; recovery_CO2_Evac=[];

end
end

