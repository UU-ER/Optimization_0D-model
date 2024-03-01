
%% write new values in data set and run the 0D model 
function [Ex_spec_CH4_all, WC_CH4_all, purity_CH4_all, recovery_CH4_all, Ex_spec_CH4_Evac, WC_CH4_Evac, purity_CH4_Evac, recovery_CH4_Evac, Ex_spec_CO2_all, WC_CO2_all, purity_CO2_all, recovery_CO2_all, Ex_spec_CO2_Evac, WC_CO2_Evac, purity_CO2_Evac, recovery_CO2_Evac, succ] = run_0D(X,c,m,file_name,file_optimization,ppm_CO2)

    %% GET ISOTHERM PARAMETERS
%     file_optimization = 'data_optim_2ppm_100C_18102023.mat';
    load(file_optimization,'data_sorbent_optim')
    cases = string(fieldnames(data_sorbent_optim));
    current_case = string(cases(c));
    names_material = string(fieldnames(data_sorbent_optim.(current_case)));
    current_material = string(names_material(m));
    
        material_name = current_material;
        % Xi_c, dH_c, alpha_c, Xi_p, dH_p, alpha_p, ns0_c, b0_c, t0_c, ns0_p, b0_p, t0_p
        p_CH4 = data_sorbent_optim.(current_case).(current_material).p;
        p_N2=data_sorbent_optim.(current_case).(current_material).p_N2;
        p_CO2=data_sorbent_optim.(current_case).(current_material).p_CO2;

        isothermModel_CH4 = data_sorbent_optim.(current_case).(current_material).model;
        isothermModel_N2=data_sorbent_optim.(current_case).(current_material).model_N2;
        isothermModel_CO2=data_sorbent_optim.(current_case).(current_material).model_CO2;

        rhoMat = data_sorbent_optim.(current_case).(current_material).rhoMat;
        rho = data_sorbent_optim.(current_case).(current_material).rho;
        cp = data_sorbent_optim.(current_case).(current_material).cp; % J/Kg/K
        
        T0 = data_sorbent_optim.(current_case).(current_material).T0;
        ppm_CH4 = data_sorbent_optim.(current_case).(current_material).ppm_CH4;
             
        sorbent_data_screen_CH4_N2_CO2;
        
        data.feed.yCH4 = ppm_CH4*1e-6; %ppm
        data.feed.yCO2 = ppm_CO2*1e-6; 
        data.feed.yN2 = 1-data.feed.yCH4-data.feed.yCO2;
        
	% chose sorbent
        m = 1;
        data.currentSorbent = m; % % 1: ex, 2: CW, 3: MIL, 4: s-shaped, 5: Ex-MCF, 6: Ex-Lew, 7: MIL-MCF, 8: MIL-Lew, 9: Lew-Lew

    %% INPUT DATA
        data.process.Tads = X(:,3); % adsorption temperature
        data.process.Tamb= 293; % ambient temperature 
        data.process.pamb = 1.001*1e5*1e-6; % total pressure, bar->Pa->MPa

    %% write decesion variables
        data.process.Tdes = X(:,1); % desoprtion temperature, K
        data.process.pvac = X(:,2); % vacuum pressure, kPa
        X(:,1)
        X(:,2)
        data.process.airVelocity = (8.9*1e-06)./(pi*(0.005/2)^2); % air velocity used in 1D-model (m/s)
        data.process.Vfeed = 8.9*1e-06; % m3/s
        
    %% run model
        [Ex_spec_CH4_all, WC_CH4_all, purity_CH4_all, recovery_CH4_all, Ex_spec_CH4_Evac, WC_CH4_Evac, purity_CH4_Evac, recovery_CH4_Evac, Ex_spec_CO2_all, WC_CO2_all, purity_CO2_all, recovery_CO2_all, Ex_spec_CO2_Evac, WC_CO2_Evac, purity_CO2_Evac, recovery_CO2_Evac] = run_0D_model(data);
        
        res = [Ex_spec_CH4_all, WC_CH4_all, purity_CH4_all, recovery_CH4_all, Ex_spec_CH4_Evac, WC_CH4_Evac, purity_CH4_Evac, recovery_CH4_Evac, Ex_spec_CO2_all, WC_CO2_all, purity_CO2_all, recovery_CO2_all, Ex_spec_CO2_Evac, WC_CO2_Evac, purity_CO2_Evac, recovery_CO2_Evac];
        
        fid = fopen(file_name);
        file_fa=file_name;
        resultsAll = [X(:,1), X(:,2), X(:,3), Ex_spec_CH4_all', WC_CH4_all', purity_CH4_all', recovery_CH4_all', Ex_spec_CH4_Evac', WC_CH4_Evac', purity_CH4_Evac', recovery_CH4_Evac', Ex_spec_CO2_all', WC_CO2_all', purity_CO2_all', recovery_CO2_all', Ex_spec_CO2_Evac', WC_CO2_Evac', purity_CO2_Evac', recovery_CO2_Evac'];
        dlmwrite(file_fa,resultsAll,'newline','pc','-append','precision',10); 
        fclose(fid);

    if isempty(res)
        succ = -1;
        Ex_spec_CH4_all = 0;
        WC_CH4_all = 0;
        purity_CH4_all = 0;
        recovery_CH4_all = 0;
        Ex_spec_CH4_Evac = 0;
        WC_CH4_Evac = 0;
        purity_CH4_Evac = 0;
        recovery_CH4_Evac = 0;
        Ex_spec_CO2_all = 0;
        WC_CO2_all  = 0;
        purity_CO2_all = 0;
        recovery_CO2_all = 0;
        Ex_spec_CO2_Evac = 0;
        WC_CO2_Evac = 0;
        purity_CO2_Evac = 0;
        recovery_CO2_Evac = 0;
      
    else
        succ = 1;
    end
%     if recovery_CH4>1
%         succ=-1;
%     end
end



    