%========================================================
%%% Sorbent isotherm parameters+ density material+ Cp solid + bed density
%%% if the material information for a specific gas does not exist, the script assumes it
%%% "constant" and input a fixed predefined value for that gas
%========================================================
i = 1;
%% CH4
switch isothermModel_CH4
    case 'toth_cp'
        % Xi_c,dH_c[J/mol],alpha_c,Xi_p,dH_p[J/mol],alpha_p,ns0_c[mol/kg],b0_c[1/MPa],t0_c,ns0_p[mol/kg],b0_p[1/MPa],t0_p
        data.sorbent(i).name = material_name;
        data.sorbent(i).CH4Isotherm.model = isothermModel_CH4;
        data.sorbent(i).CH4Isotherm.dH_ads_CH4 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CH4Isotherm.T0 = T0; % T0 [K]
        data.sorbent(i).CH4Isotherm.param(1) = p_CH4(1); % Xi_c
        data.sorbent(i).CH4Isotherm.param(2) = p_CH4(2)*1000; % dH_c[J/mol]
        data.sorbent(i).CH4Isotherm.param(3) = p_CH4(3); % alpha_c
        data.sorbent(i).CH4Isotherm.param(4) = p_CH4(4); % Xi_p
        data.sorbent(i).CH4Isotherm.param(5) = p_CH4(5)*1000; % dH_p[J/mol]
        data.sorbent(i).CH4Isotherm.param(6) = p_CH4(6); % alpha_p
        data.sorbent(i).CH4Isotherm.param(7) = p_CH4(7); % ns0_c[mol/kg]
        data.sorbent(i).CH4Isotherm.param(8) = p_CH4(8)*(1e+6); % b0_c[1/MPa]
        data.sorbent(i).CH4Isotherm.param(9) = p_CH4(9); % t0_c
        data.sorbent(i).CH4Isotherm.param(10) = p_CH4(10); % ns0_p[mol/kg]
        data.sorbent(i).CH4Isotherm.param(11) = p_CH4(11)*(1e+6); % b0_p[1/MPa]
        data.sorbent(i).CH4Isotherm.param(12) = p_CH4(12); % t0_p

    case 's_shaped'
        % q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(mol/kg),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(K), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
        data.sorbent(i).name = material_name;
        data.sorbent(i).CH4Isotherm.model = isothermModel_CH4;
        data.sorbent(i).CH4Isotherm.dH_ads_CH4 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CH4Isotherm.T0 = T0; % T0 [K]
        data.sorbent(i).CH4Isotherm.param(1) = p_CH4(1); % q_L0, [mol/kg], a1
        data.sorbent(i).CH4Isotherm.param(2) = p_CH4(2)*(1e+6); % b_L0, [1/MPa], b0
        data.sorbent(i).CH4Isotherm.param(3) = p_CH4(2)*1000; % dU_L, [J/mol], b1
        data.sorbent(i).CH4Isotherm.param(4) = p_CH4(4); % q_U0, [mol/kg], c1
        data.sorbent(i).CH4Isotherm.param(5) = p_CH4(5)*(1e+6); % b_U0, [1/MPa] d0
        data.sorbent(i).CH4Isotherm.param(6) = p_CH4(6)*1000; % dU_U, [J/mol] d1
        data.sorbent(i).CH4Isotherm.param(7) = p_CH4(7)*(1e+6); % b_H0, [mol/kg/MPa]  bb0
        data.sorbent(i).CH4Isotherm.param(8) = p_CH4(8)*1000; % dU_H, [J/mol] bb1
        data.sorbent(i).CH4Isotherm.param(9) = p_CH4(9); % xi1, [-]
        data.sorbent(i).CH4Isotherm.param(10) = p_CH4(10); % xi_2, [K] xi2
        data.sorbent(i).CH4Isotherm.param(11) = p_CH4(11)/(1e+6); % p_step0, [MPa] ps0
        data.sorbent(i).CH4Isotherm.param(12) = p_CH4(12)*1000; % dH_step, [J/mol] Hst
        data.sorbent(i).CH4Isotherm.param(13) = p_CH4(13); % gam, [-]

    case 'DSL'
        % n1 (mol/kg), b0 (m3/mol), Hb (J/mol), n2 (mol/kg), d0 (m3/mol), Hd (J/mol)
        % e.g. dataModel from Farooq
        data.sorbent(i).name = material_name;
        data.sorbent(i).CH4Isotherm.model = isothermModel_CH4;
        data.sorbent(i).CH4Isotherm.dH_ads_CH4 = 0; % Heat of adsorption [J/mol]\
        data.sorbent(i).CH4Isotherm.param(1) = p_CH4(1); % n1 (mol/kg)
        data.sorbent(i).CH4Isotherm.param(2) = p_CH4(2); % b0 (m3/Mmol)
        data.sorbent(i).CH4Isotherm.param(3) = p_CH4(3); % Hb (J/mol)
        data.sorbent(i).CH4Isotherm.param(4) = p_CH4(4); % n2 (mol/kg)
        data.sorbent(i).CH4Isotherm.param(5) = p_CH4(5); % d0 (m3/Mmol)
        data.sorbent(i).CH4Isotherm.param(6) = p_CH4(6); % Hd (J/mol)
        % sorbent specific dataModel    --> change!!!!!!!!!!!!!!
        data.sorbent(i).MaterialDensity = dataset.(name_analysis).Methane.(material_name).(doi_name).rhoMat;
        data.sorbent(i).Density = dataset.(name_analysis).Methane.(material_name).(doi_name).rhoMat*0.35; % density sorbent, kg/m3 (particle)
        data.sorbent(i).cp = 1514.2; % could also be added from Farooq data

    case 'langm_freund'
        % 1: n0 [mol/kg], 2: Xi [-], 3: t0 [-], 4: alpha [-], 5:b0 [1/Pa], 6:Q [kJ/mol]
        data.sorbent(i).name = material_name;
        data.sorbent(i).CH4Isotherm.T0 = T0; % T0 [K]
        data.sorbent(i).CH4Isotherm.model = 'langfr';
        data.sorbent(i).CH4Isotherm.dH_ads_CH4 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CH4Isotherm.param(1) = p_CH4(1); % ns0 [mol/kg]
        data.sorbent(i).CH4Isotherm.param(2) = p_CH4(2); % Xi
        data.sorbent(i).CH4Isotherm.param(3) = p_CH4(3); % t0
        data.sorbent(i).CH4Isotherm.param(4) = p_CH4(4); % alpha
        data.sorbent(i).CH4Isotherm.param(5) = p_CH4(5)*(1e+6); % b0 [1/MPa]
        data.sorbent(i).CH4Isotherm.param(6) = p_CH4(6)*1000; % dH [J/mol]

    case 'toth'
        %  Xi_p, dH_p(kJ/mol), alpha_p, ns0_p(mol/kg), b0_p (1/Pa), t0_p
        data.sorbent(i).name = material_name;
        data.sorbent(i).CH4Isotherm.model = 'toth';
        data.sorbent(i).CH4Isotherm.dH_ads_CH4 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CH4Isotherm.T0 = T0; % [K]
        data.sorbent(i).CH4Isotherm.param(1) = p_CH4(1); % Xi
        data.sorbent(i).CH4Isotherm.param(2) = p_CH4(2)*1000; % dH[J/mol]
        data.sorbent(i).CH4Isotherm.param(3) = p_CH4(3); % alpha
        data.sorbent(i).CH4Isotherm.param(4) = p_CH4(4); % ns0[mol/kg]
        data.sorbent(i).CH4Isotherm.param(5) = p_CH4(5)*(1e+6); % b0[1/MPa]
        data.sorbent(i).CH4Isotherm.param(6) = p_CH4(6); % t0

end



%% N2

switch isothermModel_N2
    case 'toth_cp'
        % Xi_c,dH_c[KJ/mol],alpha_c,Xi_p,dH_p[KJ/mol],alpha_p,ns0_c[mol/kg],b0_c[1/Pa],t0_c,ns0_p[mol/kg],b0_p[1/MPa],t0_p
        data.sorbent(i).name = material_name;
        data.sorbent(i).N2Isotherm.model = isothermModel_N2;
        data.sorbent(i).N2Isotherm.dH_ads_N2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).N2Isotherm.T0 = T0; % [K]
        data.sorbent(i).N2Isotherm.param(1) = p_N2(1); % Xi_c
        data.sorbent(i).N2Isotherm.param(2) = p_N2(2)*1000; % dH_c[J/mol]
        data.sorbent(i).N2Isotherm.param(3) = p_N2(3); % alpha_c
        data.sorbent(i).N2Isotherm.param(4) = p_N2(4); % Xi_p
        data.sorbent(i).N2Isotherm.param(5) = p_N2(5)*1000; % dH_p[J/mol]
        data.sorbent(i).N2Isotherm.param(6) = p_N2(6); % alpha_p
        data.sorbent(i).N2Isotherm.param(7) = p_N2(7); % ns0_c[mol/kg]
        data.sorbent(i).N2Isotherm.param(8) = p_N2(8)*(1e+6); % b0_c[1/MPa]
        data.sorbent(i).N2Isotherm.param(9) = p_N2(9); % t0_c
        data.sorbent(i).N2Isotherm.param(10) = p_N2(10); % ns0_p[mol/kg]
        data.sorbent(i).N2Isotherm.param(11) = p_N2(11)*(1e+6); % b0_p[1/MPa]
        data.sorbent(i).N2Isotherm.param(12) = p_N2(12); % t0_p

    case 's_shaped'
        % N2: as an example exemplary isotherm
        % q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(mol/kg),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(-), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
        data.sorbent(i).name = material_name;
        data.sorbent(i).N2Isotherm.model = isothermModel_N2;
        data.sorbent(i).N2Isotherm.dH_ads_N2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).N2Isotherm.T0 = T0; % K
        data.sorbent(i).N2Isotherm.param(1) = p_N2(1); % q_L0, [mol/kg], a1
        data.sorbent(i).N2Isotherm.param(2) = p_N2(2)*(1e+6); % b_L0, [1/MPa], b0
        data.sorbent(i).N2Isotherm.param(3) = p_N2(2)*1000; % dU_L, [J/mol], b1
        data.sorbent(i).N2Isotherm.param(4) = p_N2(4); % q_U0, [mol/kg], c1
        data.sorbent(i).N2Isotherm.param(5) = p_N2(5)*(1e+6); % b_U0, [1/MPa] d0
        data.sorbent(i).N2Isotherm.param(6) = p_N2(6)*1000; % dU_U, [J/mol] d1
        data.sorbent(i).N2Isotherm.param(7) = p_N2(7)*(1e+6); % b_H0, [mol/kg/MPa]  bb0
        data.sorbent(i).N2Isotherm.param(8) = p_N2(8)*1000; % dU_H, [J/mol] bb1
        data.sorbent(i).N2Isotherm.param(9) = p_N2(9); % xi1
        data.sorbent(i).N2Isotherm.param(10) = p_N2(10); % xi_2, [K] xi2
        data.sorbent(i).N2Isotherm.param(11) = p_N2(11)/(1e+6); % p_step0, MPa ps0
        data.sorbent(i).N2Isotherm.param(12) = p_N2(12)*1000; % dH_step, [J/mol] Hst
        data.sorbent(i).N2Isotherm.param(13) = p_N2(13); % gam

    case 'DSL'
        % n1 (mol/kg), b0 (m3/mol), Hb (J/mol), n2 (mol/kg), d0 (m3/mol), Hd (J/mol)
        % e.g. dataModel from Farooq
        data.sorbent(i).name = material_name;
        data.sorbent(i).N2Isotherm.model = isothermModel_N2;
        data.sorbent(i).N2Isotherm.dH_ads_N2 = 0; % Heat of adsorption [J/mol]\
        data.sorbent(i).N2Isotherm.param(1) = p_N2(1); % n1 (mol/kg)
        data.sorbent(i).N2Isotherm.param(2) = p_N2(2); % b0 (m3/Mmol)
        data.sorbent(i).N2Isotherm.param(3) = p_N2(3); % Hb (J/mol)
        data.sorbent(i).N2Isotherm.param(4) = p_N2(4); % n2 (mol/kg)
        data.sorbent(i).N2Isotherm.param(5) = p_N2(5); % d0 (m3/Mmol)
        data.sorbent(i).N2Isotherm.param(6) = p_N2(6); % Hd (J/mol)
        % sorbent specific dataModel    --> change!!!!!!!!!!!!!!
        data.sorbent(i).MaterialDensity = dataset.(name_analysis).Methane.(material_name).(doi_name).rhoMat;
        data.sorbent(i).Density = dataset.(name_analysis).Methane.(material_name).(doi_name).rhoMat*0.35; % density sorbent, kg/m3 (particle)
        data.sorbent(i).cp = 1514.2; % could also be added from Farooq data

    case 'langm_freund'
        % 1: n0 [mol/kg], 2: Xi [-], 3: t0 [-], 4: alpha [-], 5:b0 [1/Pa], 6:Q [kJ/mol]
        data.sorbent(i).name = material_name;
        data.sorbent(i).N2Isotherm.T0 = T0; % [K]
        data.sorbent(i).N2Isotherm.model = 'langfr';
        data.sorbent(i).N2Isotherm.dH_ads_N2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).N2Isotherm.param(1) = p_N2(1); % ns0 [mol/kg]
        data.sorbent(i).N2Isotherm.param(2) = p_N2(2); % Xi
        data.sorbent(i).N2Isotherm.param(3) = p_N2(3); % t0
        data.sorbent(i).N2Isotherm.param(4) = p_N2(4); % alpha
        data.sorbent(i).N2Isotherm.param(5) = p_N2(5)*(1e+6); % b0 [1/MPa]
        data.sorbent(i).N2Isotherm.param(6) = p_N2(6)*1000; % dH [J/mol]

    case 'toth'
        %  Xi_p, dH_p(kJ/mol), alpha_p, ns0_p(mol/kg), b0_p (1/Pa), t0_p
        data.sorbent(i).name = material_name;
        data.sorbent(i).N2Isotherm.model = 'toth';
        data.sorbent(i).N2Isotherm.dH_ads_N2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).N2Isotherm.T0 = T0; % [K]
        data.sorbent(i).N2Isotherm.param(1) = p_N2(1); % Xi
        data.sorbent(i).N2Isotherm.param(2) = p_N2(2)*1000; % dH[J/mol]
        data.sorbent(i).N2Isotherm.param(3) = p_N2(3); % alpha
        data.sorbent(i).N2Isotherm.param(4) = p_N2(4); % ns0[mol/kg]
        data.sorbent(i).N2Isotherm.param(5) = p_N2(5)*(1e+6); % b0[1/MPa]
        data.sorbent(i).N2Isotherm.param(6) = p_N2(6); % t0

    case 'no_data'
        isothermModel_N2="langm_freund";
        data.sorbent(i).N2Isotherm.model = 'langfr';
        data.sorbent(i).N2Isotherm.T0 = 293; % T0 [K]
        data.sorbent(i).N2Isotherm.param(1) = 7.777; % ns0 [mol/kg]
        data.sorbent(i).N2Isotherm.param(2) = 4.401; % Xi
        data.sorbent(i).N2Isotherm.param(3) = 1; % t0
        data.sorbent(i).N2Isotherm.param(4) = 0; % alpha
        data.sorbent(i).N2Isotherm.param(5) =  0.1343; % b0 [1/Mpa]
        data.sorbent(i).N2Isotherm.param(6) = 3.127; % dH [J/mol]

end



%% CO2

switch isothermModel_CO2
    case 'toth_cp'
        % param: Xi_c,dH_c[J/mol],alpha_c,Xi_p,dH_p[J/mol],alpha_p,ns0_c[mol/kg],b0_c[1/Pa],t0_c,ns0_p[mol/kg],b0_p[1/MPa],t0_p
        data.sorbent(i).name = material_name;
        data.sorbent(i).CO2Isotherm.model = isothermModel_CO2;
        data.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CO2Isotherm.T0 = T0; % K
        data.sorbent(i).CO2Isotherm.param(1) = p_CO2(1); % Xi_c
        data.sorbent(i).CO2Isotherm.param(2) = p_CO2(2)*1000; % dH_c[J/mol]
        data.sorbent(i).CO2Isotherm.param(3) = p_CO2(3); % alpha_c
        data.sorbent(i).CO2Isotherm.param(4) = p_CO2(4); % Xi_p
        data.sorbent(i).CO2Isotherm.param(5) = p_CO2(5)*1000; % dH_p[J/mol]
        data.sorbent(i).CO2Isotherm.param(6) = p_CO2(6); % alpha_p
        data.sorbent(i).CO2Isotherm.param(7) = p_CO2(7); % ns0_c[mol/kg]
        data.sorbent(i).CO2Isotherm.param(8) = p_CO2(8)*(1e+6); % b0_c[1/MPa]
        data.sorbent(i).CO2Isotherm.param(9) = p_CO2(9); % t0_c
        data.sorbent(i).CO2Isotherm.param(10) = p_CO2(10); % ns0_p[mol/kg]
        data.sorbent(i).CO2Isotherm.param(11) = p_CO2(11)*(1e+6); % b0_p[1/MPa]
        data.sorbent(i).CO2Isotherm.param(12) = p_CO2(12); % t0_p

    case 's_shaped'
        % q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(mol/kg),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(-), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
        data.sorbent(i).name = material_name;
        data.sorbent(i).CO2Isotherm.model = isothermModel_CO2;
        data.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CO2Isotherm.T0 = T0; % [K]
        data.sorbent(i).CO2Isotherm.param(1) = p_CO2(1); % q_L0, [mol/kg], a1
        data.sorbent(i).CO2Isotherm.param(2) = p_CO2(2)*(1e+6); % b_L0, [1/MPa], b0
        data.sorbent(i).CO2Isotherm.param(3) = p_CO2(2)*1000; % dU_L, [J/mol], b1
        data.sorbent(i).CO2Isotherm.param(4) = p_CO2(4); % q_U0, [mol/kg], c1
        data.sorbent(i).CO2Isotherm.param(5) = p_CO2(5)*(1e+6); % b_U0, [1/MPa] d0
        data.sorbent(i).CO2Isotherm.param(6) = p_CO2(6)*1000; % dU_U, J/mol d1
        data.sorbent(i).CO2Isotherm.param(7) = p_CO2(7)*(1e+6); % b_H0, [mol/kg/Pa]  bb0
        data.sorbent(i).CO2Isotherm.param(8) = p_CO2(8)*1000; % dU_H, J/mol bb1
        data.sorbent(i).CO2Isotherm.param(9) = p_CO2(9); % xi1
        data.sorbent(i).CO2Isotherm.param(10) = p_CO2(10); % xi_2, K xi2
        data.sorbent(i).CO2Isotherm.param(11) = p_CO2(11)/(1e+6); % p_step0, [MPa] ps0
        data.sorbent(i).CO2Isotherm.param(12) = p_CO2(12)*1000; % dH_step, [J/mol] Hst
        data.sorbent(i).CO2Isotherm.param(13) = p_CO2(13); % gam

    case 'DSL'
        % n1 (mol/kg), b0 (m3/mol), Hb (J/mol), n2 (mol/kg), d0 (m3/mol), Hd (J/mol)
        % e.g. dataModel from Farooq
        data.sorbent(i).name = material_name;
        data.sorbent(i).CO2Isotherm.model = isothermModel_CO2;
        data.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]\
        data.sorbent(i).CO2Isotherm.param(1) = p_CO2(1); % n1 (mol/kg)
        data.sorbent(i).CO2Isotherm.param(2) = p_CO2(2); % b0 (m3/Mmol)
        data.sorbent(i).CO2Isotherm.param(3) = p_CO2(3); % Hb (J/mol)
        data.sorbent(i).CO2Isotherm.param(4) = p_CO2(4); % CO2 (mol/kg)
        data.sorbent(i).CO2Isotherm.param(5) = p_CO2(5); % d0 (m3/Mmol)
        data.sorbent(i).CO2Isotherm.param(6) = p_CO2(6); % Hd (J/mol)
        % sorbent specific dataModel    --> change!!!!!!!!!!!!!!
        data.sorbent(i).MaterialDensity = dataset.(name_analysis).Methane.(material_name).(doi_name).rhoMat;
        data.sorbent(i).Density = dataset.(name_analysis).Methane.(material_name).(doi_name).rhoMat*0.35; % density sorbent, kg/m3 (particle)
        data.sorbent(i).cp = 1514.2; % could also be added from Farooq data

    case 'langm_freund'
        % 1: n0 [mmol/g], 2: Xi [-], 3: t0 [-], 4: alpha [-], 5:b0 [1/Pa], 6:Q [kJ/mol]
        data.sorbent(i).name = material_name;
        data.sorbent(i).CO2Isotherm.T0 = T0; % [K]
        data.sorbent(i).CO2Isotherm.model = 'langfr';
        data.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CO2Isotherm.param(1) = p_CO2(1); % ns0 (mol/kg)
        data.sorbent(i).CO2Isotherm.param(2) = p_CO2(2); % Xi
        data.sorbent(i).CO2Isotherm.param(3) = p_CO2(3); % t0
        data.sorbent(i).CO2Isotherm.param(4) = p_CO2(4); % alpha
        data.sorbent(i).CO2Isotherm.param(5) = p_CO2(5)*(1e+6); % b0 (1/MPa)
        data.sorbent(i).CO2Isotherm.param(6) = p_CO2(6)*1000; % dH (J/mol)

    case 'toth'
        %  Xi_p, dH_p(kJ/mol), alpha_p, ns0_p(mol/kg), b0_p (1/Pa), t0_p
        data.sorbent(i).name = material_name;
        data.sorbent(i).CO2Isotherm.model = 'toth';
        data.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]
        data.sorbent(i).CO2Isotherm.T0 = T0; % [K]
        data.sorbent(i).CO2Isotherm.param(1) = p_CO2(1); % Xi
        data.sorbent(i).CO2Isotherm.param(2) = p_CO2(2)*1000; % dH[J/mol]
        data.sorbent(i).CO2Isotherm.param(3) = p_CO2(3); % alpha
        data.sorbent(i).CO2Isotherm.param(4) = p_CO2(4); % ns0[mol/kg]
        data.sorbent(i).CO2Isotherm.param(5) = p_CO2(5)*(1e+6); % b0[1/MPa]
        data.sorbent(i).CO2Isotherm.param(6) = p_CO2(6); % t0

    case 'no_data'
        isothermModel_CO2="langm_freund";
        data.sorbent(i).CO2Isotherm.model = 'langfr';
        data.sorbent(i).CO2Isotherm.T0 = 293; % K
        data.sorbent(i).CO2Isotherm.param(1) = 22.23; % ns0 [mol/kg]
        data.sorbent(i).CO2Isotherm.param(2) = 5.048e-14; %  Xi
        data.sorbent(i).CO2Isotherm.param(3) = 1; % t0
        data.sorbent(i).CO2Isotherm.param(4) = 0; % alpha
        data.sorbent(i).CO2Isotherm.param(5) =  0.4728; % b0 [1/Mpa]
        data.sorbent(i).CO2Isotherm.param(6) = 2.041e+4; % dH [J/mol]

end


%% sorbent specific data    --> assumption

data.sorbent(i).Density = rho; % density sorbent, kg/m3 (particle)
data.sorbent(i).MaterialDensity = rhoMat; % density material sorbent, kg/m3
data.sorbent(i).cp = cp; %
