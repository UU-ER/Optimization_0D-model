
%% COOLING STEP

function outputCool = simulateCooling(data,outputBlowEvac)


%%% for sub-ambient cooling adsorption set Ref="yes"
Ref="yes";

%% GET ISOTHERM PARAMETERS
i=data.currentSorbent;
CO2IsothermModel = data.sorbent(i).CO2Isotherm.model;
rho = data.sorbent(i).MaterialDensity;

switch CO2IsothermModel
    case 'toth_cp'
        % CO2 isotherm
        T0 = data.sorbent(i).CO2Isotherm.T0(1);
        Xi_c = data.sorbent(i).CO2Isotherm.param(1);
        dH_c = data.sorbent(i).CO2Isotherm.param(2); % [J/mol]
        alpha_c = data.sorbent(i).CO2Isotherm.param(3);
        Xi_p = data.sorbent(i).CO2Isotherm.param(4);
        dH_p = data.sorbent(i).CO2Isotherm.param(5); % [J/mol]
        alpha_p = data.sorbent(i).CO2Isotherm.param(6);
        ns0_c = data.sorbent(i).CO2Isotherm.param(7); % [mol/kg]
        b0_c = data.sorbent(i).CO2Isotherm.param(8); % [1/MPa]
        t0_c = data.sorbent(i).CO2Isotherm.param(9);
        ns0_p = data.sorbent(i).CO2Isotherm.param(10); % [mol/kg]
        b0_p = data.sorbent(i).CO2Isotherm.param(11); % [1/MPa]
        t0_p = data.sorbent(i).CO2Isotherm.param(12);
    case 's_shaped'
        T0 = data.sorbent(i).CO2Isotherm.T0; % K
        q_L0 = data.sorbent(i).CO2Isotherm.param(1); % q_L0, mol/kg, a1
        b_L0 = data.sorbent(i).CO2Isotherm.param(2); % b_L0, 1/MPa, b0
        dU_L = data.sorbent(i).CO2Isotherm.param(3); % dU_L, J/mol, b1
        q_U0 = data.sorbent(i).CO2Isotherm.param(4); % q_U0, J/mol, c1
        b_U0 = data.sorbent(i).CO2Isotherm.param(5); % b_U0, 1/MPa d0
        dU_U = data.sorbent(i).CO2Isotherm.param(6); % dU_U, J/mol d1
        b_H0 = data.sorbent(i).CO2Isotherm.param(7); % b_H0, mol/kg/Pa  bb0
        dU_H = data.sorbent(i).CO2Isotherm.param(8); % dU_H, J/mol bb1
        xi_1 = data.sorbent(i).CO2Isotherm.param(9); % xi1
        xi_2 = data.sorbent(i).CO2Isotherm.param(10); % xi_2, 1/K xi2
        p_step0 = data.sorbent(i).CO2Isotherm.param(11); % p_step0, MPa ps0
        dH_step = data.sorbent(i).CO2Isotherm.param(12); % dH_step, J/mol Hst
        gam = data.sorbent(i).CO2Isotherm.param(13); % gam
    case 'DSL' % the sorbents by Farooq are already in J/mol -> no need for *1000
        n1 = data.sorbent(i).CO2Isotherm.param(1); % n1 (mol/kg)
        b0 = data.sorbent(i).CO2Isotherm.param(2)/(1e6); % b0 (m3/mol)
        Hb = data.sorbent(i).CO2Isotherm.param(3); % Hb (J/mol)
        n2 = data.sorbent(i).CO2Isotherm.param(4); % n2 (mol/kg)
        d0 = data.sorbent(i).CO2Isotherm.param(5)/(1e6); % d0 (m3/mol)
        Hd = data.sorbent(i).CO2Isotherm.param(6); % Hd (J/mol)
    case 'DSL2'
        T0 = data.sorbent(i).CO2Isotherm.T0; % K
        n1 = data.sorbent(i).CO2Isotherm.param(1); % n1 (mol/kg)
        b0 = data.sorbent(i).CO2Isotherm.param(2); % b0 (m3/mol)
        Hb = data.sorbent(i).CO2Isotherm.param(3); % Hb (J/mol)
        n2 = data.sorbent(i).CO2Isotherm.param(4); % n2 (mol/kg)
        d0 = data.sorbent(i).CO2Isotherm.param(5); % d0 (m3/mol)
        Hd = data.sorbent(i).CO2Isotherm.param(6); % Hd (J/mol)
    case 'toth'
        % CO2 isotherm
        T0 = data.sorbent(i).CO2Isotherm.T0(1);
        Xi = data.sorbent(i).CO2Isotherm.param(1);
        dH = data.sorbent(i).CO2Isotherm.param(2); % [J/mol]
        alpha = data.sorbent(i).CO2Isotherm.param(3);
        ns0 = data.sorbent(i).CO2Isotherm.param(4); % [mol/kg]
        b0 = data.sorbent(i).CO2Isotherm.param(5); % [1/MPa]
        t0 = data.sorbent(i).CO2Isotherm.param(6);
    case 'langfr'
        % CO2 isotherm
        T0 = data.sorbent(i).CO2Isotherm.T0(1);
        ns0 = data.sorbent(i).CO2Isotherm.param(1);
        Xi = data.sorbent(i).CO2Isotherm.param(2);
        t0 = data.sorbent(i).CO2Isotherm.param(3);
        alpha = data.sorbent(i).CO2Isotherm.param(4); % [mol/kg]
        b0 = data.sorbent(i).CO2Isotherm.param(5); % [1/MPa]
        dH = data.sorbent(i).CO2Isotherm.param(6); % [J/mol]
end


% N2
N2IsothermModel = data.sorbent(i).N2Isotherm.model;

switch N2IsothermModel
    case 'toth_cp'
        % N2 isotherm
        nT0 = data.sorbent(i).N2Isotherm.T0(1);
        nXi_c = data.sorbent(i).N2Isotherm.param(1);
        ndH_c = data.sorbent(i).N2Isotherm.param(2); % [J/mol]
        nalpha_c = data.sorbent(i).N2Isotherm.param(3);
        nXi_p = data.sorbent(i).N2Isotherm.param(4);
        ndH_p = data.sorbent(i).N2Isotherm.param(5); % [J/mol]
        nalpha_p = data.sorbent(i).N2Isotherm.param(6);
        nns0_c = data.sorbent(i).N2Isotherm.param(7); % [mol/kg]
        nb0_c = data.sorbent(i).N2Isotherm.param(8); % [1/MPa]
        nt0_c = data.sorbent(i).N2Isotherm.param(9);
        nns0_p = data.sorbent(i).N2Isotherm.param(10); % [mol/kg]
        nb0_p = data.sorbent(i).N2Isotherm.param(11); % [1/MPa]
        nt0_p = data.sorbent(i).N2Isotherm.param(12);
    case 's_shaped'
        nT0 = data.sorbent(i).N2Isotherm.T0; % K
        nq_L0 = data.sorbent(i).N2Isotherm.param(1); % q_L0, mol/kg, a1
        nb_L0 = data.sorbent(i).N2Isotherm.param(2); % b_L0, 1/MPa, b0
        ndU_L = data.sorbent(i).N2Isotherm.param(3); % dU_L, J/mol, b1
        nq_U0 = data.sorbent(i).N2Isotherm.param(4); % q_U0, J/mol, c1
        nb_U0 = data.sorbent(i).N2Isotherm.param(5); % b_U0, 1/MPa d0
        ndU_U = data.sorbent(i).N2Isotherm.param(6); % dU_U, J/mol d1
        nb_H0 = data.sorbent(i).N2Isotherm.param(7); % b_H0, mol/kg/Pa  bb0
        ndU_H = data.sorbent(i).N2Isotherm.param(8); % dU_H, J/mol bb1
        nxi_1 = data.sorbent(i).N2Isotherm.param(9); % xi1
        nxi_2 = data.sorbent(i).N2Isotherm.param(10); % xi_2, 1/K xi2
        np_step0 = data.sorbent(i).N2Isotherm.param(11); % p_step0, MPa ps0
        ndH_step = data.sorbent(i).N2Isotherm.param(12); % dH_step, J/mol Hst
        ngam = data.sorbent(i).N2Isotherm.param(13); % gam
    case 'DSL' % the sorbents by Farooq are already in J/mol -> no need for *1000
        nn1 = data.sorbent(i).N2Isotherm.param(1); % n1 (mol/kg)
        nb0 = data.sorbent(i).N2Isotherm.param(2)/(1e6); % b0 (m3/mol)
        nHb = data.sorbent(i).N2Isotherm.param(3); % Hb (J/mol)
        nn2 = data.sorbent(i).N2Isotherm.param(4); % n2 (mol/kg)
        nd0 = data.sorbent(i).N2Isotherm.param(5)/(1e6); % d0 (m3/mol)
        nHd = data.sorbent(i).N2Isotherm.param(6); % Hd (J/mol)
    case 'DSL2'
        nT0 = data.sorbent(i).N2Isotherm.T0; % K
        nn1 = data.sorbent(i).N2Isotherm.param(1); % n1 (mol/kg)
        nb0 = data.sorbent(i).N2Isotherm.param(2); % b0 (m3/mol)
        nHb = data.sorbent(i).N2Isotherm.param(3); % Hb (J/mol)
        nn2 = data.sorbent(i).N2Isotherm.param(4); % n2 (mol/kg)
        nd0 = data.sorbent(i).N2Isotherm.param(5); % d0 (m3/mol)
        nHd = data.sorbent(i).N2Isotherm.param(6); % Hd (J/mol)
    case 'toth'
        % N2 isotherm
        nT0 = data.sorbent(i).N2Isotherm.T0(1);
        nXi = data.sorbent(i).N2Isotherm.param(1);
        ndH = data.sorbent(i).N2Isotherm.param(2); % [J/mol]
        nalpha = data.sorbent(i).N2Isotherm.param(3);
        nns0 = data.sorbent(i).N2Isotherm.param(4); % [mol/kg]
        nb0 = data.sorbent(i).N2Isotherm.param(5); % [1/MPa]
        nt0 = data.sorbent(i).N2Isotherm.param(6);
    case 'langfr'
        % N2 isotherm
        nT0 = data.sorbent(i).N2Isotherm.T0(1);
        nns0 = data.sorbent(i).N2Isotherm.param(1);
        nXi = data.sorbent(i).N2Isotherm.param(2);
        nt0 = data.sorbent(i).N2Isotherm.param(3);
        nalpha = data.sorbent(i).N2Isotherm.param(4); % [mol/kg]
        nb0 = data.sorbent(i).N2Isotherm.param(5); % [1/MPa]
        ndH = data.sorbent(i).N2Isotherm.param(6); % [J/mol]
end


% CH4
CH4IsothermModel = data.sorbent(i).CH4Isotherm.model;
switch CH4IsothermModel
    case 'toth_cp'
        % CH4 isotherm
        cT0 = data.sorbent(i).CH4Isotherm.T0(1);
        cXi_c = data.sorbent(i).CH4Isotherm.param(1);
        cdH_c = data.sorbent(i).CH4Isotherm.param(2); % [J/mol]
        calpha_c = data.sorbent(i).CH4Isotherm.param(3);
        cXi_p = data.sorbent(i).CH4Isotherm.param(4);
        cdH_p = data.sorbent(i).CH4Isotherm.param(5); % [J/mol]
        calpha_p = data.sorbent(i).CH4Isotherm.param(6);
        cns0_c = data.sorbent(i).CH4Isotherm.param(7); % [mol/kg]
        cb0_c = data.sorbent(i).CH4Isotherm.param(8); % [1/MPa]
        ct0_c = data.sorbent(i).CH4Isotherm.param(9);
        cns0_p = data.sorbent(i).CH4Isotherm.param(10); % [mol/kg]
        cb0_p = data.sorbent(i).CH4Isotherm.param(11); % [1/MPa]
        ct0_p = data.sorbent(i).CH4Isotherm.param(12);
    case 's_shaped'
        cT0 = data.sorbent(i).CH4Isotherm.T0; % K
        cq_L0 = data.sorbent(i).CH4Isotherm.param(1); % q_L0, mol/kg, a1
        cb_L0 = data.sorbent(i).CH4Isotherm.param(2); % b_L0, 1/MPa, b0
        cdU_L = data.sorbent(i).CH4Isotherm.param(3); % dU_L, J/mol, b1
        cq_U0 = data.sorbent(i).CH4Isotherm.param(4); % q_U0, J/mol, c1
        cb_U0 = data.sorbent(i).CH4Isotherm.param(5); % b_U0, 1/MPa d0
        cdU_U = data.sorbent(i).CH4Isotherm.param(6); % dU_U, J/mol d1
        cb_H0 = data.sorbent(i).CH4Isotherm.param(7); % b_H0, mol/kg/Pa  bb0
        cdU_H = data.sorbent(i).CH4Isotherm.param(8); % dU_H, J/mol bb1
        cxi_1 = data.sorbent(i).CH4Isotherm.param(9); % xi1
        cxi_2 = data.sorbent(i).CH4Isotherm.param(10); % xi_2, 1/K xi2
        cp_step0 = data.sorbent(i).CH4Isotherm.param(11); % p_step0, MPa ps0
        cdH_step = data.sorbent(i).CH4Isotherm.param(12); % dH_step, J/mol Hst
        cgam = data.sorbent(i).CH4Isotherm.param(13); % gam
    case 'DSL' % the sorbents by Farooq are already in J/mol -> no need for *1000
        cn1 = data.sorbent(i).CH4Isotherm.param(1); % n1 (mol/kg)
        cb0 = data.sorbent(i).CH4Isotherm.param(2)/(1e6); % b0 (m3/mol)
        cHb = data.sorbent(i).CH4Isotherm.param(3); % Hb (J/mol)
        cn2 = data.sorbent(i).CH4Isotherm.param(4); % n2 (mol/kg)
        cd0 = data.sorbent(i).CH4Isotherm.param(5)/(1e6); % d0 (m3/mol)
        cHd = data.sorbent(i).CH4Isotherm.param(6); % Hd (J/mol)
    case 'DSL2'
        cT0 = data.sorbent(i).CH4Isotherm.T0; % K
        cn1 = data.sorbent(i).CH4Isotherm.param(1); % n1 (mol/kg)
        cb0 = data.sorbent(i).CH4Isotherm.param(2); % b0 (m3/mol)
        cHb = data.sorbent(i).CH4Isotherm.param(3); % Hb (J/mol)
        cn2 = data.sorbent(i).CH4Isotherm.param(4); % n2 (mol/kg)
        cd0 = data.sorbent(i).CH4Isotherm.param(5); % d0 (m3/mol)
        cHd = data.sorbent(i).CH4Isotherm.param(6); % Hd (J/mol)
    case 'toth'
        % CH4 isotherm
        cT0 = data.sorbent(i).CH4Isotherm.T0(1);
        cXi = data.sorbent(i).CH4Isotherm.param(1);
        cdH = data.sorbent(i).CH4Isotherm.param(2); % [J/mol]
        calpha = data.sorbent(i).CH4Isotherm.param(3);
        cns0 = data.sorbent(i).CH4Isotherm.param(4); % [mol/kg]
        cb0 = data.sorbent(i).CH4Isotherm.param(5); % [1/MPa]
        ct0 = data.sorbent(i).CH4Isotherm.param(6);
    case 'langfr'
        % CH4 isotherm
        cT0 = data.sorbent(i).CH4Isotherm.T0(1);
        cns0 = data.sorbent(i).CH4Isotherm.param(1);
        cXi = data.sorbent(i).CH4Isotherm.param(2);
        ct0 = data.sorbent(i).CH4Isotherm.param(3);
        calpha = data.sorbent(i).CH4Isotherm.param(4); % [mol/kg]
        cb0 = data.sorbent(i).CH4Isotherm.param(5); % [1/MPa]
        cdH = data.sorbent(i).CH4Isotherm.param(6); % [J/mol]
end

%% define general data
adsorbentMass = data.process.adsorbentMass;
R = data.general.gasconstant;
void = data.process.voidFraction;
V = data.process.Vol;
% boundaries
lb = [0,0,0,0];
ub = [1,1,50,1];
% solver options
options = optimoptions('lsqnonlin','Display','off');
options.Algorithm = 'trust-region-reflective';
options.OptimalityTolerance = 1e-9;
options.FunctionTolerance = 1e-9;
options.StepTolerance = 1e-9;
options.MaxFunctionEvaluations = 300;
%% RUN COOLING MODEL
for m = 1:size(data.TCoolProfileStep,1)
    % Initial condition

    % (1) yCO2, (2), yN2, (3) Nin , (4) yCH4
    x0 = [outputBlowEvac(m).yCO2(end-1),outputBlowEvac(m).yN2(end-1),0,outputBlowEvac(m).yCH4(end-1)];

    % feed stream
    p = data.process.pamb; %
    yCO2_feed = data.feed.yCO2;
    yCH4_feed=data.feed.yCH4;
    yN2_feed = data.feed.yN2;


    TempVector = data.TCoolProfileStep(m,:);
    coolingTime = data.process.noSteps;

    pressure_vector = data.PressureCoolProfileStep(m,:);

    yCO2_k = outputBlowEvac(m).yCO2(end);
    yN2_k = outputBlowEvac(m).yN2(end);
    yCH4_k=outputBlowEvac(m).yCH4(end);
    T_k = data.process.Tdes;
    p_k = p; % pressure does not change

    % create vectors for saving results
    yCO2_save = zeros(coolingTime,1);
    T_save = zeros(coolingTime,1);
    yN2_save = zeros(coolingTime,1);
    yCH4_save = zeros(coolingTime,1);
    Nin_save = zeros(coolingTime,1);
    x_save = nan(coolingTime,4);
    Nin_save_sum = zeros(coolingTime,1);
    Nin_save_CO2 = zeros(coolingTime,1);
    Nin_save_CH4 = zeros(coolingTime,1);
    Nin_save_N2 = zeros(coolingTime,1);
    Q_save_sum = zeros(coolingTime,1);
    Ntotal_CO2 = zeros(coolingTime,1);
    Ntotal_N2 = zeros(coolingTime,1);
    Ntotal_CH4 = zeros(coolingTime,1);


    for k = 1:coolingTime

        T = TempVector(k);
        p = pressure_vector(k+1);
        p_k = pressure_vector(k);


        funct = @(x) double( fcn_dry(x));

        fun = @(x) funct(x);

        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

        yCO2_save(k) = x(1);
        yN2_save(k) = x(2);
        yCH4_save(k)=x(4);
        Nin_save(k) = x(3);
        T_save(k) = T;
        x_save(k,:) = [residual];

        Nin_save_sum(k) = Nin_save_sum(k) + Nin_save(k); % Nin_save_sum(k) is the same as Nin_save(k)
        Nin_save_CO2(k) = Nin_save_CO2(k) + Nin_save(k)*x(1);
        Nin_save_N2(k) = Nin_save_N2(k) + Nin_save(k)*x(2);
        Nin_save_CH4(k) = Nin_save_CH4(k) + Nin_save(k)*x(4);

        [Q, NCO2,  NN2, NCH4, qCO2s] = Qcool(x,data.process.adsorbentMass,data.sorbent(i).cp,data.general.gasconstant);
        
        if T < data.process.Tamb
            Q_save_sum(k+1) =  Q_save_sum(k) + Q; %
        else
            Q = 0;
            Q_save_sum(k+1) =  Q_save_sum(k) + Q;
        end
        
        qCO2_save(k) = qCO2s;
        yCO2_k = x(1);
        yN2_k = x(2);
        yCH4_k= x(4);
        T_k = T;
        x0 = x;

        Ntotal_CO2(k) = NCO2;
        Ntotal_N2(k) = NN2;
        Ntotal_CH4(k) = NCH4;
    end
    if data.plot
        plot_all = 'step_2';
        separate_plot_file_new2;
    end
    delta_CO2 = yCO2_save(end) - data.feed.yCO2; %
    if delta_CO2 > 0
        sprintf('WARNING: CO2 concentration at the end of the cooling step higher than CO2 composition in feed (%0.6f)!', delta_CO2)
    end

    delta_CH4 = yCH4_save(end)-data.feed.yCH4;
    if delta_CH4 > 0
                sprintf('WARNING: CH4 concentration at the end of the cooling step higher than CH4 composition in feed (%0.6f)!', delta_CH4)
    end

    %%% calculate the work for refregration
    Tamb=data.process.Tamb;
    Tads=data.process.Tads;

    if Tads~=Tamb
        if delta_CO2<0 && delta_CH4<0 && Ref=="yes"
            E_Ref= W_Refrig(data.process.Tamb,data.process.Tads,sum(Nin_save_sum)); % E(KJ)
        else 
           E_Ref=0; 
        end
    else
        E_Ref=0;
    end
   

    %% output
%     outputCool = [yCO2_save, yN2_save, yCH4_save, Nin_save, Nin_save_CO2(2:end), Nin_save_N2(2:end), Nin_save_CH4(2:end) ,T_save, Q_save_sum(2:end), Ntotal_CO2, Ntotal_N2, Ntotal_CH4];
    %  outputCool = [yCO2_save, yN2_save, yH2O_save, yCH4_save, Nin_save, Nin_save_CO2(2:end), Nin_save_N2(2:end), Nin_save_H2O(2:end),Nin_save_CH4(2:end) ,T_save, Q_save_sum(2:end), Ntotal_CO2, Ntotal_N2, Ntotal_H2O, Ntotal_CH4];
    outputCool(m).yCO2 = yCO2_save;
    outputCool(m).yN2 = yN2_save;
    outputCool(m).yCH4 = yCH4_save;
    outputCool(m).Nin_sum = Nin_save;
    outputCool(m).Nin_CO2 = Nin_save_CO2(2:end);
    outputCool(m).Nin_N2 = Nin_save_N2(2:end);
    outputCool(m).Nin_CH4 = Nin_save_CH4(2:end);
    outputCool(m).T = T_save;
    outputCool(m).Q_sum = Q_save_sum;
    outputCool(m).Ntotal_CO2 = Ntotal_CO2;
    outputCool(m).Ntotal_N2 = Ntotal_N2;
    outputCool(m).Ntotal_CH4 = Ntotal_CH4;
    outputCool(m).E_Refrig = E_Ref;
end
%% functions
    function [Q, NCO2, NCH4, NN2, qCO2] = Qcool(x,adsorbentMass,cp,R)

        %% equation 5: energy balance
        % CO2
        % solid phase CO2, initial condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2_k = ((ns0_c*exp(Xi_c*(1-T_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./T_k))).^(1./(t0_c+alpha_c*(1-T0./T_k)))))+((ns0_p*exp(Xi_p*(1-T_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./T_k))).^(1./(t0_p+alpha_p*(1-T0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
            case 'DSL'
                qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k));
            case 'toth'
                qCO2_k = ((ns0*exp(Xi*(1-T_k/T0)).*(b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0+alpha*(1-T0./T_k))).^(1./(t0+alpha*(1-T0./T_k)))));
            case 'langfr'
                qCO2_k = (ns0.*exp(Xi.*(1-T_k./T0))).*((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k)));
        end
        NCO2_solid_k = adsorbentMass*qCO2_k;

        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-T/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T))).^(1./(t0_p+alpha_p*(1-T0./T))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(x(1)*p)./(1+(b_L0.*exp(dU_L./(R*T))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(x(1)*p)./(1+(b_U0.*exp(dU_U./(R*T))).*(x(1)*p))+(b_H0.*exp(dU_H./(R*T))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); % adsorbed amount of CO2
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(x(1)*p)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(x(1)*p)./(1e-6*R*T));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1+b0*exp(Hb/R/T).*(x(1)*p)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1+d0*exp(Hd/R/T).*(x(1)*p));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T)));
        end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid = x(1)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid CO2, final condition
        NCO2 = NCO2_solid + NCO2_fluid;


        % N2
        % solid phase N2, initial condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2_k = ((nns0_c*exp(nXi_c*(1-T_k/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_c+nalpha_c*(1-nT0./T_k))).^(1./(nt0_c+nalpha_c*(1-nT0./T_k)))))+((nns0_p*exp(nXi_p*(1-T_k/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_p+nalpha_p*(1-nT0./T_k))).^(1./(nt0_p+nalpha_p*(1-nT0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2_k = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T_k))).*(yN2_k*p_k)./(1+(nb_L0.*exp(ndU_L./(R*T_k))).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T_k))).*(yN2_k*p_k)./(1+(nb_U0.*exp(ndU_U./(R*T_k))).*(yN2_k*p_k))+(nb_H0.*exp(ndU_H./(R*T_k))).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k)))))).^ngam);
            case 'DSL'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k));
            case 'toth'
                qN2_k = ((nns0*exp(nXi*(1-T_k/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0+nalpha*(1-nT0./T_k))).^(1./(nt0+nalpha*(1-nT0./T_k)))));
            case 'langfr'
                qN2_k = (nns0.*exp(nXi.*(1-T_k./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k)));
        end
        NN2_solid_k = adsorbentMass*qN2_k;

        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p*x(2)))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p*x(2))).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-T/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T-1))).*(p*x(2)))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T-1))).*(p*x(2))).^(nt0_p+nalpha_p*(1-nT0./T))).^(1./(nt0_p+nalpha_p*(1-nT0./T))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*(x(2)*p)./(1+(nb_L0.*exp(ndU_L./(R*T))).*(x(2)*p))).*(1-((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((x(2)*p)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*(x(2)*p)./(1+(nb_U0.*exp(ndU_U./(R*T))).*(x(2)*p))+(nb_H0.*exp(ndU_H./(R*T))).*(x(2)*p)).*((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((x(2)*p)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam); % adsorbed amount of CO2
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*(x(2)*p)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*(x(2)*p)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*(x(2)*p)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*(x(2)*p)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*(x(2)*p)./(1+nb0*exp(nHb/R/T).*(x(2)*p)) + nn2*nd0*exp(nHd/R/T)*(x(2)*p)./(1+nd0*exp(nHd/R/T).*(x(2)*p));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(x(2)*p))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(x(2)*p)).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = x(2)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid N2, final condition
        NN2 = NN2_solid + NN2_fluid;
        % fluid phase N2, final condition
        %             NN2 = x(2)*p*V*void/(R*T)*1e6; % mol

        % CH4
        % solid phase CH4, initial condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4_k = ((cns0_c*exp(cXi_c*(1-T_k/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_c+calpha_c*(1-cT0./T_k))).^(1./(ct0_c+calpha_c*(1-cT0./T_k)))))+((cns0_p*exp(cXi_p*(1-T_k/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_p+calpha_p*(1-cT0./T_k))).^(1./(ct0_p+calpha_p*(1-cT0./T_k))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4_k = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k))).*(1-((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k))+(cb_H0.*exp(cdU_H./(R*T_k))).*(yCH4_k*p_k)).*((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam);
            case 'DSL'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k));
            case 'toth'
                qCH4_k = ((cns0*exp(cXi*(1-T_k/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0+calpha*(1-cT0./T_k))).^(1./(ct0+calpha*(1-cT0./T_k)))));
            case 'langfr'
                qCH4_k = (cns0.*exp(cXi.*(1-T_k./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k)));
        end
        NCH4_solid_k = adsorbentMass*qCH4_k;

        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p*x(4)))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p*x(4))).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-T/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T-1))).*(p*x(4)))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T-1))).*(p*x(4))).^(ct0_p+calpha_p*(1-cT0./T))).^(1./(ct0_p+calpha_p*(1-cT0./T))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(x(4)*p)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(x(4)*p))).*(1-((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((x(4)*p)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(x(4)*p)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(x(4)*p))+(cb_H0.*exp(cdU_H./(R*T))).*(x(4)*p)).*((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((x(4)*p)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam); % adsorbed amount of CH4
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(x(4)*p)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(x(4)*p)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(x(4)*p)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(x(4)*p)./(1e-6*R*T));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(x(4)*p)./(1+cb0*exp(cHb/R/T).*(x(4)*p)) + cn2*cd0*exp(cHd/R/T)*(x(4)*p)./(1+cd0*exp(cHd/R/T).*(x(4)*p));
            case 'toth'
                qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(x(4)*p))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(x(4)*p)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
            case 'langfr'
                qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T)));
        end
        NCH4_solid = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid = x(4)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid CH4, final condition
        NCH4 = NCH4_solid + NCH4_fluid;


        %% equation 5: energy balance
        % CO2
        switch CO2IsothermModel
            case 'toth_cp'
                dqCO2_dTequ_k = - b0_c.*ns0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)).*(((T0.*alpha_c.*log(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).*(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)))./T_k.^2 - (b0_c.*dH_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*(t0_c - alpha_c.*(T0./T_k - 1)).*(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1) - 1))./(R.*T_k.^2))./((t0_c - alpha_c.*(T0./T_k - 1)).*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)) + 1)) - (T0.*alpha_c.*log((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1))./(T_k.^2.*(t0_c - alpha_c.*(T0./T_k - 1)).^2.*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1))))) - b0_p.*ns0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)).*(((T0.*alpha_p.*log(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).*(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)))./T_k.^2 - (b0_p.*dH_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*(t0_p - alpha_p.*(T0./T_k - 1)).*(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1) - 1))./(R.*T_k.^2))./((t0_p - alpha_p.*(T0./T_k - 1)).*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1)) + 1)) - (T0.*alpha_p.*log((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1))./(T_k.^2.*(t0_p - alpha_p.*(T0./T_k - 1)).^2.*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1))))) - (Xi_c.*b0_c.*ns0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)))./(T0.*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)))) - (Xi_p.*b0_p.*ns0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)))./(T0.*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1)))) - (b0_c.*dH_c.*ns0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)))./(R.*T_k.^2.*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)))) - (b0_p.*dH_p.*ns0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)))./(R.*T_k.^2.*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1))));
                dqCO2_dpCO2_k = (b0_c.*ns0_c.*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)))./((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1))) + (b0_p.*ns0_p.*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)))./((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1))) - (b0_c.^2.*ns0_c.*(p_k.*yCO2_k).*exp((2.*dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)).*(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1) - 1))./((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)) + 1) - (b0_p.^2.*ns0_p.*(p_k.*yCO2_k).*exp((2.*dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)).*(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1) - 1))./((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1)) + 1);
            case 's_shaped'
                dqCO2_dTequ_k = gam.*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1) - (exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1).*(b_H0.*(p_k.*yCO2_k).*exp(dU_H./(R.*T_k)) + (b_U0.*(p_k.*yCO2_k).*q_U0.*exp(dU_U./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1)) - (exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam.*((b_H0.*dU_H.*(p_k.*yCO2_k).*exp(dU_H./(R.*T_k)))./(R.*T_k.^2) - (b_U0.^2.*dU_U.*(p_k.*yCO2_k).^2.*q_U0.*exp((2.*dU_U)./(R.*T_k)))./(R.*T_k.^2.*(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1).^2) + (b_U0.*dU_U.*(p_k.*yCO2_k).*q_U0.*exp(dU_U./(R.*T_k)))./(R.*T_k.^2.*(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1))) - (b_L0.*gam.*(p_k.*yCO2_k).*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1) - (exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1) - (b_L0.^2.*dU_L.*(p_k.*yCO2_k).^2.*q_L0.*exp((2.*dU_L)./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(R.*T_k.^2.*(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1).^2) + (b_L0.*dU_L.*(p_k.*yCO2_k).*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(R.*T_k.^2.*(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1));
                dqCO2_dpCO2_k =  (exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam.*(b_H0.*exp(dU_H./(R.*T_k)) + (b_U0.*q_U0.*exp(dU_U./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1) - (b_U0.^2.*(p_k.*yCO2_k).*q_U0.*exp((2.*dU_U)./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1).^2) + gam.*((exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)) - (exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2)).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1).*(b_H0.*(p_k.*yCO2_k).*exp(dU_H./(R.*T_k)) + (b_U0.*(p_k.*yCO2_k).*q_U0.*exp(dU_U./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1)) - (b_L0.*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1) + (b_L0.^2.*(p_k.*yCO2_k).*q_L0.*exp((2.*dU_L)./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1).^2 - (b_L0.*gam.*(p_k.*yCO2_k).*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)) - (exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2)).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1);
            case 'DSL'
                dqCO2_dTequ_k = (1000000.*b0.*n1.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)).*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k.^2) + (1000000.*Hb.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000.*d0.*n2.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k.^2.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*Hb.*b0.*n1.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*Hd.*d0.*n2.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*b0.*n1.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k.^2.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*d0.*n2.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)).*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k.^2) + (1000000.*Hd.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1).^2);
                dqCO2_dpCO2_k = (1000000.*b0.*n1.*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*d0.*n2.*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000000000.*b0.^2.*n1.*(yCO2_k.*p_k).*yCO2_k.^2.*exp((2.*Hb)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000000000.*d0.^2.*n2.*(yCO2_k.*p_k).*yCO2_k.^2.*exp((2.*Hd)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1).^2);
            case 'DSL2'
                dqCO2_dTequ_k = (Hb.*b0.^2.*n1.*(p_k.*yCO2_k).^2.*exp((2.*Hb)./(R.*T_k)))./(R.*T_k.^2.*(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1).^2) - (Hd.*d0.*n2.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)))./(R.*T_k.^2.*(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1)) - (Hb.*b0.*n1.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)))./(R.*T_k.^2.*(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1)) + (Hd.*d0.^2.*n2.*(p_k.*yCO2_k).^2.*exp((2.*Hd)./(R.*T_k)))./(R.*T_k.^2.*(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1).^2);
                dqCO2_dpCO2_k = (b0.*n1.*exp(Hb./(R.*T_k)))./(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1) + (d0.*n2.*exp(Hd./(R.*T_k)))./(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1) - (b0.^2.*n1.*(p_k.*yCO2_k).*exp((2.*Hb)./(R.*T_k)))./(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1).^2 - (d0.^2.*n2.*(p_k.*yCO2_k).*exp((2.*Hd)./(R.*T_k)))./(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1).^2;
            case 'toth'
                dqCO2_dTequ_k = - b0.*ns0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(((T0.*alpha.*log(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)))./T_k.^2 - (b0.*dH.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*(t0 - alpha.*(T0./T_k - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1) - 1))./(R.*T_k.^2))./((t0 - alpha.*(T0./T_k - 1)).*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1)) + 1)) - (T0.*alpha.*log((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1))./(T_k.^2.*(t0 - alpha.*(T0./T_k - 1)).^2.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1))))) - (Xi.*b0.*ns0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)))./(T0.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1)))) - (b0.*dH.*ns0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)))./(R.*T_k.^2.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1))));
                dqCO2_dpCO2_k = (b0.*ns0.*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1))) - (b0.^2.*ns0.*(p_k.*yCO2_k).*exp((2.*dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1) - 1))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1)) + 1);
            case 'langfr'
                dqCO2_dTequ_k = (ns0.*exp(-Xi.*(T_k./T0 - 1)).*((T0.*alpha.*log(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./T_k.^2 + (b0.*dH.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./(R.*T_k.^2)))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1) - (ns0.*exp(-Xi.*(T_k./T0 - 1)).*((T0.*alpha.*log(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./T_k.^2 + (b0.*dH.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./(R.*T_k.^2)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1).^2 - (Xi.*ns0.*exp(-Xi.*(T_k./T0 - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./(T0.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1));
                dqCO2_dpCO2_k = (b0.*ns0.*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1).^2 - (b0.*ns0.*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1);
        end
        % dH
        dHCO2 = R*T_k^2/(p_k*yCO2_k) *(-dqCO2_dTequ_k)/(dqCO2_dpCO2_k);
        % CO2 desorbed
        dNCO2_solid = NCO2_solid - NCO2_solid_k;


        % N2
        switch N2IsothermModel
            case 'toth_cp'
                dqN2_dTequ_k = - nb0_c.*nns0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_c.*(T_k./nT0 - 1)).*(((nT0.*nalpha_c.*log(nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).*(nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)))./T_k.^2 - (nb0_c.*ndH_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0)).*(nt0_c - nalpha_c.*(nT0./T_k - 1)).*(nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1) - 1))./(R.*T_k.^2))./((nt0_c - nalpha_c.*(nT0./T_k - 1)).*((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1).^(1./(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1)) - (nT0.*nalpha_c.*log((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1))./(T_k.^2.*(nt0_c - nalpha_c.*(nT0./T_k - 1)).^2.*((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1).^(1./(nt0_c - nalpha_c.*(nT0./T_k - 1))))) - nb0_p.*nns0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_p.*(T_k./nT0 - 1)).*(((nT0.*nalpha_p.*log(nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).*(nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)))./T_k.^2 - (nb0_p.*ndH_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0)).*(nt0_p - nalpha_p.*(nT0./T_k - 1)).*(nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1) - 1))./(R.*T_k.^2))./((nt0_p - nalpha_p.*(nT0./T_k - 1)).*((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1).^(1./(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1)) - (nT0.*nalpha_p.*log((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1))./(T_k.^2.*(nt0_p - nalpha_p.*(nT0./T_k - 1)).^2.*((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1).^(1./(nt0_p - nalpha_p.*(nT0./T_k - 1))))) - (nXi_c.*nb0_c.*nns0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_c.*(T_k./nT0 - 1)))./(nT0.*((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1).^(1./(nt0_c - nalpha_c.*(nT0./T_k - 1)))) - (nXi_p.*nb0_p.*nns0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_p.*(T_k./nT0 - 1)))./(nT0.*((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1).^(1./(nt0_p - nalpha_p.*(nT0./T_k - 1)))) - (nb0_c.*ndH_c.*nns0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_c.*(T_k./nT0 - 1)))./(R.*T_k.^2.*((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1).^(1./(nt0_c - nalpha_c.*(nT0./T_k - 1)))) - (nb0_p.*ndH_p.*nns0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_p.*(T_k./nT0 - 1)))./(R.*T_k.^2.*((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1).^(1./(nt0_p - nalpha_p.*(nT0./T_k - 1))));
                dqN2_dpN2_k = (nb0_c.*nns0_c.*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_c.*(T_k./nT0 - 1)))./((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1).^(1./(nt0_c - nalpha_c.*(nT0./T_k - 1))) + (nb0_p.*nns0_p.*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_p.*(T_k./nT0 - 1)))./((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1).^(1./(nt0_p - nalpha_p.*(nT0./T_k - 1))) - (nb0_c.^2.*nns0_c.*(p_k.*yN2_k).*exp((2.*ndH_c.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_c.*(T_k./nT0 - 1)).*(nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1) - 1))./((nb0_c.*(p_k.*yN2_k).*exp((ndH_c.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1).^(1./(nt0_c - nalpha_c.*(nT0./T_k - 1)) + 1) - (nb0_p.^2.*nns0_p.*(p_k.*yN2_k).*exp((2.*ndH_p.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi_p.*(T_k./nT0 - 1)).*(nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1) - 1))./((nb0_p.*(p_k.*yN2_k).*exp((ndH_p.*(nT0./T_k - 1))./(R.*nT0))).^(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1).^(1./(nt0_p - nalpha_p.*(nT0./T_k - 1)) + 1);
            case 's_shaped'
                dqN2_dTequ_k = ngam.*((exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1).*((ndH_step.*exp(-nxi_2.*(1./nT0 - 1./T_k)))./(R.*T_k.^2.*nxi_1) + (nxi_2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./(T_k.^2.*nxi_1)))./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1) - (exp(-(2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1).*((ndH_step.*exp(-nxi_2.*(1./nT0 - 1./T_k)))./(R.*T_k.^2.*nxi_1) + (nxi_2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./(T_k.^2.*nxi_1)))./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1).^2).*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^(ngam - 1).*(nb_H0.*(p_k.*yN2_k).*exp(ndU_H./(R.*T_k)) + (nb_U0.*(p_k.*yN2_k).*nq_U0.*exp(ndU_U./(R.*T_k)))./(nb_U0.*(p_k.*yN2_k).*exp(ndU_U./(R.*T_k)) + 1)) - (exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^ngam.*((nb_H0.*ndU_H.*(p_k.*yN2_k).*exp(ndU_H./(R.*T_k)))./(R.*T_k.^2) - (nb_U0.^2.*ndU_U.*(p_k.*yN2_k).^2.*nq_U0.*exp((2.*ndU_U)./(R.*T_k)))./(R.*T_k.^2.*(nb_U0.*(p_k.*yN2_k).*exp(ndU_U./(R.*T_k)) + 1).^2) + (nb_U0.*ndU_U.*(p_k.*yN2_k).*nq_U0.*exp(ndU_U./(R.*T_k)))./(R.*T_k.^2.*(nb_U0.*(p_k.*yN2_k).*exp(ndU_U./(R.*T_k)) + 1))) - (nb_L0.*ngam.*(p_k.*yN2_k).*nq_L0.*exp(ndU_L./(R.*T_k)).*((exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1).*((ndH_step.*exp(-nxi_2.*(1./nT0 - 1./T_k)))./(R.*T_k.^2.*nxi_1) + (nxi_2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./(T_k.^2.*nxi_1)))./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1) - (exp(-(2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1).*((ndH_step.*exp(-nxi_2.*(1./nT0 - 1./T_k)))./(R.*T_k.^2.*nxi_1) + (nxi_2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./(T_k.^2.*nxi_1)))./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1).^2).*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^(ngam - 1))./(nb_L0.*(p_k.*yN2_k).*exp(ndU_L./(R.*T_k)) + 1) - (nb_L0.^2.*ndU_L.*(p_k.*yN2_k).^2.*nq_L0.*exp((2.*ndU_L)./(R.*T_k)).*((exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^ngam - 1))./(R.*T_k.^2.*(nb_L0.*(p_k.*yN2_k).*exp(ndU_L./(R.*T_k)) + 1).^2) + (nb_L0.*ndU_L.*(p_k.*yN2_k).*nq_L0.*exp(ndU_L./(R.*T_k)).*((exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^ngam - 1))./(R.*T_k.^2.*(nb_L0.*(p_k.*yN2_k).*exp(ndU_L./(R.*T_k)) + 1));
                dqN2_dpN2_k =  (exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^ngam.*(nb_H0.*exp(ndU_H./(R.*T_k)) + (nb_U0.*nq_U0.*exp(ndU_U./(R.*T_k)))./(nb_U0.*(p_k.*yN2_k).*exp(ndU_U./(R.*T_k)) + 1) - (nb_U0.^2.*(p_k.*yN2_k).*nq_U0.*exp((2.*ndU_U)./(R.*T_k)))./(nb_U0.*(p_k.*yN2_k).*exp(ndU_U./(R.*T_k)) + 1).^2) + ngam.*((exp(-nxi_2.*(1./nT0 - 1./T_k)).*exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1))./((p_k.*yN2_k).*nxi_1.*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)) - (exp(-nxi_2.*(1./nT0 - 1./T_k)).*exp(-(2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1))./((p_k.*yN2_k).*nxi_1.*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1).^2)).*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^(ngam - 1).*(nb_H0.*(p_k.*yN2_k).*exp(ndU_H./(R.*T_k)) + (nb_U0.*(p_k.*yN2_k).*nq_U0.*exp(ndU_U./(R.*T_k)))./(nb_U0.*(p_k.*yN2_k).*exp(ndU_U./(R.*T_k)) + 1)) - (nb_L0.*nq_L0.*exp(ndU_L./(R.*T_k)).*((exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^ngam - 1))./(nb_L0.*(p_k.*yN2_k).*exp(ndU_L./(R.*T_k)) + 1) + (nb_L0.^2.*(p_k.*yN2_k).*nq_L0.*exp((2.*ndU_L)./(R.*T_k)).*((exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^ngam - 1))./(nb_L0.*(p_k.*yN2_k).*exp(ndU_L./(R.*T_k)) + 1).^2 - (nb_L0.*ngam.*(p_k.*yN2_k).*nq_L0.*exp(ndU_L./(R.*T_k)).*((exp(-nxi_2.*(1./nT0 - 1./T_k)).*exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1))./((p_k.*yN2_k).*nxi_1.*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)) - (exp(-nxi_2.*(1./nT0 - 1./T_k)).*exp(-(2.*exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1))./((p_k.*yN2_k).*nxi_1.*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1).^2)).*(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1)./(exp(-(exp(-nxi_2.*(1./nT0 - 1./T_k)).*(log(np_step0.*exp(-(ndH_step.*(1./nT0 - 1./T_k))./R)) - log((p_k.*yN2_k))))./nxi_1) + 1)).^(ngam - 1))./(nb_L0.*(p_k.*yN2_k).*exp(ndU_L./(R.*T_k)) + 1);
            case 'DSL'
                dqN2_dTequ_k = (1000000.*nb0.*nn1.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)).*((1000000.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k.^2) + (1000000.*nHb.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000.*nd0.*nn2.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k.^2.*((1000000.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*nHb.*nb0.*nn1.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*nHd.*nd0.*nn2.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*nb0.*nn1.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k.^2.*((1000000.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*nd0.*nn2.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)).*((1000000.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k.^2) + (1000000.*nHd.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k) + 1).^2);
                dqN2_dpN2_k = (1000000.*nb0.*nn1.*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k.*((1000000.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*nd0.*nn2.*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k.*((1000000.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000000000.*nb0.^2.*nn1.*(yN2_k.*p_k).*yN2_k.^2.*exp((2.*nHb)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*nb0.*(yN2_k.*p_k).*yN2_k.*exp(nHb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000000000.*nd0.^2.*nn2.*(yN2_k.*p_k).*yN2_k.^2.*exp((2.*nHd)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*nd0.*(yN2_k.*p_k).*yN2_k.*exp(nHd./(R.*T_k)))./(R.*T_k) + 1).^2);
            case 'DSL2'
                dqN2_dTequ_k = (nHb.*nb0.^2.*nn1.*(p_k.*yN2_k).^2.*exp((2.*nHb)./(R.*T_k)))./(R.*T_k.^2.*(nb0.*(p_k.*yN2_k).*exp(nHb./(R.*T_k)) + 1).^2) - (nHd.*nd0.*nn2.*(p_k.*yN2_k).*exp(nHd./(R.*T_k)))./(R.*T_k.^2.*(nd0.*(p_k.*yN2_k).*exp(nHd./(R.*T_k)) + 1)) - (nHb.*nb0.*nn1.*(p_k.*yN2_k).*exp(nHb./(R.*T_k)))./(R.*T_k.^2.*(nb0.*(p_k.*yN2_k).*exp(nHb./(R.*T_k)) + 1)) + (nHd.*nd0.^2.*nn2.*(p_k.*yN2_k).^2.*exp((2.*nHd)./(R.*T_k)))./(R.*T_k.^2.*(nd0.*(p_k.*yN2_k).*exp(nHd./(R.*T_k)) + 1).^2);
                dqN2_dpN2_k = (nb0.*nn1.*exp(nHb./(R.*T_k)))./(nb0.*(p_k.*yN2_k).*exp(nHb./(R.*T_k)) + 1) + (nd0.*nn2.*exp(nHd./(R.*T_k)))./(nd0.*(p_k.*yN2_k).*exp(nHd./(R.*T_k)) + 1) - (nb0.^2.*nn1.*(p_k.*yN2_k).*exp((2.*nHb)./(R.*T_k)))./(nb0.*(p_k.*yN2_k).*exp(nHb./(R.*T_k)) + 1).^2 - (nd0.^2.*nn2.*(p_k.*yN2_k).*exp((2.*nHd)./(R.*T_k)))./(nd0.*(p_k.*yN2_k).*exp(nHd./(R.*T_k)) + 1).^2;
            case 'toth'
                dqN2_dTequ_k = - nb0.*nns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 - (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nt0 - nalpha.*(nT0./T_k - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2))./((nt0 - nalpha.*(nT0./T_k - 1)).*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)) + 1)) - (nT0.*nalpha.*log((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1))./(T_k.^2.*(nt0 - nalpha.*(nT0./T_k - 1)).^2.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))))) - (nXi.*nb0.*nns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./(nT0.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)))) - (nb0.*ndH.*nns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./(R.*T_k.^2.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))));
                dqN2_dpN2_k = (nb0.*nns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))) - (nb0.^2.*nns0.*(p_k.*yN2_k).*exp((2.*ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)) + 1);
            case 'langfr'
                dqN2_dTequ_k = (nns0.*exp(-nXi.*(T_k./nT0 - 1)).*((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 + (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1) - (nns0.*exp(-nXi.*(T_k./nT0 - 1)).*((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 + (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1).^2 - (nXi.*nns0.*exp(-nXi.*(T_k./nT0 - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./(nT0.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1));
                dqN2_dpN2_k = (nb0.*nns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1).^2 - (nb0.*nns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1);
        end
        % dH
        dHN2 = R*T_k^2/(p_k*yN2_k) *(-dqN2_dTequ_k)/(dqN2_dpN2_k);
        % N2 desorbed
        dNN2_solid = NN2_solid - NN2_solid_k;

        % CH4
        switch CH4IsothermModel
            case 'toth_cp'
                dqCH4_dTequ_k = - cb0_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)).*(((cT0.*calpha_c.*log(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)))./T_k.^2 - (cb0_c.*cdH_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*(ct0_c - calpha_c.*(cT0./T_k - 1)).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0_c - calpha_c.*(cT0./T_k - 1)).*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha_c.*log((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0_c - calpha_c.*(cT0./T_k - 1)).^2.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1))))) - cb0_p.*cns0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)).*(((cT0.*calpha_p.*log(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)))./T_k.^2 - (cb0_p.*cdH_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*(ct0_p - calpha_p.*(cT0./T_k - 1)).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0_p - calpha_p.*(cT0./T_k - 1)).*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha_p.*log((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0_p - calpha_p.*(cT0./T_k - 1)).^2.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))))) - (cXi_c.*cb0_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./(cT0.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)))) - (cXi_p.*cb0_p.*cns0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./(cT0.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)))) - (cb0_c.*cdH_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)))) - (cb0_p.*cdH_p.*cns0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))));
                dqCH4_dpCH4_k = (cb0_c.*cns0_c.*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1))) + (cb0_p.*cns0_p.*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))) - (cb0_c.^2.*cns0_c.*(p_k.*yCH4_k).*exp((2.*cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1) - 1))./((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1) - (cb0_p.^2.*cns0_p.*(p_k.*yCH4_k).*exp((2.*cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1) - 1))./((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1);
            case 's_shaped'
                dqCH4_dTequ_k = cgam.*((exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1).*((cdH_step.*exp(-cxi_2.*(1./cT0 - 1./T_k)))./(R.*T_k.^2.*cxi_1) + (cxi_2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./(T_k.^2.*cxi_1)))./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1) - (exp(-(2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1).*((cdH_step.*exp(-cxi_2.*(1./cT0 - 1./T_k)))./(R.*T_k.^2.*cxi_1) + (cxi_2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./(T_k.^2.*cxi_1)))./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1).^2).*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^(cgam - 1).*(cb_H0.*(p_k.*yCH4_k).*exp(cdU_H./(R.*T_k)) + (cb_U0.*(p_k.*yCH4_k).*cq_U0.*exp(cdU_U./(R.*T_k)))./(cb_U0.*(p_k.*yCH4_k).*exp(cdU_U./(R.*T_k)) + 1)) - (exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^cgam.*((cb_H0.*cdU_H.*(p_k.*yCH4_k).*exp(cdU_H./(R.*T_k)))./(R.*T_k.^2) - (cb_U0.^2.*cdU_U.*(p_k.*yCH4_k).^2.*cq_U0.*exp((2.*cdU_U)./(R.*T_k)))./(R.*T_k.^2.*(cb_U0.*(p_k.*yCH4_k).*exp(cdU_U./(R.*T_k)) + 1).^2) + (cb_U0.*cdU_U.*(p_k.*yCH4_k).*cq_U0.*exp(cdU_U./(R.*T_k)))./(R.*T_k.^2.*(cb_U0.*(p_k.*yCH4_k).*exp(cdU_U./(R.*T_k)) + 1))) - (cb_L0.*cgam.*(p_k.*yCH4_k).*cq_L0.*exp(cdU_L./(R.*T_k)).*((exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1).*((cdH_step.*exp(-cxi_2.*(1./cT0 - 1./T_k)))./(R.*T_k.^2.*cxi_1) + (cxi_2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./(T_k.^2.*cxi_1)))./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1) - (exp(-(2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1).*((cdH_step.*exp(-cxi_2.*(1./cT0 - 1./T_k)))./(R.*T_k.^2.*cxi_1) + (cxi_2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./(T_k.^2.*cxi_1)))./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1).^2).*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^(cgam - 1))./(cb_L0.*(p_k.*yCH4_k).*exp(cdU_L./(R.*T_k)) + 1) - (cb_L0.^2.*cdU_L.*(p_k.*yCH4_k).^2.*cq_L0.*exp((2.*cdU_L)./(R.*T_k)).*((exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^cgam - 1))./(R.*T_k.^2.*(cb_L0.*(p_k.*yCH4_k).*exp(cdU_L./(R.*T_k)) + 1).^2) + (cb_L0.*cdU_L.*(p_k.*yCH4_k).*cq_L0.*exp(cdU_L./(R.*T_k)).*((exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^cgam - 1))./(R.*T_k.^2.*(cb_L0.*(p_k.*yCH4_k).*exp(cdU_L./(R.*T_k)) + 1));
                dqCH4_dpCH4_k =  (exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^cgam.*(cb_H0.*exp(cdU_H./(R.*T_k)) + (cb_U0.*cq_U0.*exp(cdU_U./(R.*T_k)))./(cb_U0.*(p_k.*yCH4_k).*exp(cdU_U./(R.*T_k)) + 1) - (cb_U0.^2.*(p_k.*yCH4_k).*cq_U0.*exp((2.*cdU_U)./(R.*T_k)))./(cb_U0.*(p_k.*yCH4_k).*exp(cdU_U./(R.*T_k)) + 1).^2) + cgam.*((exp(-cxi_2.*(1./cT0 - 1./T_k)).*exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1))./((p_k.*yCH4_k).*cxi_1.*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)) - (exp(-cxi_2.*(1./cT0 - 1./T_k)).*exp(-(2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1))./((p_k.*yCH4_k).*cxi_1.*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1).^2)).*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^(cgam - 1).*(cb_H0.*(p_k.*yCH4_k).*exp(cdU_H./(R.*T_k)) + (cb_U0.*(p_k.*yCH4_k).*cq_U0.*exp(cdU_U./(R.*T_k)))./(cb_U0.*(p_k.*yCH4_k).*exp(cdU_U./(R.*T_k)) + 1)) - (cb_L0.*cq_L0.*exp(cdU_L./(R.*T_k)).*((exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^cgam - 1))./(cb_L0.*(p_k.*yCH4_k).*exp(cdU_L./(R.*T_k)) + 1) + (cb_L0.^2.*(p_k.*yCH4_k).*cq_L0.*exp((2.*cdU_L)./(R.*T_k)).*((exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^cgam - 1))./(cb_L0.*(p_k.*yCH4_k).*exp(cdU_L./(R.*T_k)) + 1).^2 - (cb_L0.*cgam.*(p_k.*yCH4_k).*cq_L0.*exp(cdU_L./(R.*T_k)).*((exp(-cxi_2.*(1./cT0 - 1./T_k)).*exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1))./((p_k.*yCH4_k).*cxi_1.*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)) - (exp(-cxi_2.*(1./cT0 - 1./T_k)).*exp(-(2.*exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1))./((p_k.*yCH4_k).*cxi_1.*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1).^2)).*(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1)./(exp(-(exp(-cxi_2.*(1./cT0 - 1./T_k)).*(log(cp_step0.*exp(-(cdH_step.*(1./cT0 - 1./T_k))./R)) - log((p_k.*yCH4_k))))./cxi_1) + 1)).^(cgam - 1))./(cb_L0.*(p_k.*yCH4_k).*exp(cdU_L./(R.*T_k)) + 1);
            case 'DSL'
                dqCH4_dTequ_k = (1000000.*cb0.*cn1.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)).*((1000000.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k.^2) + (1000000.*cHb.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000.*cd0.*cn2.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k.^2.*((1000000.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*cHb.*cb0.*cn1.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*cHd.*cd0.*cn2.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*cb0.*cn1.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k.^2.*((1000000.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*cd0.*cn2.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)).*((1000000.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k.^2) + (1000000.*cHd.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k) + 1).^2);
                dqCH4_dpCH4_k = (1000000.*cb0.*cn1.*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k.*((1000000.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*cd0.*cn2.*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k.*((1000000.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000000000.*cb0.^2.*cn1.*(yCH4_k.*p_k).*yCH4_k.^2.*exp((2.*cHb)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*cb0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000000000.*cd0.^2.*cn2.*(yCH4_k.*p_k).*yCH4_k.^2.*exp((2.*cHd)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*cd0.*(yCH4_k.*p_k).*yCH4_k.*exp(cHd./(R.*T_k)))./(R.*T_k) + 1).^2);
            case 'DSL2'
                dqCH4_dTequ_k = (cHb.*cb0.^2.*cn1.*(p_k.*yCH4_k).^2.*exp((2.*cHb)./(R.*T_k)))./(R.*T_k.^2.*(cb0.*(p_k.*yCH4_k).*exp(cHb./(R.*T_k)) + 1).^2) - (cHd.*cd0.*cn2.*(p_k.*yCH4_k).*exp(cHd./(R.*T_k)))./(R.*T_k.^2.*(cd0.*(p_k.*yCH4_k).*exp(cHd./(R.*T_k)) + 1)) - (cHb.*cb0.*cn1.*(p_k.*yCH4_k).*exp(cHb./(R.*T_k)))./(R.*T_k.^2.*(cb0.*(p_k.*yCH4_k).*exp(cHb./(R.*T_k)) + 1)) + (cHd.*cd0.^2.*cn2.*(p_k.*yCH4_k).^2.*exp((2.*cHd)./(R.*T_k)))./(R.*T_k.^2.*(cd0.*(p_k.*yCH4_k).*exp(cHd./(R.*T_k)) + 1).^2);
                dqCH4_dpCH4_k = (cb0.*cn1.*exp(cHb./(R.*T_k)))./(cb0.*(p_k.*yCH4_k).*exp(cHb./(R.*T_k)) + 1) + (cd0.*cn2.*exp(cHd./(R.*T_k)))./(cd0.*(p_k.*yCH4_k).*exp(cHd./(R.*T_k)) + 1) - (cb0.^2.*cn1.*(p_k.*yCH4_k).*exp((2.*cHb)./(R.*T_k)))./(cb0.*(p_k.*yCH4_k).*exp(cHb./(R.*T_k)) + 1).^2 - (cd0.^2.*cn2.*(p_k.*yCH4_k).*exp((2.*cHd)./(R.*T_k)))./(cd0.*(p_k.*yCH4_k).*exp(cHd./(R.*T_k)) + 1).^2;
            case 'toth'
                dqCH4_dTequ_k = - cb0.*cns0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 - (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(ct0 - calpha.*(cT0./T_k - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0 - calpha.*(cT0./T_k - 1)).*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha.*log((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0 - calpha.*(cT0./T_k - 1)).^2.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))))) - (cXi.*cb0.*cns0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./(cT0.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)))) - (cb0.*cdH.*cns0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))));
                dqCH4_dpCH4_k = (cb0.*cns0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))) - (cb0.^2.*cns0.*(p_k.*yCH4_k).*exp((2.*cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)) + 1);
            case 'langfr'
                dqCH4_dTequ_k = (cns0.*exp(-cXi.*(T_k./cT0 - 1)).*((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 + (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1) - (cns0.*exp(-cXi.*(T_k./cT0 - 1)).*((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 + (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1).^2 - (cXi.*cns0.*exp(-cXi.*(T_k./cT0 - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./(cT0.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1));
                dqCH4_dpCH4_k = (cb0.*cns0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1).^2 - (cb0.*cns0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1);
        end
        % dH
        dHCH4 = R*T_k^2/(p_k*yCH4_k) *(-dqCH4_dTequ_k)/(dqCH4_dpCH4_k);
        % CH4 desorbed
        dNCH4_solid = NCH4_solid - NCH4_solid_k;





        %% overall  Energy Balance
        Q = (adsorbentMass*cp*(T-T_k) - dHCO2*dNCO2_solid  - dHN2*dNN2_solid - dHCH4*dNCH4_solid)/1000; %kJ

        out = [Q, NCO2, NN2, NCH4, qCO2];
    end

    function fvec = fcn_dry(x)
        % x(1): yCO2; x(2): yN2; x(3): Nin; x(4): yCH4

        %% equation 1: material balance CO2
        % CO2
        % solid phase CO2, initial condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2_k = ((ns0_c*exp(Xi_c*(1-T_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./T_k))).^(1./(t0_c+alpha_c*(1-T0./T_k)))))+((ns0_p*exp(Xi_p*(1-T_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./T_k))).^(1./(t0_p+alpha_p*(1-T0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
            case 'DSL'
                qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T));
            case 'DSL2'
                qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k));
            case 'toth'
                qCO2_k = ((ns0*exp(Xi*(1-T_k/T0)).*(b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0+alpha*(1-T0./T_k))).^(1./(t0+alpha*(1-T0./T_k)))));
            case 'langfr'
                qCO2_k = (ns0.*exp(Xi.*(1-T_k./T0))).*((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k)));
        end
        NCO2_solid_k = adsorbentMass*qCO2_k;
        % fluid phase CO2, initial condition
        NCO2_fluid_k = yCO2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CO2, initial condition
        NCO2_k = NCO2_solid_k + NCO2_fluid_k;

        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-T/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T))).^(1./(t0_p+alpha_p*(1-T0./T))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(x(1)*p)./(1+(b_L0.*exp(dU_L./(R*T))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(x(1)*p)./(1+(b_U0.*exp(dU_U./(R*T))).*(x(1)*p))+(b_H0.*exp(dU_H./(R*T))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); % adsorbed amount of CO2
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(x(1)*p)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(x(1)*p)./(1e-6*R*T));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1+b0*exp(Hb/R/T).*(x(1)*p)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1+d0*exp(Hd/R/T).*(x(1)*p));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T)));
        end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid = x(1)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid CO2, final condition
        NCO2 = NCO2_solid + NCO2_fluid;

        % material balance CO2
        fvec(1) = - NCO2 + yCO2_feed*x(3) + NCO2_k; % == 0



        %% CH4
        % solid phase CH4, initial condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4_k = ((cns0_c*exp(cXi_c*(1-T_k/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_c+calpha_c*(1-cT0./T_k))).^(1./(ct0_c+calpha_c*(1-cT0./T_k)))))+((cns0_p*exp(cXi_p*(1-T_k/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_p+calpha_p*(1-cT0./T_k))).^(1./(ct0_p+calpha_p*(1-cT0./T_k))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4_k = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k))).*(1-((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k))+(cb_H0.*exp(cdU_H./(R*T_k))).*(yCH4_k*p_k)).*((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam);
            case 'DSL'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T));
            case 'DSL2'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k));
            case 'toth'
                qCH4_k = ((cns0*exp(cXi*(1-T_k/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0+calpha*(1-cT0./T_k))).^(1./(ct0+calpha*(1-cT0./T_k)))));
            case 'langfr'
                qCH4_k = (cns0.*exp(cXi.*(1-T_k./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k)));
        end
        NCH4_solid_k = adsorbentMass*qCH4_k;
        % fluid phase CH4, initial condition
        NCH4_fluid_k = yCH4_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CH4, initial condition
        NCH4_k = NCH4_solid_k + NCH4_fluid_k;

        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p*x(4)))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p*x(4))).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-T/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T-1))).*(p*x(4)))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T-1))).*(p*x(4))).^(ct0_p+calpha_p*(1-cT0./T))).^(1./(ct0_p+calpha_p*(1-cT0./T))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(x(4)*p)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(x(4)*p))).*(1-((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((x(4)*p)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(x(4)*p)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(x(4)*p))+(cb_H0.*exp(cdU_H./(R*T))).*(x(4)*p)).*((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((x(4)*p)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam); % adsorbed amount of CH4
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(x(4)*p)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(x(4)*p)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(x(4)*p)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(x(4)*p)./(1e-6*R*T));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(x(4)*p)./(1+cb0*exp(cHb/R/T).*(x(4)*p)) + cn2*cd0*exp(cHd/R/T)*(x(4)*p)./(1+cd0*exp(cHd/R/T).*(x(4)*p));
            case 'toth'
                qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(x(4)*p))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(x(4)*p)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
            case 'langfr'
                qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T)));
        end
        NCH4_solid = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid = x(4)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid CH4, final condition
        NCH4 = NCH4_solid + NCH4_fluid;

        % material balance CH4
        fvec(4) = - NCH4 + yCH4_feed*x(3) + NCH4_k; % == 0

        %% equation 2 material balance N2

        %N2
        % solid phase N2, initial condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2_k = ((nns0_c*exp(nXi_c*(1-T_k/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_c+nalpha_c*(1-nT0./T_k))).^(1./(nt0_c+nalpha_c*(1-nT0./T_k)))))+((nns0_p*exp(nXi_p*(1-T_k/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_p+nalpha_p*(1-nT0./T_k))).^(1./(nt0_p+nalpha_p*(1-nT0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2_k = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T_k))).*(yN2_k*p_k)./(1+(nb_L0.*exp(ndU_L./(R*T_k))).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T_k))).*(yN2_k*p_k)./(1+(nb_U0.*exp(ndU_U./(R*T_k))).*(yN2_k*p_k))+(nb_H0.*exp(ndU_H./(R*T_k))).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T_k))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T_k)))))).^ngam);
            case 'DSL'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k)./(1e-6*R*T));
            case 'DSL2'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k));
            case 'toth'
                qN2_k = ((nns0*exp(nXi*(1-T_k/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0+nalpha*(1-nT0./T_k))).^(1./(nt0+nalpha*(1-nT0./T_k)))));
            case 'langfr'
                qN2_k = (nns0.*exp(nXi.*(1-T_k./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k)));
        end
        NN2_solid_k = adsorbentMass*qN2_k;
        % fluid phase N2, initial condition
        NN2_fluid_k = yN2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid N2, initial condition
        NN2_k = NN2_solid_k + NN2_fluid_k;

        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p*x(2)))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p*x(2))).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-T/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T-1))).*(p*x(2)))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T-1))).*(p*x(2))).^(nt0_p+nalpha_p*(1-nT0./T))).^(1./(nt0_p+nalpha_p*(1-nT0./T))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*(x(2)*p)./(1+(nb_L0.*exp(ndU_L./(R*T))).*(x(2)*p))).*(1-((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((x(2)*p)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*(x(2)*p)./(1+(nb_U0.*exp(ndU_U./(R*T))).*(x(2)*p))+(nb_H0.*exp(ndU_H./(R*T))).*(x(2)*p)).*((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((x(2)*p)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam); % adsorbed amount of CO2
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*(x(2)*p)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*(x(2)*p)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*(x(2)*p)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*(x(2)*p)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*(x(2)*p)./(1+nb0*exp(nHb/R/T).*(x(2)*p)) + nn2*nd0*exp(nHd/R/T)*(x(2)*p)./(1+nd0*exp(nHd/R/T).*(x(2)*p));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(x(2)*p))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(x(2)*p)).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = x(2)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid N2, final condition
        NN2 = NN2_solid + NN2_fluid;


        % material balance N2
%         fvec(2) =  - (NN2 + NCO2 + NCH4) + x(3) + (NCO2_k + NN2_k + NCH4_k); % == 0
                fvec(2) =  - NN2 + yN2_feed*x(3) +  NN2_k;

        %% equation 4: overall material balance
        fvec(3) = 1-x(1)-x(2)-x(4); % == 0

        fvec = real(fvec);
    end

end