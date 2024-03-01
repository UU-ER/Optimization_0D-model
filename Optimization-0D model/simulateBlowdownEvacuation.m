
%% BLOWDOWN AND REGENERATION

% output: pressure, molefraction CO2, solid loadings of CO2 and H2O,
% total moles desorbed, moles of CO2 and H2O desorbed and the energy consumption
% Matrix structure: [Pressure MoleFracCO2 SolidLoadingCO2 SolidLoadingH2O
% MolesDesorbedTotal MolesDesorbedCO2 MolesDesorbedH2O EnergyConsumption]

% from 1D model -> evacuation 1/4 and regeneration 3/4 of the time -> set
% to 1/2 and 1/2

function outputBlowEvac = simulateBlowdownEvacuation(data)

%% define global variables
% global T_reg adsorbentMass void V R T0 Xi_c dH_c alpha_c Xi_p dH_p alpha_p ns0_c b0_c t0_c ns0_p b0_p t0_p CG0 HC K0 HK Cm0 beta yCO2_k T_k p_k p yN2_k yH2O_k cp

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


%% N2 Isotherm

N2IsothermModel= data.sorbent(i).N2Isotherm.model;

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
        ns0 = data.sorbent(i).N2Isotherm.param(4); % [mol/kg]
        nb0 = data.sorbent(i).N2Isotherm.param(5); % [1/MPa]
        nt0 = data.sorbent(i).N2Isotherm.param(6);
    case 'langfr'
        % N2 isotherm
        nT0 = data.sorbent(i).N2Isotherm.T0(1);
        ns0 = data.sorbent(i).N2Isotherm.param(1);
        nXi = data.sorbent(i).N2Isotherm.param(2);
        nt0 = data.sorbent(i).N2Isotherm.param(3);
        nalpha = data.sorbent(i).N2Isotherm.param(4); % [mol/kg]
        nb0 = data.sorbent(i).N2Isotherm.param(5); % [1/MPa]
        ndH = data.sorbent(i).N2Isotherm.param(6); % [J/mol]
end

%% CH4

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
        cs0_p = data.sorbent(i).CH4Isotherm.param(10); % [mol/kg]
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
        cn2 = data.sorbent(i).CH4Isotherm.param(4); % CH4 (mol/kg)
        cd0 = data.sorbent(i).CH4Isotherm.param(5)/(1e6); % d0 (m3/mol)
        cHd = data.sorbent(i).CH4Isotherm.param(6); % Hd (J/mol)
    case 'DSL2'
        cT0 = data.sorbent(i).CH4Isotherm.T0; % K
        cn1 = data.sorbent(i).CH4Isotherm.param(1); % n1 (mol/kg)
        cb0 = data.sorbent(i).CH4Isotherm.param(2); % b0 (m3/mol)
        cHb = data.sorbent(i).CH4Isotherm.param(3); % Hb (J/mol)
        cn2 = data.sorbent(i).CH4Isotherm.param(4); % CH4 (mol/kg)
        cd0 = data.sorbent(i).CH4Isotherm.param(5); % d0 (m3/mol)
        cHd = data.sorbent(i).CH4Isotherm.param(6); % Hd (J/mol)
    case 'toth'
        % CH4 isotherm
        cT0 = data.sorbent(i).CH4Isotherm.T0(1);
        cXi = data.sorbent(i).CH4Isotherm.param(1);
        cdH = data.sorbent(i).CH4Isotherm.param(2); % [J/mol]
        calpha = data.sorbent(i).CH4Isotherm.param(3);
        cs0 = data.sorbent(i).CH4Isotherm.param(4); % [mol/kg]
        cb0 = data.sorbent(i).CH4Isotherm.param(5); % [1/MPa]
        ct0 = data.sorbent(i).CH4Isotherm.param(6);
    case 'langfr'
        % CH4 isotherm
        cT0 = data.sorbent(i).CH4Isotherm.T0(1);
        cs0 = data.sorbent(i).CH4Isotherm.param(1);
        cXi = data.sorbent(i).CH4Isotherm.param(2);
        ct0 = data.sorbent(i).CH4Isotherm.param(3);
        calpha = data.sorbent(i).CH4Isotherm.param(4); % [mol/kg]
        cb0 = data.sorbent(i).CH4Isotherm.param(5); % [1/MPa]
        cdH = data.sorbent(i).CH4Isotherm.param(6); % [J/mol]
end

%% define general dataModel
adsorbentMass = data.process.adsorbentMass;
cp = data.sorbent(i).cp; % J/Kg/K
R = data.general.gasconstant;
void = data.process.voidFraction;
V = data.process.Vol;
%% BLOWDOWN WITH PRESSURE VECTOR
startingCondBD = data.startingCondBD;

for m =1:size(data.pressureVector,1)

    % starting conditions
    % 100% saturation
    % x0 (1): yCO2, (2): yN2, (3): Nout, (4): T, (5): yCH4
    x0 = [startingCondBD(m,1),startingCondBD(m,2),0,data.process.Tamb,startingCondBD(m,3)];


    yCO2_k = x0(1);
    yN2_k = x0(2);
    yCH4_k = x0(5);
    T_k = data.process.Tads;
    p_k = data.process.pamb;

    pressureVector = data.pressureVector(m,:);
    stepsBD = data.process.noSteps;
    stepsHeat = data.process.noSteps;

    % create vectors for saving results
    factorTime = 2; % from 1D model -> evacuation 1/4 and regeneration 3/4 of the time
    length = stepsBD+stepsHeat;
    yCO2_save = nan(length,1);
    yN2_save = nan(length,1);
    yCH4_save = nan(length,1);
    Nout_save = nan(length,1);
    T_save = nan(length,1);
    P_save = nan(length,1);
    x_save = nan(length,5);
    dNCO2_solid = zeros(length,1);
    dNN2_solid=zeros(length,1);
    dNCH4_solid=zeros(length,1);
    E_total = zeros(length,1);
    Nout_save_sum = zeros(length,1);
    Nout_save_CO2 = zeros(length,1);
    Nout_save_N2 = zeros(length,1);
    Nout_save_CH4 = zeros(length,1);
    Q_save_sum = zeros(length,1);

    % Boundries (1): yCO2, (2): yN2, (3): Nout, (4): T, (5): yCH4
    % double check it
    lb = [data.startingCondBD(m,1), 1e-5, 0, data.process.Tads, data.startingCondBD(m,3)];
    ub = [1,1,50,data.process.Tamb,1];

    % solver options
    options = optimoptions('lsqnonlin','Display','off');
    options.Algorithm = 'levenberg-marquardt'; % trust-region-reflective, levenberg-marquardt
    options.OptimalityTolerance = 1e-9;
    options.FunctionTolerance = 1e-9;
    options.StepTolerance = 1e-9;

    for k = 1:stepsBD  %data.process.noSteps*(factorTime-1)
        p = pressureVector(k);

        funct = @(x) double( fcn_BD(x));

        fun = @(x) funct(x);

        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

        yCO2_k = x(1);
        yN2_k = x(2);
        T_k = x(4);
        yCH4_k = x(5);

        x0 = x;

        x_save(k,:) = [residual];
        yCO2_save(k) = x(1);
        yN2_save(k) = x(2);
        yCH4_save(k) = x(5);
        Nout_save(k) = x(3);
        T_save(k) = x(4);
        P_save(k) = p;

        Nout_save_sum(k+1) = Nout_save_sum(k) + Nout_save(k); % mol
        Nout_save_CO2(k+1) = Nout_save_CO2(k) + Nout_save(k);
        Nout_save_N2(k+1) = Nout_save_N2(k) + Nout_save(k)*x(2);
        Nout_save_CH4(k+1) = Nout_save_CH4(k) + Nout_save(k)*x(5);

        molesDesorbed_Total = Nout_save_sum(k+1)-Nout_save_sum(k); % moles desorbed during this step

        % Energy consumption till that given step

        E(1) = cmp_W_vac2(p*10,p_k*10,x(4),[0,(Nout_save(k)*x(1)),(Nout_save(k)*x(2)),(Nout_save(k)*x(5))],x); % KJ

        E_total(k+1,1) = E_total(k,1) + E(1); % J

        p_k = pressureVector(k);

    end

    %% REGENERATION WITH HEAT VECTOR-yCO2_save(stepsBD)

    % starting conditions
    x0 = [yCO2_save(stepsBD),yN2_save(stepsBD),Nout_save(stepsBD),yCH4_save(stepsBD)];

    p = pressureVector(end); %% take last pressure from pressure step
    p_k = p; %% pressure does not change
    T_k = T_save(stepsBD); %

    % new: include boundaries for solver:
    lb = [0, 0, 0, 0];
    ub = [1, 1, 50, 1];

    % solver options
    options = optimoptions('lsqnonlin','Display','off');
    options.Algorithm = 'trust-region-reflective';
    options.OptimalityTolerance = 1e-10;
    options.FunctionTolerance = 1e-10;
    options.StepTolerance = 1e-10;

    TempVector = data.THeatProfileStep(m,:); %

    for j = 1:stepsHeat %

        T_reg = TempVector(j);

        funct = @(x) double( fcn_Evacuation(x));
        fun = @(x) funct(x);

        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

        x0 = x;

        x_save(k+j,:) = [residual,0];
        yCO2_save(k+j) = x(1);
        yN2_save(k+j) = x(2);
        Nout_save(k+j) = x(3);
        yCH4_save(k+j) = x(4);

        % test
        lb = [0, 0, 0, 0];
        ub = [1, x(2), 50, 1];

        Nout_save_sum(k+j) = Nout_save_sum(k+j-1) + Nout_save(k+j);
        Nout_save_CO2(k+j) = Nout_save_CO2(k+j-1) + Nout_save(k+j)*x(1);
        Nout_save_N2(k+j) = Nout_save_N2(k+j-1) + Nout_save(k+j)*x(2);
        Nout_save_CH4(k+j) = Nout_save_CH4(k+j-1) + Nout_save(k+j)*x(4);

        molesDesorbed_Total = Nout_save_CO2(k+j)-Nout_save_CO2(k+j-1); % moles desorbed during this step

        if T_reg >= data.process.Tamb
            Q = Qex(x,data.process.adsorbentMass,data.sorbent(i).cp,data.general.gasconstant);
            Q_save_sum(k+j) =  Q_save_sum(k+j-1) + Q; % kJ
        else
            Q = 0;
            Q_save_sum(k+j) =  Q_save_sum(k+j-1) + Q; % kJ
        end


        yCO2_k = x(1);
        yN2_k = x(2);
        yCH4_k = x(4);
        T_k = TempVector(j);
        T_save(k+j) = TempVector(j);
        P_save(k+j) = p;

        % Energy consumption till that given step
        E_total(k+1+j,1) = E_total(k+1,1) + 0; % J vacuum pump CO2 prod
    end
    if data.plot
        plot_all = 'step_1';
        separate_plot_file_new2;
    end
%     outputBlowEvac = [yCO2_save, yN2_save,yCH4_save, Nout_save_sum(2:end), Nout_save_CO2(2:end), Nout_save_N2(2:end), Nout_save_CH4(2:end), E_total(2:end), T_save, Q_save_sum(2:end)];
    %  outputBlowEvac = [yCO2_save, yN2_save, yH2O_save, yCH4_save, Nout_save_sum(2:end), Nout_save_CO2(2:end), Nout_save_N2(2:end), Nout_save_H2O(2:end), Nout_save_CH4(2:end), E_total(2:end), T_save, Q_save_sum(2:end)];
    outputBlowEvac(m).yCO2 = yCO2_save;
    outputBlowEvac(m).yN2 = yN2_save;
    outputBlowEvac(m).yCH4 = yCH4_save;
    outputBlowEvac(m).Nout_sum = Nout_save_sum;
    outputBlowEvac(m).Nout_CO2 = Nout_save_CO2;
    outputBlowEvac(m).Nout_N2 = Nout_save_N2;
    outputBlowEvac(m).Nout_CH4 = Nout_save_CH4;
    outputBlowEvac(m).E_total = E_total; % vacuum pump work= Exergy (KJ)
    outputBlowEvac(m).T = T_save;
    outputBlowEvac(m).P = P_save;
    outputBlowEvac(m).Q_sum = Q_save_sum;
end

%% functions
    function fvec = fcn_BD(x)
        % x(1): yCO2; x(2): yN2; x(3): Nout; x(4): T, x(5):yCH4

        %% equation 1: material balance CO2
        % solid phase CO2, initial condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2_k = ((ns0_c*exp(Xi_c*(1-T_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./T_k))).^(1./(t0_c+alpha_c*(1-T0./T_k)))))+((ns0_p*exp(Xi_p*(1-T_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./T_k))).^(1./(t0_p+alpha_p*(1-T0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                %                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
                qCO2_k = (q_L0.*b_L0.*exp(dU_L/(R*T_k)).*(yCO2_k*p_k)./(1+b_L0.*exp(dU_L./(R*T_k)).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step/R.*(1/T0-1/T_k))))./(xi_1*exp(xi_2.*(1/T0-1/T_k))))./(1+exp((log((yCO2_k*p_k))-log(p_step0*exp(-dH_step/R.*(1/T0-1/T_k))))./(xi_1*exp(xi_2*(1/T0-1/T_k)))))).^gam))+(q_L0*b_U0.*exp(dU_U./(R*T_k)).*(yCO2_k*p_k)./(1+b_U0*exp(dU_U/(R*T_k)).*(yCO2_k*p_k))+b_H0*exp(dU_H/(R*T_k)).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0*exp(-dH_step/R*(1/T0-1/T_k))))./(xi_1*exp(xi_2*(1/T0-1/T_k))))./(1+exp((log((yCO2_k*p_k))-log(p_step0*exp(-dH_step/R*(1/T0-1/T_k))))./(xi_1*exp(xi_2*(1/T0-1/T_k)))))).^gam);
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
        % fluid phase CO2, initial condition
        NCO2_fluid_k = yCO2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CO2, initial condition
        NCO2_k = NCO2_solid_k + NCO2_fluid_k;

        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-x(4)/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./x(4)-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./x(4)-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./x(4)))).^(1./(t0_c+alpha_c*(1-T0./x(4))))))+((ns0_p*exp(Xi_p*(1-x(4)/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./x(4)-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./x(4)-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./x(4)))).^(1./(t0_p+alpha_p*(1-T0./x(4)))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*b_L0.*exp(dU_L/(R*x(4))).*(x(1)*p)./(1+b_L0.*exp(dU_L./(R*x(4))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step/R.*(1/T0-1/x(4)))))./(xi_1*exp(xi_2.*(1/T0-1/x(4)))))./(1+exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R.*(1/T0-1/x(4)))))./(xi_1*exp(xi_2*(1/T0-1/x(4))))))).^gam))+(q_L0*b_U0.*exp(dU_U./(R*x(4))).*(x(1)*p)./(1+b_U0*exp(dU_U/(R*x(4))).*(x(1)*p))+b_H0*exp(dU_H/(R*x(4))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R*(1/T0-1/x(4)))))./(xi_1*exp(xi_2*(1/T0-1/x(4)))))./(1+exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R*(1/T0-1/x(4)))))./(xi_1*exp(xi_2*(1/T0-1/x(4))))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/x(4))*(x(1)*p)./(1e-6*R*x(4))./(1+b0*exp(Hb/R/x(4)).*(x(1)*p)./(1e-6*R*x(4))) + n2*d0*exp(Hd/R/x(4))*(x(1)*p)./(1e-6*R*x(4))./(1+d0*exp(Hd/R/x(4)).*(x(1)*p)./(1e-6*R*x(4)));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/x(4))*(x(1)*p)./(1+b0*exp(Hb/R/x(4)).*(x(1)*p)) + n2*d0*exp(Hd/R/x(4))*(x(1)*p)./(1+d0*exp(Hd/R/x(4)).*(x(1)*p));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-x(4)/T0)).*(b0*exp(dH/(R*T0)*(T0./x(4)-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./x(4)-1))).*(x(1)*p)).^(t0+alpha*(1-T0./x(4)))).^(1./(t0+alpha*(1-T0./x(4))))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-x(4)./T0))).*((b0*exp(dH./(R*T0).*(T0./x(4)-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./x(4)))./(1+ ((b0*exp(dH./(R*T0).*(T0./x(4)-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./x(4))));
        end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid = x(1)*p*V*void/(R*x(4))*1e6; % mol
        % solid and liquid CO2, final condition
        NCO2 = NCO2_solid + NCO2_fluid;

        % material balance CO2
        fvec(1) = NCO2_k - x(1)*x(3) - NCO2; % == 0



        %% equation ?: material balance CH4
        % solid phase CH4, initial condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4_k = ((cns0_c*exp(cXi_c*(1-T_k/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_c+calpha_c*(1-cT0./T_k))).^(1./(ct0_c+calpha_c*(1-cT0./T_k)))))+((cs0_p*exp(cXi_p*(1-T_k/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_p+calpha_p*(1-cT0./T_k))).^(1./(ct0_p+calpha_p*(1-cT0./T_k))))); % adsorbed amount of CH4
            case 's_shaped'
                %                    qCH4_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCH4_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCH4_k*p_k))).*(1-((exp((log((yCH4_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCH4_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCH4_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCH4_k*p_k)).*((exp((log((yCH4_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
                qCH4_k = (cq_L0.*cb_L0.*exp(cdU_L/(R*T_k)).*(yCH4_k*p_k)./(1+cb_L0.*exp(cdU_L./(R*T_k)).*(yCH4_k*p_k))).*(1-((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step/R.*(1/T0-1/T_k))))./(cxi_1*exp(cxi_2.*(1/T0-1/T_k))))./(1+exp((log((yCH4_k*p_k))-log(cp_step0*exp(-cdH_step/R.*(1/T0-1/T_k))))./(cxi_1*exp(cxi_2*(1/T0-1/T_k)))))).^cgam))+(cq_L0*cb_U0.*exp(cdU_U./(R*T_k)).*(yCH4_k*p_k)./(1+cb_U0*exp(cdU_U/(R*T_k)).*(yCH4_k*p_k))+cb_H0*exp(cdU_H/(R*T_k)).*(yCH4_k*p_k)).*((exp((log((yCH4_k*p_k))-log(cp_step0*exp(-cdH_step/R*(1/T0-1/T_k))))./(cxi_1*exp(cxi_2*(1/T0-1/T_k))))./(1+exp((log((yCH4_k*p_k))-log(cp_step0*exp(-cdH_step/R*(1/T0-1/T_k))))./(cxi_1*exp(cxi_2*(1/T0-1/T_k)))))).^cgam);
            case 'DSL'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k));
            case 'toth'
                qCH4_k = ((cs0*exp(cXi*(1-T_k/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0+calpha*(1-cT0./T_k))).^(1./(ct0+calpha*(1-cT0./T_k)))));
            case 'langfr'
                qCH4_k = (cs0.*exp(cXi.*(1-T_k./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k)));
        end
        NCH4_solid_k = adsorbentMass*qCH4_k;
        % fluid phase CH4, initial condition
        NCH4_fluid_k = yCH4_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CH4, initial condition
        NCH4_k = NCH4_solid_k + NCH4_fluid_k;

        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-x(4)/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./x(4)-1))).*(p*x(5)))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./x(4)-1))).*(p*x(5))).^(ct0_c+calpha_c*(1-cT0./x(4)))).^(1./(ct0_c+calpha_c*(1-cT0./x(4))))))+((cs0_p*exp(cXi_p*(1-x(4)/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./x(4)-1))).*(p*x(5)))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./x(4)-1))).*(p*x(5))).^(ct0_p+calpha_p*(1-cT0./x(4)))).^(1./(ct0_p+calpha_p*(1-cT0./x(4)))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*cb_L0.*exp(cdU_L/(R*x(4))).*(x(5)*p)./(1+cb_L0.*exp(cdU_L./(R*x(4))).*(x(5)*p))).*(1-((exp((log((x(5)*p))-log(cp_step0.*exp(-cdH_step/R.*(1/T0-1/x(4)))))./(cxi_1*exp(cxi_2.*(1/T0-1/x(4)))))./(1+exp((log((x(5)*p))-log(cp_step0*exp(-cdH_step/R.*(1/T0-1/x(4)))))./(cxi_1*exp(cxi_2*(1/T0-1/x(4))))))).^cgam))+(cq_L0*cb_U0.*exp(cdU_U./(R*x(4))).*(x(5)*p)./(1+cb_U0*exp(cdU_U/(R*x(4))).*(x(5)*p))+cb_H0*exp(cdU_H/(R*x(4))).*(x(5)*p)).*((exp((log((x(5)*p))-log(cp_step0*exp(-cdH_step/R*(1/T0-1/x(4)))))./(cxi_1*exp(cxi_2*(1/T0-1/x(4)))))./(1+exp((log((x(5)*p))-log(cp_step0*exp(-cdH_step/R*(1/T0-1/x(4)))))./(cxi_1*exp(cxi_2*(1/T0-1/x(4))))))).^cgam);
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/x(4))*(x(5)*p)./(1e-6*R*x(4))./(1+cb0*exp(cHb/R/x(4)).*(x(5)*p)./(1e-6*R*x(4))) + cn2*cd0*exp(cHd/R/x(4))*(x(5)*p)./(1e-6*R*x(4))./(1+cd0*exp(cHd/R/x(4)).*(x(5)*p)./(1e-6*R*x(4)));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/x(4))*(x(5)*p)./(1+cb0*exp(cHb/R/x(4)).*(x(5)*p)) + cn2*cd0*exp(cHd/R/x(4))*(x(5)*p)./(1+cd0*exp(cHd/R/x(4)).*(x(5)*p));
            case 'toth'
                qCH4 = ((cs0*exp(cXi*(1-x(4)/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./x(4)-1))).*(x(5)*p))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./x(4)-1))).*(x(5)*p)).^(ct0+calpha*(1-cT0./x(4)))).^(1./(ct0+calpha*(1-cT0./x(4))))));
            case 'langfr'
                qCH4 = (cs0.*exp(cXi.*(1-x(4)./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./x(4)-1))).*(x(5)*p)).^(1./ct0+calpha.*(1-cT0./x(4)))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./x(4)-1))).*(x(5)*p)).^(1./ct0+calpha.*(1-cT0./x(4))));
        end
        NCH4_solid = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid = x(5)*p*V*void/(R*x(4))*1e6; % mol
        % solid and liquid CH4, final condition
        NCH4 = NCH4_solid + NCH4_fluid;

        % material balance CH4
        fvec(2) = NCH4_k - x(5)*x(3) - NCH4; % == 0

        % % %

        %% equation 2: material balance N2
        % solid phase N2, initial condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2_k = ((nns0_c*exp(nXi_c*(1-T_k/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_c+nalpha_c*(1-nT0./T_k))).^(1./(nt0_c+nalpha_c*(1-nT0./T_k)))))+((nns0_p*exp(nXi_p*(1-T_k/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_p+nalpha_p*(1-nT0./T_k))).^(1./(nt0_p+nalpha_p*(1-nT0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                %                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yN2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yN2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yN2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
                qN2_k = (nq_L0.*nb_L0.*exp(ndU_L/(R*T_k)).*(yN2_k*p_k)./(1+nb_L0.*exp(ndU_L./(R*T_k)).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step/R.*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2.*(1/T0-1/T_k))))./(1+exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R.*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k)))))).^ngam))+(nq_L0*nb_U0.*exp(ndU_U./(R*T_k)).*(yN2_k*p_k)./(1+nb_U0*exp(ndU_U/(R*T_k)).*(yN2_k*p_k))+nb_H0*exp(ndU_H/(R*T_k)).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k))))./(1+exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k)))))).^ngam);
            case 'DSL'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k));
            case 'toth'
                qN2_k = ((ns0*exp(nXi*(1-T_k/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0+nalpha*(1-nT0./T_k))).^(1./(nt0+nalpha*(1-nT0./T_k)))));
            case 'langfr'
                qN2_k = (ns0.*exp(nXi.*(1-T_k./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k)));
        end


        NN2_solid_k = adsorbentMass*qN2_k;
        % fluid phase N2, initial condition
        NN2_fluid_k = yN2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid N2, initial condition
        NN2_k = NN2_solid_k + NN2_fluid_k;

        % solid phase N2, final condition

        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-x(4)/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./x(4)-1))).*(p*x(2)))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./x(4)-1))).*(p*x(2))).^(nt0_c+nalpha_c*(1-nT0./x(4)))).^(1./(nt0_c+nalpha_c*(1-nT0./x(4))))))+((nns0_p*exp(nXi_p*(1-x(4)/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./x(4)-1))).*(p*x(2)))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./x(4)-1))).*(p*x(2))).^(nt0_p+nalpha_p*(1-nT0./x(4)))).^(1./(nt0_p+nalpha_p*(1-nT0./x(4)))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*nb_L0.*exp(ndU_L/(R*x(4))).*(x(2)*p)./(1+nb_L0.*exp(ndU_L./(R*x(4))).*(x(2)*p))).*(1-((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step/R.*(1/T0-1/x(4)))))./(nxi_1*exp(nxi_2.*(1/T0-1/x(4)))))./(1+exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R.*(1/T0-1/x(4)))))./(nxi_1*exp(nxi_2*(1/T0-1/x(4))))))).^ngam))+(nq_L0*nb_U0.*exp(ndU_U./(R*x(4))).*(x(2)*p)./(1+nb_U0*exp(ndU_U/(R*x(4))).*(x(2)*p))+nb_H0*exp(ndU_H/(R*x(4))).*(x(2)*p)).*((exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R*(1/T0-1/x(4)))))./(nxi_1*exp(nxi_2*(1/T0-1/x(4)))))./(1+exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R*(1/T0-1/x(4)))))./(nxi_1*exp(nxi_2*(1/T0-1/x(4))))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/x(4))*(x(2)*p)./(1e-6*R*x(4))./(1+nb0*exp(nHb/R/x(4)).*(x(2)*p)./(1e-6*R*x(4))) + nn2*nd0*exp(nHd/R/x(4))*(x(2)*p)./(1e-6*R*x(4))./(1+nd0*exp(nHd/R/x(4)).*(x(2)*p)./(1e-6*R*x(4)));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/x(4))*(x(2)*p)./(1+nb0*exp(nHb/R/x(4)).*(x(2)*p)) + nn2*nd0*exp(nHd/R/x(4))*(x(2)*p)./(1+nd0*exp(nHd/R/x(4)).*(x(2)*p));
            case 'toth'
                qN2 = ((ns0*exp(nXi*(1-x(4)/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./x(4)-1))).*(x(2)*p))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./x(4)-1))).*(x(2)*p)).^(nt0+nalpha*(1-nT0./x(4)))).^(1./(nt0+nalpha*(1-nT0./x(4))))));
            case 'langfr'
                qN2 = (ns0.*exp(nXi.*(1-x(4)./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./x(4)-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./x(4)))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./x(4)-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./x(4))));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = x(2)*p*V*void/(R*x(4))*1e6; % mol
        % solid and liquid CO2, final condition
        NN2 = NN2_solid + NN2_fluid;

        % material balance N2
        fvec(3) = NN2_k - x(2)*x(3) - NN2; % == 0
        %         fvec(3) = (NCO2_k + NN2_k + NCH4_k) - x(3) - (NN2 + NCO2+ NCH4); % == 0

        %% equation 4: overall material balance
        fvec(4) = 1-x(1)-x(2)-x(5); % == 0

        %% equation 5: energy balance
        % CO2
        % dq/dT, inital condition
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
        dNCO2_solid = NCO2_solid_k - NCO2_solid;
        % N2
        % dq/dT, inital condition

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
                dqN2_dTequ_k = - nb0.*ns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 - (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nt0 - nalpha.*(nT0./T_k - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2))./((nt0 - nalpha.*(nT0./T_k - 1)).*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)) + 1)) - (nT0.*nalpha.*log((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1))./(T_k.^2.*(nt0 - nalpha.*(nT0./T_k - 1)).^2.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))))) - (nXi.*nb0.*ns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./(nT0.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)))) - (nb0.*ndH.*ns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./(R.*T_k.^2.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))));
                dqN2_dpN2_k = (nb0.*ns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))) - (nb0.^2.*ns0.*(p_k.*yN2_k).*exp((2.*ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)) + 1);
            case 'langfr'
                dqN2_dTequ_k = (ns0.*exp(-nXi.*(T_k./nT0 - 1)).*((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 + (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1) - (ns0.*exp(-nXi.*(T_k./nT0 - 1)).*((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 + (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1).^2 - (nXi.*ns0.*exp(-nXi.*(T_k./nT0 - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./(nT0.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1));
                dqN2_dpN2_k = (nb0.*ns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1).^2 - (nb0.*ns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1);
        end

        dHN2 = R*T_k^2/(p_k*yN2_k) *(-dqN2_dTequ_k)/(dqN2_dpN2_k);
        % N2 desorbed
        dNN2_solid = NN2_solid_k - NN2_solid;


        % CH4
        % dq/dT, inital condition
        switch CH4IsothermModel
            case 'toth_cp'
                dqCH4_dTequ_k = - cb0_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)).*(((cT0.*calpha_c.*log(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)))./T_k.^2 - (cb0_c.*cdH_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*(ct0_c - calpha_c.*(cT0./T_k - 1)).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0_c - calpha_c.*(cT0./T_k - 1)).*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha_c.*log((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0_c - calpha_c.*(cT0./T_k - 1)).^2.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1))))) - cb0_p.*cs0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)).*(((cT0.*calpha_p.*log(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)))./T_k.^2 - (cb0_p.*cdH_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*(ct0_p - calpha_p.*(cT0./T_k - 1)).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0_p - calpha_p.*(cT0./T_k - 1)).*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha_p.*log((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0_p - calpha_p.*(cT0./T_k - 1)).^2.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))))) - (cXi_c.*cb0_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./(cT0.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)))) - (cXi_p.*cb0_p.*cs0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./(cT0.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)))) - (cb0_c.*cdH_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)))) - (cb0_p.*cdH_p.*cs0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))));
                dqCH4_dpCH4_k = (cb0_c.*cns0_c.*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1))) + (cb0_p.*cs0_p.*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))) - (cb0_c.^2.*cns0_c.*(p_k.*yCH4_k).*exp((2.*cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1) - 1))./((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1) - (cb0_p.^2.*cs0_p.*(p_k.*yCH4_k).*exp((2.*cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1) - 1))./((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1);
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
                dqCH4_dTequ_k = - cb0.*cs0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 - (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(ct0 - calpha.*(cT0./T_k - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0 - calpha.*(cT0./T_k - 1)).*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha.*log((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0 - calpha.*(cT0./T_k - 1)).^2.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))))) - (cXi.*cb0.*cs0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./(cT0.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)))) - (cb0.*cdH.*cs0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))));
                dqCH4_dpCH4_k = (cb0.*cs0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))) - (cb0.^2.*cs0.*(p_k.*yCH4_k).*exp((2.*cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)) + 1);
            case 'langfr'
                dqCH4_dTequ_k = (cs0.*exp(-cXi.*(T_k./cT0 - 1)).*((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 + (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1) - (cs0.*exp(-cXi.*(T_k./cT0 - 1)).*((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 + (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1).^2 - (cXi.*cs0.*exp(-cXi.*(T_k./cT0 - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./(cT0.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1));
                dqCH4_dpCH4_k = (cb0.*cs0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1).^2 - (cb0.*cs0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1);
        end
        % dH
        dHCH4 = R*T_k^2/(p_k*yCH4_k) *(-dqCH4_dTequ_k)/(dqCH4_dpCH4_k);
        % CH4 desorbed
        dNCH4_solid = NCH4_solid_k - NCH4_solid;

        %% overall energy balance

%         fvec(5) = (adsorbentMass*cp*(x(4)-T_k) - dHCO2*dNCO2_solid-dHN2*dNN2_solid-dHCH4*dNCH4_solid);
                fvec(5) = (adsorbentMass*cp*(x(4)-T_k)) ;


        fvec = real(fvec);
    end

    function fvec = fcn_Evacuation(x)
        % x(1): yCO2; x(2): yN2; x(3): Nout, x(4):yCH4

        %% equation 1: material balance CO2
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
        % fluid phase CO2, initial condition
        NCO2_fluid_k = yCO2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CO2, initial condition
        NCO2_k = NCO2_solid_k + NCO2_fluid_k;

        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T_reg/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_reg-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_reg-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T_reg))).^(1./(t0_c+alpha_c*(1-T0./T_reg)))))+((ns0_p*exp(Xi_p*(1-T_reg/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_reg-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_reg-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T_reg))).^(1./(t0_p+alpha_p*(1-T0./T_reg))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T_reg))).*(x(1)*p)./(1+(b_L0.*exp(dU_L./(R*T_reg))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_reg))))./(xi_1.*exp(xi_2.*(1./T0-1./T_reg))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_reg))))./(xi_1.*exp(xi_2.*(1./T0-1./T_reg)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_reg))).*(x(1)*p)./(1+(b_U0.*exp(dU_U./(R*T_reg))).*(x(1)*p))+(b_H0.*exp(dU_H./(R*T_reg))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_reg))))./(xi_1.*exp(xi_2.*(1./T0-1./T_reg))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_reg))))./(xi_1.*exp(xi_2.*(1./T0-1./T_reg)))))).^gam); % adsorbed amount of CO2
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T_reg)*(x(1)*p)./(1e-6*R*T_reg)./(1+b0*exp(Hb/R/T_reg).*(x(1)*p)./(1e-6*R*T_reg)) + n2*d0*exp(Hd/R/T_reg)*(x(1)*p)./(1e-6*R*T_reg)./(1+d0*exp(Hd/R/T_reg).*(x(1)*p)./(1e-6*R*T_reg));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T_reg)*(x(1)*p)./(1+b0*exp(Hb/R/T_reg).*(x(1)*p)) + n2*d0*exp(Hd/R/T_reg)*(x(1)*p)./(1+d0*exp(Hd/R/T_reg).*(x(1)*p));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T_reg/T0)).*(b0*exp(dH/(R*T0)*(T0./T_reg-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T_reg-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T_reg))).^(1./(t0+alpha*(1-T0./T_reg)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T_reg./T0))).*((b0*exp(dH./(R*T0).*(T0./T_reg-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T_reg))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_reg-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T_reg)));
        end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid = x(1)*p*V*void/(R*T_reg)*1e6; % mol
        % solid and liquid CO2, final condition
        NCO2 = NCO2_solid + NCO2_fluid;

        % material balance CO2
        fvec(1) = NCO2_k - x(1)*x(3) - NCO2; % == 0


        %% equation 1: material balance CH4
        % solid phase CH4, initial condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4_k = ((cns0_c*exp(cXi_c*(1-T_k/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_c+calpha_c*(1-cT0./T_k))).^(1./(ct0_c+calpha_c*(1-cT0./T_k)))))+((cs0_p*exp(cXi_p*(1-T_k/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_p+calpha_p*(1-cT0./T_k))).^(1./(ct0_p+calpha_p*(1-cT0./T_k))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4_k = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k))).*(1-((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k))+(cb_H0.*exp(cdU_H./(R*T_k))).*(yCH4_k*p_k)).*((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam);
            case 'DSL'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k));
            case 'toth'
                qCH4_k = ((cs0*exp(cXi*(1-T_k/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0+calpha*(1-cT0./T_k))).^(1./(ct0+calpha*(1-cT0./T_k)))));
            case 'langfr'
                qCH4_k = (cs0.*exp(cXi.*(1-T_k./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k)));
        end
        NCH4_solid_k = adsorbentMass*qCH4_k;
        % fluid phase CH4, initial condition
        NCH4_fluid_k = yCH4_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CH4, initial condition
        NCH4_k = NCH4_solid_k + NCH4_fluid_k;

        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T_reg/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_reg-1))).*(p*x(4)))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_reg-1))).*(p*x(4))).^(ct0_c+calpha_c*(1-cT0./T_reg))).^(1./(ct0_c+calpha_c*(1-cT0./T_reg)))))+((cs0_p*exp(cXi_p*(1-T_reg/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_reg-1))).*(p*x(4)))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_reg-1))).*(p*x(4))).^(ct0_p+calpha_p*(1-cT0./T_reg))).^(1./(ct0_p+calpha_p*(1-cT0./T_reg))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T_reg))).*(x(4)*p)./(1+(cb_L0.*exp(cdU_L./(R*T_reg))).*(x(4)*p))).*(1-((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_reg))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_reg))))./(1+exp(((log((x(4)*p)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_reg))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_reg)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T_reg))).*(x(4)*p)./(1+(cb_U0.*exp(cdU_U./(R*T_reg))).*(x(4)*p))+(cb_H0.*exp(cdU_H./(R*T_reg))).*(x(4)*p)).*((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_reg))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_reg))))./(1+exp(((log((x(4)*p)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_reg))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_reg)))))).^cgam); % adsorbed amount of CH4
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T_reg)*(x(4)*p)./(1e-6*R*T_reg)./(1+cb0*exp(cHb/R/T_reg).*(x(4)*p)./(1e-6*R*T_reg)) + cn2*cd0*exp(cHd/R/T_reg)*(x(4)*p)./(1e-6*R*T_reg)./(1+cd0*exp(cHd/R/T_reg).*(x(4)*p)./(1e-6*R*T_reg));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T_reg)*(x(4)*p)./(1+cb0*exp(cHb/R/T_reg).*(x(4)*p)) + cn2*cd0*exp(cHd/R/T_reg)*(x(4)*p)./(1+cd0*exp(cHd/R/T_reg).*(x(4)*p));
            case 'toth'
                qCH4 = ((cs0*exp(cXi*(1-T_reg/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_reg-1))).*(x(4)*p))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_reg-1))).*(x(4)*p)).^(ct0+calpha*(1-cT0./T_reg))).^(1./(ct0+calpha*(1-cT0./T_reg)))));
            case 'langfr'
                qCH4 = (cs0.*exp(cXi.*(1-T_reg./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_reg-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T_reg))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_reg-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T_reg)));
        end
        NCH4_solid = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid = x(4)*p*V*void/(R*T_reg)*1e6; % mol
        % solid and liquid CH4, final condition
        NCH4 = NCH4_solid + NCH4_fluid;

        % material balance CH4
        fvec(4) = NCH4_k - x(4)*x(3) - NCH4; % == 0


        %% equation 2: material balance N2
        % solid phase N2, initial condition

        switch N2IsothermModel
            case 'toth_cp'
                qN2_k = ((nns0_c*exp(nXi_c*(1-T_k/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_c+nalpha_c*(1-nT0./T_k))).^(1./(nt0_c+nalpha_c*(1-nT0./T_k)))))+((nns0_p*exp(nXi_p*(1-T_k/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_p+nalpha_p*(1-nT0./T_k))).^(1./(nt0_p+nalpha_p*(1-nT0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                %                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yN2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yN2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yN2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
                qN2_k = (nq_L0.*nb_L0.*exp(ndU_L/(R*T_k)).*(yN2_k*p_k)./(1+nb_L0.*exp(ndU_L./(R*T_k)).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step/R.*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2.*(1/T0-1/T_k))))./(1+exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R.*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k)))))).^ngam))+(nq_L0*nb_U0.*exp(ndU_U./(R*T_k)).*(yN2_k*p_k)./(1+nb_U0*exp(ndU_U/(R*T_k)).*(yN2_k*p_k))+nb_H0*exp(ndU_H/(R*T_k)).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k))))./(1+exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k)))))).^ngam);
            case 'DSL'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k));
            case 'toth'
                qN2_k = ((ns0*exp(nXi*(1-T_k/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0+nalpha*(1-nT0./T_k))).^(1./(nt0+nalpha*(1-nT0./T_k)))));
            case 'langfr'
                qN2_k = (ns0.*exp(nXi.*(1-T_k./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k)));
        end

        NN2_solid_k = adsorbentMass*qN2_k;
        % fluid phase N2, initial condition
        NN2_fluid_k = yN2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid N2, initial condition
        NN2_k = NN2_solid_k + NN2_fluid_k;


        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T_reg/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_reg-1))).*(p*x(2)))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_reg-1))).*(p*x(2))).^(nt0_c+nalpha_c*(1-nT0./T_reg))).^(1./(nt0_c+nalpha_c*(1-nT0./T_reg)))))+((nns0_p*exp(nXi_p*(1-T_reg/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_reg-1))).*(p*x(2)))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_reg-1))).*(p*x(2))).^(nt0_p+nalpha_p*(1-nT0./T_reg))).^(1./(nt0_p+nalpha_p*(1-nT0./T_reg))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*nb_L0.*exp(ndU_L/(R*T_reg)).*(x(2)*p)./(1+nb_L0.*exp(ndU_L./(R*T_reg)).*(x(2)*p))).*(1-((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step/R.*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2.*(1/T0-1/T_reg))))./(1+exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R.*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2*(1/T0-1/T_reg)))))).^ngam))+(nq_L0*nb_U0.*exp(ndU_U./(R*T_reg)).*(x(2)*p)./(1+nb_U0*exp(ndU_U/(R*T_reg)).*(x(2)*p))+nb_H0*exp(ndU_H/(R*T_reg)).*(x(2)*p)).*((exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2*(1/T0-1/T_reg))))./(1+exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2*(1/T0-1/T_reg)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T_reg)*(x(2)*p)./(1e-6*R*T_reg)./(1+nb0*exp(nHb/R/T_reg).*(x(2)*p)./(1e-6*R*T_reg)) + nn2*nd0*exp(nHd/R/T_reg)*(x(2)*p)./(1e-6*R*T_reg)./(1+nd0*exp(nHd/R/T_reg).*(x(2)*p)./(1e-6*R*T_reg));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T_reg)*(x(2)*p)./(1+nb0*exp(nHb/R/T_reg).*(x(2)*p)) + nn2*nd0*exp(nHd/R/T_reg)*(x(2)*p)./(1+nd0*exp(nHd/R/T_reg).*(x(2)*p));
            case 'toth'
                qN2 = ((ns0*exp(nXi*(1-T_reg/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_reg-1))).*(x(2)*p))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_reg-1))).*(x(2)*p)).^(nt0+nalpha*(1-nT0./T_reg))).^(1./(nt0+nalpha*(1-nT0./T_reg)))));
            case 'langfr'
                qN2 = (ns0.*exp(nXi.*(1-T_reg./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_reg-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T_reg))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_reg-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T_reg)));
        end

        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = x(2)*p*V*void/(R*T_reg)*1e6; % mol
        % solid and liquid N2, final condition
        NN2 = NN2_solid + NN2_fluid;

        % material balance N2
        fvec(2) = (NCO2_k + NN2_k + NCH4_k) - x(3) - (NN2 + NCO2 + NCH4); % == 0



        %% equation 4: overall material balance
        fvec(3) = 1-x(1)-x(2)-x(4); % == 0

        fvec = real(fvec);
    end

    function [Q, qCO2] = Qex(x,adsorbentMass,cp,R)

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
                qCO2 = ((ns0_c*exp(Xi_c*(1-T_reg/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_reg-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_reg-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T_reg))).^(1./(t0_c+alpha_c*(1-T0./T_reg)))))+((ns0_p*exp(Xi_p*(1-T_reg/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_reg-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_reg-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T_reg))).^(1./(t0_p+alpha_p*(1-T0./T_reg))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*b_L0.*exp(dU_L/(R*T_reg)).*(x(1)*p)./(1+b_L0.*exp(dU_L./(R*T_reg)).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step/R.*(1/T0-1/T_reg))))./(xi_1*exp(xi_2.*(1/T0-1/T_reg))))./(1+exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R.*(1/T0-1/T_reg))))./(xi_1*exp(xi_2*(1/T0-1/T_reg)))))).^gam))+(q_L0*b_U0.*exp(dU_U./(R*T_reg)).*(x(1)*p)./(1+b_U0*exp(dU_U/(R*T_reg)).*(x(1)*p))+b_H0*exp(dU_H/(R*T_reg)).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R*(1/T0-1/T_reg))))./(xi_1*exp(xi_2*(1/T0-1/T_reg))))./(1+exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R*(1/T0-1/T_reg))))./(xi_1*exp(xi_2*(1/T0-1/T_reg)))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T_reg)*(x(1)*p)./(1e-6*R*T_reg)./(1+b0*exp(Hb/R/T_reg).*(x(1)*p)./(1e-6*R*T_reg)) + n2*d0*exp(Hd/R/T_reg)*(x(1)*p)./(1e-6*R*T_reg)./(1+d0*exp(Hd/R/T_reg).*(x(1)*p)./(1e-6*R*T_reg));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T_reg)*(x(1)*p)./(1+b0*exp(Hb/R/T_reg).*(x(1)*p)) + n2*d0*exp(Hd/R/T_reg)*(x(1)*p)./(1+d0*exp(Hd/R/T_reg).*(x(1)*p));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T_reg/T0)).*(b0*exp(dH/(R*T0)*(T0./T_reg-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T_reg-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T_reg))).^(1./(t0+alpha*(1-T0./T_reg)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T_reg./T0))).*((b0*exp(dH./(R*T0).*(T0./T_reg-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T_reg))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_reg-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T_reg)));
        end
        NCO2_solid = adsorbentMass*qCO2;

        % N2
        switch N2IsothermModel
            case 'toth_cp'
                qN2_k = ((nns0_c*exp(nXi_c*(1-T_k/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_c+nalpha_c*(1-nT0./T_k))).^(1./(nt0_c+nalpha_c*(1-nT0./T_k)))))+((nns0_p*exp(nXi_p*(1-T_k/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0_p+nalpha_p*(1-nT0./T_k))).^(1./(nt0_p+nalpha_p*(1-nT0./T_k))))); % adsorbed amount of CO2
            case 's_shaped'
                %                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yN2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yN2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yN2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yN2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam);
                qN2_k = (nq_L0.*nb_L0.*exp(ndU_L/(R*T_k)).*(yN2_k*p_k)./(1+nb_L0.*exp(ndU_L./(R*T_k)).*(yN2_k*p_k))).*(1-((exp((log((yN2_k*p_k))-log(np_step0.*exp(-ndH_step/R.*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2.*(1/T0-1/T_k))))./(1+exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R.*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k)))))).^ngam))+(nq_L0*nb_U0.*exp(ndU_U./(R*T_k)).*(yN2_k*p_k)./(1+nb_U0*exp(ndU_U/(R*T_k)).*(yN2_k*p_k))+nb_H0*exp(ndU_H/(R*T_k)).*(yN2_k*p_k)).*((exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k))))./(1+exp((log((yN2_k*p_k))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_k))))./(nxi_1*exp(nxi_2*(1/T0-1/T_k)))))).^ngam);
            case 'DSL'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1e-6*R*T_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qN2_k = nn1*nb0*exp(nHb/R/T_k)*(yN2_k*p_k)./(1+nb0*exp(nHb/R/T_k).*(yN2_k*p_k)) + nn2*nd0*exp(nHd/R/T_k)*(yN2_k*p_k)./(1+nd0*exp(nHd/R/T_k).*(yN2_k*p_k));
            case 'toth'
                qN2_k = ((ns0*exp(nXi*(1-T_k/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_k-1))).*(p_k*yN2_k)).^(nt0+nalpha*(1-nT0./T_k))).^(1./(nt0+nalpha*(1-nT0./T_k)))));
            case 'langfr'
                qN2_k = (ns0.*exp(nXi.*(1-T_k./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_k-1))).*(p_k*yN2_k)).^(1./nt0+nalpha.*(1-nT0./T_k)));
        end

        NN2_solid_k = adsorbentMass*qN2_k;

        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T_reg/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_reg-1))).*(p*x(2)))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T_reg-1))).*(p*x(2))).^(nt0_c+nalpha_c*(1-nT0./T_reg))).^(1./(nt0_c+nalpha_c*(1-nT0./T_reg)))))+((nns0_p*exp(nXi_p*(1-T_reg/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_reg-1))).*(p*x(2)))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./T_reg-1))).*(p*x(2))).^(nt0_p+nalpha_p*(1-nT0./T_reg))).^(1./(nt0_p+nalpha_p*(1-nT0./T_reg))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*nb_L0.*exp(ndU_L/(R*T_reg)).*(x(2)*p)./(1+nb_L0.*exp(ndU_L./(R*T_reg)).*(x(2)*p))).*(1-((exp((log((x(2)*p))-log(np_step0.*exp(-ndH_step/R.*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2.*(1/T0-1/T_reg))))./(1+exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R.*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2*(1/T0-1/T_reg)))))).^ngam))+(nq_L0*nb_U0.*exp(ndU_U./(R*T_reg)).*(x(2)*p)./(1+nb_U0*exp(ndU_U/(R*T_reg)).*(x(2)*p))+nb_H0*exp(ndU_H/(R*T_reg)).*(x(2)*p)).*((exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2*(1/T0-1/T_reg))))./(1+exp((log((x(2)*p))-log(np_step0*exp(-ndH_step/R*(1/T0-1/T_reg))))./(nxi_1*exp(nxi_2*(1/T0-1/T_reg)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T_reg)*(x(2)*p)./(1e-6*R*T_reg)./(1+nb0*exp(nHb/R/T_reg).*(x(2)*p)./(1e-6*R*T_reg)) + nn2*nd0*exp(nHd/R/T_reg)*(x(2)*p)./(1e-6*R*T_reg)./(1+nd0*exp(nHd/R/T_reg).*(x(2)*p)./(1e-6*R*T_reg));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T_reg)*(x(2)*p)./(1+nb0*exp(nHb/R/T_reg).*(x(2)*p)) + nn2*nd0*exp(nHd/R/T_reg)*(x(2)*p)./(1+nd0*exp(nHd/R/T_reg).*(x(2)*p));
            case 'toth'
                qN2 = ((ns0*exp(nXi*(1-T_reg/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T_reg-1))).*(x(2)*p))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T_reg-1))).*(x(2)*p)).^(nt0+nalpha*(1-nT0./T_reg))).^(1./(nt0+nalpha*(1-nT0./T_reg)))));
            case 'langfr'
                qN2 = (ns0.*exp(nXi.*(1-T_reg./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T_reg-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T_reg))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T_reg-1))).*(x(2)*p)).^(1./nt0+nalpha.*(1-nT0./T_reg)));
        end
        NN2_solid = adsorbentMass*qN2;

        %%
        % CH4
        % solid phase CH4, initial condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4_k = ((cns0_c*exp(cXi_c*(1-T_k/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_c+calpha_c*(1-cT0./T_k))).^(1./(ct0_c+calpha_c*(1-cT0./T_k)))))+((cs0_p*exp(cXi_p*(1-T_k/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0_p+calpha_p*(1-cT0./T_k))).^(1./(ct0_p+calpha_p*(1-cT0./T_k))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4_k = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_L0.*exp(cdU_L./(R*T_k))).*(yCH4_k*p_k))).*(1-((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k)./(1+(cb_U0.*exp(cdU_U./(R*T_k))).*(yCH4_k*p_k))+(cb_H0.*exp(cdU_H./(R*T_k))).*(yCH4_k*p_k)).*((exp((log((yCH4_k*p_k))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k))))./(1+exp(((log((yCH4_k*p_k)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T_k))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T_k)))))).^cgam);
            case 'DSL'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1e-6*R*T_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k)./(1e-6*R*T_k));
            case 'DSL2'
                qCH4_k = cn1*cb0*exp(cHb/R/T_k)*(yCH4_k*p_k)./(1+cb0*exp(cHb/R/T_k).*(yCH4_k*p_k)) + cn2*cd0*exp(cHd/R/T_k)*(yCH4_k*p_k)./(1+cd0*exp(cHd/R/T_k).*(yCH4_k*p_k));
            case 'toth'
                qCH4_k = ((cs0*exp(cXi*(1-T_k/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_k-1))).*(p_k*yCH4_k)).^(ct0+calpha*(1-cT0./T_k))).^(1./(ct0+calpha*(1-cT0./T_k)))));
            case 'langfr'
                qCH4_k = (cs0.*exp(cXi.*(1-T_k./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_k-1))).*(p_k*yCH4_k)).^(1./ct0+calpha.*(1-cT0./T_k)));
        end

        NCH4_solid_k = adsorbentMass*qCH4_k;

        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T_reg/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_reg-1))).*(p*x(4)))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T_reg-1))).*(p*x(4))).^(ct0_c+calpha_c*(1-cT0./T_reg))).^(1./(ct0_c+calpha_c*(1-cT0./T_reg)))))+((cs0_p*exp(cXi_p*(1-T_reg/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_reg-1))).*(p*x(4)))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./T_reg-1))).*(p*x(4))).^(ct0_p+calpha_p*(1-cT0./T_reg))).^(1./(ct0_p+calpha_p*(1-cT0./T_reg))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*cb_L0.*exp(cdU_L/(R*T_reg)).*(x(4)*p)./(1+cb_L0.*exp(cdU_L./(R*T_reg)).*(x(4)*p))).*(1-((exp((log((x(4)*p))-log(cp_step0.*exp(-cdH_step/R.*(1/T0-1/T_reg))))./(cxi_1*exp(cxi_2.*(1/T0-1/T_reg))))./(1+exp((log((x(4)*p))-log(cp_step0*exp(-cdH_step/R.*(1/T0-1/T_reg))))./(cxi_1*exp(cxi_2*(1/T0-1/T_reg)))))).^cgam))+(cq_L0*cb_U0.*exp(cdU_U./(R*T_reg)).*(x(4)*p)./(1+cb_U0*exp(cdU_U/(R*T_reg)).*(x(4)*p))+cb_H0*exp(cdU_H/(R*T_reg)).*(x(4)*p)).*((exp((log((x(4)*p))-log(cp_step0*exp(-cdH_step/R*(1/T0-1/T_reg))))./(cxi_1*exp(cxi_2*(1/T0-1/T_reg))))./(1+exp((log((x(4)*p))-log(cp_step0*exp(-cdH_step/R*(1/T0-1/T_reg))))./(cxi_1*exp(cxi_2*(1/T0-1/T_reg)))))).^cgam);
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T_reg)*(x(4)*p)./(1e-6*R*T_reg)./(1+cb0*exp(cHb/R/T_reg).*(x(4)*p)./(1e-6*R*T_reg)) + cn2*cd0*exp(cHd/R/T_reg)*(x(4)*p)./(1e-6*R*T_reg)./(1+cd0*exp(cHd/R/T_reg).*(x(4)*p)./(1e-6*R*T_reg));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T_reg)*(x(4)*p)./(1+cb0*exp(cHb/R/T_reg).*(x(4)*p)) + cn2*cd0*exp(cHd/R/T_reg)*(x(4)*p)./(1+cd0*exp(cHd/R/T_reg).*(x(4)*p));
            case 'toth'
                qCH4 = ((cs0*exp(cXi*(1-T_reg/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T_reg-1))).*(x(4)*p))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T_reg-1))).*(x(4)*p)).^(ct0+calpha*(1-cT0./T_reg))).^(1./(ct0+calpha*(1-cT0./T_reg)))));
            case 'langfr'
                qCH4 = (cs0.*exp(cXi.*(1-T_reg./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T_reg-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T_reg))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T_reg-1))).*(x(4)*p)).^(1./ct0+calpha.*(1-cT0./T_reg)));
        end
        NCH4_solid = adsorbentMass*qCH4;






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
                dqN2_dTequ_k = - nb0.*ns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 - (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nt0 - nalpha.*(nT0./T_k - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2))./((nt0 - nalpha.*(nT0./T_k - 1)).*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)) + 1)) - (nT0.*nalpha.*log((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1))./(T_k.^2.*(nt0 - nalpha.*(nT0./T_k - 1)).^2.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))))) - (nXi.*nb0.*ns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./(nT0.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)))) - (nb0.*ndH.*ns0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./(R.*T_k.^2.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))));
                dqN2_dpN2_k = (nb0.*ns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1))) - (nb0.^2.*ns0.*(p_k.*yN2_k).*exp((2.*ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(nt0 - nalpha.*(nT0./T_k - 1)) + 1).^(1./(nt0 - nalpha.*(nT0./T_k - 1)) + 1);
            case 'langfr'
                dqN2_dTequ_k = (ns0.*exp(-nXi.*(T_k./nT0 - 1)).*((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 + (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1) - (ns0.*exp(-nXi.*(T_k./nT0 - 1)).*((nT0.*nalpha.*log(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./T_k.^2 + (nb0.*ndH.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./(R.*T_k.^2)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1).^2 - (nXi.*ns0.*exp(-nXi.*(T_k./nT0 - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)))./(nT0.*((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1));
                dqN2_dpN2_k = (nb0.*ns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1).^2 - (nb0.*ns0.*exp((ndH.*(nT0./T_k - 1))./(R.*nT0)).*exp(-nXi.*(T_k./nT0 - 1)).*(nalpha.*(nT0./T_k - 1) - 1./nt0).*(nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1) - 1))./((nb0.*(p_k.*yN2_k).*exp((ndH.*(nT0./T_k - 1))./(R.*nT0))).^(1./nt0 - nalpha.*(nT0./T_k - 1)) + 1);
        end
        % dH
        dHN2 = R*T_k^2/(p_k*yN2_k) *(-dqN2_dTequ_k)/(dqN2_dpN2_k);

        % N2 desorbed
        dNN2_solid = NN2_solid - NN2_solid_k;


        % CH4
        switch CH4IsothermModel
            case 'toth_cp'
                dqCH4_dTequ_k = - cb0_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)).*(((cT0.*calpha_c.*log(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)))./T_k.^2 - (cb0_c.*cdH_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*(ct0_c - calpha_c.*(cT0./T_k - 1)).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0_c - calpha_c.*(cT0./T_k - 1)).*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha_c.*log((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0_c - calpha_c.*(cT0./T_k - 1)).^2.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1))))) - cb0_p.*cs0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)).*(((cT0.*calpha_p.*log(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)))./T_k.^2 - (cb0_p.*cdH_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*(ct0_p - calpha_p.*(cT0./T_k - 1)).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0_p - calpha_p.*(cT0./T_k - 1)).*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha_p.*log((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0_p - calpha_p.*(cT0./T_k - 1)).^2.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))))) - (cXi_c.*cb0_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./(cT0.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)))) - (cXi_p.*cb0_p.*cs0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./(cT0.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)))) - (cb0_c.*cdH_c.*cns0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)))) - (cb0_p.*cdH_p.*cs0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))));
                dqCH4_dpCH4_k = (cb0_c.*cns0_c.*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)))./((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1))) + (cb0_p.*cs0_p.*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)))./((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1))) - (cb0_c.^2.*cns0_c.*(p_k.*yCH4_k).*exp((2.*cdH_c.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_c.*(T_k./cT0 - 1)).*(cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1) - 1))./((cb0_c.*(p_k.*yCH4_k).*exp((cdH_c.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1).^(1./(ct0_c - calpha_c.*(cT0./T_k - 1)) + 1) - (cb0_p.^2.*cs0_p.*(p_k.*yCH4_k).*exp((2.*cdH_p.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi_p.*(T_k./cT0 - 1)).*(cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1) - 1))./((cb0_p.*(p_k.*yCH4_k).*exp((cdH_p.*(cT0./T_k - 1))./(R.*cT0))).^(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1).^(1./(ct0_p - calpha_p.*(cT0./T_k - 1)) + 1);
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
                dqCH4_dTequ_k = - cb0.*cs0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 - (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(ct0 - calpha.*(cT0./T_k - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2))./((ct0 - calpha.*(cT0./T_k - 1)).*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)) + 1)) - (cT0.*calpha.*log((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1))./(T_k.^2.*(ct0 - calpha.*(cT0./T_k - 1)).^2.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))))) - (cXi.*cb0.*cs0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./(cT0.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)))) - (cb0.*cdH.*cs0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./(R.*T_k.^2.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))));
                dqCH4_dpCH4_k = (cb0.*cs0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1))) - (cb0.^2.*cs0.*(p_k.*yCH4_k).*exp((2.*cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(ct0 - calpha.*(cT0./T_k - 1)) + 1).^(1./(ct0 - calpha.*(cT0./T_k - 1)) + 1);
            case 'langfr'
                dqCH4_dTequ_k = (cs0.*exp(-cXi.*(T_k./cT0 - 1)).*((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 + (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1) - (cs0.*exp(-cXi.*(T_k./cT0 - 1)).*((cT0.*calpha.*log(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./T_k.^2 + (cb0.*cdH.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./(R.*T_k.^2)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1).^2 - (cXi.*cs0.*exp(-cXi.*(T_k./cT0 - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)))./(cT0.*((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1));
                dqCH4_dpCH4_k = (cb0.*cs0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1).^2 - (cb0.*cs0.*exp((cdH.*(cT0./T_k - 1))./(R.*cT0)).*exp(-cXi.*(T_k./cT0 - 1)).*(calpha.*(cT0./T_k - 1) - 1./ct0).*(cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1) - 1))./((cb0.*(p_k.*yCH4_k).*exp((cdH.*(cT0./T_k - 1))./(R.*cT0))).^(1./ct0 - calpha.*(cT0./T_k - 1)) + 1);
        end
        % dH
        dHCH4 = R*T_k^2/(p_k*yCH4_k) *(-dqCH4_dTequ_k)/(dqCH4_dpCH4_k);

        % CH4 desorbed
        dNCH4_solid = NCH4_solid - NCH4_solid_k;




        %% overall
        Q = (adsorbentMass*cp*(T_reg-T_k) - dHCO2*dNCO2_solid-dHN2*dNN2_solid-dHCH4*dNCH4_solid)/1000; % kJ
    end



    function NCO2_solid = check_vac(x,adsorbentMass,cp,R)

        Tequ_k = T_k;
        Tequ = x(4);
        %% equation 5: energy balance
        % CO2
        % solid phase CO2, initial condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2_k = ((ns0_c*exp(Xi_c*(1-Tequ_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./Tequ_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./Tequ_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./Tequ_k))).^(1./(t0_c+alpha_c*(1-T0./Tequ_k)))))+((ns0_p*exp(Xi_p*(1-Tequ_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tequ_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tequ_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./Tequ_k))).^(1./(t0_p+alpha_p*(1-T0./Tequ_k))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*Tequ_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*Tequ_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/Tequ_k))))./(xi_1.*exp(xi_2.*(1./T0-1./Tequ_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/Tequ_k))))./(xi_1.*exp(xi_2.*(1./T0-1./Tequ_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*Tequ_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*Tequ_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*Tequ_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/Tequ_k))))./(xi_1.*exp(xi_2.*(1./T0-1./Tequ_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/Tequ_k))))./(xi_1.*exp(xi_2.*(1./T0-1./Tequ_k)))))).^gam);
            case 'DSL'
                qCO2_k = n1*b0*exp(Hb/R/Tequ_k)*(yCO2_k*p_k)./(1e-6*R*Tequ_k)./(1+b0*exp(Hb/R/Tequ_k).*(yCO2_k*p_k)./(1e-6*R*Tequ_k)) + n2*d0*exp(Hd/R/Tequ_k)*(yCO2_k*p_k)./(1e-6*R*Tequ_k)./(1+d0*exp(Hd/R/Tequ_k).*(yCO2_k*p_k)./(1e-6*R*Tequ_k));
            case 'DSL2'
                qCO2_k = n1*b0*exp(Hb/R/Tequ_k)*(yCO2_k*p_k)./(1+b0*exp(Hb/R/Tequ_k).*(yCO2_k*p_k)) + n2*d0*exp(Hd/R/Tequ_k)*(yCO2_k*p_k)./(1+d0*exp(Hd/R/Tequ_k).*(yCO2_k*p_k));
            case 'toth'
                qCO2_k = ((ns0*exp(Xi*(1-Tequ_k/T0)).*(b0*exp(dH/(R*T0)*(T0./Tequ_k-1))).*(p_k*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./Tequ_k-1))).*(p_k*yCO2_k)).^(t0+alpha*(1-T0./Tequ_k))).^(1./(t0+alpha*(1-T0./Tequ_k)))));
            case 'langfr'
                qCO2_k = (ns0.*exp(Xi.*(1-Tequ_k./T0))).*((b0*exp(dH./(R*T0).*(T0./Tequ_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./Tequ_k))./(1+ ((b0*exp(dH./(R*T0).*(T0./Tequ_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./Tequ_k)));
        end

        NCO2_solid_k = adsorbentMass*qCO2_k;

        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-Tequ/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./Tequ-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./Tequ-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./Tequ))).^(1./(t0_c+alpha_c*(1-T0./Tequ)))))+((ns0_p*exp(Xi_p*(1-Tequ/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tequ-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tequ-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./Tequ))).^(1./(t0_p+alpha_p*(1-T0./Tequ))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*b_L0.*exp(dU_L/(R*Tequ)).*(x(1)*p)./(1+b_L0.*exp(dU_L./(R*Tequ)).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step/R.*(1/T0-1/Tequ))))./(xi_1*exp(xi_2.*(1/T0-1/Tequ))))./(1+exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R.*(1/T0-1/Tequ))))./(xi_1*exp(xi_2*(1/T0-1/Tequ)))))).^gam))+(q_L0*b_U0.*exp(dU_U./(R*Tequ)).*(x(1)*p)./(1+b_U0*exp(dU_U/(R*Tequ)).*(x(1)*p))+b_H0*exp(dU_H/(R*Tequ)).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R*(1/T0-1/Tequ))))./(xi_1*exp(xi_2*(1/T0-1/Tequ))))./(1+exp((log((x(1)*p))-log(p_step0*exp(-dH_step/R*(1/T0-1/Tequ))))./(xi_1*exp(xi_2*(1/T0-1/Tequ)))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/Tequ)*(x(1)*p)./(1e-6*R*Tequ)./(1+b0*exp(Hb/R/Tequ).*(x(1)*p)./(1e-6*R*Tequ)) + n2*d0*exp(Hd/R/Tequ)*(x(1)*p)./(1e-6*R*Tequ)./(1+d0*exp(Hd/R/Tequ).*(x(1)*p)./(1e-6*R*Tequ));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/Tequ)*(x(1)*p)./(1+b0*exp(Hb/R/Tequ).*(x(1)*p)) + n2*d0*exp(Hd/R/Tequ)*(x(1)*p)./(1+d0*exp(Hd/R/Tequ).*(x(1)*p));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-Tequ/T0)).*(b0*exp(dH/(R*T0)*(T0./Tequ-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./Tequ-1))).*(x(1)*p)).^(t0+alpha*(1-T0./Tequ))).^(1./(t0+alpha*(1-T0./Tequ)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-Tequ./T0))).*((b0*exp(dH./(R*T0).*(T0./Tequ-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./Tequ))./(1+ ((b0*exp(dH./(R*T0).*(T0./Tequ-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./Tequ)));
        end
        NCO2_solid = adsorbentMass*qCO2;
    end
end
