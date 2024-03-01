
%% Simulate adsorption step
%% GET ISOTHERM PARAMETERS

function outAdsorption = simulateAdsorption(data,outputCool)


%%% for sub-ambient cooling adsorption set Ref="yes"
Ref="yes";

i=data.currentSorbent;

% CO2
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

% N2 Isotherm
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
        % CO2 isotherm
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


%% solver options
options = optimoptions('lsqnonlin','Display','off');
options.Algorithm = 'trust-region-reflective'; % levenberg-marquardt, trust-region-reflective
options.OptimalityTolerance = 1e-16;
options.FunctionTolerance = 1e-16;
options.StepTolerance = 1e-16;
options.MaxFunctionEvaluations = 6000;

%% RUN ADSORPTION MODEL
for m = 1:size(data.process.Tdes,1)

    Tads = data.process.Tads;
    Tamb = data.process.Tamb;
    p_amb = data.process.pamb;
    T = data.process.Tads;
    yCO2_feed = data.feed.yCO2;
    yN2_feed = data.feed.yN2;
    yCH4_feed = data.feed.yCH4;

    yCO2_k = outputCool(m).yCO2(end);
    yN2_k = outputCool(m).yN2(end);
    yCH4_k = outputCool(m).yCH4(end);


    %% Adsorption with assigned saturation level for CO2
    % calculate composition at given saturation
    yCO2_sat = data.startingCondBD(m,1); % If ful_saturation=true equal to feed composition
    yN2_sat=yN2_feed;
    yCH4_sat=yCH4_feed;

    %% Compare (inlet/capacity) of the components to identify the saturation order (CH4/CO2)
    %saturation_speed: 1) CO2, 2)CH4, 3)N2
    saturation_speed = adsorptioncompare(data.process.adsorbentMass,data.process.voidFraction,data.process.Vol,data.general.gasconstant,yCO2_k,yN2_k,yCH4_k,yCO2_sat,yN2_sat,yCH4_sat);
    if saturation_speed(2)<0
        sprintf('WARNING: CH4 Concentration at the end of the cooling step is higher than the Feed concentration')
    else
        %         sprintf('CH4 concentration at the end of cooling is not higher than feed concentration')

    end

    if saturation_speed(2)>saturation_speed(1) % saturates with methane first
        CH4_saturation=1;
        CO2_saturation=0;
        ysat_1=yCH4_sat;
        ysat_2=yCO2_sat;
        yi_k=yCO2_k;
    else
        CO2_saturation=1;
        CH4_saturation=0;
        ysat_1=yCO2_sat;
        ysat_2=yCH4_sat;
        yi_k=yCH4_k;
    end

    % Adsorption step (saturate with CH4)
    % Initial condition
    % x(1): N_feed, x(2): yN2, x(3):CO2/CH4, x(4): N_waste, saturates with CH4
    x0 = [1000, yN2_k, yi_k, 1000 ];

    % boundaries
    lbDim = [0, 0, 0, 0];
    ubDim = [5500, 1, 1, 5500];

    x0 = (x0-lbDim)./(ubDim-lbDim); % normalized

    lb = [0, 0, 0, 0];
    ub = [1, 1, 1, 1];

    funct = @(x) double(adsorption1(x,data.process.adsorbentMass,data.process.voidFraction,data.process.Vol,data.general.gasconstant,CH4_saturation,ysat_1,yCO2_k,yN2_k,yCH4_k,lbDim,ubDim));
    fun = @(x) funct(x);

    [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

    % dimensioned
    xopt = x;
    x = xopt.*(ubDim-lbDim)+lbDim;

    N_feed_1 = x(1);
    yN2_save_1 = x(2);

    if CH4_saturation==1
        yCO2_save_1 = x(3);
        yCH4_save_1=yCH4_sat;
    else
        yCH4_save_1=x(3);
        yCO2_save_1=yCO2_sat;
    end

    Nwaste_save_1 = x(4);
    x_save(:) = residual;



    Nwaste_save_sum(1) = Nwaste_save_1;
    Nfeed_save_sum(1) = N_feed_1;
    Nwaste_save_CO2(1) = Nwaste_save_1*yCO2_save_1;
    Nwaste_save_N2(1) = Nwaste_save_1*yN2_save_1;
    Nwaste_save_CH4(1) = Nwaste_save_1*yCH4_save_1;



    % Energy consumption till that given step
    E_total(1) = cmp_air(data.process.pamb,data.general.densityAir,x(1),data.general.MMCO2); % J, air blower

    %% Adsorption step 2(second component saturation)
    % Initial condition
    % x(1): N_feed, x(2): yN2, x(3):yCH4/CO2, x(4): N_waste, saturates with
    % the second component




    % Nin, yN2, Nout
    x0 = [N_feed_1, yN2_k, Nwaste_save_1];
    % boundaries
    lbDim = [0, 0, 0];
    ubDim = [5500, 1, 5500];

    x0 = (x0-lbDim)./(ubDim-lbDim); % normalized

    lb = [0, 0, 0];
    ub = [1, 1, 1];

    yCO2_k = yCO2_save_1;
    yN2_k = yN2_save_1;
    yCH4_k = yCH4_save_1;

    if CH4_saturation==1 % saturates with methane first
        second_saturation="CO2";
        ysat_2=yCO2_sat;
    else
        second_saturation="CH4";
        ysat_2=yCH4_sat;
    end





    funct = @(x) double(adsorption2(x,data.process.adsorbentMass,data.process.voidFraction,data.process.Vol,data.general.gasconstant,second_saturation,ysat_2,yCO2_k,yN2_k,yCH4_k,lbDim,ubDim));
    fun = @(x) funct(x);

    [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

    % dimensioned
    xopt = x;
    x = xopt.*(ubDim-lbDim)+lbDim;

    N_feed_2 = x(1);
    yCO2_save_2 = yCO2_sat;
    yN2_save_2 = x(2); %
    yCH4_save_2 = yCH4_sat;
    Nwaste_save_2 = x(3);
    x_save2(:) = residual;


    Nwaste_save_sum(2) = Nwaste_save_2;
    Nfeed_save_sum(2) = N_feed_2;
    Nwaste_save_CO2(2) = Nwaste_save_2*yCO2_save_2;
    Nwaste_save_N2(2) = Nwaste_save_2*yN2_save_2;
    Nwaste_save_CH4(2) = Nwaste_save_2*yCH4_save_2;




    % Energy consumption till that given step
    E_total(2) = cmp_air(data.process.pamb,data.general.densityAir,x(1),data.general.MMCO2); % J, air blower


    %% Combine both steps
    yCO2_save = [yCO2_save_1; yCO2_save_2];
    yN2_save = [yN2_save_1; yN2_save_2];
    yCH4_save = [yCH4_save_1; yCH4_save_2];

    % Energy consumption Refrigeration
    if Tads~=Tamb
        if Ref=="yes"
            E_Ref= W_Refrig(data.process.Tamb,data.process.Tads,sum(Nfeed_save_sum)); % E(KJ)
        end
    else
        E_Ref=0;
    end



    %outAdsorption = [yCO2_save  yN2_save   yCH4_save  Nwaste_save_sum'  Nwaste_save_CO2'  Nwaste_save_N2'   Nwaste_save_CH4'  E_total'  Nfeed_save_sum'];
    % outAdsorption = [yCO2_save  yN2_save  yH2O_save yCH4_save  Nwaste_save_sum'  Nwaste_save_CO2'  Nwaste_save_N2'  Nwaste_save_H2O' Nwaste_save_CH4'  E_total'  Nfeed_save_sum'];
    if data.plot
        plot_all = 'step_3';
        separate_plot_file_new3;
    end
    outAdsorption(m).yCO2 = yCO2_save; % saturation in bed
    outAdsorption(m).yN2 = yN2_save;
    outAdsorption(m).yCH4 = yCH4_save;
    outAdsorption(m).Nwaste_sum = Nwaste_save_sum;
    outAdsorption(m).Nwaste_CO2 = Nwaste_save_CO2;
    outAdsorption(m).Nwaste_N2 = Nwaste_save_N2;
    outAdsorption(m).Nwaste_CH4 = Nwaste_save_CH4;
    outAdsorption(m).E_total = E_total;
    outAdsorption(m).Nfeed_sum = Nfeed_save_sum;
    outAdsorption(m).E_Refrig= E_Ref;

end

%% functions
    function [yi_dq] = adsorptioncompare(adsorbentMass,void,V,R,yCO2_cool,yN2_cool, yCH4_cool,yCO2_sat, yN2_sat,yCH4_sat)
        % x(1): N_feed, x(2): yCO2, x(3): yN2, x(4): N_waste, x(5): yCH4

        %% equation 1: material balance CO2
        % CO2 solid and fluid from cooling
        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tads/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_cool))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_cool)).^(t0_p+alpha_p*(1-T0./Tads))).^(1./(t0_p+alpha_p*(1-T0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb))).*(1-((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_cool*p_amb)).*((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T)));
        end
        NCO2_solid_cool = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid_cool = yCO2_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition
        Ntot_CO2_cool = NCO2_solid_cool + NCO2_fluid_cool;

        % solid and fluid at end
        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tads/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_sat))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_sat)).^(t0_p+alpha_p*(1-T0./Tads))).^(1./(t0_p+alpha_p*(1-T0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb))).*(1-((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_sat*p_amb)).*((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T)));
        end
        NCO2_solid_end = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid_end = yCO2_sat*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition
        Ntot_CO2_end = NCO2_solid_end + NCO2_fluid_end;

        delqCO2=Ntot_CO2_end-Ntot_CO2_cool;



        % material balance CH4
        % CH4 solid and fluid from cooling
        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool)).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-Tads/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_cool))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_cool)).^(ct0_p+calpha_p*(1-cT0./Tads))).^(1./(ct0_p+calpha_p*(1-cT0./Tads))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_cool*p_amb)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_cool*p_amb))).*(1-((exp((log((yCH4_cool*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_cool*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_cool*p_amb)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_cool*p_amb))+(cb_H0.*exp(cdU_H./(R*T))).*(yCH4_cool*p_amb)).*((exp((log((yCH4_cool*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_cool*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam);
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_cool*p_amb)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(yCH4_cool*p_amb)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(yCH4_cool*p_amb)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(yCH4_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_cool*p_amb)./(1+cb0*exp(cHb/R/T).*(yCH4_cool*p_amb)) + cn2*cd0*exp(cHd/R/T)*(yCH4_cool*p_amb)./(1+cd0*exp(cHd/R/T).*(yCH4_cool*p_amb));
            case 'toth'
                qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
            case 'langfr'
                qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_cool)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_cool)).^(1./ct0+calpha.*(1-cT0./T)));
        end
        NCH4_solid_cool = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid_cool = yCH4_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CH4, final condition
        Ntot_CH4_cool = NCH4_solid_cool + NCH4_fluid_cool;

        % solid and fluid at end
        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat)).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-Tads/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_sat))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_sat)).^(ct0_p+calpha_p*(1-cT0./Tads))).^(1./(ct0_p+calpha_p*(1-cT0./Tads))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_sat*p_amb)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_sat*p_amb))).*(1-((exp((log((yCH4_sat*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_sat*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_sat*p_amb)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_sat*p_amb))+(cb_H0.*exp(cdU_H./(R*T))).*(yCH4_sat*p_amb)).*((exp((log((yCH4_sat*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_sat*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam);
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_sat*p_amb)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(yCH4_sat*p_amb)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(yCH4_sat*p_amb)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(yCH4_sat*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_sat*p_amb)./(1+cb0*exp(cHb/R/T).*(yCH4_sat*p_amb)) + cn2*cd0*exp(cHd/R/T)*(yCH4_sat*p_amb)./(1+cd0*exp(cHd/R/T).*(yCH4_sat*p_amb));
            case 'toth'
                qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
            case 'langfr'
                qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_sat)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_sat)).^(1./ct0+calpha.*(1-cT0./T)));
        end
        NCH4_solid_end = adsorbentMass*qCH4;
        %         fluid phase CH4, final condition
        NCH4_fluid_end = yCH4_sat*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CH4, final condition
        Ntot_CH4_end = NCH4_solid_end + NCH4_fluid_end;
        delqCH4=Ntot_CH4_end-Ntot_CH4_cool;
        %     fvec(1) = N_CH4_cool + yCH4_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CH4_end - yCH4_cool * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0

        %% equation 3: material balance N2

        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool)).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-Tads/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_cool))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_cool)).^(nt0_p+nalpha_p*(1-nT0./Tads))).^(1./(nt0_p+nalpha_p*(1-nT0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*(yN2_cool*p_amb)./(1+(nb_L0.*exp(ndU_L./(R*T))).*(yN2_cool*p_amb))).*(1-((exp((log((yN2_cool*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_cool*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*(yN2_cool*p_amb)./(1+(nb_U0.*exp(ndU_U./(R*T))).*(yN2_cool*p_amb))+(nb_H0.*exp(ndU_H./(R*T))).*(yN2_cool*p_amb)).*((exp((log((yN2_cool*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_cool*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_cool*p_amb)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*(yN2_cool*p_amb)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*(yN2_cool*p_amb)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*(yN2_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_cool*p_amb)./(1+nb0*exp(nHb/R/T).*(yN2_cool*p_amb)) + nn2*nd0*exp(nHd/R/T)*(yN2_cool*p_amb)./(1+nd0*exp(nHd/R/T).*(yN2_cool*p_amb));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool)).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_cool)).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_cool)).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid_cool = adsorbentMass*qN2;

        % solid and fluid at end
        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_sat))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_sat)).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-Tads/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_sat))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_sat)).^(nt0_p+nalpha_p*(1-nT0./Tads))).^(1./(nt0_p+nalpha_p*(1-nT0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*(yN2_sat*p_amb)./(1+(nb_L0.*exp(ndU_L./(R*T))).*(yN2_sat*p_amb))).*(1-((exp((log((yN2_sat*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_sat*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*(yN2_sat*p_amb)./(1+(nb_U0.*exp(ndU_U./(R*T))).*(yN2_sat*p_amb))+(nb_H0.*exp(ndU_H./(R*T))).*(yN2_sat*p_amb)).*((exp((log((yN2_sat*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_sat*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_sat*p_amb)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*(yN2_sat*p_amb)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*(yN2_sat*p_amb)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*(yN2_sat*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_sat*p_amb)./(1+nb0*exp(nHb/R/T).*(yN2_sat*p_amb)) + nn2*nd0*exp(nHd/R/T)*(yN2_sat*p_amb)./(1+nd0*exp(nHd/R/T).*(yN2_sat*p_amb));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_sat))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_sat)).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_sat)).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_sat)).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid_end = adsorbentMass*qN2;

        delqN2=NN2_solid_end-NN2_solid_cool;


        yi_dq = [yCO2_feed/delqCO2, yCH4_feed/delqCH4, yN2_feed/delqN2 ];

    end


    function fvec = adsorption1(x,adsorbentMass,void,V,R,CH4_saturation,ysat_1,yCO2_cool,yN2_cool,yCH4_cool,lb,ub)
        % x(1): N_feed, x(2): yN2, x(3):CO2/CH4, x(4): N_waste, saturates with CH4

        if CH4_saturation==1 % saturates with methane
            yCH4_out=ysat_1;
            yCO2_out=(x(3)*(ub(3)-lb(3))+lb(3));
        else
            yCO2_out=ysat_1;
            yCH4_out=(x(3)*(ub(3)-lb(3))+lb(3));
        end




        %% equation 1: material balance CO2
        % CO2 solid and fluid from cooling
        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tads/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_cool))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_cool)).^(t0_p+alpha_p*(1-T0./Tads))).^(1./(t0_p+alpha_p*(1-T0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb))).*(1-((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_cool*p_amb)).*((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T)));
        end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid = yCO2_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition
        N_CO2_cool = NCO2_solid + NCO2_fluid;

        % solid and fluid at end
        % solid phase CO2, final condition
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_out))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_out))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T))))+((ns0_p*exp(Xi_p*(1-Tads/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_out))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_out)).^(t0_p+alpha_p*(1-T0./Tads))).^(1./(t0_p+alpha_p*(1-T0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_out*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_out*p_amb))).*(1-((exp((log((yCO2_out*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_out*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_out*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_out*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_out*p_amb)).*((exp((log((yCO2_out*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_out*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam);
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_out*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_out*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_out*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_out*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_out*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_out*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_out*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_out*p_amb));
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_out))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_out)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_out)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_out)).^(1./t0+alpha.*(1-T0./T)));
        end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition
        NCO2_fluid = yCO2_out*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition
        N_CO2_end = NCO2_solid + NCO2_fluid;

        %         fvec(1) = N_CO2_cool + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end - yCO2_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

        fvec(1) = N_CO2_cool + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end; % = 0
        %% equation 1: material balance CH4
        % CH4 solid and fluid from cooling
        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool)).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-Tads/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_cool))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_cool)).^(ct0_p+calpha_p*(1-cT0./Tads))).^(1./(ct0_p+calpha_p*(1-cT0./Tads))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_cool*p_amb)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_cool*p_amb))).*(1-((exp((log((yCH4_cool*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_cool*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_cool*p_amb)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_cool*p_amb))+(cb_H0.*exp(cdU_H./(R*T))).*(yCH4_cool*p_amb)).*((exp((log((yCH4_cool*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_cool*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam);
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_cool*p_amb)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(yCH4_cool*p_amb)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(yCH4_cool*p_amb)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(yCH4_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_cool*p_amb)./(1+cb0*exp(cHb/R/T).*(yCH4_cool*p_amb)) + cn2*cd0*exp(cHd/R/T)*(yCH4_cool*p_amb)./(1+cd0*exp(cHd/R/T).*(yCH4_cool*p_amb));
            case 'toth'
                qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
            case 'langfr'
                qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_cool)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_cool)).^(1./ct0+calpha.*(1-cT0./T)));
        end
        NCH4_solid = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid = yCH4_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CH4, final condition
        N_CH4_cool = NCH4_solid + NCH4_fluid;

        % solid and fluid at end
        % solid phase CH4, final condition
        switch CH4IsothermModel
            case 'toth_cp'
                qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_out))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_out)).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-Tads/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_out))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_out)).^(ct0_p+calpha_p*(1-cT0./Tads))).^(1./(ct0_p+calpha_p*(1-cT0./Tads))))); % adsorbed amount of CH4
            case 's_shaped'
                qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_out*p_amb)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_out*p_amb))).*(1-((exp((log((yCH4_out*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_out*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_out*p_amb)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_out*p_amb))+(cb_H0.*exp(cdU_H./(R*T))).*(yCH4_out*p_amb)).*((exp((log((yCH4_out*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_out*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam);
            case 'DSL'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_out*p_amb)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(yCH4_out*p_amb)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(yCH4_out*p_amb)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(yCH4_out*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_out*p_amb)./(1+cb0*exp(cHb/R/T).*(yCH4_out*p_amb)) + cn2*cd0*exp(cHd/R/T)*(yCH4_out*p_amb)./(1+cd0*exp(cHd/R/T).*(yCH4_out*p_amb));
            case 'toth'
                qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_out))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_out)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
            case 'langfr'
                qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_out)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_out)).^(1./ct0+calpha.*(1-cT0./T)));
        end
        NCH4_solid = adsorbentMass*qCH4;
        % fluid phase CH4, final condition
        NCH4_fluid = yCH4_out*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CH4, final condition
        N_CH4_end = NCH4_solid + NCH4_fluid;

        %         fvec(3) = N_CH4_cool + yCH4_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CH4_end - yCH4_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

        fvec(3) = N_CH4_cool + yCH4_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CH4_end; % = 0


        %% equation 3: material balance N2

        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool)).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-Tads/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_cool))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_cool)).^(nt0_p+nalpha_p*(1-nT0./Tads))).^(1./(nt0_p+nalpha_p*(1-nT0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*(yN2_cool*p_amb)./(1+(nb_L0.*exp(ndU_L./(R*T))).*(yN2_cool*p_amb))).*(1-((exp((log((yN2_cool*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_cool*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*(yN2_cool*p_amb)./(1+(nb_U0.*exp(ndU_U./(R*T))).*(yN2_cool*p_amb))+(nb_H0.*exp(ndU_H./(R*T))).*(yN2_cool*p_amb)).*((exp((log((yN2_cool*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_cool*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_cool*p_amb)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*(yN2_cool*p_amb)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*(yN2_cool*p_amb)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*(yN2_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_cool*p_amb)./(1+nb0*exp(nHb/R/T).*(yN2_cool*p_amb)) + nn2*nd0*exp(nHd/R/T)*(yN2_cool*p_amb)./(1+nd0*exp(nHd/R/T).*(yN2_cool*p_amb));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool)).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_cool)).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_cool)).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = yN2_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid N2, final condition
        N_N2_cool = NN2_solid + NN2_fluid;

        % solid and fluid at end
        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-Tads/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(nt0_p+nalpha_p*(1-nT0./Tads))).^(1./(nt0_p+nalpha_p*(1-nT0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+(nb_L0.*exp(ndU_L./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))).*(1-((exp((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+(nb_U0.*exp(ndU_U./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))+(nb_H0.*exp(ndU_H./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)).*((exp((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+nb0*exp(nHb/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)) + nn2*nd0*exp(nHd/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+nd0*exp(nHd/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = (x(2)*(ub(2)-lb(2))+lb(2))*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid N2, final condition
        N_N2_end = NN2_solid + NN2_fluid;


        %%%%%%%%% (x(2)*(ub(2)-lb(2))+lb(2)) = yN2_cool

        %         fvec(2) = N_N2_cool + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - yN2_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

        fvec(2) = N_N2_cool + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

        %% equation 4: overall concentration balance
        fvec(4) =  yCH4_out + (x(2)*(ub(2)-lb(2))+lb(2)) + yCO2_out  - 1; % = 0

        fvec = real(fvec);
    end

    function fvec = adsorption2(x,adsorbentMass,void,V,R,second_saturation,ysat_2,yCO2_cool,yN2_cool,yCH4_cool,lb,ub)
        % x(1): N_feed, x(2): yN2, x(3): N_waste, saturates with second component
        %% equation 1: material balance CO2
        % CO2 solid and fluid from cooling
        % solid phase CO2, final condition

        if second_saturation=="CH4"
            yCH4_sat=ysat_2;
            ysat_out=yCO2_sat;

            switch CH4IsothermModel
                case 'toth_cp'
                    qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool)).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-Tads/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_cool))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_cool)).^(ct0_p+calpha_p*(1-cT0./Tads))).^(1./(ct0_p+calpha_p*(1-cT0./Tads))))); % adsorbed amount of CH4
                case 's_shaped'
                    qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_cool*p_amb)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_cool*p_amb))).*(1-((exp((log((yCH4_cool*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_cool*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_cool*p_amb)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_cool*p_amb))+(cb_H0.*exp(cdU_H./(R*T))).*(yCH4_cool*p_amb)).*((exp((log((yCH4_cool*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_cool*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam);
                case 'DSL'
                    qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_cool*p_amb)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(yCH4_cool*p_amb)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(yCH4_cool*p_amb)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(yCH4_cool*p_amb)./(1e-6*R*T));
                case 'DSL2'
                    qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_cool*p_amb)./(1+cb0*exp(cHb/R/T).*(yCH4_cool*p_amb)) + cn2*cd0*exp(cHd/R/T)*(yCH4_cool*p_amb)./(1+cd0*exp(cHd/R/T).*(yCH4_cool*p_amb));
                case 'toth'
                    qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_cool)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
                case 'langfr'
                    qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_cool)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_cool)).^(1./ct0+calpha.*(1-cT0./T)));
            end
            NCH4_solid = adsorbentMass*qCH4;
            % fluid phase CH4, final condition
            NCH4_fluid = yCH4_cool*p_amb*V*void/(R*T)*1e6; % mol
            % solid and fluid CH4, final condition
            N_CH4_cool = NCH4_solid + NCH4_fluid;

            % solid and fluid at end
            % solid phase CH4, final condition
            switch CH4IsothermModel
                case 'toth_cp'
                    qCH4 = ((cns0_c*exp(cXi_c*(1-T/cT0)).*(cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat))./((1+((cb0_c*exp(cdH_c/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat)).^(ct0_c+calpha_c*(1-cT0./T))).^(1./(ct0_c+calpha_c*(1-cT0./T)))))+((cns0_p*exp(cXi_p*(1-Tads/cT0)).*(cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_sat))./((1+((cb0_p*exp(cdH_p/(R*cT0)*(cT0./Tads-1))).*(p_amb*yCH4_sat)).^(ct0_p+calpha_p*(1-cT0./Tads))).^(1./(ct0_p+calpha_p*(1-cT0./Tads))))); % adsorbed amount of CH4
                case 's_shaped'
                    qCH4 = (cq_L0.*(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_sat*p_amb)./(1+(cb_L0.*exp(cdU_L./(R*T))).*(yCH4_sat*p_amb))).*(1-((exp((log((yCH4_sat*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_sat*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam))+(cq_U0.*(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_sat*p_amb)./(1+(cb_U0.*exp(cdU_U./(R*T))).*(yCH4_sat*p_amb))+(cb_H0.*exp(cdU_H./(R*T))).*(yCH4_sat*p_amb)).*((exp((log((yCH4_sat*p_amb))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T))))./(1+exp(((log((yCH4_sat*p_amb)))-log(cp_step0.*exp(-cdH_step./R.*(1/cT0-1/T))))./(cxi_1.*exp(cxi_2.*(1./cT0-1./T)))))).^cgam);
                case 'DSL'
                    qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_sat*p_amb)./(1e-6*R*T)./(1+cb0*exp(cHb/R/T).*(yCH4_sat*p_amb)./(1e-6*R*T)) + cn2*cd0*exp(cHd/R/T)*(yCH4_sat*p_amb)./(1e-6*R*T)./(1+cd0*exp(cHd/R/T).*(yCH4_sat*p_amb)./(1e-6*R*T));
                case 'DSL2'
                    qCH4 = cn1*cb0*exp(cHb/R/T)*(yCH4_sat*p_amb)./(1+cb0*exp(cHb/R/T).*(yCH4_sat*p_amb)) + cn2*cd0*exp(cHd/R/T)*(yCH4_sat*p_amb)./(1+cd0*exp(cHd/R/T).*(yCH4_sat*p_amb));
                case 'toth'
                    qCH4 = ((cns0*exp(cXi*(1-T/cT0)).*(cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat))./((1+((cb0*exp(cdH/(R*cT0)*(cT0./T-1))).*(p_amb*yCH4_sat)).^(ct0+calpha*(1-cT0./T))).^(1./(ct0+calpha*(1-cT0./T)))));
                case 'langfr'
                    qCH4 = (cns0.*exp(cXi.*(1-T./cT0))).*((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_sat)).^(1./ct0+calpha.*(1-cT0./T))./(1+ ((cb0*exp(cdH./(R*cT0).*(cT0./T-1))).*(p_amb*yCH4_sat)).^(1./ct0+calpha.*(1-cT0./T)));
            end
            NCH4_solid = adsorbentMass*qCH4;
            % fluid phase CH4, final condition
            NCH4_fluid = yCH4_sat*p_amb*V*void/(R*T)*1e6; % mol
            % solid and fluid CH4, final condition
            N_CH4_end = NCH4_solid + NCH4_fluid;

            %             fvec(1) = N_CH4_cool + yCH4_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CH4_end - yCH4_cool * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0


            fvec(1) = N_CH4_cool + yCH4_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CH4_end ; % = 0



        elseif second_saturation=="CO2"
            yCO2_sat=ysat_2;
            ysat_out=yCH4_sat;



            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tads/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_cool))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_cool)).^(t0_p+alpha_p*(1-T0./Tads))).^(1./(t0_p+alpha_p*(1-T0./Tads))))); % adsorbed amount of CO2
                case 's_shaped'
                    qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb))).*(1-((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_cool*p_amb)).*((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam);
                case 'DSL'
                    qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T));
                case 'DSL2'
                    qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb));
                case 'toth'
                    qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
                case 'langfr'
                    qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T)));
            end
            NCO2_solid = adsorbentMass*qCO2;
            % fluid phase CO2, final condition
            NCO2_fluid = yCO2_cool*p_amb*V*void/(R*T)*1e6; % mol
            % solid and fluid CO2, final condition
            N_CO2_cool = NCO2_solid + NCO2_fluid;

            % solid and fluid at end
            % solid phase CO2, final condition
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tads/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_sat))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tads-1))).*(p_amb*yCO2_sat)).^(t0_p+alpha_p*(1-T0./Tads))).^(1./(t0_p+alpha_p*(1-T0./Tads))))); % adsorbed amount of CO2
                case 's_shaped'
                    qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb))).*(1-((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_sat*p_amb)).*((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam);
                case 'DSL'
                    qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T));
                case 'DSL2'
                    qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb));
                case 'toth'
                    qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));
                case 'langfr'
                    qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T)));
            end
            NCO2_solid = adsorbentMass*qCO2;
            % fluid phase CO2, final condition
            NCO2_fluid = yCO2_sat*p_amb*V*void/(R*T)*1e6; % mol
            % solid and fluid CO2, final condition
            N_CO2_end = NCO2_solid + NCO2_fluid;

            %             fvec(1) = N_CO2_cool + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end - yCO2_cool * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0
            fvec(1) = N_CO2_cool + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end; % = 0


        end




        %% equation 2: material balance N2

        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool)).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-Tads/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_cool))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*yN2_cool)).^(nt0_p+nalpha_p*(1-nT0./Tads))).^(1./(nt0_p+nalpha_p*(1-nT0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*(yN2_cool*p_amb)./(1+(nb_L0.*exp(ndU_L./(R*T))).*(yN2_cool*p_amb))).*(1-((exp((log((yN2_cool*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_cool*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*(yN2_cool*p_amb)./(1+(nb_U0.*exp(ndU_U./(R*T))).*(yN2_cool*p_amb))+(nb_H0.*exp(ndU_H./(R*T))).*(yN2_cool*p_amb)).*((exp((log((yN2_cool*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log((yN2_cool*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_cool*p_amb)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*(yN2_cool*p_amb)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*(yN2_cool*p_amb)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*(yN2_cool*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*(yN2_cool*p_amb)./(1+nb0*exp(nHb/R/T).*(yN2_cool*p_amb)) + nn2*nd0*exp(nHd/R/T)*(yN2_cool*p_amb)./(1+nd0*exp(nHd/R/T).*(yN2_cool*p_amb));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*yN2_cool)).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_cool)).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*yN2_cool)).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = yN2_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid N2, final condition
        N_N2_cool = NN2_solid + NN2_fluid;

        % solid and fluid at end
        % solid phase N2, final condition
        switch N2IsothermModel
            case 'toth_cp'
                qN2 = ((nns0_c*exp(nXi_c*(1-T/nT0)).*(nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((nb0_c*exp(ndH_c/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(nt0_c+nalpha_c*(1-nT0./T))).^(1./(nt0_c+nalpha_c*(1-nT0./T)))))+((nns0_p*exp(nXi_p*(1-Tads/nT0)).*(nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((nb0_p*exp(ndH_p/(R*nT0)*(nT0./Tads-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(nt0_p+nalpha_p*(1-nT0./Tads))).^(1./(nt0_p+nalpha_p*(1-nT0./Tads))))); % adsorbed amount of CO2
            case 's_shaped'
                qN2 = (nq_L0.*(nb_L0.*exp(ndU_L./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+(nb_L0.*exp(ndU_L./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))).*(1-((exp((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam))+(nq_U0.*(nb_U0.*exp(ndU_U./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+(nb_U0.*exp(ndU_U./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))+(nb_H0.*exp(ndU_H./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)).*((exp((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T))))./(1+exp(((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)))-log(np_step0.*exp(-ndH_step./R.*(1/nT0-1/T))))./(nxi_1.*exp(nxi_2.*(1./nT0-1./T)))))).^ngam);
            case 'DSL'
                qN2 = nn1*nb0*exp(nHb/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)./(1+nb0*exp(nHb/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)) + nn2*nd0*exp(nHd/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)./(1+nd0*exp(nHd/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T));
            case 'DSL2'
                qN2 = nn1*nb0*exp(nHb/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+nb0*exp(nHb/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)) + nn2*nd0*exp(nHd/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+nd0*exp(nHd/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb));
            case 'toth'
                qN2 = ((nns0*exp(nXi*(1-T/nT0)).*(nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((nb0*exp(ndH/(R*nT0)*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(nt0+nalpha*(1-nT0./T))).^(1./(nt0+nalpha*(1-nT0./T)))));
            case 'langfr'
                qN2 = (nns0.*exp(nXi.*(1-T./nT0))).*((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(1./nt0+nalpha.*(1-nT0./T))./(1+ ((nb0*exp(ndH./(R*nT0).*(nT0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(1./nt0+nalpha.*(1-nT0./T)));
        end
        NN2_solid = adsorbentMass*qN2;
        % fluid phase N2, final condition
        NN2_fluid = (x(2)*(ub(2)-lb(2))+lb(2))*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid N2, final condition
        N_N2_end = NN2_solid + NN2_fluid;


        %%%%%%%%% (x(2)*(ub(2)-lb(2))+lb(2)) = yN2_cool

        fvec(2) = N_N2_cool + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - (1-ysat_out) * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0


        %% equation 4: overall concentration balance
        fvec(3) =  yCO2_sat + (x(2)*(ub(2)-lb(2))+lb(2)) + yCH4_sat  - 1; % = 0

        fvec = real(fvec);
    end

end