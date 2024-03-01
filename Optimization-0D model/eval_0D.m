
%% fitness function for multiobjective optimization of 0D model

function res = eval_0D(X,c,m,file_name,file_optimization,ppm_CO2)

[Ex_spec_CH4_all, WC_CH4_all, purity_CH4_all, recovery_CH4_all, Ex_spec_CH4_Evac, WC_CH4_Evac, purity_CH4_Evac, recovery_CH4_Evac, Ex_spec_CO2_all, WC_CO2_all, purity_CO2_all, recovery_CO2_all, Ex_spec_CO2_Evac, WC_CO2_Evac, purity_CO2_Evac, recovery_CO2_Evac, succ] = run_0D(X,c,m,file_name,file_optimization,ppm_CO2);

spec_PurCO2 = 0.7;

if succ==1
    % constraint in purity of CO2 (actual purity - purity constraint)
    %         C = min(0., (purity_dry - spec_PurCO2)*10); % C>=0 to meet constraints
    %         violation = sum((min(0,C)).^2);

%     ObjF = [-purity_CH4', -recovery_CH4' Q_spec'];
%     ObjF = [WC_CH4_all',Ex_spec_CH4_all'];
    ObjF = [purity_CH4_all',recovery_CH4_all'];

    %         Perf1 = ObjF(:,1)+violation;
    %         Perf2 = ObjF(:,2)+violation;
    Perf1 = ObjF(:,1); % purity CH4
    Perf2 = ObjF(:,2); % recovery CH4
%     Perf3 = ObjF(:,3); % Thermal Energy
    %         Perf3 = ObjF(:,3);
%     if recovery_CH4>1
%         Perf2 = 100;
%     end
    if isnan(ObjF)
        Perf1 = +inf;
        Perf2 = +inf;
%         Perf3 = +inf;
    end
else
    Perf1 = +inf;
    Perf2 = +inf;
%     Perf3 = +inf;
end
% res = [Perf1, Perf2, Perf3];
    res = [Perf1, Perf2];
end

