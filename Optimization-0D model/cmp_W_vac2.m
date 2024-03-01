function E= cmp_W_vac2(pin,pout,T,Nout,x) % kJ vacuum pump CO2 prod
% x: (1): yCO2, (2): yN2, (3): Nout, (4): T, (5): yCH4
% varying efficiency with vacuum level
% 	poleffHH = 0.72; % efficiency vacuum pump p = 1.01 bar
% 	poleffH = 0.7; % efficiency vacuum pump p = 0.8 bar
% 	poleffM = 0.6; % efficiency vaccum pump p = 0.1 bar
% 	poleffL = 0.3; % efficiency vacuum pump p = 0.01 bar
% 	poleff=interp1([0.01 0.1 0.8 1.01],[poleffL poleffM poleffH poleffHH],pin);

R0 = 8.314; % J/mol/K

%%% Molar mass (g/mol)
MW_CO2 = 44.01; MW_N2 = 28.013; MW_CH4 = 16.043;
%%% Cp (J/g/K)*MW(g/mol)= Cp(J/mol/k), @ 300 K
cp_CO2 = 0.846*MW_CO2; cp_N2 = 1.039*MW_N2; cp_CH4 = 2.2537*MW_CH4;
cp_mix=x(1)*cp_CO2+x(2)*cp_N2+x(5)*cp_CH4; % cp_mix=sum(xi*Cpi)
% 	cp_in = 37/44;
% 	gamma = cp_in/(cp_in - R0/44);
gamma = cp_mix/(cp_mix - R0);
%     theta = (gamma-1)/gamma;

% 	E = 1/poleff*1/theta*R0*T*((pout/pin)^(theta)-1)*sum(Nout)/1000; %
% 	Arvind with varying effciency

% Webley
poleff = 0.75;
k = gamma;
R = R0; % J/mol/K
Tin = T;
% 	E = 1/poleff*k*R*Tin/(k-1)*((pout/pin).^((k-1)/k)-1)/1000; % kJ???	*sum(Nout) original
E = sum(Nout)/poleff*k*R*Tin/(k-1)*((pout/pin).^((k-1)/k)-1)/1000; % kJ

if isnan(E)
    E=0;
else
end

end