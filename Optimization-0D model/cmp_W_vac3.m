function E= cmp_W_vac3(pin,pout,Tin,mol_in)

poleffHH = 0.72; % efficiency vacuum pump p = 1.01 bar
poleffH = 0.7; % efficiency vacuum pump p = 0.8 bar
poleffM = 0.6; % efficiency vaccum pump p = 0.1 bar
poleffL = 0.3; % efficiency vacuum pump p = 0.01 bar
poleff=interp1([0.01 0.1 0.8 1.01],[poleffL poleffM poleffH poleffHH],pin);

H2O=1;
CO2=2;
N2=3;
MM=1;
I(H2O,MM)=18; % kg/kmol
I(CO2,MM)=44.009799999999998; % kg/kmol
I(N2,MM)=28.013480000000001; % kg/kmol
molFrac=2;
I(:,molFrac)=(mol_in')./sum(mol_in); % H2O, CO2, N2

mass=mol_in*I(:,MM)/1000; % kg of air

beta=pout/pin;

R = 8.314; % universal gas constant, J/K/mol
cp_in=37.11; % J/K/mol
cp_in=cp_in/I(CO2,MM); % J/K/kg
gamma=cp_in/(cp_in - R/I(CO2,MM));
theta=(gamma-1)/gamma;

E=mass*cp_in*Tin/poleff*(beta^theta-1)/1000; % kJ

end