

%% calculate temperature profile and heating time

function [TCoolProfileStep, TimeCoolProfilePlot, pressureVectorCool] = coolingProfile_1D(data,Tdes,pvac)

	for k = 1:length(Tdes)
		pvac2 = pvac(k);
		Tdes2 = Tdes(k);
		
		p = [0.174853800090986,-0.0262921350934970];

		timeCooling = data.process.coolingTime-1; % directly from 1D model
		
		xm = linspace(0,44,44);
		pp = [0.000221560995722903,0.00170366982722347,0.141834668956426,0.000322209875008416,0.00168615364119241];
		pressureVector = pvac2.*10.*exp(((pp(1))-(pp(4)).*pvac2.*10).*xm)+((pp(2))-(pp(5)).*pvac2.*10).*exp((pp(3)).*xm); % bar
		pressureVector = [pressureVector./10, ones(1,(timeCooling-30)).*data.process.pamb]; % MPa

		% for the calculation, reduce/increase no of steps to 100
		max_steps = data.process.noSteps;

		timeEnd = 0;
		timeEnd2 = 0;
		TCoolProfileStep(k,1) = Tdes2;
		a = (timeCooling-1)/max_steps;
		b = (timeCooling-1)/max_steps; % for time profile
		pressureVectorCool(k,1) = pvac2;
		for m = (1:max_steps)
			pressureVectorCool(k,m+1) = pressureVector(round(m*a));
			TCoolProfileStep(k,m+1) = 293 + (p(1)).*Tdes2*exp((p(2)).*timeEnd); % 
			timeEnd = timeEnd+a;
			TimeCoolProfilePlot(k,m) = m*b;
			timeEnd2 = timeEnd2+b;
		end
	end
end