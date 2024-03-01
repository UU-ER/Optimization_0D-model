

%% calculate temperature profile and heating time

function [THeatProfileStep, timeHeatingNN, TimeHeatProfilePlot] = temperatureProfile(data,Tdes)

	% 1. calculate time dependent on the specific sorbent (rho) and time
    % 2. make the profile the same length as the time
	
	% for the calculation, reduce/increase no of steps to 100
    max_steps = data.process.noSteps;
	pvac = data.process.pvac;
	
	% starting point
	THeatProfileStep = ones(length(Tdes),max_steps);   
	TimeHeatProfilePlot = ones(length(Tdes),max_steps);

    for k = 1:length(Tdes)
        Tdes2 = Tdes(k);
		pvac2 = pvac(k);
		% calculate time using NN
		load('net_new2','netNN_new');
		ParamNN = netNN_new([data.sorbent(data.currentSorbent).MaterialDensity;...
			data.sorbent(data.currentSorbent).Density;...
			data.sorbent(data.currentSorbent).cp; ...
			Tdes2;pvac2]);
		timeHeatingNN(k) = ParamNN; % 
		
		% profile
			p = [-687.316091760085,661.535735641533,42.0052174647580,-2791.97322593205,-4.29623797695230,-0.345664545572582,-34162.5562758762,0.00867252521860531]; % from fitting
			a = (timeHeatingNN(k)-1)/(max_steps-1); % 
			b = (2606-1)/(max_steps-1); % for time profile
			timeEnd = b;
			THeatProfileStep(k,1) = data.process.Tamb;
			TimeHeatProfilePlot(k,1) = 0;
				for m = (2:(max_steps-1))
					THeatProfileStep(k,m) = (p(1)+Tdes2-p(4).*Tdes2.^p(6)).*atan((timeEnd-p(2).*1e9.*Tdes2.^p(5))./p(3))+(Tdes2-p(7)).*p(8);
					TimeHeatProfilePlot(k,m) = m*a;
					timeEnd = timeEnd+b;
				end
			
		% end point
			THeatProfileStep(k,max_steps) =  Tdes2;
			TimeHeatProfilePlot(k,max_steps) = timeHeatingNN(k);
    end
end