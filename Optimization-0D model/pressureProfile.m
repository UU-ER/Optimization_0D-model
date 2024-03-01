

%% calculate temperature profile and heating time

function [pressureVector, timeBD, TimeBDProfilePlot] = pressureProfile(data,pvac)

    for k = 1:length(pvac)

        pvac2 = pvac(k)*10; % bar
        x = 1;
        y = 0;
        deltaP = abs(pvac2-y);
		pp = [0.186277700879350,-0.123808829948893,-7.13805471361795e-05];
        while deltaP>1e-3 && x<10000
            y = (pp(1))./(pvac2.^0.85).*exp((pp(2)).*x) +...
            pvac2.*exp((pp(3))./(pvac2.^0.85).*x);
            x = x+1;
            deltaP = abs(pvac2-y);
        end
		
%	timeEnd1 = (1:(x-1)); % 2000/100.*x
%    pressureVector = ((pp(1))./(pvac.^0.85).*exp((pp(2)).*x) +...
%            pvac.*exp((pp(3))./(pvac.^0.85).*x))./10; % MPa

    % pvac = -0.1*x+2.07;;
    timeBD(k) = floor((-10*pvac2+20.7)./100*200);

    % for the calculation, reduce/increase no of steps to 100
    max_steps = data.process.noSteps;
    timeEnd1 = 0;
    a = (x-1)/(max_steps-2);
    b = (timeBD(k)-1)/(max_steps-1); % for time profile
    timeEnd2 = b;
    pressureVector(k,1) = data.process.pamb;
    TimeBDProfilePlot(k,1) = 0;
    for m = (2:(max_steps-1))
        pressureVector(k,m) = ((pp(1))./(pvac2.^0.85).*exp((pp(2)).*timeEnd1) +...
            pvac2.*exp((pp(3))./(pvac2.^0.85).*timeEnd1))./10; % MPa
%         pressureVector(m) = (pvac+0.170255738564685./pvac.*exp(-0.536951393269455.*timeEnd1))./10; % MPa
        timeEnd1 = timeEnd1+a;
        TimeBDProfilePlot(k,m) = m*b;
        timeEnd2 = timeEnd2+b;
    end
    pressureVector(k,max_steps) = pvac2/10;
    TimeBDProfilePlot(k,max_steps) = timeBD(k);	
    end
end