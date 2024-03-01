%% Optimization by Thompson sampling efficient multiobjective optimization” (TSEMO)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2020-23-05.
% Implemented for 0D model optimization by Kian Karimi, Utrecht University,
% Group of Dr. Matteo Gazzani
%% Step 1: Set path to folder
% Matlab Home, Set Path, Add folder with subfolders, select "TS-EMO-master"
% folder

%% Step 3: Specify problem
no_outputs = 2;               % number of objectives [Purity, Recovery]
no_inputs  = 3;               % number of decision variables [Tdes, Vacuum Pressure, Feed Velocity]
lb = [303, 0.1/10 253];    % K, MPa, K % define lower bound on decision variables, [lb1,lb2,...]
ub = [373, 0.8/10 293];    % K, MPa, K % define upper bound on decision variables, [ub1,ub2,...]

file_optimization = 'data_optim_500ppm_100C_15122023.mat';
load(file_optimization,'data_sorbent_optim') % If this names is changed, it also needs to be updated in 0D model

ppm_CO2 = 1800; % CO2 concentration in feed 1800 ppm
Tamb=293;

cases = string(fieldnames(data_sorbent_optim));
l=0;

c=1;% case 1
current_case = string(cases(c));
names_materials = string(fieldnames(data_sorbent_optim.(current_case)));
error_opt=0;
for m = 1:size(names_materials,1)
    lb = [303, 0.1/10 253];    % K, MPa, K % define lower bound on decision variables, [lb1,lb2,...]
    ub = [373, 0.8/10 293];    % K, MPa, K % define upper bound on decision variables, [ub1,ub2,...]

    l=l+1;
    current_material = string(names_materials(m));
    name_real = convert_string_name_to_normal (current_material,'name')
    Tdes = data_sorbent_optim.(current_case).(current_material).Tdes;
    Concentration_CH4 = data_sorbent_optim.(current_case).(current_material).ppm_CH4;
    file_name = string(sprintf('results_case_%2.2d_ppm_%2.2d_C_%s.txt',Concentration_CH4, Tdes,current_material));
    fid = fopen(file_name,'w');
    fclose(fid);

    % Step 2: Create function file with multiobjective function to be optimized
    f = @(X)eval_0D(X,c,m,file_name,file_optimization,ppm_CO2);
   
    TSEMO_find_limit_Temp;
    if lb(1)<ub(1)
    TSEMO_find_limit_Vacuum;
  
    if ismember(1,Y(:,1:2)~=inf)

        if -1<=Y(:,1:2)
            %% Step 4: Generate initial dataset
            dataset_size = 5*no_inputs;             % initial dataset size
            X = lhsdesign(dataset_size,no_inputs);  % Latin hypercube design

            Y = zeros(dataset_size,no_outputs);     % corresponding matrix of response data

            for k = 1:size(X,1)
                X(k,:) = X(k,:).*(ub-lb)+lb;        % adjustment of bounds
                Y(k,:) = f(X(k,:));                 % calculation of response data
            end

            if Y(:,1:2)<=0
                if -1<=Y(:,1:2)

                    opt = TSEMO_options;             % call options for solver, see TSEMO_options file to adjust
                    opt.maxeval = 40;                % number of function evaluations before termination
                    opt.NoOfBachSequential = 1;      % number of function evaluations per iteration
                    % Total number of iterations = opt.maxeval/opt.NoOfBachSequential

                    %% Step 5: Start algorithm to find Pareto front
                    try
                    [Xpareto,Ypareto,X,Y,XparetoGP,YparetoGP,YparetoGPstd,hypf] = TSEMO_V4(f,X,Y,lb,ub,opt);

                    % INPUTS
                    % f denotes the function to be optimized
                    % X and Y are the initial datasets to create a surrogate model
                    % lb and ub are the lower and upper bound of the decision variables
                    % opt is the option structure of the algorithm

                    % OUTPUTS
                    %   Xpareto and Ypareto correspond to the current best Pareto set and Pareto
                    %   front respectively.
                    %   X and Y are the complete dataset of the decision variables and
                    %   the objectives respectively.
                    %   XparetoGP and YparetoGP represent the Pareto set and Pareto front of the
                    %   final Gaussian process model within the algorithm. It is recommended to
                    %   use these as final result for problems with measurement noise.
                    %   YparetoGPstd denotes the standard deviations of the predictions of
                    %   the GP pareto front YparetoGP
                    %   hypf represents the final hyperparameters found for analysis

                    % For each iteration the current iteration number is displayed, the
                    % predicted hypervolume improvement and the time taken.

                    % TS-EMO creates a log file named "TSEMO_log.txt" that contains all relevant information
                    % over the entire algorithm run.

                    %% Step 6: Visualise results
                    figure
                    % hold on
                    plot(-Y(1:dataset_size,1),-Y(1:dataset_size,2),'.','MarkerSize',14)
%                     scatter3(-Y(1:dataset_size,1),-Y(1:dataset_size,2),Y(1:dataset_size,3),'.','blue')
                    hold on
                    plot(-Y(dataset_size+1:end,1),-Y(dataset_size+1:end,2),'x','MarkerSize',8,'LineWidth',2)
%                     scatter3(-Y(dataset_size+1:end,1),-Y(dataset_size+1:end,2),Y(dataset_size+1:end,3),'x',"red")
                    hold on
                    plot(-Ypareto(:,1),-Ypareto(:,2),'O','MarkerSize',8,'LineWidth',2)
%                     scatter3(-Ypareto(:,1),-Ypareto(:,2),Ypareto(:,3),'O','green')
                    hold on
                    % plot(-YparetoGP(:,1),-YparetoGP(:,2),'x','MarkerSize',8,'LineWidth',2)
                    % errorbar(-YparetoGP(:,1),-YparetoGP(:,2),-YparetoGPstd(:,2),-YparetoGPstd(:,2),-YparetoGPstd(:,1),-YparetoGPstd(:,1),'O','MarkerSize',8,'LineWidth',0.5)
                    legend('Initial LHC','Algorithm','Pareto front','GP Pareto front','Location','Northeast')
                    title(name_real)
                    %                     xlabel('Working Capacity(molCH4/KgAdsorbent)')
                    %                     ylabel('Exergy Consumption (MJ/KgCH4)')
                    xlabel('Purity')
                    ylabel('Recovery')
                    %                     zlabel('Energy Consumption (MJ/kgCO2eq)')
                    %                     zlabel('Exergy (MJ/kgCO2eq)')
                    pause(0.2)
                    
                    catch
                        error_opt=error_opt+1;
                    end
                else
                end
            else
            end
        else
        end
    else
    end
    end
end