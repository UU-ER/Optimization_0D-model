
clear all
close all
rng(1,'twister')   
    
% Case 1
k=1;

%% optimize multi objective: Productivity and energy
% X: Tdes, pvac
MultiObj.var_min = [363, 0.1/10, 7e-10]; % K, MPa
MultiObj.var_max = [373, 0.8/10,  8.9e-06]; % K, MPa
file_optimization = 'data_sorbent_optim_21092023.mat';
load(file_optimization,'data_sorbent_optim') % If this names is changed, it also needs to be updated in 0D model

cases = string(fieldnames(data_sorbent_optim));
l=0;
current_case = string(cases(k));
names_mat = string(fieldnames(data_sorbent_optim.(current_case)));

for m = 1:size(names_mat,1)
    l=l+1;
    current_mat = string(names_mat(m));
    Tdes = data_sorbent_optim.(current_case).(current_mat).Tdes;
    yCO2 = data_sorbent_optim.(current_case).(current_mat).pvac;
    file_name = string(sprintf('results_case_%2.2d_Pa_%2.2d_C_%s.txt',yCO2, Tdes,current_mat));
    fid = fopen(file_name,'w');
    fclose(fid);

    MultiObj.fun = @(X) eval_0D(X,k,m,file_name);
    MultiObj.nVar = 3;

    % Parameters
    params.Np = 30;        % Population size 30
    params.Nr = 30;        % Repository size 30
    params.maxgen = 30;    % Maximum number of generations 30
    params.W = 0.4;         % Inertia weight
    params.C1 = 2;          % Individual confidence factor
    params.C2 = 2;          % Swarm confidence factor
    params.ngrid = 20;      % Number of grids in each dimension
    params.maxvel = 5;      % Maxmium vel in percentage
    params.u_mut = 0.5;     % Uniform mutation percentage

    % MOPSO
    REP = MOPSO(params,MultiObj);
end

