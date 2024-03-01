
%% Plot results of several Pareto curves for comparison
clear all
close all
optimization_file='data_optim_10ppm_100C_08122023';
% Figures
bar_charts_pur_rec = true;
desVar = false;
paretoCharts = false;
isosteric_heat = false;

graph_colored_AB = false;
graph_colored_AB_2 = false;
graph_colored_AB_3 = true;

legendPlot = true;

%% input
load(optimization_file,'data_sorbent_optim')

cases = string(fieldnames(data_sorbent_optim));
% define which case:
    k=1;
    current_case = string(cases(k));
materials = string(fieldnames(data_sorbent_optim.(current_case))); % list material names
materials_2 = materials;

color_list = {[0.03,0.28,0.45],[0.11,0.43,0.66],[0.18,0.67,0.72],[0.62,0.83,0.29],[0.86,0.85,0.34],[0.89,0.60,0.15],[0.97,0.46,0.33],     [0.03,0.28,0.45],[0.11,0.43,0.66],[0.18,0.67,0.72],[0.62,0.83,0.29],[0.86,0.85,0.34],[0.89,0.60,0.15],[0.97,0.46,0.33]};
color_shade = {[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00],[0.80,0.80,0.80],[1.00,1.00,1.00]};

m = 1;
for k = 1:(length(materials))
    % 0D model
   filename(m) = strcat('results_',current_case,'_',materials(m),'.txt'); 
   m=m+1;
end

% real names materials
materials = convert_string_name_to_normal ( materials, "name");

%     materials(3) = "APDES-NFC";
%     materials(7) = "Lewatit";
%     materials(8) = "MIL-101(Cr)-PEI-800"; 

%% plot lines to get nice legend
if legendPlot
      
    
    if bar_charts_pur_rec
        figure(3)
        subplot(1,2,2) % choose where the legend should be
        hold on
        for k = 1:(length(materials))
            plot([-100 -101],[-100 -101],'-','linewidth',1,'Color',color_list{k})
        end
    end
    if desVar
        figure(4)
        subplot(2,3,6) % choose where the legend should be
        hold on
        for k = 1:(length(materials))
            plot([-100 -101],[-100 -101],'-','linewidth',1,'Color',color_list{k})
        end
    end
    if paretoCharts
        figure(5)
        hold on
        for k = 1:(length(materials))
            plot([-100 -101],[-100 -101],'-','linewidth',1,'Color',color_list{k})
        end
    end
    
    if graph_colored_AB
       figure(8)
       hold on
       set(gcf, 'Position', [35.666666666666664,301.6666666666666,1256.666666666667,262.6666666666666])
        set(gcf,'color','w');
    end
end

%% calc data
j=0;
m=1;
for k = 1:(length(materials))
    out_0D = read_res(filename(m));
    Pr_0D = out_0D.FrontPts_ProEn(:,2);
    Q_0D = out_0D.FrontPts_ProEn(:,1);
    
    save_Pr_low_0D(k) = Pr_0D(1);
    save_Pr_high_0D(k) = Pr_0D(end);
    
    save_Q_low_0D(k) = Q_0D(1);
    save_Q_high_0D(k) = Q_0D(end);
    
    purity_0D = out_0D.purity_FrontPts;
    recovery_0D = out_0D.recovery_FrontPts;
    
    desVar_0D = out_0D.FrontPts_DesVarEnProd;
    
    tAds_0D = out_0D.tAds;
    tBD_0D = out_0D.tBD;
    tHeat_0D = out_0D.tHeat;
    
    % delete capture rate > 1
            indices2 = find(recovery_0D>1.01);
            indices2 = find(purity_0D>1.01);
    
    if size(Pr_0D,1) == 1
        Pr_0D(2) = out_0D.FrontPts_ProEn(:,1)+out_0D.FrontPts_ProEn(:,1)/100;
        Q_0D(2) = out_0D.FrontPts_ProEn(:,2)+out_0D.FrontPts_ProEn(:,2)/100;

        purity_0D(2) = out_0D.purity_FrontPts+out_0D.purity_FrontPts/100;
        recovery_0D(2) = out_0D.recovery_FrontPts+out_0D.recovery_FrontPts/100;

        desVar_0D(2,:) = out_0D.FrontPts_DesVarEnProd+out_0D.FrontPts_DesVarEnProd/100;

        tAds_0D(2) = out_0D.tAds+out_0D.tAds/100;
        tBD_0D(2) = out_0D.tBD+out_0D.tBD/100;
        tHeat_0D(2) = out_0D.tHeat+out_0D.tHeat/100;
        
    end
       
    %% Plot       
    if bar_charts_pur_rec
        figure(3)
        subplot(1,2,1)
        hold on
        sz = 25;
        errorbar(m,(min(purity_0D)+max(purity_0D))/2,(max(purity_0D)-min(purity_0D))/2,'linewidth',1,'Color',color_list{m})
        scatter(m,min(purity_0D),sz,'^','MarkerEdgeColor', [color_list{m}],'MarkerFaceColor',[color_list{m}])
        scatter(m,max(purity_0D),sz,'^','MarkerEdgeColor', [color_list{m}],'MarkerFaceColor',[color_list{m}])

        subplot(1,2,2)
        hold on
        errorbar(m,(min(recovery_0D)+max(recovery_0D))/2,(max(recovery_0D)-min(recovery_0D))/2,'linewidth',1,'Color',color_list{m})
        scatter(m,min(recovery_0D),sz,'^','MarkerEdgeColor', [color_list{m}],'MarkerFaceColor',[color_list{m}])
        scatter(m,max(recovery_0D),sz,'^','MarkerEdgeColor', [color_list{m}],'MarkerFaceColor',[color_list{m}])
    end
    
    if desVar
        figure(4)
         subplot(2,3,1)
         hold on 
         plot(Pr_0D,desVar_0D(:,1),'-o','Color', [color_list{m}])
         xlabel('CO_2 Productivity (kg/(m^3 h))','FontSize',10)
         ylabel('T_{des} (K)','fontSize',10)
         hold on
         grid on
         box on
         
         subplot(2,3,2)
         hold on 
         plot(Pr_0D,desVar_0D(:,2).*10,'-o','Color', [color_list{m}])
         xlabel('CO_2 Productivity (kg/(m^3 h))','FontSize',10)
         ylabel('p_{vac} (bar)','fontSize',10)
         hold on
         grid on
         box on
         
         subplot(2,3,3)
         hold on 
         plot(Pr_0D,desVar_0D(:,3),'-o','Color', [color_list{m}])
         xlabel('CO_2 Productivity (kg/(m^3 h))','FontSize',10)
         ylabel('V_{feed} (m3/s)','fontSize',10)
         hold on
         grid on
         box on
         
         subplot(2,3,4)
         hold on 
         plot(Pr_0D,tAds_0D,'-o','Color', [color_list{m}])
         xlabel('CO_2 Productivity (kg/(m^3 h))','FontSize',10)
         ylabel('t_{adsorption} (s)','fontSize',10)
         hold on
         grid on
         box on
         
         subplot(2,3,5)
         hold on 
         plot(Pr_0D,tBD_0D,'-o','Color', [color_list{m}])
         xlabel('CO_2 Productivity (kg/(m^3 h))','FontSize',10)
         ylabel('t_{BD} (s)','fontSize',10)
         hold on
         grid on
         box on
         
         subplot(2,3,6)
         hold on 
         plot(Pr_0D,tHeat_0D,'-o','Color', [color_list{m}])
         xlabel('CO_2 Productivity (kg/(m^3 h))','FontSize',10)
         ylabel('t_{heat} (s)','fontSize',10)
         hold on
         grid on
         box on
    end
    
    if paretoCharts
        figure(5)
        hold on
        plot(Pr_0D,Q_0D,'-','linewidth',2,'Color',color_list{m})
    end
       
    vector_lowPr(:,k) = [j j j+1 j+1];
    j=j+2;
    vector_highPr(:,k) = [j-1 j-1 j j];
    vector_ticks(k) = j-1;
    
    vector_Pr2(:,k) = [m-1 m-1 m m]; 
    vector_Pr2_Y(:,k) = [save_Pr_high_0D(k) save_Pr_low_0D(k) save_Pr_low_0D(k) save_Pr_high_0D(k)];
    vector_E2(:,k) = [save_Q_high_0D(k) save_Q_low_0D(k) save_Q_low_0D(k) save_Q_high_0D(k)];
       
    m=m+1;
end

if graph_colored_AB
    figure(8)
    xlim([0 vector_highPr(end,end)+3])
    
    yticks([0.5 1.5])
    yticklabels({'Low Q_{th}','High Pr'})
    
    xticks(vector_ticks)
    xticklabels(materials)
    
    grid off
    
    hold on
    % Productivity
        % low produtivity/energy
            Pr_X = vector_lowPr;
            Pr_Y = [ones(1,length(materials)); zeros(1,length(materials)); zeros(1,length(materials)); ones(1,length(materials))];
            C_Pr = save_Pr_low_0D;
        p_Pr1 = patch(Pr_X,Pr_Y,C_Pr,'LineStyle','none');
        % high productivity/energy
            Pr_X = vector_lowPr;
            Pr_Y = [ones(1,length(materials)).*2; ones(1,length(materials)); ones(1,length(materials)); ones(1,length(materials)).*2];
            C_Pr = save_Pr_high_0D;
        p_Pr2 = patch(Pr_X,Pr_Y,C_Pr,'LineStyle','none');
            
    colormap(gca,brewermap(256,'OrRd'))
    cb1 = colorbar('eastoutside');
        x1=get(gca,'position');
        x=get(cb1,'Position');

    caxis([0 25])
    xlabel(cb1,'Produtvitvity (kg/(m^3 h)') 

    % Energy
    cb2 = newcolorbar('east');
        x2=get(gca,'position');
        x22=get(cb2,'Position');
        x22(4)=x(4); % same length
        x22(2)=x(2); % same y
        set(cb2,'Position',x22)
        set(gca,'position',x2)
    % low produtivity/energy
            E_X = vector_highPr;
            E_Y = [ones(1,length(materials)); zeros(1,length(materials)); zeros(1,length(materials)); ones(1,length(materials))];
            C_E = save_Q_low_0D;
        p_E1 = patch(E_X,E_Y,C_E,'LineStyle','none');
        % high productivity/energy
            E_X = vector_highPr;
            E_Y = [ones(1,length(materials)).*2; ones(1,length(materials)); ones(1,length(materials)); ones(1,length(materials)).*2];
            C_E = save_Q_high_0D;
        p_E2 = patch(E_X,E_Y,C_E,'LineStyle','none');
    colormap(gca,brewermap(256,'*OrRd')) %BuPu, GnBu
    caxis([0 18])
    xlabel(cb2,'Thermal energy (MJ/kg_{CO2})')   
   
    j=0;
    for k=1:length(materials)-1
        j=j+2;
       plot([j j],[0 2],'--','LineWidth',1,'Color',[0,0,0]) % 1.00,1.00,1.00
    end
    j=j+2;
    plot([0 j],[0 0],'LineWidth',0.5,'Color',[0,0,0]) % 1.00,1.00,1.00 
    plot([0 j],[2 2],'LineWidth',0.5,'Color',[0,0,0]) % 1.00,1.00,1.00 
    plot([0 0],[0 2],'LineWidth',0.5,'Color',[0,0,0]) % 1.00,1.00,1.00 
    plot([j j],[0 2],'LineWidth',0.5,'Color',[0,0,0]) % 1.00,1.00,1.00 
    plot([j j+3],[0 0],'LineWidth',0.5,'Color',[1,1,1]) % 1.00,1.00,1.00
    
    plot([0 vector_highPr(end,end)],[1 1],'LineWidth',0.5,'Color',[0,0,0])
end

if graph_colored_AB_2
    figure(9)
    hold on
        Pr_X = vector_Pr2;
        Pr_Y = vector_Pr2_Y;
        C_Pr = vector_E2;
    p_Pr1 = patch(Pr_X,Pr_Y,C_Pr);
    colormap(gca,brewermap(256,'OrRd')) %
    cb1 = colorbar('eastoutside');
    caxis([0 30])
    xlabel(cb1,'Q_{th}')
    ylabel('Productivity (kg/(m^3 h)')
    
    xticks(linspace(0.5,length(materials)-0.5,length(materials)))
    xticklabels(materials)
    grid off
    box on
    set(gcf, 'Position', [36,181,1213,383])
    set(gcf,'color','w');
end

if graph_colored_AB_3
    figure(11)
    set(gcf, 'Position', [272,258,799,293])
    hold on
        Pr_Y = vector_Pr2;
        Pr_X = vector_Pr2_Y;
        C_Pr = vector_E2;
    p_Pr1 = patch(Pr_X,Pr_Y,C_Pr);
    colormap(gca,brewermap(256,'OrRd')) %
    cb1 = colorbar('eastoutside');
    caxis([4 22])
    xlabel(cb1,'Q_{th} (MJ/kg_{CO2})')

    yticks(linspace(0.5,length(materials)-0.5,length(materials)))
    yticklabels(materials)
%     yticklabels(["Ca-X","Carbon","APDES-NFC","Cr-MIL(101)","CuBTC","Exemplary","Lewatit","MIL-101(Cr)-PEI-800","MIL-101","MIL-53(Al)","Zeolite Na-LSX","Zn-DABCO"])
    
    xlabel('Purity CH4')
    grid on
    box on
    
    set(gcf,'color','w');
    
    ylim([-0.5 length(materials)+0.5])  
end



%% Figure properties
if bar_charts_pur_rec
    figure(3)
    subplot(1,2,1)
        xlim([0 (length(materials)+1)])
        ylim([0 1])        
        grid on
        box on
        xlabel('materials')
        ylabel('Purity')
    subplot(1,2,2)
        xlim([0 (length(materials)+1)])
        ylim([0 1])
        grid on
        box on
        xlabel('materials')
        ylabel('Recovery')
        if legendPlot
            legend(materials);
        end
    set(gcf,'color','w');
end

if desVar
    figure(4)
    set(gcf, 'Position', [50, 45, 1200, 600])
    set(gcf,'color','w');
        if legendPlot
            subplot(2,3,6)
            xlim([0 12])
            ylim([0 6000])
            legend(materials);
        end
    set(gcf,'color','w');
end

if paretoCharts
    figure(5)
    set(gcf, 'Position', [335 49 769 582])
    grid on
    box on
    xlabel('Productivity, Pr (kg/(m^3 h))')
    ylabel('Thermal energy (MJ/kg_{CO2})')
    set(gcf,'color','w');
    if legendPlot
        xlim([0 25])
        ylim([0 45])
        legend(materials);
    end
end

%% reading results
function out = read_res(filename)
    fileID = fopen(filename, 'r');
    results_temp = cell2mat(textscan(fileID, repmat('%f ',[1 12]),'Delimiter',','));
    fclose(fileID);

    % results_temp:
    % 1: temperature [K],   2: vacuum pressure [MPa], 3: Q_spec [MJ/kgCO2eq],
    % 4: E_spec [MJ/kgCO2eq],  5: Pr [kg/(m^3 h)], 6: purity CO2,
    % 7: Purity CH4, 8: recovery CH4, 9: WC_CH4 [mol_CO2_recovered/kg_sorbent]
    % 10: t Ads, 11: t BD, 12: t Heating

%     % delte capture rate >1
%     [RowNrs,ColNrs] = find(results_temp(:,9)>1.02);
%     results_temp(RowNrs,:)=[];
    
    % delte Pr < 0.05
%     [RowNrs,ColNrs] = find(results_temp(:,6)<0.05);
%     results_temp(RowNrs,:)=[];

    % ordering results for purity CH4
    [results_temp,idx] = sortrows(results_temp,7);
    
    % ordering results for recovery CH4
%     [results_temp,idx] = sortrows(results_temp,8);
    
    % ordering results for Energy Consumption (MJ/kgCO2eq)
%     [results_temp,idx] = sortrows(results_temp,3);
  

    %% PARETO DATA
    Objectives_EnProd = [results_temp(:,3) results_temp(:,7) results_temp(:,8)]; % thermalenergy + purity + recovery

    front_ProEn= paretoGroup([Objectives_EnProd(:,1),-Objectives_EnProd(:,2), -Objectives_EnProd(:,3)]);
    out.FrontPts_ProEn = [Objectives_EnProd(front_ProEn,1) Objectives_EnProd(front_ProEn,2) Objectives_EnProd(front_ProEn,3)];
    out.FrontPts_DesVarEnProd = [results_temp(front_ProEn,1) results_temp(front_ProEn,2)];

    out.purity_FrontPts = results_temp(front_ProEn,7);
    out.recovery_FrontPts = results_temp(front_ProEn,8);
    
    out.tAds = results_temp(front_ProEn,10);
    out.tBD = results_temp(front_ProEn,11);
    out.tHeat = results_temp(front_ProEn,12);
end

