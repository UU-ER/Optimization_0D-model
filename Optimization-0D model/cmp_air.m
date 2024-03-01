% updated 21.05.2020 Alexa
function E= cmp_air(p_amb,densityAir,FeedAdsorptionStep,MMCO2) % MPa, kg/m3, mol, g/mol 
                       
    eff = 0.6;
    
    dp = 0.00001; % bar

    V_air = FeedAdsorptionStep*MMCO2/1000/densityAir; % m3
       
    E = 1/eff*1/(p_amb*1e6-dp*1e5)*V_air; % (J) energy per unit
%%% constant pressure drop, and using the MMAir instead of CO2
end