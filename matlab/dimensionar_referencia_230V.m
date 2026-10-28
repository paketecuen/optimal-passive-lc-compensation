function r = dimensionar_referencia_230V(PNominal_W)
%DIMENSIONAR_REFERENCIA_230V Escala el ejemplo elegido a 230 V RMS y potencia.
% Referencia teorica, NO seleccion de componentes/valvulas comerciales.
% PNominal_W es obligatorio: la tension por si sola no dimensiona el banco.
arguments
    PNominal_W (1,1) double {mustBePositive,mustBeFinite}
end
base=caso_m3; f=familia_optima_m3(base.u,base.i,base.D,base.omega);
b=banco_familia_m3(f,27/13);
U0=sqrt(mean(base.u.^2)); P0=mean(base.u.*base.i);
lambda=230/U0; k=lambda^2*P0/PNominal_W;
bank=escalar_banco_m3(b,1/k);
row=table(PNominal_W,230,k,bank.Cp_F*1e6,bank.Ls_H*1e3,bank.Cs_F*1e6, ...
    bank.f_serie_Hz,PNominal_W/230,'VariableNames', ...
    {'P_carga_W','U_RMS_V','Escala_impedancia','Cp_uF','Ls_mH','Cs_uF','f_serie_Hz','I_ideal_A'});
r=struct('tabla',row,'banco',bank,'escala_tension',lambda,'escala_impedancia',k, ...
    'parametros_simulink',struct('Amplitud',lambda,'CargaEscala',1/k,'Banco',bank), ...
    'alcance',"Escalado de referencia ideal; no dimensiona un montaje comercial ni sus protecciones");
end
