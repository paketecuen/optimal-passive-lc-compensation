function cfg = configurar_fuente_tsc(cfg,opciones)
%CONFIGURAR_FUENTE_TSC Anade e(t) detras de Rg+Lg*D al caso existente.
% Los 230 V RMS corresponden a e, no a la tension en bornes u.
% La carga sigue siendo una inyeccion sintetica independiente de u.
arguments
    cfg (1,1) struct
    opciones.R_ohm (1,1) double {mustBeNonnegative,mustBeFinite} = .4
    opciones.L_H (1,1) double {mustBeNonnegative,mustBeFinite} = .8e-3
end
cfg.Fuente=struct('R_ohm',opciones.R_ohm,'L_H',opciones.L_H);
cfg.NombreModelo="tsc_fuente_impedante";
% Arranque periodico de la red y Cp, con ambas LC abiertas. No se eliminan
% los transitorios posteriores ni se precargan las ramas conmutables.
reference=resolver_periodico_fuente_tsc(cfg,cfg.opciones.EscalaCarga(1),[0 0]);
cfg.vCp0=reference.vCp(1); cfg.Fuente.i0_A=reference.iFuente(1);
cfg.Fuente.referenciaInicial=reference;
cfg.alcance="Solo simulacion: fuente Thevenin R-L, carga exogena sintetica, valvula idealizada y componentes fijos. El supervisor conserva su prediccion de tension medida repetida; no se afirma optimalidad con impedancia de fuente.";
end
