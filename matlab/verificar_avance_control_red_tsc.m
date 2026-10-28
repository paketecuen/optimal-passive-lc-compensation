function report = verificar_avance_control_red_tsc(campaign)
%VERIFICAR_AVANCE_CONTROL_RED_TSC Refinamiento, causalidad y guardia real.
folder=fullfile(fileparts(mfilename('fullpath')),'resultados','control_red','pruebas_auxiliares');
if ~isfolder(folder), mkdir(folder); end
base=campaign.registros{2,5}; cfg=base.config;
% Kernel de evaluacion temporal y contraccion en la metrica de energia.
c=cfg.ControlRed; kernelError=max(abs(sum(c.Kcota.^2,2)-1));
assert(kernelError<1e-11,'TSCRed:Kernel','La cota de derivada requiere norma cardinal unitaria.');
growth=0;
for k=1:4
    active=true(1,6); active([3 5])=logical(cfg.p.bits(k,:));
    for j=1:size(c.topologias{k}.F,3)
        F=c.topologias{k}.F(active,active,j);
        growth=max(growth,norm(F,2)-1);
    end
end
assert(growth<1e-9,'TSCRed:Contraccion','La propagacion crea energia libre.');
% Zg=0: coste acoplado y coste con tension congelada deben coincidir.
ideal=simular_tsc_lazo_cerrado(configurar_control_red_tsc(configurar_fuente_tsc( ...
    configurar_tsc_lazo_cerrado(Duracion_s=.12),R_ohm=0,L_H=0)));
verifyIdeal=verificar_control_red_tsc(ideal);
idealError=max(abs(ideal.decisiones.CostesRed-ideal.decisiones.CostesBornes),[],'all');
assert(idealError<1e-8,'TSCRed:LimiteIdeal','No se recupera el coste de tension impuesta.');
save(fullfile(folder,'limite_ideal.mat'),'ideal','verifyIdeal','-v7.3');
% Dos refinamientos distintos: solver de planta y malla de planificacion.
fineCfg=cfg; fineCfg.PasoMaximo_s=cfg.PasoMaximo_s/2;
fineCfg.opciones.PasoMaximo_s=fineCfg.PasoMaximo_s;
fine=simular_tsc_lazo_cerrado(fineCfg); verifyFine=verificar_control_red_tsc(fine);
assert(isequal(base.eventos.Tipo,fine.eventos.Tipo) && isequal(base.eventos.Rama,fine.eventos.Rama));
solverError=abs(fine.resumen.I_fuente_RMS_A/base.resumen.I_fuente_RMS_A-1);
solverPeak=abs(fine.resumen.VC_pico_V-base.resumen.VC_pico_V);
assert(solverError<1e-4 && solverPeak<.1,'TSCRed:Refinamiento','Diferencia material al refinar el solver.');
save(fullfile(folder,'refinamiento_solver.mat'),'fine','verifyFine','-v7.3');
fineCfg=configurar_control_red_tsc(cfg,MaxVC_V=c.MaxVC_V,MaxI_A=c.MaxI_A, ...
    ModeloCoste=c.ModeloCoste,Margen=c.Margen,PasosPlan=2*c.PasosPlan,MuestrasCota=c.MuestrasCota);
grid=simular_tsc_lazo_cerrado(fineCfg); verifyGrid=verificar_control_red_tsc(grid);
assert(verifyGrid.control.Limites_cumplidos_en_traza,'TSCRed:Malla','Se incumplen limites al refinar la malla.');
gridError=abs(grid.resumen.I_fuente_RMS_A/base.resumen.I_fuente_RMS_A-1);
save(fullfile(folder,'refinamiento_malla.mat'),'grid','verifyGrid','-v7.3');
% No debe haber anticipacion ni de las rampas ni del valor futuro de carga.
args=cfg.opciones; args.EscalaCarga(3:4)=[1.8 .3]; ac=namedargs2cell(args);
futureCfg=configurar_control_red_tsc(configurar_fuente_tsc(configurar_tsc_lazo_cerrado(ac{:}), ...
    R_ohm=cfg.Fuente.R_ohm,L_H=cfg.Fuente.L_H),MaxVC_V=c.MaxVC_V,MaxI_A=c.MaxI_A);
future=simular_tsc_lazo_cerrado(futureCfg); verifyFuture=verificar_control_red_tsc(future);
mask=base.decisiones.Tiempo_s<.4-1e-10;
prefixError=max(abs(base.decisiones.CostesRed(mask,:)-future.decisiones.CostesRed(mask,:)),[],'all');
assert(prefixError<1e-9 && isequaln(base.decisiones.EstadoSolicitado(mask),future.decisiones.EstadoSolicitado(mask)), ...
    'TSCRed:Futuro','La carga futura altera el pasado.');
save(fullfile(folder,'futuro.mat'),'future','verifyFuture','-v7.3');
% Corriente como unica restriccion finita: comprueba que no es una etiqueta
% inactiva oculta tras la cota de tension.
currentCfg=configurar_control_red_tsc(cfg,MaxVC_V=[Inf Inf],MaxI_A=[8 8]);
currentOnly=simular_tsc_lazo_cerrado(currentCfg); verifyCurrent=verificar_control_red_tsc(currentOnly);
assert(verifyCurrent.control.ON_admitidos>0 && verifyCurrent.control.Rechazos_plan>0 && ...
    verifyCurrent.control.Limites_cumplidos_en_traza,'TSCRed:LimiteI','La prueba de limite I no es efectiva.');
save(fullfile(folder,'solo_limite_corriente.mat'),'currentOnly','verifyCurrent','-v7.3');
% Perturbacion deliberada despues de decidir y antes del ON programado.
% No representa un escenario nominal ni una especificacion de red.
abruptCfg=configurar_control_red_tsc(configurar_fuente_tsc( ...
    configurar_tsc_lazo_cerrado(Duracion_s=.12,EscalaCarga=[.45 8], ...
    InstantesCarga_s=[0 .021],Transicion_s=.004),R_ohm=1.5,L_H=3e-3));
abrupt=simular_tsc_lazo_cerrado(abruptCfg); verifyAbrupt=verificar_control_red_tsc(abrupt);
assert(any(~verifyAbrupt.admisiones.Admitido & verifyAbrupt.admisiones.Razon==1), ...
    'TSCRed:Revalidacion','La perturbacion no ejercita el rechazo por estado actual.');
save(fullfile(folder,'revalidacion_ON.mat'),'abrupt','verifyAbrupt','-v7.3');
% Alteraciones de evidencia: coste, estado del evento y adquisicion futura.
changed=base; changed.decisiones.CostesRed(1,1)=changed.decisiones.CostesRed(1,1)+.1;
mustFail(changed,'TSCRed:Costes');
changed=base; changed.ventanas_decision{1}.t(end)=changed.decisiones.Tiempo_s(1)+cfg.TsSupervision;
mustFail(changed,'TSCRed:Causalidad');
changed=base; original=changed.guardias{1};
ts=timeseries(original.Data,original.Time); v=ts.Data;
if size(v,2)==15, v(:,10)=v(:,10)+1; else, v(10,:,:)=v(10,:,:)+1; end
ts.Data=v; changed.guardias{1}=ts;
mustFail(changed,'TSCRed:EstadoEvento');
report=struct('passed',true,'Kernel_error',kernelError,'Crecimiento_operador',growth, ...
    'Coste_limite_ideal_error_A2',idealError,'Solver_error_RMS_relativo',solverError, ...
    'Solver_error_pico_V',solverPeak,'Malla_error_RMS_relativo',gridError, ...
    'Malla_pico_V',grid.resumen.VC_pico_V,'Decisiones_causales',nnz(mask), ...
    'Error_costes_previos_A2',prefixError,'I_solo_pico_A',currentOnly.resumen.I_rama_pico_A, ...
    'I_solo_bloqueos',verifyCurrent.control.Rechazos_plan, ...
    'Rechazos_ON_por_estado',nnz(~verifyAbrupt.admisiones.Admitido & verifyAbrupt.admisiones.Razon==1), ...
    'Pruebas_corrupcion',3);
fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
crear_modelo_tsc_lazo_cerrado(campaign.registros{1,5}.config);
disp(report);
end
function mustFail(r,id)
try
    verificar_control_red_tsc(r,InspeccionarModelo=false);
catch err
    assert(strcmp(err.identifier,id),'TSCRed:Defensa','Esperado %s; observado %s.',id,err.identifier); return
end
error('TSCRed:Defensa','Corrupcion no detectada: %s.',id);
end
