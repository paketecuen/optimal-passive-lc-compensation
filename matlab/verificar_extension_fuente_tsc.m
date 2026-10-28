function informe = verificar_extension_fuente_tsc(campana)
%VERIFICAR_EXTENSION_FUENTE_TSC Limites R/L, estado periodico y causalidad.
arguments
    campana (1,1) struct
end
folder=fullfile(fileparts(mfilename('fullpath')),'resultados','fuente_impedante','pruebas_auxiliares');
if ~isfolder(folder), mkdir(folder); end
% Zg=0 debe recuperar la campaña anterior sin modificarla ni regenerarla.
old=load(fullfile(fileparts(mfilename('fullpath')),'resultados', ...
    'tsc_lazo_cerrado','moderada','lazo_cerrado.mat'),'saved');
old=old.saved; ideal=campana.registros(1,:);
errRMS=zeros(1,6); errVC=errRMS;
for k=1:6
    errRMS(k)=abs(ideal{k}.resumen.I_fuente_RMS_A/old.registros{k}.resumen.I_fuente_RMS_A-1);
    errVC(k)=abs(ideal{k}.resumen.VC_pico_V-old.registros{k}.resumen.VC_pico_V);
    assert(isequal(ideal{k}.eventos.Tipo,old.registros{k}.eventos.Tipo) && ...
        isequal(ideal{k}.eventos.Rama,old.registros{k}.eventos.Rama), ...
        'TSCFuente:LimiteIdeal','Cambia la secuencia de maniobras para Zg=0.');
end
assert(max(errRMS)<1e-5 && max(errVC)<.02,'TSCFuente:LimiteIdeal','No se recupera la referencia ideal.');
% La matriz temporal se contrasta con la planta continua sin conmutaciones.
resistances=[.4,0,.4]; inductances=[0,.8e-3,.8e-3]; stationary=cell(3,1); rows=struct([]);
for k=1:3
    cfg=configurar_fuente_tsc(configurar_tsc_lazo_cerrado(Modo="fijo",BitsFijos=[0 0], ...
        Duracion_s=.08,EscalaCarga=.45,InstantesCarga_s=0), ...
        R_ohm=resistances(k),L_H=inductances(k));
    stationary{k}=simular_tsc_lazo_cerrado(cfg); verificar_fuente_tsc(stationary{k});
    a=stationary{k}.trazas; ref=resolver_periodico_fuente_tsc(cfg,.45,[0 0]);
    predicted=interpolar_periodico([ref.u,ref.iFuente,ref.vCp],cfg.p.omega*a.t);
    delta=max(abs([a.u,a.iFuente,a.vCp]-predicted),[],1);
    assert(delta(1)<1e-3 && delta(2)<1e-4 && delta(3)<1e-3, ...
        'TSCFuente:Estacionario','La planta no reproduce la solucion temporal periodica.');
    assert(ref.Residuo_KVL_V<1e-8 && ref.Residuo_KCL_A<1e-9 && abs(ref.Balance_W)<1e-8, ...
        'TSCFuente:Operador','La solucion periodica no satisface las leyes del circuito.');
    rows=[rows;struct('R_ohm',resistances(k),'L_mH',1e3*inductances(k), ...
        'Error_U_V',delta(1),'Error_Is_A',delta(2),'Error_VCp_V',delta(3))]; %#ok<AGROW>
end
% Mayor impedancia: refinamiento temporal sin cambiar ni fuente ni control.
base=campana.registros{3,1}; cfg=base.config;
cfg.PasoMaximo_s=cfg.PasoMaximo_s/2; cfg.opciones.PasoMaximo_s=cfg.PasoMaximo_s;
refined=simular_tsc_lazo_cerrado(cfg); verificar_fuente_tsc(refined);
assert(isequal(base.eventos.Tipo,refined.eventos.Tipo) && isequal(base.eventos.Rama,refined.eventos.Rama), ...
    'TSCFuente:Refinamiento','La secuencia de maniobras cambia al refinar.');
fine=struct('Error_RMS_relativo',abs(refined.resumen.I_fuente_RMS_A/base.resumen.I_fuente_RMS_A-1), ...
    'Error_pico_VC_V',abs(refined.resumen.VC_pico_V-base.resumen.VC_pico_V), ...
    'Error_eventos_s',max(abs(base.eventos.Tiempo_s-refined.eventos.Tiempo_s)));
assert(fine.Error_RMS_relativo<1e-4 && fine.Error_pico_VC_V<.1 && fine.Error_eventos_s<1e-5, ...
    'TSCFuente:Refinamiento','El refinamiento produce diferencias materiales.');
args=base.config.opciones; args.EscalaCarga(3:4)=[1.8 .3]; ac=namedargs2cell(args);
cfg=configurar_fuente_tsc(configurar_tsc_lazo_cerrado(ac{:}), ...
    R_ohm=base.config.Fuente.R_ohm,L_H=base.config.Fuente.L_H);
future=simular_tsc_lazo_cerrado(cfg); verificar_fuente_tsc(future);
mask=base.decisiones.Tiempo_s<.4-1e-10;
d1=base.decisiones(mask,:); d2=future.decisiones(mask,:);
costNames={'Coste00','Coste01','Coste10','Coste11'};
errCosts=max(abs(d1{:,costNames}-d2{:,costNames}),[],'all');
e1=base.eventos(base.eventos.Tiempo_s<.4-1e-10,:); e2=future.eventos(future.eventos.Tiempo_s<.4-1e-10,:);
assert(errCosts<1e-9 && isequal(d1.EstadoSolicitado,d2.EstadoSolicitado) && ...
    isequal(e1.Tipo,e2.Tipo) && isequal(e1.Rama,e2.Rama) && ...
    max(abs(e1.Tiempo_s-e2.Tiempo_s))<1e-9,'TSCFuente:Causalidad','El futuro altera decisiones previas.');
% La envolvente de estado incluye Lg/Cp; no se usa la cota de una LC aislada.
startup=cell(3,1);
for k=1:3, startup{k}=auditar_arranque_fuente_tsc(campana.registros{k,1}); end
% Una perdida de energia de red y una caida Rg falsa deben detectarse.
corrupt=base; corrupt.trazas.integralFuente_J(end)=corrupt.trazas.integralFuente_J(end)+1;
mustFail(corrupt,'TSCFuente:Energia');
corrupt=base; corrupt.trazas.e=corrupt.trazas.e+1;
mustFail(corrupt,'TSCFuente:KVL');
informe=struct('passed',true,'Limite_ideal_error_RMS',max(errRMS), ...
    'Limite_ideal_error_VC_V',max(errVC),'Estacionarios',struct2table(rows), ...
    'Refinamiento',fine,'Decisiones_causales',height(d1),'Error_costes_previos_A2',errCosts, ...
    'Arranque',{startup},'Corrupciones_detectadas',2);
save(fullfile(folder,'extension_verificada.mat'),'informe','stationary','refined','future','-v7.3');
writetable(informe.Estacionarios,fullfile(folder,'estacionarios.csv'));
small=informe;
for k=1:3, small.Arranque{k}=rmfield(small.Arranque{k},{'t','energiaDesajuste_J','perdidasDesajuste_W'}); end
fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(small,PrettyPrint=true));
crear_modelo_tsc_lazo_cerrado(campana.registros{2,1}.config);
disp(small);
end

function mustFail(r,identifier)
try
    verificar_fuente_tsc(r,InspeccionarModelo=false);
catch err
    assert(strcmp(err.identifier,identifier),'TSCFuente:Defensa','Identificador inesperado: %s.',err.identifier);
    return
end
error('TSCFuente:Defensa','No se detecta la corrupcion %s.',identifier);
end
