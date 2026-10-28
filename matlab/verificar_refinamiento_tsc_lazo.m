function informe = verificar_refinamiento_tsc_lazo(base,opciones)
%VERIFICAR_REFINAMIENTO_TSC_LAZO Refinamiento del solver y causalidad dinamica.
% La segunda carga coincide hasta 0.4s y cambia despues. No se comunica el
% cambio futuro al supervisor; se exige el mismo prefijo de decisiones.
arguments
    base (1,1) struct
    opciones.Exportar (1,1) logical = true
end
assert(base.config.Modo=="ga" && base.config.Duracion_s==.8, ...
    'TSCLazo:BaseAuxiliar','Se requiere la corrida GA de 0.8 s.');
args=base.config.opciones; args.PasoMaximo_s=args.PasoMaximo_s/2;
ac=namedargs2cell(args); cfg=configurar_tsc_lazo_cerrado(ac{:});
refinado=simular_tsc_lazo_cerrado(cfg); auditarRef=verificar_tsc_lazo_cerrado(refinado);
e=base.eventos; f=refinado.eventos;
assert(isequal(e.Tipo,f.Tipo) && isequal(e.Rama,f.Rama), ...
    'TSCLazo:ConvergenciaEventos','El refinamiento cambia la secuencia de maniobras.');
errorT=max([0;abs(e.Tiempo_s-f.Tiempo_s)]);
errorRMS=abs(refinado.resumen.I_fuente_RMS_A/base.resumen.I_fuente_RMS_A-1);
assert(errorT<1e-5 && errorRMS<2e-4,'TSCLazo:Convergencia', ...
    'La trayectoria no converge al reducir el paso del solver.');
args=base.config.opciones; args.EscalaCarga=[.45 1.35 .8 .9];
ac=namedargs2cell(args); cfg=configurar_tsc_lazo_cerrado(ac{:});
futuro=simular_tsc_lazo_cerrado(cfg); auditarFuturo=verificar_tsc_lazo_cerrado(futuro);
d=base.decisiones(base.decisiones.Tiempo_s<.4-1e-10,:);
df=futuro.decisiones(futuro.decisiones.Tiempo_s<.4-1e-10,:);
assert(isequal(d.EstadoSolicitado,df.EstadoSolicitado) && isequal(d.Cambio,df.Cambio), ...
    'TSCLazo:CargaFutura','La carga futura altera decisiones previas.');
v1=[d.Coste00,d.Coste01,d.Coste10,d.Coste11];
v2=[df.Coste00,df.Coste01,df.Coste10,df.Coste11];
errorCoste=max(abs(v1-v2),[],'all');
assert(errorCoste<1e-7,'TSCLazo:CosteFuturo','La carga futura altera costes previos.');
eb=e(e.Tiempo_s<.39,:); ef=futuro.eventos(futuro.eventos.Tiempo_s<.39,:);
assert(isequal(eb.Tipo,ef.Tipo) && isequal(eb.Rama,ef.Rama) && ...
    max([0;abs(eb.Tiempo_s-ef.Tiempo_s)])<1e-8, ...
    'TSCLazo:ManiobraFutura','La carga futura altera maniobras previas.');
informe=struct('passed',true,'paso_base_s',base.config.PasoMaximo_s, ...
    'paso_refinado_s',refinado.config.PasoMaximo_s,'error_instantes_s',errorT, ...
    'error_RMS_relativo',errorRMS,'error_pico_VC_V', ...
    abs(refinado.resumen.VC_pico_V-base.resumen.VC_pico_V), ...
    'decisiones_prefijo_comprobadas',height(d),'error_costes_prefijo_A2',errorCoste, ...
    'auditoria_refinada',auditarRef,'auditoria_futuro',auditarFuturo);
if opciones.Exportar
    folder=fullfile(fileparts(mfilename('fullpath')),'resultados','tsc_lazo_cerrado', ...
        base.config.CasoTension,'pruebas_auxiliares');
    if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,'refinamiento_y_causalidad.mat'),'informe','refinado','futuro','-v7.3');
    fid=fopen(fullfile(folder,'refinamiento_y_causalidad.json'),'w'); assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(informe,PrettyPrint=true));
end
disp(informe);
end
