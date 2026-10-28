function r = experimento_tsc_lazo_cerrado(opciones)
%EXPERIMENTO_TSC_LAZO_CERRADO Seis estrategias calculadas dentro de Simulink.
% La referencia MATLAB anterior solo se lee despues para comparar resultados.
arguments
    opciones.CasoTension (1,1) string {mustBeMember(opciones.CasoTension,["original","moderada"])} = "moderada"
    opciones.PNominal_W (1,1) double {mustBePositive,mustBeFinite} = 1000
    opciones.PasoMaximo_s (1,1) double {mustBePositive,mustBeFinite} = 2.5e-5
    opciones.Visible (1,1) string {mustBeMember(opciones.Visible,["on","off"])} = "off"
    opciones.Exportar (1,1) logical = true
end
labels=["Supervisor GA","Reactiva fundamental","Fijo 00","Fijo media","Fijo unidad","Fijo ambas"];
fixed=[0 0;1 0;0 1;1 1]; runs=cell(6,1); rows=struct([]); audits=cell(6,1);
root=fileparts(mfilename('fullpath'));
folder=fullfile(root,'resultados','tsc_lazo_cerrado',opciones.CasoTension);
if opciones.Exportar && ~isfolder(folder), mkdir(folder); end
for k=1:6
    args=struct('CasoTension',opciones.CasoTension,'PNominal_W',opciones.PNominal_W, ...
        'PasoMaximo_s',opciones.PasoMaximo_s);
    if k==1, args.Modo="ga";
    elseif k==2, args.Modo="convencional";
    else, args.Modo="fijo"; args.BitsFijos=fixed(k-2,:);
    end
    ac=namedargs2cell(args); cfg=configurar_tsc_lazo_cerrado(ac{:});
    runs{k}=simular_tsc_lazo_cerrado(cfg);
    audits{k}=verificar_tsc_lazo_cerrado(runs{k});
    a=runs{k}.trazas; t=a.t; average=@(x)trapz(t,x)/t(end);
    rows=[rows;struct('Caso',labels(k),'I_fuente_RMS_A',sqrt(average(a.iFuente.^2)), ...
        'I_sin_banco_RMS_A',sqrt(average(a.iCarga.^2)), ...
        'P_carga_media_W',average(a.u.*a.iCarga),'P_fuente_media_W',average(a.u.*a.iFuente), ...
        'P_perdidas_media_W',a.integralPerdidas_J(end)/t(end), ...
        'Cambio_energia_J',a.energia_J(end)-a.energia_J(1), ...
        'Pico_I_rama_A',max(abs(a.iRamas),[],'all'),'Pico_VC_V',max(abs(a.vCs),[],'all'), ...
        'Cambios_solicitados',nnz(diff(a.estadoSolicitado)), ...
        'Maniobras_aplicadas',runs{k}.resumen.Maniobras, ...
        'Rechazos_ON',runs{k}.resumen.Rechazos, ...
        'Tiempo_pendiente_s',trapz(t,double(any(a.pendiente~=0,2))))]; %#ok<AGROW>
    fprintf('%s, %s: I RMS %.8f A; %d maniobras\n',opciones.CasoTension,labels(k), ...
        rows(end).I_fuente_RMS_A,rows(end).Maniobras_aplicadas);
    if opciones.Exportar
        writetable(runs{k}.eventos,fullfile(folder,sprintf('eventos_%d.csv',k)));
        writetable(runs{k}.decisiones,fullfile(folder,sprintf('decisiones_%d.csv',k)));
    end
end
summary=struct2table(rows); [~,best]=min(summary.I_fuente_RMS_A(3:6));
r=struct('registros',{runs},'resumen',summary,'auditorias',{audits}, ...
    'mejor_fijo_indice',best+2,'opciones',opciones, ...
    'alcance',"Control y planta calculados dentro de Simulink; carga sintetica, fuente impuesta y valvula idealizada.");
reference=fullfile(root,'resultados','control_tsc');
if opciones.CasoTension=="moderada", reference=fullfile(reference,'tension_moderada'); end
file=fullfile(reference,'control_tsc.mat');
if isfile(file) && opciones.PNominal_W==1000
    old=load(file,'saved'); old=old.saved;
    rmsOld=old.resumen.I_fuente_RMS_A;
    r.comparacion_matlab=table(labels',rmsOld,summary.I_fuente_RMS_A, ...
        100*(summary.I_fuente_RMS_A./rmsOld-1), ...
        'VariableNames',{'Caso','RMS_MATLAB_A','RMS_Simulink_A','Diferencia_porcentaje'});
end
r.figuras=graficas_control_tsc(r,Visible=opciones.Visible);
if opciones.Exportar
    writetable(summary,fullfile(folder,'comparacion.csv'));
    if isfield(r,'comparacion_matlab'), writetable(r.comparacion_matlab,fullfile(folder,'comparacion_matlab.csv')); end
    names=["01_corriente_comparada","02_orden_y_aplicacion","03_maniobras_y_estados","04_energia_y_coste"];
    for k=1:numel(r.figuras)
        exportgraphics(r.figuras(k),fullfile(folder,names(k)+".png"),'Resolution',160);
        savefig(r.figuras(k),fullfile(folder,names(k)+".fig"));
    end
    saved=rmfield(r,'figuras'); save(fullfile(folder,'lazo_cerrado.mat'),'saved','-v7.3');
    fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(struct('casos',labels,'auditorias',{audits}),PrettyPrint=true));
end
disp(summary);
end
