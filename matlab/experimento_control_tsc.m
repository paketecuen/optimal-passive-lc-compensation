function r = experimento_control_tsc(opciones)
%EXPERIMENTO_CONTROL_TSC Comparacion dinamica con el mismo banco y ejecutor.
% 230 V RMS; la potencia es una hipotesis configurable de escala nominal.
% El mejor fijo se elige a posteriori y no representa una politica causal.
arguments
    opciones.PNominal_W (1,1) double {mustBePositive,mustBeFinite} = 1000
    opciones.PasoSalida_s (1,1) double {mustBePositive} = 1e-4
    opciones.CasoTension (1,1) string {mustBeMember(opciones.CasoTension,["original","moderada"])} = "original"
    opciones.Visible (1,1) string {mustBeMember(opciones.Visible,["on","off"])} = "on"
    opciones.Exportar (1,1) logical = true
end
labels=["Supervisor GA","Reactiva fundamental","Fijo 00", ...
    "Fijo media","Fijo unidad","Fijo ambas"];
fixed=[0,0;1,0;0,1;1,1]; runs=cell(6,1); rows=struct([]);
coef=[1,1,1]; if opciones.CasoTension=="moderada", coef=[1,.03,.03]; end
for k=1:6
    args=struct('PNominal_W',opciones.PNominal_W,'PasoSalida_s',opciones.PasoSalida_s, ...
        'CoefTension',coef);
    if k==1, args.Modo="ga";
    elseif k==2, args.Modo="convencional";
    else, args.Modo="fijo"; args.BitsFijos=fixed(k-2,:);
    end
    argumentsCell=namedargs2cell(args);
    runs{k}=simular_control_tsc(argumentsCell{:});
    fprintf('%s: %s, I RMS %.6f A\n',opciones.CasoTension,labels(k),runs{k}.resumen.I_fuente_RMS_A);
    a=runs{k}.trazas; t=a.t; duration=t(end)-t(1);
    average=@(x)trapz(t,x)/duration;
    energyDelta=a.energia_J(end)-a.energia_J(1);
    rows=[rows;struct('Caso',labels(k),'I_fuente_RMS_A',sqrt(average(a.iFuente.^2)), ...
        'I_sin_banco_RMS_A',sqrt(average(a.iCarga.^2)), ...
        'P_carga_media_W',average(a.u.*a.iCarga), ...
        'P_fuente_media_W',average(a.u.*a.iFuente), ...
        'P_perdidas_media_W',average(a.perdidas_W),'Cambio_energia_J',energyDelta, ...
        'Pico_I_rama_A',max(abs(a.iRamas),[],'all'), ...
        'Pico_VC_V',max(abs(a.vCs),[],'all'), ...
        'Cambios_solicitados',nnz(startsWith(string(runs{k}.eventos.Tipo),"solicitud_")), ...
        'Maniobras_aplicadas',nnz(ismember(string(runs{k}.eventos.Tipo),["ON","OFF"])), ...
        'Rechazos_ON',nnz(string(runs{k}.eventos.Tipo)=="rechazo_ON"), ...
        'Tiempo_pendiente_s',trapz(t,double(any(a.pendiente~=0,2))))]; %#ok<AGROW>
end
summary=struct2table(rows);
[~,j]=min(summary.I_fuente_RMS_A(3:6)); bestFixed=j+2;
r=struct('registros',{runs},'resumen',summary,'mejor_fijo_indice',bestFixed, ...
    'opciones',opciones,'alcance', ...
    "Supervisor y ejecutor hibrido integrados; fuente fijada y carga de corriente sintetica. Valvula idealizada, no hardware.");
r.figuras=graficas_control_tsc(r,Visible=opciones.Visible);
if opciones.Exportar
    folder=fullfile(fileparts(mfilename('fullpath')),'resultados','control_tsc');
    if opciones.CasoTension=="moderada", folder=fullfile(folder,'tension_moderada'); end
    if ~isfolder(folder), mkdir(folder); end
    writetable(summary,fullfile(folder,'comparacion.csv'));
    for k=1:6
        writetable(runs{k}.eventos,fullfile(folder,sprintf('eventos_%d.csv',k)));
        writetable(runs{k}.decisiones,fullfile(folder,sprintf('decisiones_%d.csv',k)));
    end
    names=["01_corriente_comparada","02_orden_y_aplicacion", ...
        "03_maniobras_y_estados","04_energia_y_coste"];
    for k=1:numel(r.figuras)
        exportgraphics(r.figuras(k),fullfile(folder,names(k)+".png"),'Resolution',150);
        savefig(r.figuras(k),fullfile(folder,names(k)+".fig"));
    end
    saved=rmfield(r,'figuras'); save(fullfile(folder,'control_tsc.mat'),'saved');
end
disp(summary); fprintf('Mejor fijo de esta trayectoria: %s\n',labels(bestFixed));
end
