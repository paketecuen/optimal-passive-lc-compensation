function r = simular_tsc_lazo_cerrado(cfg,opciones)
%SIMULAR_TSC_LAZO_CERRADO Ejecuta el lazo dentro de Simulink y recoge evidencia.
% MATLAB prepara excitacion/parametros y recoge resultados; no calcula un
% cronograma previo de decisiones, disparos o extinciones.
arguments
    cfg (1,1) struct
    opciones.Exportar (1,1) logical = false
    opciones.CarpetaSalida (1,1) string = ""
end
path=crear_modelo_tsc_lazo_cerrado(cfg); [~,name]=fileparts(path);
load_system(path); cleanup=onCleanup(@()close_system(name,0));
out=sim(name); raw=out.registro;
[t,idx]=unique(raw.Time,'last'); y=datos_timeseries(raw,19); y=y(idx,:); p=cfg.p;
ir=y(:,[6 8]); vc=y(:,[7 9]);
energy=.5*p.Cp_F*y(:,5).^2+.5*sum(ir.^2.*p.L_H+vc.^2.*p.C_F,2);
loss=p.Rp_ohm*y(:,4).^2+sum(ir.^2.*p.R_ohm+vc.^2./p.Rb_ohm,2);
z=struct('t',t,'u',y(:,1),'iCarga',y(:,2),'iFuente',y(:,3), ...
    'iCp',y(:,4),'vCp',y(:,5),'iRamas',ir,'vCs',vc, ...
    'estadoAplicado',y(:,10:11),'estadoSolicitado',y(:,14:15), ...
    'disparoProgramado_s',y(:,16:17),'pendiente',y(:,18:19), ...
    'energia_J',energy,'integralEntrada_J',y(:,12), ...
    'integralPerdidas_J',y(:,13),'perdidas_W',loss, ...
    'iCargaReconstruida',y(:,3)-y(:,4)-sum(ir,2));
executors={out.ejecutor1,out.ejecutor2}; events=extraer_eventos(executors);
if isfield(cfg,'ControlRed')
    [decisions,windows]=decodificar_supervisor_red_tsc(out,cfg);
else
    [decisions,windows]=decodificar_registro_supervisor(out,cfg);
end
for k=1:numel(windows)
    actualLoad=interp1(t,z.iCarga,windows{k}.t,'linear');
    decisions.ErrorReconstruccion_A(k)=max(abs(windows{k}.iCargaReconstruida-actualLoad));
end
scale=repmat(cfg.opciones.EscalaCarga(1),size(t));
for j=2:numel(cfg.opciones.EscalaCarga)
    phase=max(0,min(1,(t-cfg.opciones.InstantesCarga_s(j))/cfg.opciones.Transicion_s));
    scale=scale+(cfg.opciones.EscalaCarga(j)-cfg.opciones.EscalaCarga(j-1))*(1-cos(pi*phase))/2;
end
z.escalaCarga=scale;
balance=energy-energy(1)-y(:,12)+y(:,13);
summary=struct('Modo',cfg.Modo,'CasoTension',cfg.CasoTension, ...
    'Duracion_s',t(end),'U_RMS_V',sqrt(trapz(t,z.u.^2)/t(end)), ...
    'I_fuente_RMS_A',sqrt(trapz(t,z.iFuente.^2)/t(end)), ...
    'I_sin_banco_RMS_A',sqrt(trapz(t,z.iCarga.^2)/t(end)), ...
    'I_rama_pico_A',max(abs(ir),[],'all'),'VC_pico_V',max(abs(vc),[],'all'), ...
    'P_perdidas_media_W',y(end,13)/t(end),'Balance_max_J',max(abs(balance)), ...
    'Maniobras',sum(events.Tipo=="ON" | events.Tipo=="OFF"), ...
    'Rechazos',sum(events.Tipo=="rechazo_ON"), ...
    'Pendientes_finales',nnz(z.pendiente(end,:)), ...
    'Error_reconstruccion_A',max(abs(z.iCargaReconstruida-z.iCarga)), ...
    'Ramas_conmutables',2,'Ramas_fijas',1,'Componentes_LC',5);
r=struct('trazas',z,'eventos',events,'decisiones',decisions, ...
    'ventanas_decision',{windows},'parametros',p,'config',cfg,'modelo',path, ...
    'registro_original',raw,'ejecutores',{executors}, ...
    'supervisor_diag',out.supervisor_diag,'resumen',summary,'alcance',cfg.alcance);
if isfield(cfg,'Fuente')
    rawSource=out.registro_fuente;
    source=datos_timeseries(rawSource,6); source=source(idx,:);
    assert(isequal(rawSource.Time,raw.Time),'TSCFuente:Registro','Tiempos de registro distintos.');
    z=r.trazas; z.e=source(:,1); z.iRed=source(:,2);
    z.energiaRed_J=.5*cfg.Fuente.L_H*z.iRed.^2;
    z.integralFuente_J=source(:,4); z.integralCarga_J=source(:,5);
    z.integralPerdidasRed_J=source(:,6);
    z.energiaTotal_J=z.energia_J+z.energiaRed_J;
    z.balanceTotal_J=z.energiaTotal_J-z.energiaTotal_J(1) ...
        -z.integralFuente_J+z.integralCarga_J+z.integralPerdidas_J+z.integralPerdidasRed_J;
    r.trazas=z; r.registro_fuente=rawSource;
    r.resumen.R_fuente_ohm=cfg.Fuente.R_ohm; r.resumen.L_fuente_mH=1e3*cfg.Fuente.L_H;
    r.resumen.E_RMS_V=sqrt(trapz(t,z.e.^2)/t(end));
    r.resumen.Caida_RMS_V=sqrt(trapz(t,(z.e-z.u).^2)/t(end));
    r.resumen.U_pico_V=max(abs(z.u));
    r.resumen.P_fuente_W=z.integralFuente_J(end)/t(end);
    r.resumen.P_carga_W=z.integralCarga_J(end)/t(end);
    r.resumen.P_perdidas_red_W=z.integralPerdidasRed_J(end)/t(end);
    r.resumen.Balance_total_max_J=max(abs(z.balanceTotal_J));
    summary=r.resumen;
end
if isfield(cfg,'ControlRed')
    r.contextos_admision=out.contextos_admision;
    r.guardias={out.guardia1,out.guardia2};
    r.resumen.Rechazos_plan=sum(decisions.RechazosPlan);
    r.resumen.MaxVC_config_V=cfg.ControlRed.MaxVC_V;
    r.resumen.MaxI_config_A=cfg.ControlRed.MaxI_A;
    summary=r.resumen;
end
if opciones.Exportar
    folder=opciones.CarpetaSalida;
    if strlength(folder)==0
        folder=fullfile(fileparts(mfilename('fullpath')),'resultados','tsc_lazo_cerrado',cfg.CasoTension,cfg.Modo);
    end
    if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,'simulacion.mat'),'r');
    writetable(events,fullfile(folder,'eventos.csv'));
    writetable(decisions,fullfile(folder,'decisiones.csv'));
    writetable(struct2table(summary),fullfile(folder,'resumen.csv'));
end
end

function ev=extraer_eventos(executors)
rows=struct([]);
for j=1:2
    ts=executors{j}; data=datos_timeseries(ts,10);
    counts=data(:,4:6); jumps=[zeros(1,3);diff(counts)];
    indices=find(any(jumps>0,2));
    for k=indices'
        kind="ON";
        if jumps(k,2)>0, kind="OFF"; end
        if jumps(k,3)>0, kind="rechazo_ON"; end
        rows=[rows;struct('Tipo',kind,'Rama',j,'Tiempo_s',ts.Time(k), ...
            'I_pre_A',data(k,7),'VC_pre_V',data(k,8), ...
            'DiDt_A_s',data(k,9),'Admisible',kind~="rechazo_ON")]; %#ok<AGROW>
    end
end
if isempty(rows)
    ev=table('Size',[0 7],'VariableTypes',{'string','double','double','double','double','double','logical'}, ...
        'VariableNames',{'Tipo','Rama','Tiempo_s','I_pre_A','VC_pre_V','DiDt_A_s','Admisible'});
else
    ev=sortrows(struct2table(rows),{'Tiempo_s','Rama'});
end
end

function values=datos_timeseries(ts,width)
data=ts.Data; count=numel(ts.Time);
if size(data,1)==count && size(data,2)==width
    values=reshape(data,count,width);
elseif size(data,1)==width && size(data,ndims(data))==count
    values=reshape(data,width,count)';
else
    error('TSCLazo:FormatoRegistro','Registro incompatible: %s.',mat2str(size(data)));
end
end
