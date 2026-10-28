function [decisiones,ventanas] = decodificar_registro_supervisor(out,cfg)
%DECODIFICAR_REGISTRO_SUPERVISOR Convierte salidas registradas, sin redecidir.
% out.supervisor_diag: Timeseries de los 30 diagnosticos del supervisor.
% out.supervisor_ventanas: Timeseries de 693 valores, muestreado cada T.
% No se reconstruye la trayectoria a partir del motor MATLAB anterior.
% ErrorReconstruccion_A queda NaN: el runner debe contrastarlo con la carga
% exogena registrada, que deliberadamente no es una entrada del supervisor.
dts=out.supervisor_diag; wts=out.supervisor_ventanas;
D=rows(dts,30); W=rows(wts,693); times=wts.Time(:);
eligible=find(D(:,1)>0); [~,first]=unique(D(eligible,1),'stable');
D=D(eligible(first),:); n=size(D,1);
if n==0, decisiones=table; ventanas=cell(0,1); return; end
decisiones=table(D(:,2),D(:,3),D(:,2),D(:,4), ...
    D(:,18),D(:,19),D(:,20),D(:,21),D(:,6),D(:,7),D(:,8), ...
    logical(D(:,9)),D(:,10),logical(D(:,12)),D(:,11),NaN(n,1),D(:,22),D(:,5), ...
    'VariableNames',{'Tiempo_s','InicioVentana_s','FinVentana_s','UltimoPuntoInterpolado_s', ...
    'Coste00','Coste01','Coste10','Coste11','EstadoActual','EstadoCandidato', ...
    'EstadoSolicitado','Pendiente','Permanencia_s','Cambio','Ganancia_A2', ...
    'ErrorReconstruccion_A','ErrorGram_A2','MaxMuestraUsada_s'});
decisiones.UltimaAdquisicionVentana_s=D(:,4);
decisiones.RechazosPlan=D(:,13);
decisiones.UltimaConfirmacion_s=D(:,23);
decisiones.DiDtPredicha1_A_s=D(:,24); decisiones.DiDtPredicha2_A_s=D(:,25);
decisiones.DesajustePredicho1_J=D(:,26); decisiones.DesajustePredicho2_J=D(:,27);
decisiones.ProgramacionON1_s=D(:,28); decisiones.ProgramacionON2_s=D(:,29);
decisiones.NumeroMuestras=D(:,30);
ventanas=cell(n,1);
for k=1:n
    [distance,j]=min(abs(times-D(k,2)));
    assert(distance<1e-10,'TSCLazo:RegistroVentana', ...
        'No consta la ventana de la decision t=%.12g s.',D(k,2));
    M=reshape(W(j,:),63,11); t=D(k,3)+(0:62)'*cfg.TsSupervision;
    Phi=M(:,8:9); base=M(:,6)+M(:,7);
    cs=cos(cfg.p.omega*t); sn=sin(cfg.p.omega*t);
    ac=2*mean(M(:,1).*cs); as=2*mean(M(:,1).*sn);
    quadrature=sqrt(2)*(ac*sn-as*cs)/hypot(ac,as);
    states={M(:,[8 10]),M(:,[9 11])};
    ventanas{k}=struct('t',t,'u',M(:,1),'iFuente',M(:,2), ...
        'iCp',M(:,3),'iRamas',M(:,4:5),'iCargaReconstruida',M(:,6), ...
        'iCpReferencia',M(:,7),'phiRamas',Phi, ...
        'cuadraturaFundamental',quadrature,'G',Phi'*Phi/63, ...
        'h',Phi'*base/63,'c',mean(base.^2),'costesRMS',D(k,14:17), ...
        'costesControl',D(k,18:21),'referenciaEstado',{states}, ...
        'tiempos_adquisicion_s',t,'procedencia', ...
        "Muestras directas del supervisor Simulink; periodo anterior sin interpolacion.");
end
end

function values=rows(ts,width)
data=ts.Data; count=numel(ts.Time);
if size(data,1)==count && size(data,2)==width
    values=reshape(data,count,width);
elseif size(data,1)==width && size(data,ndims(data))==count
    values=reshape(data,width,count)';
else
    error('TSCLazo:FormatoRegistro','Registro incompatible: %s, esperado %dx%d.', ...
        mat2str(size(data)),count,width);
end
end
