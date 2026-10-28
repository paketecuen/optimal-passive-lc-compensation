function [d,windows] = decodificar_supervisor_red_tsc(out,cfg)
%DECODIFICAR_SUPERVISOR_RED_TSC Recupera costes y predicciones realmente usadas.
D=rows(out.supervisor_diag,66); W=rows(out.supervisor_ventanas,2583);
eligible=find(D(:,1)>0); [~,first]=unique(D(eligible,1),'stable'); D=D(eligible(first),:);
n=size(D,1); windows=cell(n,1);
d=table(D(:,2),D(:,3),D(:,4),D(:,6),D(:,7),D(:,8),logical(D(:,9)),D(:,10), ...
    logical(D(:,12)),D(:,11),D(:,13),D(:,18),D(:,19),D(:,20),D(:,21),D(:,31), ...
    'VariableNames',{'Tiempo_s','InicioVentana_s','UltimoPuntoInterpolado_s', ...
    'EstadoActual','EstadoCandidato','EstadoSolicitado','Pendiente','Permanencia_s', ...
    'Cambio','Ganancia_A2','RechazosPlan','Coste00','Coste01','Coste10','Coste11','EstadoSeleccionado'});
d.DiDtPredicha1_A_s=D(:,24); d.DiDtPredicha2_A_s=D(:,25);
d.ProgramacionON1_s=D(:,28); d.ProgramacionON2_s=D(:,29);
d.ErrorReconstruccion_A=NaN(n,1); d.MaxMuestraUsada_s=D(:,5);
d.Factibles=D(:,32:35); d.Alcanzables=D(:,36:39);
d.CotasVC_V=D(:,40:43); d.CotasI_A=D(:,44:47);
d.InstantesCandidatos_s=D(:,48:51); d.EnergiasCandidatas_J=D(:,52:55);
d.EstadoMedido=D(:,56:61); d.CostesBornes=D(:,62:65); d.CostesRed=D(:,14:17);
for k=1:n
    [distance,idx]=min(abs(out.supervisor_ventanas.Time-D(k,2)));
    assert(distance<1e-10,'TSCRed:Registro','Falta una ventana de decision.');
    m=reshape(W(idx,:),63,41); refs=cell(4,1);
    for j=1:4, refs{j}=struct('u',m(:,9+j),'is',m(:,13+j),'x',m(:,18+6*(j-1):23+6*(j-1))); end
    windows{k}=struct('t',D(k,3)+(0:62)'*cfg.TsSupervision, ...
        'medidas',m(:,1:8),'u',m(:,1),'iFuente',m(:,2),'iCp',m(:,3), ...
        'iRamas',m(:,4:5),'e',m(:,8),'iCargaReconstruida',m(:,9),'referencias',{refs});
end
end
function v=rows(ts,width)
data=ts.Data; count=numel(ts.Time);
if size(data,1)==count && size(data,2)==width, v=reshape(data,count,width);
else, v=reshape(data,width,count)'; end
end
