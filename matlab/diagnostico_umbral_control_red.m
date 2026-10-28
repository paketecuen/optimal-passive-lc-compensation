function report = diagnostico_umbral_control_red
%DIAGNOSTICO_UMBRAL_CONTROL_RED Busca donde congelar u cambia la decision.
% Barrido declarado de carga constante; no selecciona casos por RMS dinamico.
% Primero se calcula el diccionario periodico, despues se contrastan dos
% decisiones y trayectorias en Simulink desde el mismo estado descargado.
cfg=configurar_control_red_tsc(configurar_fuente_tsc(configurar_tsc_lazo_cerrado, ...
    R_ohm=1.5,L_H=3e-3),MaxVC_V=[Inf Inf],MaxI_A=[Inf Inf]);
c=cfg.ControlRed; il=cfg.p.iCargaReferencia; e=cfg.p.uReferencia;
alpha=(.05:.002:2.5)'; m=numel(alpha); red=zeros(m,4); frozen=red;
base=c.topologias{1}; u0=base.Ue*e; u1=base.Ul*il;
for k=1:4
    op=c.topologias{k}; a=op.Y*(op.Ue*e); b=il+op.Y*(op.Ul*il);
    red(:,k)=mean((a+b.*alpha').^2,1)';
    a=op.Y*u0; b=il+op.Y*u1; frozen(:,k)=mean((a+b.*alpha').^2,1)';
end
[~,bestRed]=min(red(:,1:3),[],2); [~,bestFrozen]=min(frozen(:,1:3),[],2);
for k=1:m
    if red(k,1)-red(k,bestRed(k))<=cfg.p.histeresis_A2, bestRed(k)=1; end
    if frozen(k,1)-frozen(k,bestFrozen(k))<=cfg.p.histeresis_A2, bestFrozen(k)=1; end
end
penalty=red(sub2ind(size(red),(1:m)',bestFrozen))-red(sub2ind(size(red),(1:m)',bestRed));
different=bestRed~=bestFrozen; candidates=find(different);
tableScan=table(alpha,bestRed-1,bestFrozen-1,penalty, ...
    'VariableNames',{'EscalaCarga','Estado_red','Estado_bornes','Penalizacion_coste_real_A2'});
folder=fullfile(fileparts(mfilename('fullpath')),'resultados','control_red','pruebas_auxiliares');
if ~isfolder(folder), mkdir(folder); end
writetable(tableScan,fullfile(folder,'barrido_umbral.csv'));
report=struct('Puntos',m,'Decisiones_distintas',nnz(different),'R_ohm',1.5,'L_H',3e-3);
if isempty(candidates), disp(report); return; end
[~,ii]=max(penalty(candidates)); idx=candidates(ii); chosen=alpha(idx);
runs=cell(2,1); audits=cell(2,1); models=["bornes","red"];
for k=1:2
    cfg=configurar_control_red_tsc(configurar_fuente_tsc(configurar_tsc_lazo_cerrado( ...
        EscalaCarga=chosen,InstantesCarga_s=0,Duracion_s=.12),R_ohm=1.5,L_H=3e-3), ...
        MaxVC_V=[Inf Inf],MaxI_A=[Inf Inf],ModeloCoste=models(k));
    runs{k}=simular_tsc_lazo_cerrado(cfg); audits{k}=verificar_control_red_tsc(runs{k});
end
assert(runs{1}.decisiones.EstadoSolicitado(1)==bestFrozen(idx)-1 && ...
    runs{2}.decisiones.EstadoSolicitado(1)==bestRed(idx)-1, ...
    'TSCRed:Umbral','La simulacion no reproduce la decision prevista en el barrido.');
report.Escala_elegida=chosen; report.Estado_bornes=bestFrozen(idx)-1;
report.Estado_red=bestRed(idx)-1; report.Diferencia_coste_periodico_A2=penalty(idx);
report.RMS_periodico_bornes_A=sqrt(red(idx,bestFrozen(idx)));
report.RMS_periodico_red_A=sqrt(red(idx,bestRed(idx)));
report.RMS_transitorio_bornes_A=runs{1}.resumen.I_fuente_RMS_A;
report.RMS_transitorio_red_A=runs{2}.resumen.I_fuente_RMS_A;
report.Alcance="Caso encontrado por barrido de amplitud para mostrar diferencia de decision estacionaria. No es evidencia de superioridad dinamica general.";
save(fullfile(folder,'umbral.mat'),'runs','audits','report','-v7.3');
fid=fopen(fullfile(folder,'umbral.json'),'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true)); disp(report);
end
