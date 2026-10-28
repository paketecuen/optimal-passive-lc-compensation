function figs = graficas_control_tsc(r,opciones)
%GRAFICAS_CONTROL_TSC Corrientes reales, solicitudes, energia y maniobras.
arguments
    r (1,1) struct
    opciones.Visible (1,1) string = "on"
end
ga=r.registros{1}; a=ga.trazas; p=ga.parametros;
ids=[1,2,r.mejor_fijo_indice]; labels=r.resumen.Caso(ids);
figs=gobjects(4,1);
figs(1)=figure('Visible',opciones.Visible,'Color','w','Position',[60,60,1200,760]);
tiledlayout(2,1,'TileSpacing','compact'); nexttile; hold on
for k=1:numel(ids)
    b=r.registros{ids(k)}.trazas;
    [times,values]=rms_periodos(b.t,b.iFuente,p.T);
    plot(times,values,'LineWidth',1.4,'DisplayName',labels(k));
end
[times,values]=rms_periodos(a.t,a.iCarga,p.T);
plot(times,values,'k--','DisplayName','Carga sin banco');
ylabel('I fuente RMS por periodo (A)'); legend('Location','best'); grid on
title('Corriente dinamica observada; cada punto resume el periodo anterior');
nexttile; plot(a.t,a.escalaCarga,'LineWidth',1.4); grid on
xlabel('Tiempo (s)'); ylabel('Multiplicador de corriente de carga');
title('Carga sintetica suave, sin anticipacion de sus cambios por el supervisor');
figs(2)=figure('Visible',opciones.Visible,'Color','w','Position',[90,90,1200,800]);
tiledlayout(3,1,'TileSpacing','compact');
for k=1:2
    nexttile;
    stairs(a.t,double(a.estadoSolicitado(:,k)),'--','LineWidth',1.5); hold on
    stairs(a.t,double(a.estadoAplicado(:,k)),'LineWidth',1.1);
    area(a.t,.2*double(a.pendiente(:,k)~=0),'FaceAlpha',.2,'EdgeColor','none');
    ylabel(sprintf('Rama LC %d',k)); ylim([-.1,1.25]); grid on
    legend('Solicitado','Aplicado','Pendiente (banda)','Location','best');
end
age=ga.decisiones.Permanencia_s; age(~isfinite(age))=NaN;
% Antes de la primera confirmacion, algunos ejecutores usan un centinela
% finito para habilitar la primera orden. No representa una edad fisica.
age(age>ga.trazas.t(end))=NaN;
nexttile; stairs(ga.decisiones.Tiempo_s,age,'LineWidth',1.3);
xlabel('Tiempo (s)'); ylabel('Edad de estado aplicado (s)'); grid on
title('La permanencia se cuenta desde la ultima maniobra confirmada');
figs(3)=figure('Visible',opciones.Visible,'Color','w','Position',[120,120,1200,800]);
tiledlayout(2,2,'TileSpacing','compact');
[~,peakrow]=max(max(abs(a.vCs),[],2)); center=a.t(peakrow);
window=[max(0,center-.03),min(a.t(end),center+.03)];
nexttile; plot(a.t,a.iRamas,'LineWidth',1.2); xlim(window); grid on
xlabel('Tiempo (s)'); ylabel('Corriente LC (A)'); legend('Media','Unidad');
title('Corrientes en la ventana del mayor pico de tension C');
nexttile; plot(a.t,[a.u,a.vCs],'LineWidth',1.1); xlim(window); grid on
xlabel('Tiempo (s)'); ylabel('Tension (V)'); legend('Bornes','C media','C unidad');
title('La tension de C se conserva al cambiar de estado');
nexttile; hold on
for k=1:2
    plot(sqrt(p.L_H(k))*a.iRamas(:,k),sqrt(p.C_F(k))*a.vCs(:,k),'LineWidth',1.1);
end
xlabel('sqrt(L) i'); ylabel('sqrt(C) v_C'); grid on; axis equal
legend('Media','Unidad'); title('Trayectorias en el marco {P_+,N_-}');
nexttile; plot(a.t,a.energia_J,'LineWidth',1.2); grid on
xlabel('Tiempo (s)'); ylabel('Energia total almacenada (J)');
title('Incluye el condensador fijo y ambas LC');
figs(4)=figure('Visible',opciones.Visible,'Color','w','Position',[150,150,1250,800]);
tiledlayout(2,2,'TileSpacing','compact'); nexttile;
plot(a.t,a.energia_J-a.energia_J(1),'LineWidth',1.3); hold on
plot(a.t,a.integralEntrada_J-a.integralPerdidas_J,'--','LineWidth',1.2);
xlabel('Tiempo (s)'); ylabel('Energia (J)'); grid on
legend('E(t)-E(0)','Integral entrada - perdidas','Location','best');
title('Balance del banco durante maniobras');
nexttile; plot(a.t,a.perdidas_W); grid on
xlabel('Tiempo (s)'); ylabel('Disipacion instantanea (W)');
nexttile; bar(r.resumen.I_fuente_RMS_A); grid on
xticks(1:height(r.resumen)); xticklabels(r.resumen.Caso); xtickangle(20);
ylabel('I RMS global (A)'); title('Misma carga y condiciones iniciales');
nexttile; bar([r.resumen.Cambios_solicitados,r.resumen.Maniobras_aplicadas]); grid on
xticks(1:height(r.resumen)); xticklabels(r.resumen.Caso); xtickangle(20);
ylabel('Cambios de bit'); legend('Solicitados','Aplicados','Location','best');
title('Mejor fijo seleccionado a posteriori entre cuatro');
for k=1:numel(figs)
    sgtitle(figs(k),sprintf('230 V RMS, escala nominal %.0f W; tension %s', ...
        r.opciones.PNominal_W,r.opciones.CasoTension));
end
end

function [times,values]=rms_periodos(t,y,T)
[t,indices]=unique(t,'last'); y=y(indices);
count=floor((t(end)+1e-10)/T); times=(1:count)'*T; values=zeros(count,1);
for k=1:count
    lo=(k-1)*T; hi=k*T; take=t>lo & t<hi;
    tx=[lo;t(take);hi]; yy=interp1(t,y,tx,'linear');
    values(k)=sqrt(trapz(tx,yy.^2)/T);
end
end
