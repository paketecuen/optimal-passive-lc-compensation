function figs = graficas_control_red_tsc(r,opciones)
%GRAFICAS_CONTROL_RED_TSC Corriente observada, cotas y tensiones predichas.
arguments
    r (1,1) struct
    opciones.Visible (1,1) string = "on"
    opciones.CarpetaSalida (1,1) string = ""
end
figs=gobjects(4,1);
for k=1:4, figs(k)=figure('Visible',opciones.Visible,'Color','w','Position',[60+20*k,70+20*k,1250,800]); end
ids=[1 3 5 6];
figure(figs(1)); tiledlayout(2,1,'TileSpacing','compact');
for net=1:2
    nexttile; hold on
    for k=ids
        a=r.registros{net,k}.trazas; [tt,yy]=rmsPeriods(a.t,a.iFuente,.02);
        plot(tt,yy,'LineWidth',1.3,'DisplayName',r.controles(k));
    end
    grid on; xlabel('Tiempo (s)'); ylabel('Is RMS por periodo (A)');
    title(r.redes(net),'Interpreter','none'); legend('Location','best');
end
sgtitle('Corriente observada al cambiar el predictor y los limites de arranque');
figure(figs(2)); tiledlayout(2,2,'TileSpacing','compact');
new=r.registros{1,5}; a=new.trazas; old=r.registros{1,1}.trazas;
nexttile; plot(a.t,a.vCs,'LineWidth',1.2); hold on
yline(new.config.ControlRed.MaxVC_V(1),'k--'); yline(-new.config.ControlRed.MaxVC_V(1),'k--');
grid on; xlabel('Tiempo (s)'); ylabel('VC (V)'); title('Red intermedia; limite 1100 V');
legend('Rama media','Rama unidad','Limite','Location','best');
nexttile; plot(old.t,old.iRamas(:,2),':','LineWidth',1.1); hold on
plot(a.t,a.iRamas,'LineWidth',1.2); yline(8,'k--'); yline(-8,'k--');
grid on; xlabel('Tiempo (s)'); ylabel('Corriente de rama (A)');
title('Limite 8 A: comparar con rama unidad del control anterior');
legend('Unidad anterior','Media limitada','Unidad limitada','Location','best');
nexttile; stairs(a.t,2*a.estadoAplicado(:,1)+a.estadoAplicado(:,2),'LineWidth',1.3); hold on
stairs(old.t,2*old.estadoAplicado(:,1)+old.estadoAplicado(:,2),'--','LineWidth',1.1);
yticks(0:3); yticklabels({'00','01','10','11'}); grid on; xlabel('Tiempo (s)'); ylabel('Estado aplicado');
legend('Limitado','Anterior','Location','best'); title('Una rama por decision, OFF al cero natural');
nexttile; stem(new.decisiones.Tiempo_s,new.decisiones.RechazosPlan,'filled'); grid on
xlabel('Tiempo (s)'); ylabel('Candidatos bloqueados'); title('Mejorarian el coste, pero su ON no es admisible');
sgtitle('Restricciones de arranque; los limites no son datos de un dispositivo comercial');
figure(figs(3)); tiledlayout(2,2,'TileSpacing','compact');
fields={'Is_RMS_A','PicoVC_V','PicoI_A','Maniobras'};
yl={'Is RMS global (A)','Pico VC (V)','Pico corriente rama (A)','Maniobras'};
for k=1:4
    nexttile; bar(reshape(r.resumen.(fields{k}),6,2)); grid on
    xticks(1:6); xticklabels(r.controles); xtickangle(20); ylabel(yl{k});
end
lg=legend(r.redes,'Interpreter','none','Orientation','horizontal'); lg.Layout.Tile='south';
sgtitle('Misma fuente, carga y arranque descargado dentro de cada red');
figure(figs(4)); tiledlayout(2,2,'TileSpacing','compact');
rr=r.registros{2,3}; w=rr.ventanas_decision{1}; labels=["00","01","10","11"];
nexttile; hold on
for j=1:4, plot(1e3*(w.t-w.t(1)),w.referencias{j}.u,'LineWidth',1.2,'DisplayName',labels(j)); end
grid on; xlabel('Tiempo dentro del periodo (ms)'); ylabel('u predicha (V)'); legend('Location','best');
title('Cada topologia cambia la tension de nodo');
nexttile; hold on
for j=1:4, plot(1e3*(w.t-w.t(1)),w.referencias{j}.is,'LineWidth',1.2,'DisplayName',labels(j)); end
grid on; xlabel('Tiempo dentro del periodo (ms)'); ylabel('is predicha (A)'); legend('Location','best');
title('Cuatro vectores temporales de corriente, con red incluida');
nexttile; hold on
for net=1:2
    tab=r.auditorias{net,5}.admisiones;
    if ~isempty(tab)
        scatter(max(tab{:,{'CotaVC1_V','CotaVC2_V'}},[],2), ...
            max(tab{:,{'CotaI1_A','CotaI2_A'}},[],2),55,'filled','DisplayName',r.redes(net));
    end
end
xline(1100,'k--','1100 V','HandleVisibility','off');
yline(8,'k--','8 A','HandleVisibility','off'); grid on
ax=gca; ax.XLim=ax.XLim+[-2 2]; ax.YLim=ax.YLim+[-.05 .05];
xlabel('Cota de tension en el evento ON (V)'); ylabel('Cota de corriente (A)');
title('Revalidacion con el estado medido al conectar'); legend('Location','best','Interpreter','none');
nexttile; bar([sqrt(rr.decisiones.CostesBornes(1,:));sqrt(rr.decisiones.CostesRed(1,:))]');
xticks(1:4); xticklabels(labels); grid on; ylabel('RMS previsto (A)');
legend('Tension congelada','Red acoplada','Location','best'); title('Mismas medidas; cambia el modelo del coste');
sgtitle('Prediccion estacionaria y admision energetica son comprobaciones distintas');
if strlength(opciones.CarpetaSalida)>0
    names=["01_corriente","02_limites_y_estados","03_comparacion","04_prediccion_y_admision"];
    for k=1:4
        savefig(figs(k),fullfile(opciones.CarpetaSalida,names(k)+".fig"));
        exportgraphics(figs(k),fullfile(opciones.CarpetaSalida,names(k)+".png"),'Resolution',150);
    end
end
end
function [tt,yy]=rmsPeriods(t,y,T)
tt=(1:floor((t(end)+1e-10)/T))'*T; yy=zeros(size(tt));
for k=1:numel(tt)
    tx=[tt(k)-T;t(t>tt(k)-T & t<tt(k));tt(k)];
    yy(k)=sqrt(trapz(tx,interp1(t,y,tx,'linear').^2)/T);
end
end
