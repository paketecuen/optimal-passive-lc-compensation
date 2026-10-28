function figs = graficas_fuente_tsc(r,opciones)
%GRAFICAS_FUENTE_TSC Efecto de la red sobre bornes, arranque y consumo.
arguments
    r (1,1) struct
    opciones.Visible (1,1) string = "on"
    opciones.CarpetaSalida (1,1) string = ""
end
labels=["R_g=0, L_g=0","R_g=0.4 ohm, L_g=0.8 mH","R_g=1.5 ohm, L_g=3 mH"];
figs=gobjects(4,1); colors=lines(3);
for k=1:4
    figs(k)=figure('Visible',opciones.Visible,'Color','w','Position',[60+25*k,60+25*k,1220,790]);
end
figure(figs(1)); tiledlayout(3,1,'TileSpacing','compact');
for panel=1:3
    nexttile; hold on
    for net=1:3
        a=r.registros{net,1}.trazas;
        if panel==1, y=a.u; elseif panel==2, y=a.iFuente; else, y=a.e-a.u; end
        [t,v]=rmsPeriodos(a.t,y,.02); plot(t,v,'LineWidth',1.4,'Color',colors(net,:));
    end
    grid on; legend(labels,'Location','best');
    if panel==1
        ylabel('U bornes RMS (V)'); title('Supervisor GA; cada punto resume el periodo anterior');
    elseif panel==2
        ylabel('I fuente RMS (A)');
    else
        ylabel('RMS(e-u) (V)'); xlabel('Tiempo (s)');
    end
end
sgtitle('La fuente mantiene e RMS=230 V; la tension u responde al banco y a la red');
figure(figs(2)); tiledlayout(2,2,'TileSpacing','compact');
for panel=1:4
    nexttile; hold on
    for net=1:3
        a=r.registros{net,1}.trazas; take=a.t<=.12;
        if panel==1, y=a.vCs(:,1); elseif panel==2, y=a.u-a.e;
        elseif panel==3, y=a.iRamas(:,1); else, y=a.energiaRed_J; end
        plot(1e3*a.t(take),y(take),'LineWidth',1.2,'Color',colors(net,:));
    end
    grid on; xlabel('Tiempo (ms)'); legend(labels,'Location','best');
    if panel==1
        ylabel('VC rama media (V)'); title('Arranque con condensadores descargados');
    elseif panel==2
        ylabel('u-e (V)'); title('Perturbacion de la tension en bornes');
    elseif panel==3
        ylabel('Corriente rama media (A)'); title('Oscilacion libre y respuesta forzada');
    else
        ylabel('Energia en L_g (J)'); title('Nuevo estado dinamico de la fuente');
    end
end
sgtitle('La impedancia de fuente no elimina por si sola la sobretension de arranque');
figure(figs(3)); tiledlayout(2,2,'TileSpacing','compact');
fields={'I_fuente_RMS_A','VC_pico_V','P_carga_W','P_perdidas_red_W'};
ys={'I fuente RMS (A)','Pico VC (V)','Potencia media de carga (W)','Perdidas medias de red (W)'};
titles={'Trayectoria completa, incluidos arranques','No son limites de dispositivo', ...
    'La misma corriente exogena no implica la misma potencia','R_g multiplica el RMS cuadrado de la corriente'};
for panel=1:4
    nexttile; values=reshape(r.resumen.(fields{panel}),6,3);
    bar(values); grid on; xticks(1:6); xticklabels(r.estrategias); xtickangle(20);
    ylabel(ys{panel}); title(titles{panel});
end
lg=legend(labels,'Orientation','horizontal'); lg.Layout.Tile='south';
sgtitle('GA, reactiva fundamental y cuatro estados fijos; mismas condiciones iniciales');
figure(figs(4)); tiledlayout(3,1,'TileSpacing','compact');
for net=1:3
    nexttile; a=r.registros{net,1}.trazas;
    yyaxis left
    plot(a.t,a.energiaTotal_J-a.energiaTotal_J(1),'LineWidth',1.3); hold on
    supply=a.integralFuente_J-a.integralCarga_J-a.integralPerdidasRed_J-a.integralPerdidas_J;
    plot(a.t,supply,'--','LineWidth',1.3); ylabel('Cambio de energia (J)');
    yyaxis right; plot(a.t,1e6*a.balanceTotal_J,':','LineWidth',1.1); ylabel('Residuo (microJ)');
    grid on; xlabel('Tiempo (s)'); title(labels(net));
    legend('E(t)-E(0)','Entrada - carga - perdidas','Residuo','Location','best');
end
sgtitle('Balance total: inductancia de fuente, banco, carga y perdidas');
if strlength(opciones.CarpetaSalida)>0
    names=["01_bornes_y_corriente","02_arranque","03_estrategias","04_balance_total"];
    for k=1:4
        savefig(figs(k),fullfile(opciones.CarpetaSalida,names(k)+".fig"));
        exportgraphics(figs(k),fullfile(opciones.CarpetaSalida,names(k)+".png"),'Resolution',150);
    end
end
end

function [times,values]=rmsPeriodos(t,y,T)
count=floor((t(end)+1e-10)/T); times=(1:count)'*T; values=zeros(count,1);
for k=1:count
    lo=(k-1)*T; hi=k*T; take=t>lo & t<hi;
    tx=[lo;t(take);hi]; yy=interp1(t,y,tx,'linear');
    values(k)=sqrt(trapz(tx,yy.^2)/T);
end
end
