function r = experimento_banco_discreto_m3(opciones)
%EXPERIMENTO_BANCO_DISCRETO_M3 Objetivos discretos de un supervisor estacionario.
% Ramas fisicas fijas; NO simula maniobras ni autoriza disparos de tiristores.
% Las C puras son una referencia ideal, no modulos TSC industriales terminados.
arguments
    opciones.Visible (1,1) string = "on"
    opciones.Exportar (1,1) logical = true
end
root=fileparts(mfilename('fullpath')); folder=fullfile(root,'resultados','banco_discreto_m3');
if opciones.Exportar && ~isfolder(folder), mkdir(folder); end
N=63; D=derivador_periodico(N); tau=(0:N-1)'*2*pi/N;
omega0=2*pi*50; Cp=1/(26*omega0); Ls=169/(175*omega0); Cs=175/(351*omega0);
% Dividing branch admittance by two doubles L and halves C, preserving pole.
branches=table(["C media";"C unidad";"LC media";"LC unidad"], ...
    ["C";"C";"LC";"LC"],[.5;1;.5;1], ...
    [0;0;2*Ls;Ls],[.5*Cp;Cp;.5*Cs;Cs], ...
    'VariableNames',{'Rama','Tipo','Escala','L_H','C_F'});
% Illustrative fixed detuned LC modules: p=X_L/X_C at 50 Hz. Unit branch
% susceptance is 1 S at 50 Hz, without selecting commercial components.
p=[.07;.07;.14;.14]; weight=[.5;1;.5;1];
detuned=table(["LC 7% media";"LC 7% unidad";"LC 14% media";"LC 14% unidad"], ...
    repmat("LC",4,1),weight,p./(omega0*(1-p).*weight),(1-p).*weight/omega0, ...
    'VariableNames',branches.Properties.VariableNames);
branches=[branches;detuned]; branches.p_50Hz=omega0^2*branches.L_H.*branches.C_F;
designs={4,[2,4],[2,3,4],1:4,5:8}; nd=numel(designs);
labels=["Academico 1 rama","Academico 2 ramas","Academico 3 ramas", ...
    "Academico 4 ramas","LC limitadas a p7%/14%"];
K=301; dt=.2; t=(0:K-1)'*dt;
loadScale=.9+.5*sin(2*pi*t/60-pi/2);
freq=50+sin(2*pi*t/45); shape=1+.1*sin(2*pi*t/36);
noiseFraction=.01; rng(7411,'twister'); noise=randn(N,K);
threshold=.0025; dwell=10; records=cell(nd,1); summary=struct([]);
for d=1:nd
    ids=designs{d}; B=numel(ids); bits=dec2bin(0:2^B-1,B)-'0'; S=size(bits,1);
    costs=zeros(K,S); exactCosts=costs; lower=zeros(K,1); original=lower;
    pcomp=zeros(K,S); PhiAll=zeros(N,B,K); measured=zeros(N,K); actual=measured;
    qstar=zeros(N,K); voltage=measured; branchRms=zeros(K,B); coneLower=lower;
    for k=1:K
        u=sqrt(2)*(4*cos(tau).^3-2*cos(tau)+shape(k)*(2*cos(tau).^2-1));
        u=u-mean(u); a=freq(k)/50; Dt=omega0*a*D; du=Dt*u;
        aug=[Dt,ones(N,1);ones(1,N),0]; temp=aug\[u;0]; ju=temp(1:N);
        loadCurrent=loadScale(k)*(.8*u+11/(20*omega0)*du+17*omega0/10*ju ...
            +3*omega0/4*((Dt*Dt+6*omega0^2*eye(N))\du));
        Phi=zeros(N,B);
        for b=1:B
            comp=branches(ids(b),:);
            if comp.Tipo=="C"
                Phi(:,b)=comp.C_F*du;
            else
                Phi(:,b)=(Dt*Dt+1/(comp.L_H*comp.C_F)*eye(N))\du/comp.L_H;
            end
        end
        m=modelo_temporal_m3(u,D); qs=-m.Q*(m.Q'*loadCurrent/N);
        im=loadCurrent+noiseFraction*sqrt(mean(loadCurrent.^2))*noise(:,k);
        % RMS Gram matrix is the coordinate expression of inner products of
        % temporal vectors. Evaluating its vertices enumerates the dictionary.
        G=Phi'*Phi/N; h=Phi'*im/N;
        costs(k,:)=mean(im.^2)+2*(bits*h)'+sum((bits*G).*bits,2)';
        exactCosts(k,:)=mean((loadCurrent+Phi*bits').^2,1);
        pcomp(k,:)=u'*(Phi*bits')/N;
        lower(k)=mean((loadCurrent+qs).^2); original(k)=mean(loadCurrent.^2);
        cone=cono_temporal(Phi,loadCurrent); coneLower(k)=cone.I2;
        PhiAll(:,:,k)=Phi; actual(:,k)=loadCurrent; measured(:,k)=im;
        qstar(:,k)=qs; voltage(:,k)=u; branchRms(k,:)=sqrt(mean(Phi.^2,1));
    end
    [~,raw]=min(costs,[],2); [~,oracle]=min(exactCosts,[],2);
    timer=tic; [selected,history]=supervisor(costs,bits,threshold,dwell); elapsedHost=toc(timer);
    take=@(indices) exactCosts(sub2ind([K,S],(1:K)',indices));
    iRaw=sqrt(take(raw)); iSelected=sqrt(take(selected)); iOracle=sqrt(take(oracle));
    swaps=@(indices) nnz(any(diff(bits(indices,:),1,1),2));
    requests=@(indices) sum(abs(diff(bits(indices,:),1,1)),'all');
    T=table(t,freq,loadScale,shape,sqrt(original),sqrt(lower),sqrt(coneLower),iOracle,iRaw,iSelected, ...
        raw-1,selected-1,oracle-1,'VariableNames',{'Tiempo_hipotetico_s','f_Hz','Carga','Forma', ...
        'I_sin_banco_A','Cota_proyeccion_A','Cota_cono_topologia_A','I_oraculo_diccionario_A','I_directo_A', ...
        'I_supervisor_A','Estado_directo','Estado_supervisor','Estado_oraculo'});
    summary=[summary;struct('Diseno',labels(d),'Ramas_instaladas',B,'Componentes_LC',sum(1+(branches.Tipo(ids)=="LC")), ...
        'Estados',S,'Media_I_oraculo_A',mean(iOracle),'Media_I_supervisor_A',mean(iSelected), ...
        'Peor_I_supervisor_A',max(iSelected),'Exceso_medio_I2_sobre_proyeccion_A2',mean(iSelected.^2-lower), ...
        'Exceso_medio_I2_sobre_diccionario_A2',mean(iSelected.^2-iOracle.^2), ...
        'Brecha_cono_proyeccion_A2',mean(coneLower-lower),'Ramas_activas_max',max(sum(bits(selected,:),2)), ...
        'Tiempo_host_supervisor_us_por_decision',elapsedHost/K*1e6, ...
        'Solicitudes_estado_directo',swaps(raw),'Solicitudes_estado_supervisor',swaps(selected), ...
        'Cambios_bits_directo',requests(raw),'Cambios_bits_supervisor',requests(selected))]; %#ok<AGROW>
    records{d}=struct('ramas',ids,'bits',bits,'costes_medidos',costs,'costes_reales',exactCosts, ...
        'seleccion',selected,'directo',raw,'oraculo',oracle,'historia',history,'tabla',T, ...
        'Phi',PhiAll,'i',actual,'i_medido',measured,'q_ideal',qstar,'u',voltage, ...
        'potencia_compensador',pcomp,'I_ramas_conectadas_A',branchRms.*bits(selected,:));
end
% Adversarial cost fixture: raw best state alternates, while the improvement
% remains below the prescribed hysteresis. Not a thyristor simulation.
fixture=[1+.001*(-1).^(1:120)',ones(120,1)];
[steady,fixtureHistory]=supervisor(fixture,[0;1],threshold,dwell);
[~,fixtureRaw]=min(fixture,[],2);
report=struct2table(summary); disp(branches); disp(report);
r=struct('ramas',branches,'resumen',report,'registros',{records},'dt_s',dt, ...
    'histeresis_A2',threshold,'permanencia_muestras',dwell,'ruido_i_relativo',noiseFraction, ...
    'semilla',7411,'anti_chatter',struct('costes',fixture,'seleccion',steady, ...
    'directo',fixtureRaw,'historia',fixtureHistory), ...
    'alcance',"Objetivos estacionarios para un ejecutor TSC separado; sin maniobras fisicas ni transitorios");
f1=figure('Visible',opciones.Visible,'Color','w','Name','Banco discreto: corriente y estados');
tiledlayout(3,1); T=records{4}.tabla;
nexttile; plot(t,[T.I_sin_banco_A,T.Cota_proyeccion_A,T.I_oraculo_diccionario_A,T.I_supervisor_A],'LineWidth',1.2);
ylabel('I fuente RMS (A)'); grid on; legend('Sin banco','Proyeccion ideal','Mejor de 16','Supervisor','Location','best');
title('Secuencia de estados estacionarios: las transiciones no estan simuladas');
nexttile; stairs(t,records{4}.bits(records{4}.seleccion,:)+(0:3),'LineWidth',1.2);
ylabel('Orden binaria + desplazamiento'); yticks(.5:1:3.5); yticklabels(branches.Rama(1:4)); grid on;
nexttile; yyaxis left; plot(t,loadScale); ylabel('Escala carga');
yyaxis right; plot(t,freq); ylabel('f (Hz)'); xlabel('Tiempo supuesto entre evaluaciones (s)'); grid on;
f2=figure('Visible',opciones.Visible,'Color','w','Name','Presupuesto de ramas'); tiledlayout(1,3);
nexttile; bar([report.Media_I_oraculo_A,report.Media_I_supervisor_A]);
xticks(1:nd); xticklabels(["1","2","3","4","4 p7/14"]); xlabel('Ramas instaladas');
ylabel('Media de I RMS (A)'); legend('Mejor estado real','Supervisor'); grid on;
nexttile; bar([report.Cambios_bits_directo,report.Cambios_bits_supervisor]);
xticks(1:nd); xticklabels(["1","2","3","4","4 p7/14"]); xlabel('Ramas instaladas');
ylabel('Cambios binarios solicitados'); legend('Sin retencion','Con retencion'); grid on;
nexttile; plot(t,records{4}.I_ramas_conectadas_A,'LineWidth',1.2); ylabel('I de rama RMS (A)');
xlabel('Tiempo supuesto (s)'); legend(branches.Rama(1:4),'Location','best'); grid on;
f3=figure('Visible',opciones.Visible,'Color','w','Name','Estados realizables del diccionario');
c=caso_m3(N); m=modelo_temporal_m3(c.u,c.D); Phi=records{4}.Phi(:,:,1);
points=m.Q'*Phi*records{4}.bits'/N;
plot3(points(1,:),points(2,:),points(3,:),'o','LineWidth',1.3); hold on;
qs=m.Q'*records{4}.q_ideal(:,1)/N; plot3(qs(1),qs(2),qs(3),'rx','MarkerSize',12,'LineWidth',2);
for b=1:4
    for s=1:16
        if records{4}.bits(s,b)==0
            target=records{4}.bits(s,:); target(b)=1;
            j=find(all(records{4}.bits==target,2));
            plot3(points(1,[s,j]),points(2,[s,j]),points(3,[s,j]),'Color',[.7,.7,.7]);
        end
    end
end
xlabel('q en e_1 RMS'); ylabel('q en e_2 RMS'); zlabel('q en e_3 RMS'); grid on; view(35,25);
title('16 estados binarios; segmentos solo como guia geometrica');
f4=figure('Visible',opciones.Visible,'Color','w','Name','Coste de restringir la topologia');
T=records{5}.tabla;
plot(t,[T.I_sin_banco_A,T.Cota_proyeccion_A,T.Cota_cono_topologia_A,T.I_oraculo_diccionario_A,T.I_supervisor_A],'LineWidth',1.3);
legend('Sin banco','Proyeccion ideal','Cono positivo p7/14','16 estados p7/14','Supervisor p7/14','Location','best');
xlabel('Tiempo supuesto (s)'); ylabel('I fuente RMS (A)'); grid on;
title('Cuatro LC fijas p=7%/14%: limites de topologia, discretizacion y retencion');
r.figuras=[f1,f2,f3,f4];
if opciones.Exportar
    writetable(branches,fullfile(folder,'ramas.csv')); writetable(report,fullfile(folder,'comparacion.csv'));
    for d=1:nd, writetable(records{d}.tabla,fullfile(folder,sprintf('trayectoria_diseno_%d.csv',d))); end
    names=["01_corriente_estados","02_presupuesto","03_geometria","04_limite_topologia"];
    for k=1:4
        exportgraphics(r.figuras(k),fullfile(folder,names(k)+".png"),'Resolution',150);
        savefig(r.figuras(k),fullfile(folder,names(k)+".fig"));
    end
    saved=rmfield(r,'figuras'); save(fullfile(folder,'experimento_banco_discreto.mat'),'saved');
end
end

function [selected,h] = supervisor(costs,bits,threshold,dwell)
K=size(costs,1); selected=zeros(K,1); gain=zeros(K,1); changed=false(K,1);
[~,selected(1)]=min(costs(1,:)); last=1; elapsed=zeros(K,1);
for k=2:K
    old=selected(k-1); best=min(costs(k,:));
    ties=find(abs(costs(k,:)-best)<1e-12); distances=sum(abs(bits(ties,:)-bits(old,:)),2);
    [~,j]=min(distances); candidate=ties(j);
    gain(k)=costs(k,old)-costs(k,candidate); elapsed(k)=k-last;
    selected(k)=old;
    if gain(k)>threshold && elapsed(k)>=dwell
        selected(k)=candidate; changed(k)=candidate~=old; last=k;
    end
end
h=struct('ganancia_A2',gain,'edad_muestras',elapsed,'cambio',changed);
end
