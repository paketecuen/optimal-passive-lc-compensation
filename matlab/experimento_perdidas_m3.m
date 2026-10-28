function r = experimento_perdidas_m3(opciones)
%EXPERIMENTO_PERDIDAS_M3 Mismo banco fijo con R de bobina y ESR de condensadores.
% Q a 50 Hz define resistencias constantes; no modelo comercial ni termico.
% Respuesta periodica forzada con tension impuesta. Sin FFT ni reajuste LC.
arguments
    opciones.Visible (1,1) string = "on"
    opciones.Exportar (1,1) logical = true
    opciones.N (1,1) double {mustBeInteger,mustBePositive} = 31
    opciones.Barrido (1,1) logical = true
end
N=opciones.N; D=derivador_periodico(N); w0=2*pi*50;
b=struct('Cp_F',1/(26*w0),'Ls_H',169/(175*w0),'Cs_F',175/(351*w0));
QL=[Inf,200,100,50,20]; QC=[Inf,1000,500,200,100];
rows=struct([]); records=cell(5,1);
for k=1:numel(QL)
    rec=resolver(b,D,50,1,1,QL(k),QC(k)); records{k}=rec;
    rows=[rows;resumen(rec)]; %#ok<AGROW>
end
nominal=struct2table(rows); rows=struct([]);
if opciones.Barrido
    for ql=QL
        for qc=QC
            for f=49:.25:51
                for loadFactor=[.8,1,1.2]
                    for shape=[.8,1,1.2]
                        rec=resolver(b,D,f,loadFactor,shape,ql,qc);
                        rows=[rows;resumen(rec)]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
campaign=struct2table(rows);
r=struct('banco',b,'N',N,'referencia_Hz',50,'QL',QL,'QC',QC, ...
    'nominal',nominal,'registros_nominales',{records},'barrido',campaign, ...
    'alcance',"R y ESR constantes, tension impuesta y respuesta periodica; sin reajustar L/C");
figs=grafica(r,opciones.Visible); r.figuras=figs; disp(nominal);
if opciones.Exportar
    out=fullfile(fileparts(mfilename('fullpath')),'resultados','perdidas_m3');
    if ~isfolder(out), mkdir(out); end
    names=["01_nominal","02_formas_onda","03_dominio"];
    for k=1:numel(figs)
        exportgraphics(figs(k),fullfile(out,names(k)+".png"),'Resolution',150);
        savefig(figs(k),fullfile(out,names(k)+".fig"));
    end
    writetable(nominal,fullfile(out,'nominal.csv'));
    if opciones.Barrido, writetable(campaign,fullfile(out,'barrido.csv')); end
    saved=rmfield(r,'figuras'); save(fullfile(out,'experimento_perdidas.mat'),'saved');
end
end

function rec=resolver(b,D,f,loadFactor,shape,QL,QC)
N=size(D,1); tau=(0:N-1)'*2*pi/N; omega=2*pi*f; w0=2*pi*50;
Dt=omega*D; E=eye(N); one=ones(N,1);
u=sqrt(2)*(4*cos(tau).^3-2*cos(tau)+shape*(2*cos(tau).^2-1));
u=u-mean(u); du=Dt*u;
aug=[Dt,one;one',0]; temp=aug\[u;0]; ju=temp(1:N);
% Load is the same passive ideal physical circuit as caso_m3, at actual f.
Ll=1/(w0*(3/4)); Cl=(3/4)/(6*w0);
i=loadFactor*(.8*u+(11/20)/w0*du+(17/10)*w0*ju+ ...
    (Ll*Dt*Dt+E/Cl)\du);
RL=w0*b.Ls_H/QL; Rp=1/(w0*b.Cp_F*QC); Rs=1/(w0*b.Cs_F*QC);
% RC current and RLC current equations, both in physical seconds.
ip=(E+Rp*b.Cp_F*Dt)\(b.Cp_F*du);
is=(b.Ls_H*Dt*Dt+(RL+Rs)*Dt+E/b.Cs_F)\du;
vp=u-Rp*ip; vl=b.Ls_H*Dt*is; vs=u-(RL+Rs)*is-vl;
q=ip+is; source=i+q; loss=Rp*ip.^2+(RL+Rs)*is.^2;
energy=(b.Cp_F*vp.^2+b.Ls_H*is.^2+b.Cs_F*vs.^2)/2;
activeConductance=mean(u.*q)/mean(u.^2);
orthogonalSource=source-(.8*loadFactor+activeConductance)*u;
rec=struct('f',f,'carga',loadFactor,'forma',shape,'QL',QL,'QC',QC, ...
    'RL',RL,'Rp',Rp,'Rs',Rs,'banco',b,'t',tau/omega,'u',u,'i',i,'q',q, ...
    'fuente',source,'ip',ip,'is',is,'vp',vp,'vs',vs,'vl',vl, ...
    'perdidas',loss,'energia',energy,'G_perdidas',activeConductance, ...
    'fuente_ortogonal',orthogonalSource, ...
    'error_kcl',max(abs(q-ip-is)), ...
    'error_rc',max(abs(ip-b.Cp_F*Dt*vp)), ...
    'error_rlc',max(abs(is-b.Cs_F*Dt*vs)), ...
    'error_energia',max(abs(Dt*energy-u.*q+loss)));
end

function s=resumen(r)
rmsv=@(x) sqrt(mean(x.^2));
s=struct('f_Hz',r.f,'Carga',r.carga,'Forma',r.forma,'QL_ref',r.QL,'QC_ref',r.QC, ...
    'R_bobina_Ohm',r.RL,'ESR_Cp_Ohm',r.Rp,'ESR_Cs_Ohm',r.Rs, ...
    'U_RMS_V',rmsv(r.u),'I_sin_banco_A',rmsv(r.i),'I_fuente_A',rmsv(r.fuente), ...
    'P_carga_W',mean(r.u.*r.i),'P_fuente_W',mean(r.u.*r.fuente), ...
    'P_perdidas_W',mean(r.perdidas),'G_perdidas_S',r.G_perdidas, ...
    'I_ortogonal_A',rmsv(r.fuente_ortogonal), ...
    'I_Cp_A',rmsv(r.ip),'I_serie_A',rmsv(r.is), ...
    'V_Cp_V',rmsv(r.vp),'V_Cs_V',rmsv(r.vs),'V_Ls_V',rmsv(r.vl), ...
    'V_Cs_max_nodos_V',max(abs(r.vs)),'I_serie_max_nodos_A',max(abs(r.is)), ...
    'Energia_media_J',mean(r.energia), ...
    'Error_balance_W',abs(mean(r.u.*r.fuente)-mean(r.u.*r.i)-mean(r.perdidas)), ...
    'Error_energia_W',r.error_energia, ...
    'Error_RC_A',r.error_rc,'Error_RLC_A',r.error_rlc);
end

function figs=grafica(r,visible)
T=r.nominal; names=["Ideal","Q_L=200 / Q_C=1000","100 / 500","50 / 200","20 / 100"];
x=categorical(names,names,'Ordinal',true);
figs=gobjects(0); figs(1)=figure('Visible',visible,'Color','w','Position',[100,100,1200,680]);
tiledlayout(2,2,'TileSpacing','compact');
nexttile; bar(x,T.I_fuente_A); yline(T.I_sin_banco_A(1),'--','Sin banco');
ylabel('I fuente RMS (A)'); grid on; title('Mismos L y C, distintas resistencias');
nexttile; bar(x,1000*T.P_perdidas_W); ylabel('Potencia disipada (mW)'); grid on;
nexttile; plot(1:5,T.I_Cp_A,'o-',1:5,T.I_serie_A,'s-'); xticks(1:5); xticklabels(names);
ylabel('Corriente RMS (A)'); legend('Rama C_p','Rama L_s-C_s','Location','best'); grid on;
nexttile; plot(1:5,T.V_Cs_V,'o-',1:5,T.V_Ls_V,'s-'); xticks(1:5); xticklabels(names);
ylabel('Tension interna RMS (V)'); legend('C_s','L_s','Location','best'); grid on;
sgtitle('Perdidas nominales a 50 Hz: resistencias exploratorias constantes');
figs(2)=figure('Visible',visible,'Color','w','Position',[100,100,1200,720]);
tiledlayout(3,1,'TileSpacing','compact'); colors=lines(3); which=[1,3,5];
tau=linspace(0,2*pi,1001)'; waves=cell(3,1);
for k=1:3
    a=r.registros_nominales{which(k)};
    waves{k}=interpolar_periodico([a.fuente,a.q,a.perdidas],tau);
end
nexttile; hold on;
for k=1:3, plot(tau/(2*pi*50)*1000,waves{k}(:,1),'Color',colors(k,:),'LineWidth',1.5); end
ylabel('i fuente (A)'); grid on; legend(names(which),'Location','eastoutside');
nexttile; hold on;
for k=1:3, plot(tau/(2*pi*50)*1000,waves{k}(:,2),'Color',colors(k,:),'LineWidth',1.5); end
ylabel('i banco (A)'); grid on;
nexttile; hold on;
for k=1:3, plot(tau/(2*pi*50)*1000,1000*waves{k}(:,3),'Color',colors(k,:),'LineWidth',1.5); end
ylabel('p disipada (mW)'); xlabel('Tiempo (ms)'); grid on;
sgtitle('Respuesta periodica: el banco con perdidas absorbe potencia activa');
if isempty(r.barrido), return; end
figs(3)=figure('Visible',visible,'Color','w','Position',[100,100,1250,550]);
tiledlayout(1,2,'TileSpacing','compact');
A=zeros(5); B=A;
for j=1:5
    for k=1:5
        mask=r.barrido.QL_ref==r.QL(j) & r.barrido.QC_ref==r.QC(k);
        A(j,k)=max(r.barrido.I_fuente_A(mask)); B(j,k)=1000*max(r.barrido.P_perdidas_W(mask));
    end
end
nexttile; imagesc(A); colorbar; title('Mayor I fuente en los 81 estados (A)');
xticks(1:5); xticklabels(string(r.QC)); yticks(1:5); yticklabels(string(r.QL));
xlabel('Q_C de referencia'); ylabel('Q_L de referencia');
for j=1:5, for k=1:5, text(k,j,sprintf('%.3f',A(j,k)),'HorizontalAlignment','center','Color','k'); end, end
nexttile; imagesc(B); colorbar; title('Mayor potencia disipada (mW)');
xticks(1:5); xticklabels(string(r.QC)); yticks(1:5); yticklabels(string(r.QL));
xlabel('Q_C de referencia'); ylabel('Q_L de referencia');
for j=1:5, for k=1:5, text(k,j,sprintf('%.1f',B(j,k)),'HorizontalAlignment','center','Color','k'); end, end
sgtitle('2025 casos: 49-51 Hz; carga y forma 0.8-1.2; sin tolerancias adicionales');
end
