function r = comparacion_estimacion_paper
%COMPARACION_ESTIMACION_PAPER Mismos registros y ocho periodos de entrenamiento.
% mSOGI: Hackl/Landerer, TPEL2020, (11)-(13), frecuencia conocida,
% sin FLL. place como Remark II.2. Polos elegidos aqui, no los del ensayo original.
root=fileparts(fileparts(mfilename('fullpath'))); folder=fullfile(root,'paper','estimacion');
if ~isfolder(folder), mkdir(folder); end
data=load(fullfile(root,'matlab','resultados','medidas_simulink','experimento_medidas.mat'),'saved');
base=data.saved; bank=base.banco_nominal; rows=table; audit=0; lastTrace=[];
for a=1:size(base.registros_ruido,1)
    for b=1:size(base.registros_ruido,2)
        rec=base.registros_ruido{a,b}; N=rec.opciones.N; P=rec.opciones.Periodos/2;
        train=1:N*P; U=reshape(rec.u_medida(train),N,P); I=reshape(rec.i_medida(train),N,P);
        labels=["un_periodo","promedio","regularizado","mSOGI a=40","mSOGI a=100"];
        for j=1:5
            if j<=3
                fit=ajustar_banco_medido(U,I,50,bank.z,labels(j));
            else
                decay=[40 100]; [est,err]=observer([rec.u_medida(train),rec.i_medida(train)], ...
                    rec.t(train),N,decay(j-3)); audit=max(audit,err);
                fit=ajustar_banco_medido(est(:,1),est(:,2),50,bank.z,"promedio");
                if a==4 && b==1 && j==5, lastTrace=table(rec.t(1:N),rec.u(1:N),est(:,1),rec.i(1:N),est(:,2), ...
                        'VariableNames',{'t_s','u_real_V','u_estimada_V','i_real_A','i_estimada_A'}); end
            end
            D=derivador_periodico(N); u=rec.u(1:N); il=rec.i(1:N);
            J=[D,ones(N,1);ones(1,N),0]\[u;0];
            Phi=[D*u,J(1:N),(D*D+bank.z*eye(N))\(D*u)];
            current=il+Phi*fit.coef; ideal=mean(u.*il)/mean(u.^2)*u;
            row=table(100*rec.opciones.Ruido,rec.opciones.Semilla,labels(j),P, ...
                norm(fit.coef-bank.coef)/norm(bank.coef),sqrt(mean(current.^2)), ...
                mean(current.^2)-mean(ideal.^2),sqrt(mean((fit.u_ajustada-u).^2)), ...
                'VariableNames',{'Ruido_pct','Semilla','Metodo','Periodos_disponibles', ...
                'Error_coef_rel','Is_RMS_A','Exceso_A2','Error_u_V'});
            rows=[rows;row]; %#ok<AGROW>
        end
    end
end
% Las tres primeras columnas reproducen el ensayo anterior de circuito Simulink.
delta=0;
for k=1:height(base.ruido)
    old=base.ruido(k,:); idx=rows.Ruido_pct==old.Ruido_pct & rows.Semilla==old.Semilla & rows.Metodo==old.Metodo;
    delta=max(delta,abs(rows.Is_RMS_A(idx)-old.I_fuente_A));
end
assert(delta<1e-6 && audit<1e-8);
r=struct('tabla',rows,'traza',lastTrace,'Error_reproduccion_A',delta,'Error_observador',audit);
writetable(rows,fullfile(folder,'comparacion.csv')); writetable(lastTrace,fullfile(folder,'traza.csv'));
save(fullfile(folder,'estimacion.mat'),'r');
report=struct('passed',true,'Filas',height(rows),'Registros',numel(base.registros_ruido), ...
    'Error_reproduccion_A',delta,'Error_observador',audit, ...
    'fuente',"https://arxiv.org/html/1902.04653v5#S2.SS4", ...
    'alcance',"Ocho registros sinteticos originales, dos semillas, 0/0.1/1/3% de ruido. Misma ventana pasada de ocho periodos. mSOGI conoce frecuencias 50/100/150 Hz y usa polos -a+j*n*omega, a=40/100 s^-1; sin FLL. Propagacion FOH exacta entre muestras pasadas y prediccion periodica al final del entrenamiento. Auditoria lsim del mismo observador. No es ensayo de frecuencia desconocida ni replica del benchmark publicado.");
fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0); fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true)); fclose(fid);
disp(rows);
end

function [period,delta]=observer(y,t,N,decay)
w=2*pi*50; J=[0 -1;1 0]; A=w*blkdiag(J,2*J,3*J); C=[1 0 1 0 1 0];
poles=reshape([-decay+1i*w*(1:3);-decay-1i*w*(1:3)],[],1);
L=place(A',C',poles)'; F=A-L*C; dt=t(2)-t(1);
H=expm([F,L,zeros(6,1);zeros(1,7),1;zeros(1,8)]*dt);
x=zeros(6,2); states=zeros(numel(t),6,2);
for k=2:numel(t)
    x=H(1:6,1:6)*x+H(1:6,7)*y(k-1,:)+H(1:6,8)*(y(k,:)-y(k-1,:))/dt;
    states(k,:,:)=reshape(x,1,6,2);
end
delta=0;
for j=1:2
    ref=lsim(ss(F,L,eye(6),zeros(6,1)),y(:,j),t-t(1),zeros(6,1),'foh');
    delta=max(delta,max(abs(ref-states(:,:,j)),[],'all'));
end
period=zeros(N,2); phase=(0:N-1)'*.02/N;
% Propagar desde el ultimo dato hasta la siguiente frontera de periodo.
for k=1:N, period(k,:)=C*expm(A*(phase(k)+dt))*x; end
end
