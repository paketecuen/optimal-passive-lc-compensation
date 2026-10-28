function r = comparacion_sintonias_paper
%COMPARACION_SINTONIAS_PAPER Redes fijas, perdidas y arranque desde LC descargada.
% P/N: x=[sqrt(L)i,sqrt(C)vC], E=|x|^2/2. La solucion es
% x=x_per-exp(A*(t-ton))*x_per(ton), con x(ton)=0.
root=fileparts(fileparts(mfilename('fullpath'))); folder=fullfile(root,'paper','sintonias');
if ~isfolder(folder), mkdir(folder); end
cases=paper_casos; c=cases{5}; scale=dimensionar_referencia_230V(1000); kZ=scale.escala_impedancia;
original=c; amp=230/sqrt(3); original.e=amp*c.e; original.u=amp*c.u;
original.i=amp*c.i/kZ; original.qstar=amp*c.qstar/kZ;
original.er=amp*c.er; original.Y=c.Y/kZ; original.Yr=c.Yr/kZ;
names=["Sin banco","Familia z=1.45","Familia z=1.8","Familia z=2.04115", ...
    "Familia tres elementos","LC media fija","Jeon tres LC","p=14% ajustado","p=7% ajustado", ...
    "p=14% ajuste 3%","p=7% ajuste 3%"];
banks=cell(size(names)); banks{1}=struct('ramas',zeros(0,3));
zvalues=[1.45 1.8 2.041148796 27/13];
for j=1:4
    z=zvalues(j); cp=(7-3*z)/20; gamma=max(0,(27-13*z)/(5*z)); alpha=3*(z-1)*(z-4)*(z-9)/(20*z);
    branches=[0 cp/c.omega 0;1/(c.omega*alpha) alpha/(c.omega*z) 0];
    if gamma>1e-12, branches(end+1,:)=[1/(c.omega*gamma) 0 0]; end %#ok<AGROW>
    branches(:,1)=branches(:,1)*kZ; branches(:,2)=branches(:,2)/kZ;
    banks{j+1}=struct('ramas',branches);
end
banks{6}=banks{5}; banks{6}.ramas(2,1)=2*banks{6}.ramas(2,1); banks{6}.ramas(2,2)=banks{6}.ramas(2,2)/2;
banks{7}=paper_disenar(original,"Jeon polos previos");
for j=1:4
    z=1./[.14 .07 .14 .07]; z=z(j); fit=original; if j>2, fit=cases{6}; end
    D=fit.Dr; v=D*fit.er;
    Phi=[v,(z*eye(size(D))+D*D)\v]; coef=lsqnonneg(Phi,-fit.Yr*fit.er);
    b=zeros(0,3);
    if coef(1)>1e-12, b(end+1,:)=[0 coef(1)/c.omega 0]; end %#ok<AGROW>
    if coef(2)>1e-12, b(end+1,:)=[1/(c.omega*coef(2)) coef(2)/(c.omega*z) 0]; end %#ok<AGROW>
    banks{7+j}=struct('ramas',b);
end
rows=table; records=cell(2,numel(names),2); components=table; worstProjection=0;
for j=1:numel(names)
    b=banks{j};
    for a=1:size(b.ramas,1)
        v=b.ramas(a,:); p=c.omega^2*v(1)*v(2);
        components=[components;table(names(j),a,v(1),v(2),p, ...
            'VariableNames',{'Diseno','Rama','L_H','C_F','p'})]; %#ok<AGROW>
    end
    for wave=1:2
        test=original; label="Original 230 V"; if wave==2, test=cases{6}; label="Moderada 3%"; end
        for loss=1:2
            b=banks{j}; QL=Inf; QC=Inf;
            if loss==2
                QL=100; QC=500;
                for a=1:size(b.ramas,1)
                    L=b.ramas(a,1); C=b.ramas(a,2);
                    b.ramas(a,3)=test.omega*L/QL;
                    if C>0, b.ramas(a,3)=b.ramas(a,3)+1/(test.omega*C*QC); end
                end
            end
            value=paper_evaluar_red(test,b); startup=arranque(test,b,4096);
            if loss==1 && ismember(j,[2:5 7])
                error=norm(value.q_nodos-test.qstar)/sqrt(test.N);
                worstProjection=max(worstProjection,error); assert(error<1e-8);
            end
            row=table(label,names(j),QL,QC,value.Elementos,value.Is_RMS_A,value.PF_fuente, ...
                value.P_banco_W,value.E_media_J,value.V_elemento_pico_V,value.I_rama_pico_A, ...
                startup.VC_pico_V,startup.I_pico_A,startup.Espera_max_ms, ...
                'VariableNames',{'Tension','Diseno','QL','QC','Elementos','Is_RMS_A','PF', ...
                'Perdidas_W','E_media_J','V_elemento_periodico_V','I_rama_periodico_A', ...
                'VC_arranque_V','I_arranque_A','Espera_max_ms'});
            rows=[rows;row]; records{wave,j,loss}=struct('banco',b,'periodico',value,'arranque',startup); %#ok<AGROW>
        end
    end
end
% Auditoria Simulink de las LC que afrontan el pico; matriz de estados P/N.
model='paper_arranque_lc'; if bdIsLoaded(model), close_system(model,0); end
load_system('simulink'); new_system(model);
set_param(model,'Solver','ode4','SolverType','Fixed-step','FixedStep','h', ...
    'StopTime','tend','SaveOutput','off','SaveTime','off','ReturnWorkspaceOutputs','on');
add_block('simulink/Sources/Clock',[model '/Tiempo']);
add_block('simulink/User-Defined Functions/Fcn',[model '/Tension'], ...
    'Expr','sqrt(2)*(U1*cos(w*(u+ton))+U2*cos(2*w*(u+ton))+U3*cos(3*w*(u+ton)))');
add_block('simulink/Continuous/State-Space',[model '/Estados_PN'], ...
    'A','A','B','B','C','eye(2)','D','zeros(2,1)','X0','zeros(2,1)');
add_block('simulink/Sinks/To Workspace',[model '/Registro'], ...
    'VariableName','xPN','SaveFormat','Timeseries');
add_line(model,'Tiempo/1','Tension/1'); add_line(model,'Tension/1','Estados_PN/1'); add_line(model,'Estados_PN/1','Registro/1');
simErrors=table; peakDelta=0;
for wave=1:2
    for j=[5 6 8 10 11]
        test=original; if wave==2, test=cases{6}; end
        record=records{wave,j,2};
        for a=1:numel(record.arranque.ramas)
            branch=record.arranque.ramas{a}; if isempty(branch), continue; end
            U=test.U; if wave==1, U=amp*U; end
            par=struct('A',branch.A,'B',[1/sqrt(branch.L);0],'h',.02/4096, ...
                'tend',.1,'w',test.omega,'ton',branch.ton,'U1',U(1),'U2',U(2),'U3',U(3));
            input=Simulink.SimulationInput(model); fields=fieldnames(par);
            ws=get_param(model,'ModelWorkspace');
            for a2=1:numel(fields), input=input.setVariable(fields{a2},par.(fields{a2})); assignin(ws,fields{a2},par.(fields{a2})); end
            out=sim(input); ref=trayectoria(test,branch,out.xPN.Time);
            ei=max(abs(out.xPN.Data(:,1)-ref(:,1)))/sqrt(branch.L);
            ev=max(abs(out.xPN.Data(:,2)-ref(:,2)))/sqrt(branch.C);
            assert(ei<1e-5 && ev<1e-3);
            simErrors=[simErrors;table(wave,names(j),a,ei,ev, ...
                'VariableNames',{'Tension_id','Diseno','Rama','Error_A','Error_V'})]; %#ok<AGROW>
        end
        fine=arranque(test,record.banco,8192);
        peakDelta=max(peakDelta,abs(fine.VC_pico_V-record.arranque.VC_pico_V));
    end
end
save_system(model,fullfile(folder,[model '.slx'])); close_system(model,0);
r=struct('tabla',rows,'componentes',components,'registros',{records},'auditoria',simErrors);
writetable(rows,fullfile(folder,'comparacion.csv')); writetable(components,fullfile(folder,'componentes.csv'));
writetable(simErrors,fullfile(folder,'simulink.csv')); save(fullfile(folder,'sintonias.mat'),'r','-v7.3');
report=struct('passed',true,'Casos',height(rows),'Simulaciones_Simulink',height(simErrors), ...
    'Error_proyeccion_A',worstProjection,'Error_Simulink_A',max(simErrors.Error_A), ...
    'Error_Simulink_V',max(simErrors.Error_V),'Refinamiento_pico_V',peakDelta, ...
    'p_min_familia',13/27,'QL',100,'QC',500, ...
    'alcance',"Fuente ideal fijada, carga LTI, escala 230 V/1000 W. Redes fijas. R/ESR elegidos a 50 Hz. Arranque de cada LC descargada con fase de minimo desajuste en un periodo, observacion 5 periodos; ramas simples en regimen. No son picos globales bajo todas las maniobras ni simulacion de semiconductores.");
fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0); fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true)); fclose(fid);
disp(rows);
end

function r=arranque(c,b,N)
D=c.omega*c.Dr; I=eye(size(D)); phase=(0:N-1)'*2*pi/N;
K=interpolar_periodico(c.V,phase); records=cell(size(b.ramas,1),1); vp=0; ip=0; delay=0;
for a=1:size(b.ramas,1)
    L=b.ramas(a,1); C=b.ramas(a,2); R=b.ramas(a,3); if L==0 || C==0, continue; end
    q=(I+R*C*D+L*C*D*D)\(C*D*c.er); vc=c.er-L*D*q-R*q;
    X=K*[sqrt(L)*q,sqrt(C)*vc]; [~,idx]=min(sum(X.^2,2));
    branch=struct('L',L,'C',C,'A',[-R/L,-1/sqrt(L*C);1/sqrt(L*C),0], ...
        'coef',[sqrt(L)*q,sqrt(C)*vc],'ton',phase(idx)/c.omega);
    t=(0:5*N)'*.02/N; x=trayectoria(c,branch,t);
    vp=max(vp,max(abs(x(:,2)))/sqrt(C)); ip=max(ip,max(abs(x(:,1)))/sqrt(L));
    delay=max(delay,branch.ton*1000); branch.t=t; branch.x=x; records{a}=branch;
end
r=struct('VC_pico_V',vp,'I_pico_A',ip,'Espera_max_ms',delay,'ramas',{records});
end

function x=trayectoria(c,b,t)
K=interpolar_periodico(c.V,c.omega*(t+b.ton)); per=K*b.coef;
initial=interpolar_periodico(c.V,c.omega*b.ton)*b.coef;
[V,E]=eig(b.A); eigs=diag(E); transient=real((V*((V\initial').*exp(eigs*t')))');
x=per-transient;
end
