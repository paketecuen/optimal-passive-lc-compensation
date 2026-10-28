function r = comparacion_maniobras_paper
%COMPARACION_MANIOBRAS_PAPER Misma planta, solicitudes y condiciones por criterio.
root=fileparts(fileparts(mfilename('fullpath'))); folder=fullfile(root,'paper','maniobras');
if ~isfolder(folder), mkdir(folder); end
data=load(fullfile(root,'matlab','resultados','conmutacion_tsc','conmutacion_tsc.mat'),'saved');
base=data.saved; p=base.parametros; rows=table; events=table; records=cell(5,3);
policies=["inmediato","valvula","estado"];
for k=1:5
    pk=p; pk.R=base.resumen.R_ohm(k); pk.Rb=base.resumen.Rb_ohm(k); pk.v0=base.resumen.VC_inicial_V(k);
    for j=1:3
        label="Caso "+k+" / "+policies(j); a=paper_maniobra(pk,policies(j),label); records{k,j}=a;
        on=a.events.Accion=="conectar";
        row=table(k,policies(j),pk.R,pk.Rb,pk.v0,max(abs(a.x(:,1))),max(abs(a.x(:,2))), ...
            max(a.events.Retardo_ms(on)),mean(a.events.Retardo_ms(on)), ...
            max(abs(a.events.dI_dt_inicial_A_s(on))),max(a.events.Desajuste_mJ(on)), ...
            max(a.events.Valvula_min_posible_V(on)),a.balance_error,a.off_error, ...
            'VariableNames',{'Caso','Criterio','R_ohm','Rb_ohm','VC0_V','I_pico_A','VC_pico_V', ...
            'Espera_max_ms','Espera_media_ms','Pendiente_max_A_s','Desajuste_max_mJ', ...
            'Minimo_valvula_max_V','Residuo_energia_J','Residuo_OFF_A'});
        rows=[rows;row]; events=[events;a.events]; %#ok<AGROW>
    end
end
assert(max(rows.Residuo_energia_J)<2e-8 && max(rows.Residuo_OFF_A)<1e-8);
% Misma politica que la corrida historica: reproduce sus magnitudes.
oldPolicy=[1 2 3 3 2]; oldError=0;
for k=1:5
    a=records{k,oldPolicy(k)};
    oldError=max(oldError,abs(max(abs(a.x(:,2)))-base.resumen.VC_pico_V(k)));
end
assert(oldError<1e-8);
pk=p; pk.R=.05; pk.Rb=1000; pk.v0=0; pk.h=p.h/2;
fine=paper_maniobra(pk,"estado","Refinamiento"); a=records{3,3};
deltaV=abs(max(abs(fine.x(:,2)))-max(abs(a.x(:,2))));
deltaI=abs(max(abs(fine.x(:,1)))-max(abs(a.x(:,1))));
assert(deltaV<.01 && deltaI<.01);
r=struct('tabla',rows,'eventos',events,'registros',{records},'parametros',p);
writetable(rows,fullfile(folder,'comparacion.csv')); writetable(events,fullfile(folder,'eventos.csv'));
save(fullfile(folder,'maniobras.mat'),'r','-v7.3');
report=struct('passed',true,'Simulaciones',15,'Refinamientos',1,'Error_base_V',oldError, ...
    'Refinamiento_V',deltaV,'Refinamiento_A',deltaI,'Residuo_energia_J',max(rows.Residuo_energia_J), ...
    'alcance',"Los cinco casos historicos se cruzan con tres criterios; los casos 1-3 tienen las mismas condiciones fisicas. Tension de prueba no escalada a 230 V. Planta ODE, no semiconductor.");
fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true)); fclose(fid); disp(rows);
end
