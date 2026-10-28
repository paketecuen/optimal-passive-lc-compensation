function r = experimento_fuente_tsc(opciones)
%EXPERIMENTO_FUENTE_TSC Campana reproducible de red y arranque en Simulink.
% Rg/Lg son hipotesis de ensayo, no limites normativos ni datos de una red.
arguments
    opciones.Visible (1,1) string {mustBeMember(opciones.Visible,["on","off"])} = "off"
    opciones.Reanudar (1,1) logical = true
end
root=fileparts(mfilename('fullpath'));
folder=fullfile(root,'resultados','fuente_impedante');
if ~isfolder(folder), mkdir(folder); end
names=["Ideal","R04_L08","R15_L30"]; resistances=[0 .4 1.5]; inductances=[0 .8e-3 3e-3];
labels=["GA","Fundamental","Fijo 00","Fijo 10","Fijo 01","Fijo 11"];
fixed=[0 0;1 0;0 1;1 1]; runs=cell(3,6); audits=cell(3,6); rows=struct([]);
for net=1:3
    for k=1:6
        args=struct('CasoTension',"moderada");
        if k==1, args.Modo="ga";
        elseif k==2, args.Modo="convencional";
        else, args.Modo="fijo"; args.BitsFijos=fixed(k-2,:);
        end
        ac=namedargs2cell(args); cfg=configurar_fuente_tsc(configurar_tsc_lazo_cerrado(ac{:}), ...
            R_ohm=resistances(net),L_H=inductances(net));
        file=fullfile(folder,sprintf('%s_%d.mat',names(net),k));
        if opciones.Reanudar && isfile(file)
            stored=load(file,'run','audit');
            assert(isequaln(stored.run.config,cfg),'TSCFuente:Cache','La configuracion guardada ha cambiado.');
            run=stored.run; audit=stored.audit;
        else
            fprintf('Red %s, %s: simulando...\n',names(net),labels(k));
            run=simular_tsc_lazo_cerrado(cfg); audit=verificar_fuente_tsc(run);
            save(file,'run','audit','-v7.3');
        end
        runs{net,k}=run; audits{net,k}=audit;
        row=run.resumen; row.Red=names(net); row.Estrategia=labels(k);
        row.Rechazos_plan=sum(run.decisiones.RechazosPlan);
        rows=[rows;row]; %#ok<AGROW>
        fprintf('  U %.4f V RMS, Is %.6f A RMS, VC pico %.3f V, %d maniobras\n', ...
            row.U_RMS_V,row.I_fuente_RMS_A,row.VC_pico_V,row.Maniobras);
        writetable(struct2table(rows),fullfile(folder,'comparacion.csv'));
        writetable(run.eventos,fullfile(folder,sprintf('%s_%d_eventos.csv',names(net),k)));
        writetable(run.decisiones,fullfile(folder,sprintf('%s_%d_decisiones.csv',names(net),k)));
    end
end
summary=struct2table(rows);
comparisons=table;
for net=1:3
    subset=summary(summary.Red==names(net),:); [minimum,j]=min(subset.I_fuente_RMS_A(3:6));
    line=table(names(net),subset.I_fuente_RMS_A(1),subset.I_fuente_RMS_A(2), ...
        subset.Estrategia(j+2),minimum,100*(subset.I_fuente_RMS_A(1)/minimum-1), ...
        'VariableNames',{'Red','GA_RMS_A','Fundamental_RMS_A','Mejor_fijo','Mejor_fijo_RMS_A','Exceso_GA_porcentaje'});
    comparisons=[comparisons;line]; %#ok<AGROW>
end
r=struct('registros',{runs},'auditorias',{audits},'resumen',summary, ...
    'comparacion',comparisons,'redes',names,'estrategias',labels, ...
    'alcance',"Simulacion con carga exogena y fuente R-L. Arranque descargado. Comparacion de la trayectoria completa, incluidos transitorios y perdidas.");
writetable(comparisons,fullfile(folder,'comparacion_estrategias.csv'));
% El indice evita duplicar los registros continuos ya guardados por corrida.
index=rmfield(r,'registros'); index.archivos=strings(3,6);
for net=1:3
    for k=1:6, index.archivos(net,k)=sprintf('%s_%d.mat',names(net),k); end
end
save(fullfile(folder,'campana.mat'),'index');
r.figuras=graficas_fuente_tsc(r,Visible=opciones.Visible,CarpetaSalida=folder);
% Deja el modelo inspeccionable con el caso intermedio y supervisor GA.
crear_modelo_tsc_lazo_cerrado(runs{2,1}.config);
disp(comparisons);
end
