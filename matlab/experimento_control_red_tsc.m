function r = experimento_control_red_tsc(opciones)
%EXPERIMENTO_CONTROL_RED_TSC Aisla efecto de predictor y limites de arranque.
arguments
    opciones.Reanudar (1,1) logical = true
    opciones.Visible (1,1) string = "off"
end
root=fileparts(mfilename('fullpath')); folder=fullfile(root,'resultados','control_red');
if ~isfolder(folder), mkdir(folder); end
nets=["R04_L08","R15_L30"]; R=[.4 1.5]; L=[.8e-3 3e-3];
labels=["Anterior","Coste bornes","Coste red","Red 1100V 16A","Red 1100V 8A","Red 900V 16A"];
vlim=[Inf Inf 1100 1100 900]; ilim=[Inf Inf 16 8 16];
runs=cell(2,6); audits=cell(2,6); files=strings(2,6); summary=table;
for net=1:2
    for k=1:6
        if k==1
            files(net,k)=fullfile('..','fuente_impedante',nets(net)+"_1.mat");
            data=load(fullfile(folder,files(net,k)),'run','audit'); run=data.run; audit=data.audit;
            limV=Inf; limI=Inf; rejected=NaN; admitted=nnz(run.eventos.Tipo=="ON");
        else
            model="red"; if k==2, model="bornes"; end
            limV=vlim(k-1); limI=ilim(k-1);
            cfg=configurar_control_red_tsc(configurar_fuente_tsc(configurar_tsc_lazo_cerrado, ...
                R_ohm=R(net),L_H=L(net)),ModeloCoste=model, ...
                MaxVC_V=[limV limV],MaxI_A=[limI limI]);
            files(net,k)=sprintf('%s_%d.mat',nets(net),k);
            file=fullfile(folder,files(net,k));
            if opciones.Reanudar && isfile(file)
                data=load(file,'run','audit');
                assert(isequaln(data.run.config,cfg),'TSCRed:Cache','Configuracion distinta de la guardada.');
                run=data.run; audit=data.audit;
            else
                fprintf('%s, %s: simulando...\n',nets(net),labels(k));
                run=simular_tsc_lazo_cerrado(cfg); audit=verificar_control_red_tsc(run);
                save(file,'run','audit','-v7.3');
            end
            rejected=audit.control.Rechazos_plan; admitted=audit.control.ON_admitidos;
        end
        runs{net,k}=run; audits{net,k}=audit; z=run.resumen;
        meets=all(max(abs(run.trazas.vCs),[],1)<=[limV limV]) && ...
            all(max(abs(run.trazas.iRamas),[],1)<=[limI limI]);
        states=unique(2*run.trazas.estadoAplicado(:,1)+run.trazas.estadoAplicado(:,2))';
        row=table(nets(net),labels(k),limV,limI,z.U_RMS_V,z.I_fuente_RMS_A,z.VC_pico_V, ...
            z.I_rama_pico_A,z.P_carga_W,z.P_perdidas_red_W,z.P_perdidas_media_W, ...
            z.Maniobras,rejected,z.Rechazos,admitted,meets,join(string(states),"/"), ...
            'VariableNames',{'Red','Control','LimiteVC_V','LimiteI_A','U_RMS_V','Is_RMS_A', ...
            'PicoVC_V','PicoI_A','P_carga_W','P_red_W','P_banco_W','Maniobras', ...
            'Candidatos_bloqueados','ON_rechazados','ON_admitidos','Limites_observados','Estados'});
        summary=[summary;row]; %#ok<AGROW>
        fprintf('  Is %.6f A; VC %.2f V; I rama %.3f A; %d bloqueos, %d ON\n', ...
            z.I_fuente_RMS_A,z.VC_pico_V,z.I_rama_pico_A,rejected,admitted);
        writetable(summary,fullfile(folder,'comparacion.csv'));
        writetable(run.decisiones,fullfile(folder,sprintf('%s_%d_decisiones.csv',nets(net),k)));
        if k>1 && ~isempty(audit.admisiones)
            writetable(audit.admisiones,fullfile(folder,sprintf('%s_%d_admisiones.csv',nets(net),k)));
        end
    end
end
r=struct('registros',{runs},'auditorias',{audits},'resumen',summary,'archivos',files, ...
    'redes',nets,'controles',labels,'alcance', ...
    "Comparacion causal en Simulink; limites de arranque condicionales y carga exogena. La tension interna e se mide. No se usa el futuro de carga.");
index=rmfield(r,'registros'); save(fullfile(folder,'campana.mat'),'index');
r.figuras=graficas_control_red_tsc(r,Visible=opciones.Visible,CarpetaSalida=folder);
crear_modelo_tsc_lazo_cerrado(runs{1,5}.config);
disp(summary);
end
