function r = comparacion_estado_del_arte(opciones)
%COMPARACION_ESTADO_DEL_ARTE Casos comunes y reimplementaciones transparentes.
arguments
    opciones.Reanudar (1,1) logical = true
    opciones.Recalcular (1,:) string = strings(1,0)
end
root=fileparts(fileparts(mfilename('fullpath'))); folder=fullfile(root,'paper','comparacion');
if ~isfolder(folder), mkdir(folder); end
cases=paper_casos; methods=["Sin compensar","C Kusters-Moore","Paralelo Page", ...
    "Formula m2","Busqueda directa LC","Busqueda LC restringida","Czarnecki L fija","Bezout Cauer", ...
    "Bezout Foster","Foster polos libres","Jeon polos previos","LP polos fijos","C optimizada circuito"];
previous=[];
if opciones.Reanudar && isfile(fullfile(folder,'comparacion.mat'))
    cached=load(fullfile(folder,'comparacion.mat'),'r'); previous=cached.r;
end
rows=table; records=cell(numel(cases),numel(methods));
for k=1:numel(cases)
    c=cases{k}; fprintf('Comparacion: %s\n',c.id);
    for j=1:numel(methods)
        method=methods(j);
        version=1;
        if ismember(method,["Busqueda directa LC","Czarnecki L fija"]), version=2; end
        if ~isempty(previous) && ~ismember(method,opciones.Recalcular)
            oldCase=find(cellfun(@(x)x.id==c.id,previous.casos),1);
            oldMethod=find(previous.metodos==method,1);
            oldVersion=1;
            if ~isempty(oldCase) && ~isempty(oldMethod) && isfield(previous.registros{oldCase,oldMethod},'version')
                oldVersion=previous.registros{oldCase,oldMethod}.version;
            end
            if ~isempty(oldCase) && ~isempty(oldMethod) && isequaln(c,previous.casos{oldCase}) && oldVersion==version
                records{k,j}=previous.registros{oldCase,oldMethod};
                row=previous.tabla(previous.tabla.Caso==c.id & previous.tabla.Metodo==method,:);
                if c.Rg>0 && ismember(method,["C Kusters-Moore","Paralelo Page"])
                    row.Estado="formula de tension fijada aplicada a u inicial; no es optimo con red";
                    records{k,j}.estado=row.Estado;
                end
                rows=[rows;row]; %#ok<AGROW>
                continue
            end
        end
        tic; [bank,status,extra]=paper_disenar(c,method); first=toc;
        measured=NaN;
        unavailable=startsWith(status,"no aplica") || startsWith(status,"q*") || contains(status,"no dan");
        if ~unavailable
            times=zeros(1,3);
            for rep=1:3, tic; paper_disenar(c,method); times(rep)=toc; end
            measured=median(times); value=paper_evaluar_red(c,bank);
            metrics=[value.Is_RMS_A,value.PF_fuente,value.PF_bornes,value.Elementos, ...
                value.E_media_J,value.E_max_J,value.V_elemento_pico_V,value.I_rama_pico_A, ...
                value.P_banco_W,value.P_red_W,value.P_carga_W];
        else, value=struct; metrics=NaN(1,11);
        end
        row=table(c.id,method,status,metrics(1),metrics(2),metrics(3),metrics(4), ...
            metrics(5),metrics(6),metrics(7),metrics(8),metrics(9),metrics(10),metrics(11), ...
            measured,first,'VariableNames',{'Caso','Metodo','Estado','Is_RMS_A','PF_fuente','PF_bornes', ...
            'Elementos','E_media_J','E_max_J','V_elemento_pico_V','I_rama_pico_A', ...
            'P_banco_W','P_red_W','P_carga_W','Tiempo_mediana_s','Tiempo_primera_s'});
        rows=[rows;row]; %#ok<AGROW>
        records{k,j}=struct('banco',bank,'estado',status,'extra',extra,'resultado',value,'version',version);
        fprintf('  %-23s %.6g A, %.3g ms\n',method,metrics(1),1e3*measured);
    end
    writetable(rows,fullfile(folder,'comparacion.csv'));
end
r=struct('casos',{cases},'metodos',methods,'registros',{records},'tabla',rows, ...
    'alcance',"Circuitos periodicos. Picos en 4096 puntos, no transitorios. Tiempo de diseno sin generacion del caso ni metricas. Reimplementaciones especificadas; no replicas de todos los algoritmos originales.");
save(fullfile(folder,'comparacion.mat'),'r','-v7.3');
disp(rows(:,{'Caso','Metodo','Is_RMS_A','Elementos'}));
end
