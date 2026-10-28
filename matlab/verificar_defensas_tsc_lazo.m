function informe = verificar_defensas_tsc_lazo(r, opciones)
%VERIFICAR_DEFENSAS_TSC_LAZO Comprueba que la auditoria detecta corrupciones.
% Modifica copias locales de un registro ya simulado; no cambia el SLX,
% sus resultados guardados ni el manuscrito. No vuelve a simular.
arguments
    r (1,1) struct
    opciones.Exportar (1,1) logical = false
    opciones.CarpetaSalida (1,1) string = ""
end
baseline=verificar_tsc_lazo_cerrado(r,InspeccionarModelo=false);
assert(baseline.passed,'TSCLazo:DefensasBase','La traza original debe superar la auditoria.');
first=find(~cellfun(@isempty,r.ventanas_decision),1);
event=find(ismember(r.eventos.Tipo,["ON","OFF"]),1);
assert(~isempty(first) && ~isempty(event), ...
    'TSCLazo:DefensasDatos','Se necesita al menos una decision y una maniobra confirmada.');
names=["coste_RMS_alterado";"muestra_futura";"maniobra_omitida";"energia_entrada_alterada"];
expected=["TSCLazo:Gram";"TSCLazo:VentanaCausal"; ...
    "TSCLazo:RegistroManiobras";"TSCLazo:Energia"];
observed=strings(4,1); captured=false(4,1);
deltaEnergy=max(.01,1e-3*max(r.trazas.energia_J));
for j=1:4
    changed=r;
    switch j
        case 1
            changed.ventanas_decision{first}.costesRMS(1)= ...
                changed.ventanas_decision{first}.costesRMS(1)+.1;
        case 2
            changed.ventanas_decision{first}.t(end)= ...
                changed.decisiones.Tiempo_s(first)+changed.config.TsSupervision;
        case 3
            changed.eventos(event,:)=[];
        case 4
            changed.trazas.integralEntrada_J(end)= ...
                changed.trazas.integralEntrada_J(end)+deltaEnergy;
    end
    try
        verificar_tsc_lazo_cerrado(changed,InspeccionarModelo=false);
    catch failure
        observed(j)=string(failure.identifier);
        captured(j)=observed(j)==expected(j);
    end
    assert(captured(j),'TSCLazo:DefensaNoDetectada', ...
        'Corrupcion %s: esperado %s; observado %s.',names(j),expected(j),observed(j));
end
tests=table(names,expected,observed,captured, ...
    'VariableNames',{'Alteracion','ErrorEsperado','ErrorObservado','Detectada'});
informe=struct('passed',all(captured),'original_valido',baseline.passed, ...
    'pruebas',tests,'incremento_energia_J',deltaEnergy, ...
    'alcance','Ensayos negativos de la auditoria sobre copias locales; ninguna nueva simulacion.');
if opciones.Exportar
    folder=opciones.CarpetaSalida;
    if strlength(folder)==0
        folder=fullfile(fileparts(mfilename('fullpath')),'resultados', ...
            'tsc_lazo_cerrado','pruebas_auxiliares');
    end
    if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,'defensas_auditoria.mat'),'informe');
    writetable(tests,fullfile(folder,'defensas_auditoria.csv'));
    copy=informe; copy.pruebas=table2struct(tests);
    fid=fopen(fullfile(folder,'defensas_auditoria.json'),'w'); assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(copy,PrettyPrint=true));
end
disp(tests);
end
