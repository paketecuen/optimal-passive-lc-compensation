function agregar_control_red_tsc(name,cfg)
%AGREGAR_CONTROL_RED_TSC Cablea medidas, estados y revalidacion en cada ON.
set_param([name,'/Supervisor_temporal'],'FunctionName','tsc_supervisor_red_sfun');
set_param([name,'/Medidas'],'Inputs','8');
add_line(name,'Tension_continua/1','Medidas/8','autorouting','on');
delete_line(name,'VC_actuales/1','Supervisor_temporal/3');
delete_line(name,'LC_1/2','VC_actuales/1');
delete_line(name,'LC_2/2','VC_actuales/2');
delete_block([name,'/VC_actuales']);
add_block('simulink/Signal Routing/Mux',[name,'/Estado_actual'], ...
    'Inputs','6','Position',[1085,240,1090,350]);
is='Estado_I_red/1';
if cfg.Fuente.L_H==0
    add_block('simulink/Sources/Constant',[name,'/Coordenada_red_nula'], ...
        'Value','0','Position',[1000,245,1040,275]);
    is='Coordenada_red_nula/1';
end
inputs={is,'Estado_VCp/1','LC_1/3','LC_1/2','LC_2/3','LC_2/2'};
for k=1:6, add_line(name,inputs{k},sprintf('Estado_actual/%d',k),'autorouting','on'); end
add_line(name,'Estado_actual/1','Supervisor_temporal/3');
add_block('simulink/Sinks/To Workspace',[name,'/Guardar_contextos'], ...
    'VariableName','contextos_admision','SaveFormat','Timeseries','MaxDataPoints','inf', ...
    'SampleTime','tsc_cfg.p.T','Position',[1435,335,1595,375]);
add_line(name,'Supervisor_temporal/4','Guardar_contextos/1','autorouting','on');
for j=1:2
    ex=[name,sprintf('/Ejecutor_%d',j)]; lat=[ex,'/Latch_evento'];
    set_param([ex,'/Entrada_latch'],'Inputs','8');
    inputs={'Estado_red','Contexto'};
    for k=1:2
        add_block('simulink/Ports & Subsystems/In1',[ex,'/',inputs{k}], ...
            'Port',num2str(k+6),'Position',[20,440+60*k,50,460+60*k]);
        add_line(ex,[inputs{k},'/1'],sprintf('Entrada_latch/%d',k+6),'autorouting','on');
    end
    set_param([lat,'/Actualizar'],'FunctionName','tsc_ejecutor_red_sfun', ...
        'Parameters',sprintf('tsc_cfg,%d',j));
    add_block('simulink/Ports & Subsystems/Out1',[lat,'/Guardia'], ...
        'Port','2','InitialOutput','zeros(15,1)','Position',[320,165,350,185]);
    add_line(lat,'Actualizar/2','Guardia/1');
    add_block('simulink/Ports & Subsystems/Out1',[ex,'/Guardia'], ...
        'Port','2','Position',[720,440,750,460]);
    add_line(ex,'Latch_evento/2','Guardia/1','autorouting','on');
    add_line(name,'Estado_actual/1',sprintf('Ejecutor_%d/7',j),'autorouting','on');
    add_line(name,'Supervisor_temporal/4',sprintf('Ejecutor_%d/8',j),'autorouting','on');
    % El registro conserva tiempo y estados del evento dentro de la señal;
    % su lectura T/63 no altera la deteccion continua de la maniobra.
    add_block('simulink/Sinks/To Workspace',[name,sprintf('/Guardar_guardia_%d',j)], ...
        'VariableName',sprintf('guardia%d',j),'SaveFormat','Timeseries','MaxDataPoints','inf', ...
        'SampleTime','tsc_cfg.TsSupervision','Position',[850,460+200*(j-1),1010,495+200*(j-1)]);
    add_line(name,sprintf('Ejecutor_%d/2',j),sprintf('Guardar_guardia_%d/1',j),'autorouting','on');
end
set_param(name,'Description',char(cfg.alcance));
end
