function agregar_fuente_impedante_tsc(name,cfg)
%AGREGAR_FUENTE_IMPEDANTE_TSC Red R-L continua, sin retrasos artificiales.
% Lg>0: Lg*dis/dt=e-Rg*is-u; u=vCp+Rp*(is-il-i1-i2).
% Lg=0: u=[Rp*e+Rg*vCp-Rg*Rp*(il+i1+i2)]/(Rp+Rg).
% El puerto de estado previo al reset evita feedthrough del reset de iLC.
f=cfg.Fuente; ws=get_param(name,'ModelWorkspace');
assignin(ws,'tm_Rg',f.R_ohm); assignin(ws,'tm_Lg',f.L_H); assignin(ws,'tm_iG0',f.i0_A);
set_param(name,'Description',char(cfg.alcance));
% Conserva e(t) como excitacion, y sustituye todos sus destinos anteriores
% por la nueva tension u(t), incluida la medida que recibe el supervisor.
ports=get_param([name,'/Tension_continua'],'PortHandles');
oldLine=get_param(ports.Outport,'Line'); destinations=get_param(oldLine,'DstPortHandle');
delete_line(oldLine);
add_block('simulink/Signal Routing/Mux',[name,'/Variables_red'], ...
    'Inputs','6','Position',[120,990,125,1160]);
add_block('simulink/User-Defined Functions/Fcn',[name,'/Tension_bornes'], ...
    'Expr','u(2)+tm_Rp*(u(6)-u(3)-u(4)-u(5))', ...
    'Position',[185,995,485,1045]);
if f.L_H>0
    add_block('simulink/Continuous/Integrator',[name,'/Estado_I_red'], ...
        'InitialCondition','tm_iG0','Position',[720,1070,770,1110]);
    add_block('simulink/Signal Routing/Mux',[name,'/Variables_KVL_red'], ...
        'Inputs','3','Position',[390,1080,395,1160]);
    add_block('simulink/User-Defined Functions/Fcn',[name,'/Derivada_I_red'], ...
        'Expr','(u(1)-tm_Rg*u(2)-u(3))/tm_Lg','Position',[460,1090,650,1130]);
    is='Estado_I_red/1';
    add_line(name,'Tension_continua/1','Variables_KVL_red/1','autorouting','on');
    add_line(name,is,'Variables_KVL_red/2','autorouting','on');
    add_line(name,'Tension_bornes/1','Variables_KVL_red/3','autorouting','on');
    add_line(name,'Variables_KVL_red/1','Derivada_I_red/1');
    add_line(name,'Derivada_I_red/1','Estado_I_red/1');
else
    set_param([name,'/Tension_bornes'],'Expr', ...
        '(tm_Rp*u(1)+tm_Rg*u(2)-tm_Rg*tm_Rp*(u(3)+u(4)+u(5)))/(tm_Rp+tm_Rg)');
    % La sexta entrada no se usa para Lg=0; conectarla a is crearia un
    % ciclo estructural innecesario aunque su coeficiente fuera cero.
    is='Corriente_fuente/1';
end
inputs={'Tension_continua/1','Estado_VCp/1','Carga_exogena/1','LC_1/3','LC_2/3',is};
if f.L_H==0, inputs{6}='Tension_continua/1'; end
for k=1:6, add_line(name,inputs{k},sprintf('Variables_red/%d',k),'autorouting','on'); end
add_line(name,'Variables_red/1','Tension_bornes/1');
newPorts=get_param([name,'/Tension_bornes'],'PortHandles');
for dst=reshape(destinations,1,[]), add_line(name,newPorts.Outport,dst,'autorouting','on'); end
if f.L_H>0
    % Se mide el estado de corriente de la fuente. La suma de ramas queda
    % como comprobacion independiente de la ecuacion de nodo.
    for dst={'Medidas/2','Registro/3'}
        delete_line(name,'Corriente_fuente/1',dst{1});
        add_line(name,is,dst{1},'autorouting','on');
    end
end
add_block('simulink/Signal Routing/Mux',[name,'/Variables_energia_red'], ...
    'Inputs','4','Position',[850,990,855,1090]);
inputs={'Tension_continua/1',is,'Tension_bornes/1','Carga_exogena/1'};
for k=1:4, add_line(name,inputs{k},sprintf('Variables_energia_red/%d',k),'autorouting','on'); end
names={'Fuente','Carga','Red'}; expr={'u(1)*u(2)','u(3)*u(4)','tm_Rg*u(2)^2'};
for k=1:3
    y=990+70*(k-1);
    add_block('simulink/User-Defined Functions/Fcn',[name,'/Potencia_',names{k}], ...
        'Expr',expr{k},'Position',[900,y,1060,y+40]);
    add_block('simulink/Continuous/Integrator',[name,'/Energia_',names{k}], ...
        'Position',[1110,y,1150,y+40]);
    add_line(name,'Variables_energia_red/1',['Potencia_',names{k},'/1'],'autorouting','on');
    add_line(name,['Potencia_',names{k},'/1'],['Energia_',names{k},'/1']);
end
add_block('simulink/Signal Routing/Mux',[name,'/Registro_red'],'Inputs','6', ...
    'Position',[1210,990,1215,1150]);
inputs={'Tension_continua/1',is,'Tension_bornes/1','Energia_Fuente/1','Energia_Carga/1','Energia_Red/1'};
for k=1:6, add_line(name,inputs{k},sprintf('Registro_red/%d',k),'autorouting','on'); end
add_block('simulink/Sinks/To Workspace',[name,'/Guardar_red'], ...
    'VariableName','registro_fuente','SaveFormat','Timeseries','MaxDataPoints','inf', ...
    'Position',[1280,1040,1390,1080]);
add_line(name,'Registro_red/1','Guardar_red/1');
end
