function r = resolver_lc_m2(u,i,D,omega,opciones)
%RESOLVER_LC_M2 Proyeccion GA temporal y realizacion LC serie de dimension dos.
% D es D_tau, tau=omega*t. J elige la primitiva periodica de media cero.
% No usa FFT, autovalores, busqueda de polos ni un ajuste de error de tension.
% La comprobacion de cierre pertenece al modelo finito (D,u) suministrado.
% Las bandas numericas de signo no prueban ceros exactos ni pasividad fisica.
arguments
    u (:,1) double {mustBeReal,mustBeFinite}
    i (:,1) double {mustBeReal,mustBeFinite}
    D (:,:) double {mustBeReal,mustBeFinite}
    omega (1,1) double {mustBePositive,mustBeFinite} = 2*pi*50
    opciones.ToleranciaCierre (1,1) double {mustBePositive} = 1e-10
end
N = numel(u);
assert(numel(i)==N && isequal(size(D),[N,N]),'LCm2:Dimensiones','Dimensiones incompatibles.');
dotr = @(x,y) sum(x.*y)/N;
normr = @(x) sqrt(dotr(x,x));
U = normr(u);
assert(U>0,'LCm2:TensionNula','Se requiere tension alterna no nula.');
assert(abs(mean(u))<=1e-12*U,'LCm2:DC','La tension debe tener media cero.');
assert(norm(D+D','fro')<=1e-12*max(norm(D,'fro'),realmin), ...
    'LCm2:Derivador','D debe satisfacer D*=-D en la metrica RMS.');
assert(norm(D*ones(N,1))<=1e-12*max(norm(D,'fro')*sqrt(N),realmin), ...
    'LCm2:Derivador','D debe anular las constantes.');
% El sistema aumentado fija la media sin elegir una constante arbitraria.
M = [D,ones(N,1);ones(1,N),0];
assert(rcond(M)>128*eps,'LCm2:Primitiva','Primitiva no identificable con este D.');
J = @(x) primitiva(M,x,N);
v = D*u; j = J(u); d2u = D*v; A = -D*D;
assert(normr(D*j-u)<=1e-11*U,'LCm2:Primitiva','DJu no reproduce u.');
M0 = dotr(u,u); M1 = dotr(v,v); Mm = dotr(j,j); M2 = dotr(d2u,d2u);
delta = Mm*M1-M0^2;
assert(delta>1e-12*max(Mm*M1,realmin),'LCm2:DimensionNoDos', ...
    'Du y Ju no definen un plano bien resuelto: dimension uno o mal condicionamiento.');
s = (Mm*M2-M0*M1)/delta;
p = (M0*M2-M1^2)/delta;
d4u = D*D*d2u;
closure = d4u+s*d2u+p*u;
closureRelative = normr(closure)/(normr(d4u)+abs(s)*normr(d2u)+abs(p)*U);
assert(closureRelative<=opciones.ToleranciaCierre,'LCm2:CierreNoM2', ...
    'Los cuatro momentos no bastan: el residuo del cierre no es compatible con m2.');
disc = s*s-4*p;
assert(p>0 && s>0 && disc>128*eps*(s*s+4*p),'LCm2:CierreDegenerado', ...
    'Se requieren dos direcciones ciclicas positivas distintas.');
rho1 = dotr(i,v); rhom = dotr(i,j);
c = -(Mm*rho1+M0*rhom)/delta;
gamma = -(M0*rho1+M1*rhom)/delta;
qstar = c*v+gamma*j; % -(i contraction B_Q)*inverse(B_Q), B_Q=Du wedge Ju.
K = p*c*c-s*c*gamma+gamma*gamma;
nv = sqrt(M1); nj = sqrt(Mm);
Gn = [1,dotr(v,j)/(nv*nj);dotr(v,j)/(nv*nj),1];
eta = 128*eps*cond(Gn,inf);
amplitudeTolerance = eta*max([normr(i),normr(qstar),realmin]);
ctol = amplitudeTolerance/nv; gtol = amplitudeTolerance/nj;
Ktol = eta*max(p*c*c+abs(s*c*gamma)+gamma*gamma,realmin) ...
    +abs(2*p*c-s*gamma)*ctol+abs(2*gamma-s*c)*gtol;
czero = abs(c)<=ctol; gzero = abs(gamma)<=gtol; Kzero = abs(K)<=Ktol;
finite = ~czero && ~gzero && ~Kzero && sign(c)==sign(gamma) && sign(c)==sign(K);
if finite
    tipo = "LC serie positiva";
elseif czero && gzero
    tipo = "Sin componente compensable resuelta";
elseif gzero && c>ctol
    tipo = "Frontera C dentro de tolerancia";
elseif czero && gamma>gtol
    tipo = "Frontera L dentro de tolerancia";
elseif Kzero
    tipo = "Cerca de una recta resonante";
else
    tipo = "Proyeccion no realizable con una LC finita";
end

% L, C y paralelo: resultados previos de los autores. Las plantillas son oblicuas;
% se comparan caras, nunca se recortan por separado los dos coeficientes.
cC = max(-rho1/M1,0); gL = max(-rhom/Mm,0);
% En cada semirrecta se omite una amplitud no resuelta. Esto no equivale
% a recortar por separado los coeficientes de la proyeccion del paralelo.
if cC<=ctol, cC=0; end
if gL<=gtol, gL=0; end
qC = cC*v; qL = gL*j;
baseCoefficients = [0,0;0,gL;cC,0];
baseQ = [zeros(N,1),qL,qC];
baseI = sqrt(sum((i+baseQ).^2,1)/N);
[~,idx] = min(baseI);
parallel = baseCoefficients(idx,:);
if c>ctol && gamma>gtol
    parallel = [c,gamma];
end
qParallel = parallel(1)*v+parallel(2)*j;

% Fronteras: radicales del cierre, SOLO para los proyectores temporales
% Pi_lambda=(A-lambda_otro I)/(lambda-lambda_otro), no para la formula LC.
lambda = [(s-sqrt(disc))/2,(s+sqrt(disc))/2];
v1 = (A*v-lambda(2)*v)/(lambda(1)-lambda(2));
v2 = (A*v-lambda(1)*v)/(lambda(2)-lambda(1));
boundaryGain = [cC^2*M1,gL^2*Mm,dotr(i,v1)^2/dotr(v1,v1), ...
    dotr(i,v2)^2/dotr(v2,v2),0];
boundaryNames = ["C sola","L sola","Limite resonante 1","Limite resonante 2","Sin rama"];
if czero && gzero
    boundaryGain(:)=0;
end
[gain,ib] = max(boundaryGain);
if czero && gzero, ib=5; end
if finite, gain=dotr(qstar,qstar); end
infimum = sqrt(max(0,dotr(i,i)-gain));
L = NaN; C = NaN; qseries = NaN(N,1); circuitError = NaN;
voltageError = NaN; bladeNorm2 = NaN; reciprocalL = NaN; reciprocalS = NaN;
if finite
    L = gamma/K; S = p*c/K; C = 1/S;
    qseries = ((D*D+(S/L)*eye(N))\v)/L;
    circuitError = normr(qseries-qstar);
    x = D*qstar; y = J(qstar);
    % Expansion del marco reciproco de B_E=Dq* wedge Jq*.
    Gxy = [dotr(x,x),dotr(x,y);dotr(x,y),dotr(y,y)];
    bladeNorm2 = det(Gxy); % ||B_E||^2; NO es B_E^2.
    LS = Gxy\[dotr(u,x);dotr(u,y)];
    reciprocalL = LS(1); reciprocalS = LS(2);
    voltageError = normr(u-L*x-S*y);
end

% Tabla: las cotas no se presentan como circuitos fisicos instalados.
labels = ["Sin compensacion";"L sola";"C sola";"LC paralelo";"LC serie";"Proyeccion ideal";"Infimo LC serie"];
qq = [zeros(N,1),qL,qC,qParallel,qseries,qstar];
currents = [sqrt(sum((i+qq).^2,1)/N),infimum]';
P = dotr(u,i);
power = [sum(u.*(i+qq),1)/N,NaN]';
pf = P./(U*currents);
Lvalues = [NaN;inverseOrNan(gL)/omega;NaN;inverseOrNan(parallel(2))/omega;L/omega;NaN;NaN];
Cvalues = [NaN;NaN;positiveOrNan(cC)/omega;positiveOrNan(parallel(1))/omega;C/omega;NaN;NaN];
states = ["Referencia";"Rama opcional";"Rama opcional";"Optimo del cono paralelo";tipo;"Cota de proyeccion GA";"Limite; no afirma realizacion finita"];
if gL==0, states(2)="Optimo: omitir la bobina"; end
if cC==0, states(3)="Optimo: omitir el condensador"; end
if czero || gzero, states(4)="Cono paralelo: cara dentro de tolerancia"; end
if finite, states(7)="Alcanzado por la LC serie"; end
comparison = table(labels,1e3*Lvalues,1e6*Cvalues,currents,power,pf,states, ...
    'VariableNames',{'Alternativa','L_mH','C_uF','I_RMS_A','P_W','FP','Interpretacion'});
r = struct('u',u,'i',i,'D',D,'Ju',j,'Du',v,'omega',omega,'c',c,'gamma',gamma,'K',K, ...
    's',s,'p',p,'delta',delta,'closure_relative',closureRelative,'lambda',lambda, ...
    'q_ideal',qstar,'q_serie',qseries,'q_paralelo',qParallel,'serie_admisible',finite, ...
    'tipo',tipo,'cerca_frontera_numerica',czero||gzero||Kzero, ...
    'tolerancias',struct('c',ctol,'gamma',gtol,'K',Ktol), ...
    'L_normalizada',L,'C_normalizada',C,'L_H',L/omega,'C_F',C/omega, ...
    'error_corriente_rms',circuitError,'error_tension_rms',voltageError, ...
    'blade_norma_cuadrada',bladeNorm2,'L_marco_reciproco',reciprocalL, ...
    'S_marco_reciproco',reciprocalS,'tabla',comparison,'infimo_rms',infimum, ...
    'mejor_frontera',boundaryNames(ib),'ganancias_fronteras',boundaryGain, ...
    'coeficientes_paralelo',parallel,'coeficientes_mejor_L_o_C',baseCoefficients(idx,:), ...
    'mejor_L_o_C',labels(idx),'norma_tension',U,'potencia',P);
end

function j = primitiva(M,x,N)
y = M\[x;0]; j = y(1:N);
end

function y = inverseOrNan(x)
if x>0, y=1/x; else, y=NaN; end
end

function y = positiveOrNan(x)
if x>0, y=x; else, y=NaN; end
end
