function cases = paper_casos(N)
%PAPER_CASOS Fuentes comunes; cargas LTI se construyen mediante D_tau.
% Los datos armonicos solo identifican las excitaciones y auditan las redes.
arguments, N (1,1) double = 31, end
D=derivador_periodico(N); t=(0:N-1)'*2*pi/N; I=eye(N);
J=[D,ones(N,1);ones(1,N),0]\[I;zeros(1,N)]; J=J(1:N,:);
ids=["Jeltsema2014","Czarnecki1987_ideal","Czarnecki1987_red", ...
    "Czarnecki1981","Proyecto_m3","Proyecto_3pct","Proyecto_m2","Czarnecki1981_tabla"];
cases=cell(numel(ids),1);
for k=1:numel(ids)
    omega=1; Rg=0; Lg=0;
    switch k
        case 1
            n=[1 5]; U=[10 5]; Y=(I+2*D)\I;
            scope="RL 1 ohm/2 H, omega=1 rad/s, valores publicados.";
        case {2,3}
            n=[1 5 7]; U=[100 3.5 2]; Y=(.5*I+sqrt(3)/2*D)\I;
            if k==3, Rg=.05/sqrt(26); Lg=5*Rg; end
            scope="RL |Z1|=1 ohm, PF1=.5, omega=1 rad/s; Zg=.05 ohm y Xg/Rg=5 en variante red.";
        case 4
            n=[1 2 3]; U=[100 50 30]; Yp=D/12+(I+6*D)\I;
            Y=(I+(.5*I+.5*D)*Yp)\Yp;
            scope="Carga exacta de Fig.1 Czarnecki1981: (.5H+.5ohm) serie [C=1/12F paralelo (6H+1ohm)].";
        case {5,6}
            omega=2*pi*50; n=[1 2 3]; U=[1 1 1];
            Y=.8*I+.55*D+1.7*J+.75*((D*D+6*I)\D);
            scope="m3 original, U_rms=sqrt(3) V, componentes SI a 50 Hz.";
            if k==6
                U=230*[1 .03 .03]/norm([1 .03 .03]);
                scale=dimensionar_referencia_230V(1000); Y=Y/scale.escala_impedancia;
                scope="Misma carga LTI escalada del proyecto, u=230 V RMS, omega=2*pi*50, P=1000 W; sin perdidas del compensador.";
            end
        case 7
            omega=2*pi*50; n=[1 2]; U=[1 1]; Y=.8*I+.3*D+.7*J;
            scope="Ejemplo m2 positivo original, tension de prueba sin escalar, 50 Hz.";
        case 8
            n=[1 2 3]; U=[100 50 30]; values=[.02749-.07419i,.01231+.09060i,.04264+.26810i];
            Y=.01*I;
            for a=1:3
                co=cos(n(a)*t); si=sin(n(a)*t);
                Y=Y+2/N*((real(values(a))-.01)*(co*co'+si*si')+ ...
                    imag(values(a))*(co*si'-si*co'));
            end
            scope="Admitancias impresas de Czarnecki1981, (12), modelo finito en 1/2/3 rad/s; no se equiparan a los componentes exactos de Fig.1.";
    end
    e=sqrt(2)*cos(t*n)*U'; Z=Rg*I+omega*Lg*D;
    if k==7, e=interpolar_periodico([1;.2;-.7;-.6;.1],t); end
    u=(I+Z*Y)\e; i=Y*u; m=numel(n);
    % Marco ciclico real de tensiones y cuadraturas, ortonormal RMS.
    E=[]; v=e; A=-D*D;
    for j=1:m, E=[E v D*v]; v=A*v; end %#ok<AGROW>
    [V,~]=qr(E,0); V=V*sqrt(N);
    Dr=V'*D*V/N; Yr=V'*Y*V/N; er=V'*e/N;
    Q=[]; v=D*u;
    for j=1:m, Q=[Q v]; v=A*v; end %#ok<AGROW>
    [Q,~]=qr(Q,0); Q=Q*sqrt(N); qstar=-Q*(Q'*i/N);
    cases{k}=struct('id',ids(k),'N',N,'D',D,'J',J,'tau',t,'omega',omega, ...
        'n',n,'U',U,'e',e,'u',u,'i',i,'Y',Y,'Rg',Rg,'Lg',Lg, ...
        'V',V,'Dr',Dr,'Yr',Yr,'er',er,'qstar',qstar,'scope',scope);
end
end
