function r = comparacion_coste_paper
%COMPARACION_COSTE_PAPER Mismo registro y costes: plantillas, FFT y enumeracion.
root=fileparts(fileparts(mfilename('fullpath'))); folder=fullfile(root,'paper','coste');
if ~isfolder(folder), mkdir(folder); end
raw=readtable(fullfile(root,'matlab','resultados','medidas_simulink','registro_nominal.csv'));
period=.02; endTime=raw.t_s(end); start=endTime-period;
base=caso_m3(31); D=base.D; rows=table; maxError=0;
for N=[63 255 1023]
    tau=(0:N-1)'*2*pi/N; i=interp1(raw.t_s,raw.i_medida_A,start+(0:N-1)'*period/N,'linear');
    for M=[2 4 8]
        poles=linspace(1.3,8.7,M); poles=poles+.017; Phi0=zeros(31,M);
        for k=1:M, Phi0(:,k)=(D*D+poles(k)*eye(31))\(D*base.u); end
        Phi=interpolar_periodico(Phi0,tau); states=double(dec2bin(0:2^M-1,M)-'0')';
        G=Phi'*Phi/N; stateEnergy=sum(states.*(G*states),1);
        F=fft(Phi); active=find(max(abs(F),[],2)>1e-10*max(abs(F),[],'all'));
        H=F(active,:); seeds=fft(i);
        kernels={@()temporal(i,Phi,states,stateEnergy), ...
            @()fourier(i,F,1:N,states,stateEnergy), ...
            @()fourier(i,H,active,states,stateEnergy), ...
            @()enumerate(i,Phi,states)};
        names=["Temporal plantillas","FFT completa","FFT soporte conocido","Enumeracion directa"];
        reference=kernels{1}();
        for j=1:4
            actual=kernels{j}(); error=max(abs(actual-reference)); maxError=max(maxError,error);
            assert(error<1e-8*max(reference));
            seconds=timeit(kernels{j});
            if j==1
                bytes=8*(numel(Phi)+numel(G)+numel(states)+numel(stateEnergy)+N);
            elseif j==2
                bytes=16*(numel(F)+numel(seeds))+8*(numel(states)+numel(stateEnergy));
            elseif j==3
                bytes=16*(numel(H)+numel(seeds))+8*(numel(states)+numel(stateEnergy));
            else
                bytes=8*(numel(Phi)+numel(states)+N+N*size(states,2));
            end
            row=table(N,M,2^M,names(j),seconds,bytes,error, ...
                'VariableNames',{'N','M','Estados','Metodo','Tiempo_s','Memoria_estimacion_bytes','Error_coste_A2'});
            rows=[rows;row]; %#ok<AGROW>
        end
    end
end
[~,cpu]=system('sysctl -n machdep.cpu.brand_string');
r=struct('passed',true,'tabla',rows,'CPU',strtrim(cpu),'MATLAB',version,'Error_max_A2',maxError, ...
    'alcance',"Timeit en el mismo host y registro remuestreado. Plantillas/G/FFT de plantillas preparados antes de medir. FFT soporte conocido usa exactamente el soporte de las plantillas, no estima frecuencias. Memoria de arrays, sin overhead del runtime ni temporales de biblioteca.");
writetable(rows,fullfile(folder,'comparacion.csv')); save(fullfile(folder,'coste.mat'),'r');
fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(rmfield(r,'tabla'),PrettyPrint=true)); fclose(fid); disp(rows);
end
function value=temporal(i,Phi,B,e)
value=mean(i.^2)+2*(i'*Phi/numel(i))*B+e;
end
function value=fourier(i,H,active,B,e)
N=numel(i); f=fft(i); v=real(f(active)'*H)/N^2;
value=sum(abs(f).^2)/N^2+2*v*B+e;
end
function value=enumerate(i,Phi,B)
value=mean((i+Phi*B).^2,1);
end
