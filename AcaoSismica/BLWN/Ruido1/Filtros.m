%Ruido Filtrado

clear 
s=tf('s');

%filtro 1 - passa-banda entre fnf1 e fnf1
%propriedades
fn1f=0.1;wn1f=2*pi*fn1f; %frequencia do passa-alto
fn2f=3;wn2f=2*pi*fn2f;  %frequncia do passa baixo
zetaf=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp=minreal(s^2/(s^2+2*zetaf*wn1f*s+wn1f^2)/(s^2+2*zetaf*wn2f*s+wn2f^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag,phase] = bode(FTfp,(wn2f+wn1f)/2);
FTf=1/mag*FTfp; %ajuste do ganho

%filtro 2
%propriedades
fn1f2=3;wn1f2=2*pi*fn1f2; %frequencia do passa-alto
fn2f2=30;wn2f2=2*pi*fn2f2;  %frequncia do passa baixo
zetaf2=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp2=minreal(s^2/(s^2+2*zetaf2*wn1f2*s+wn1f2^2)/(s^2+2*zetaf2*wn2f2*s+wn2f2^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag2,phase2] = bode(FTfp2,(wn2f2+wn1f2)/2);
FTf2=1/mag2*FTfp2; %ajuste do ganho



figure,bode(FTf,FTf2),grid

%dados para Simulink
stime=10^-3; %sample time
Intime=0;    %tempo de inicio da simulação
Sttime=20;   %tempo de finalização da simulação


%dados para Simulink
stime=10^-3; %sample time
Intime=0;    %tempo de inicio da simulação
Sttime=20;   %tempo de finalização da simulação
tsub=2;      %tempo de transição entre 0 e o valor maximo das accoes impostas
             %no inicio e fim da serie temporal


%Ficheiro Simulink
SerieRuidoFiltrado


%% projecto de dois filtros de segunda ordem para construcao duma serie 
%temporal nao estacionaria c/3 andamentos

clear 
s=tf('s');

%filtro 1 - passa-banda entre fnf1 e fnf1
%propriedades
fn1f=0.1;wn1f=2*pi*fn1f; %frequencia do passa-alto
fn2f=0.5;wn2f=2*pi*fn2f;  %frequncia do passa baixo
zetaf=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp=minreal(s^2/(s^2+2*zetaf*wn1f*s+wn1f^2)/(s^2+2*zetaf*wn2f*s+wn2f^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag,phase] = bode(FTfp,(wn2f+wn1f)/2);
FTf=1/mag*FTfp; %ajuste do ganho

%filtro 2
%propriedades
fn1f2=0.5;wn1f2=2*pi*fn1f2; %frequencia do passa-alto
fn2f2=1.5;wn2f2=2*pi*fn2f2;  %frequncia do passa baixo
zetaf2=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp2=minreal(s^2/(s^2+2*zetaf2*wn1f2*s+wn1f2^2)/(s^2+2*zetaf2*wn2f2*s+wn2f2^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag2,phase2] = bode(FTfp2,(wn2f2+wn1f2)/2);
FTf2=1/mag2*FTfp2; %ajuste do ganho


%filtro 3
%propriedades
fn1f3=1.5;wn1f3=2*pi*fn1f3; %frequencia do passa-alto
fn2f3=30;wn2f3=2*pi*fn2f3;  %frequncia do passa baixo
zetaf3=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp3=minreal(s^2/(s^2+2*zetaf3*wn1f3*s+wn1f3^2)/(s^2+2*zetaf3*wn2f3*s+wn2f3^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag3,phase3] = bode(FTfp3,(wn2f3+wn1f3)/2);
FTf3=1/mag3*FTfp3; %ajuste do ganho

figure,bode(FTf,FTf2,FTf3),grid

%dados para Simulink
stime=10^-3; %sample time
Intime=0;    %tempo de inicio da simulação
t2=15;       %tempo de para transição andamento 1 para 2
t3=30;       %tempo de para transição andamento 2 para 3
Sttime=40;   %tempo de finalização da simulação

%Ficheiro Simulink
SerieRuidoEntrada2


%% projecto de dois filtros de segunda ordem para construcao duma serie 
%temporal nao estacionaria c/4andamentos

clear 
s=tf('s');

%filtro 1 - passa-banda entre fnf1 e fnf1
%propriedades
fn1f=0.1;wn1f=2*pi*fn1f; %frequencia do passa-alto
fn2f=0.5;wn2f=2*pi*fn2f;  %frequncia do passa baixo
zetaf=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp=minreal(s^2/(s^2+2*zetaf*wn1f*s+wn1f^2)/(s^2+2*zetaf*wn2f*s+wn2f^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag,phase] = bode(FTfp,(wn2f+wn1f)/2);
FTf=1/mag*FTfp; %ajuste do ganho

%filtro 2
%propriedades
fn1f2=0.5;wn1f2=2*pi*fn1f2; %frequencia do passa-alto
fn2f2=1.5;wn2f2=2*pi*fn2f2;  %frequncia do passa baixo
zetaf2=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp2=minreal(s^2/(s^2+2*zetaf2*wn1f2*s+wn1f2^2)/(s^2+2*zetaf2*wn2f2*s+wn2f2^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag2,phase2] = bode(FTfp2,(wn2f2+wn1f2)/2);
FTf2=1/mag2*FTfp2; %ajuste do ganho


%filtro 3
%propriedades
fn1f3=1.5;wn1f3=2*pi*fn1f3; %frequencia do passa-alto
fn2f3=5;wn2f3=2*pi*fn2f3;  %frequncia do passa baixo
zetaf3=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp3=minreal(s^2/(s^2+2*zetaf3*wn1f3*s+wn1f3^2)/(s^2+2*zetaf3*wn2f3*s+wn2f3^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag3,phase3] = bode(FTfp3,(wn2f3+wn1f3)/2);
FTf3=1/mag3*FTfp3; %ajuste do ganho


%filtro 3
%propriedades
fn1f4=5;wn1f4=2*pi*fn1f4; %frequencia do passa-alto
fn2f4=30;wn2f4=2*pi*fn2f4;  %frequncia do passa baixo
zetaf4=1;  %amortecimento: considera-se sem pico na ressancia

%funcao de transferencia
FTfp4=minreal(s^2/(s^2+2*zetaf4*wn1f4*s+wn1f4^2)/(s^2+2*zetaf4*wn2f4*s+wn2f4^2));

%ganho na frequencia intermedia entre fn1 e fn2
[mag4,phase4] = bode(FTfp4,(wn2f4+wn1f4)/2);
FTf4=1/mag4*FTfp4; %ajuste do ganho


figure,bode(FTf,FTf2,FTf3,FTf4),grid

%dados para Simulink
stime=10^-3; %sample time
Intime=0;    %tempo de inicio da simulação
t2=15;       %tempo de para transição andamento 1 para 2
t3=30;       %tempo de para transição andamento 2 para 3
t4=35;       %tempo de para transição andamento 3 para 4
Sttime=40;   %tempo de finalização da simulação
tsub=2;      %tempo de transição entre 0 e o valor maximo das accoes impostas
             %no inicio e fim da serie temporal
%Ficheiro Simulink
SerieRuidoEntrada3