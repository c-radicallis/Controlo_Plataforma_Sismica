function F = RayleighDamp2DOF_fo(x)

%funcao objetivo do problema de optimizacao para determinacao da
%matriz de amortecimento do sistema 2DOF com base nos fatores de amortecimento
%identificados utilizando o modelo de Rayleigh

global M K fnzId


alpha=x(1);  %fator de amortecimento do isolamento
beta=x(2);  %fator de amortecimento da estrutura


%matriz de amortecimento
C=alpha*M+beta*K;

%matriz da dinamica
A=[zeros(2,2) eye(2);-M^-1*K -M^-1*C];


%det dos polos
[WnEp,ZEp] = damp(A);
fnzEp=sortrows([WnEp([1 3])/2/pi ZEp([1 3])]);


%funcao objetivo: eqm
F=sum((fnzEp(:,2)-fnzId(:,2)).^2)*10^20;