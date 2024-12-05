function F = Stiffness2DOF_fo(x)

%funcao objetivo do problema de optimizacao para determinacao da rigidez do
%sistema 2DOF com base nas frequencias naturais identificadas

global M fnzId

Kfo=[x(1)+x(2) -x(2);-x(2) x(2)];

%Det dos valores e vetores proprios
[Vs Ds] = eig(Kfo,M);
fns=sortrows(diag(sqrt(Ds)))/2/pi;


%funcao objetivo
F=sum((fns-fnzId(:,1)).^2)*10^20;