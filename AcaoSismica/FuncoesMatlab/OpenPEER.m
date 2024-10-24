
clear all

input='LomaPrieta.txt';
fid = fopen(input);

%fid2 = fopen('LomaPrieta.txt');

%colocar o ficheiro numa matriz de caracteres
s = '';
while ~feof(fid) 
   line = fgetl(fid);
   if isempty(line), break, end
   s = strvcat(s,line);   
end
%disp(s)


%procedimento para encontrar a linha e coluna com informacao do nº de
%pontos e discretizacao

ln=0;       %linha
tl=[];      %cluna
while isempty(tl)
    ln=ln+1;
    tl=findstr('DT', s(ln,:));
end

%nº de pontos
posn1 = findstr('=', s(ln,1:tl));
posn2 = findstr(',', s(ln,1:tl));
nps=s(ln,posn1+1:posn2-1);
np= str2double(nps);   

%discretizacao: intervalo de tempo entre amostras
posdt1 = findstr('=', s(ln,tl:end));
posdt2 = findstr('SEC', s(ln,:));
dts=s(ln,tl+posdt1:posdt2-1);
dt = str2double(dts);       


num1=1;  %indice do numero

time=0:dt:(np-1)*dt;    %vetor de tempos
ac=zeros(size(time));
for ln2=ln+1:length(s)

%determinacao dos indices em cada linha    
clear bign bign2 bign3 ind inc fin    
bign=isspace(s(ln2,:));
bign2=find(bign>0);
bign3=diff(bign2);
ind=find(bign3>1);
inc=bign2(ind);         %indice de inicio do numero
fin=bign2(ind+1);       %indice de fim do numero

%colocacao dos dados numa so coluna
for i=1:length(ind)
    ac(num1)=str2double(s(ln2,inc(i)+1:fin(i)-1))*9.86;
    num1=num1+1;
end

end
    
%ver resultados    
figure,plot(time,ac),xlabel('t (s)'),ylabel('a (m/s^2)')







