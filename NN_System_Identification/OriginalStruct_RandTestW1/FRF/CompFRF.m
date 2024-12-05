% Visualizar FRFs
%%

clear all 
close all

%ficheiro com FRFs
%gain={'1';'09';'08';'07';'06';'05';'04'};
%gain={'1';'09';'08';'07';'06';'05'};
gain={'1';'09';'08';'07';'06'};

%gama de frequencias
frange=[0.2 5];

%funcao a visualizar
FRFv='A'; %'all'-todos;'AgS2Xg';'XigS2Xg';'XsiS2Xg';'XiXg';'XsXg';'XigAg';'XsiAg';'AiAg';'AsAg';


%if strcmp(FRFv,'none')==0

if strcmp(FRFv,'all'),
FRFn ={'XigS2Xg';'XsiS2Xg';'XiXg';'XsXg';'XigAg';'XsiAg';'AiAg';'AsAg'};
elseif strcmp(FRFv,'D'),
FRFn ={'XigS2Xg';'XsiS2Xg';'XiXg';'XsXg'};
elseif strcmp(FRFv,'A'),
FRFn ={'XigAg';'XsiAg';'AiAg';'AsAg'}; 
else FRFn=cellstr(FRFv); end


for j=1:length(FRFn)
figure('Position',[300 10 600 700]),
for i=1:length(gain)

clear data
data=load(strcat('2DOF_FreqResp_G',gain{i}));

hold all
subplot(311),semilogy(data.FRF.f,abs(eval(['data.FRF.',FRFn{j}])))
xlim(frange),ylabel('Mag (m/m/s^2)'),title(FRFn{j})
hold off

hold all

%ajuste dos picos da fase
if strcmp(FRFn{j},'AiAg')
    dph=60;
else
    dph=40;
end

phase0=angle(eval(['data.FRF.',FRFn{j}]));
for kk=1:length(phase0)-1
    if phase0(kk)>0 && phase0(kk+1)<0
        if phase0(kk+1)+2*pi<(180+dph)*pi/180
           phase0(kk+1)=phase0(kk+1)+2*pi;
        end
    elseif phase0(kk)<0 && phase0(kk+1)>0
        if phase0(kk+1)-2*pi>-(180+dph)*pi/180
            phase0(kk+1)=phase0(kk+1)-2*pi;
        end
    end
end


subplot(312),plot(data.FRF.f,phase0*180/pi),xlim(frange),ylim([-180 180]),ylabel('Phase (º)')
set(gca,'YTick',-180:90:180),set(gca,'YTickLabel',{'-180','-90','0','90','180'})
hold off

hold all
subplot(313),plot(data.FRF.f,eval(['data.CO.',FRFn{j}])),xlim(frange),ylabel('Coherence'),ylim([0 1])
xlabel('f (Hz)')
hold off

end
leg=gain;legend(leg)
end

