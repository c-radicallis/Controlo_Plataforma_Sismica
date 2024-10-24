function F = CorrecaoIntegral_fo(x)

%funcao objetivo do problema de optimizacao para determinacao dos
%parametros da curva para ajuste da media do deslocamento em zero

global i xc1 In

        m1=x(1);
        m2=x(2);
        m3=x(3);
            
       
       xc=xc1(1:end,i-1)+m1*(In(1:end,1)-In(1,1))+m2*(In(1:end,1)-In(1,1)).^2+m3*(In(1:end,1)-In(1,1)).^3;

%funcao objetivo: media^2
F=mean(xc)^2;