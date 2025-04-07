function [ picos_ddx_m , picos_x_m ] = ResponseSpectrum( t_vector , accel , f_vector , do_displacement )

    m=1;%1kg
    zeta=0.05; %damping ratio (%)s
    
    picos_ddx_m = zeros( length(f_vector) ,1 );
    picos_x_m = picos_ddx_m;
    
    for i=1:length(f_vector)    
      
        % from the desired natural frequency, determining stiffness and damping  
        k = m*(2*pi*f_vector(i))^2; %N/m
        c = zeta*2*m*2*pi*f_vector(i); %N/m/s
          
        % simulacao do modelo 
        s=tf('s');
        ddx_m = lsim( (c*s+k)/(m*s^2+c*s+k) , accel , t_vector(1:size(accel,1)) ,'zoh'); % Funçao de tranferencia de mola massa amortecedor
        picos_ddx_m(i)=max(abs( ddx_m(:,1) ));  
        
        if do_displacement==1
            x_m = lsim( (c*s+k)/(m*s^2+c*s+k)*1/s^2 , accel  , t_vector(1:size(accel,1)) ,'zoh');
            picos_x_m(i)=max(abs( x_m(:,1) ));  
        end

    end
end
