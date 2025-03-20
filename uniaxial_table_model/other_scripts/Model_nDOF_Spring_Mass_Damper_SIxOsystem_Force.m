% Model_nDOF_Spring_Mass_Damper_SIxOsystem_Force(n,M,K,C)
% 
%  Modified:
%
%    03 August 2015
%
%  Author:
%
%    Chandramouli Gnansambandham    
%    Technical University of Kaiserslautern
%
%  Reference:
%
%    http://web.itu.edu.tr/~gundes/2dof.pdf    
%
% INPUT
%==================
%           K   [N/m]   stiffness value of the spring connected to each body 
%           M   [kg]    mass of each body
%           C   [Ns/m]  damping coefficent of the damper connected to each body
%           n   [#]     is the no. of masses in the system
%
% OUTPUT
%==================
%           State space matrices A and B
%
% BEISPIEL
%==================
%
%       [A,B]= Model_nDOF_Spring_Mass_Damper_SISOsystem_Force(n,M,C,K);
%
function [A,B]=Model_nDOF_Spring_Mass_Damper_SISOsystem_Force(n,M,C,K)
%%
if n > 1
    % initialize state space matrices
    B= zeros(n*2,1);
    A_K= zeros(n);
    A_C= zeros(n);
    
    %% calculating state space matrix A
    body=1;
    for row= 1:n
        for col= 1:n
            m= M(body);
            if row==1
                if col==body
                    A_K(row,col)= -K(body)/m;
                    A_C(row,col)= -C(body)/m;
                elseif col==body+1
                    A_K(row,col)= K(body)/m;
                    A_C(row,col)= C(body)/m;
                end
            elseif row==n
                if col==body
                    A_K(row,col)= -(K(body)+K(body-1))/m;
                    A_C(row,col)= -(C(body)+C(body-1))/m;
                elseif col==body-1
                    A_K(row,col)= K(body-1)/m;
                    A_C(row,col)= C(body-1)/m;
                end
            else
                if col==body
                    A_K(row,col)= -(K(body)+K(body-1))/m;
                    A_C(row,col)= -(C(body)+C(body-1))/m;
                elseif col== body-1
                    A_K(row,col)= K(body-1)/m;
                    A_C(row,col)= C(body-1)/m;
                elseif col== body+1
                    A_K(row,col)= K(body)/m;
                    A_C(row,col)= C(body)/m;
                end
            end
        end
        body= body+1;
    end
    A= [zeros(n) eye(n); A_K A_C] ;
else
    A=[0 1;
        -K/M -C/M];
end

%% calculating state space matrix B
B(n+1,1)= 1/M(1) ;

end