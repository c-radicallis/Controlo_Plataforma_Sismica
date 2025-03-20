% %%  White Noise
% 
% P = 1e-4;  
% noise = sqrt(P)*randn(size(x_ref,1),1);
% ddx_ref_noise = ddx_ref+noise;
% 
% %%
% % displacements in milimeters
% axes(ax3); % Activate the existing axes
% x_T_noise = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref_noise ,t_vector,'foh');
% erro = x_T_noise-x_ref;
% mse = mean(erro.^2);
% plot(t_vector,x_T_noise,"DisplayName","Tuned MSE="+string(mse))
% 
% % 5th plot
% axes(ax5); % Activate the existing axes
% plot(t_vector,erro,"DisplayName","Tuned")
% 
% % Fourth plot
% axes(ax4); % Activate the existing axes
% ddx_T_noise = lsim(G_xT_xref, ddx_ref_noise ,t_vector,'foh');
% erro = ddx_T_noise-ddx_ref;
% mse = mean(erro.^2);
% plot(t_vector,ddx_T_noise,"DisplayName","Tuned MSE="+string(mse))
% 
% % 6th plot
% axes(ax6); % Activate the existing axes
% hold on
% plot(t_vector,erro,"DisplayName","Tuned")
% 
% 
% % Second plot
% axes(ax2); % Activate the existing axes
% hold on
% i_sv = lsim(G_c ,   (x_ref-x_T_noise)*1e-3  ,t_vector,'foh');
% plot(t_vector,i_sv,"DisplayName","Tuned")
% %  plot
% axes(ax7); % Activate the existing axes
% hold on
% F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
% plot(t_vector,F_p_isv*1e-3,"DisplayName","Tuned")
% 
% 
% 
% % Finding Response Spectre for table tuned with noise
% 
% [picos_ddx_table_noise , picos_x_table_noise ] = ResponseSpectrum( t_vector , ddx_T_noise, f_vector , 1 );
% 
% 
% figure(fig8);
% subplot(121)
% hold on
% mse = mean((picos_ddx_table_noise-picos_ddx_ground).^2);
% plot(f_vector, picos_ddx_table_noise(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName', sprintf('Tuned Platform - MSE= %.2e', mse(1)));
% %semilogx(f_vector, picos_ddx_table_noise(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', 'Tuned Platform - Parallel');
% 
% subplot(122)
% hold on
% mse = mean((picos_x_table_noise-picos_x_ground).^2);
% plot(f_vector, picos_x_table_noise(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName',  sprintf('Tuned Platform - MSE= %.2e', mse(1)));
% %semilogx(f_vector, picos_x_table_noise(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', 'Tuned Platform - Parallel');
% 