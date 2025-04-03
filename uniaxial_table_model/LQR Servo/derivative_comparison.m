[x_T_LQG, t_out, x] = lsim(clsys, x_ref, t_vector,'foh'); % Simulate the closed-loop response using lsim:
diff_x_T_LQG=diff(x_T_LQG)./diff(t_out);
diff_x_T_LQG=diff([diff_x_T_LQG;0])./diff(t_out);
[ddx_T_LQG, t_out_ddx, x_ddx] = lsim(clsys, ddx_ref ,t_vector,'foh');

%%
figure(2); hold on; grid on; legend; % Plot the response:
plot(t_vector,x_ref,'-.')
plot(t_out, x_T_LQG,'.')
s2_ddx = lsim(1/s^2 , ddx_T_LQG,t_out);
plot(t_out, s2_ddx,'.') % #######################################################
plot(t_out , s2_ddx-x_T_LQG)
% plot(t_out,x(:,1:nx))
% plot(t_out,x(:,nx+1:end-1),'--')
% plot(t_out,x(:,end),':')
xlabel('Time (s)')
ylabel('Output y')
title('Closed-Loop Response with LQG Tracking Controller')

figure(3)
hold on
plot(t_vector,ddx_ref,'-.')
plot(t_out, ddx_T_LQG,'.')
plot(t_out, [diff_x_T_LQG;0],'.')
