%% Mass spring damper - 2 DOF
% Simulation and animation of a mass spring damper system with 2 degrees of
% freedom subjected to base excitation.
%
%%

clear ;  close all ; clc

%% Parameters

% System
M   = 1500;                     % Sprung mass                  [kg]
m   = 150;                      % Unsprung mass                [kg]
Ks  = 48000;                    % Spring constant suspension   [N/m]
Kt  = 200000;                   % Spring constant tire         [N/m]
Cs  = 4000;                     % Damping constant suspension  [N.s/m]

% Animation model
L0_1  = 0.4;                    % Tire height                   [m]
L0_2  = 0.7;                    % Spring relaxed length         [m]
h1  = 0.2;                      % Height of the unsprung mass   [m]
h2  = 0.4;                      % Height of the sprung mass     [m]
a   = 0.8;                      % Width of the block            [m]

% Video
tF      = 15;                   % Final time                    [s]
fR      = 60;                   % Frame rate                    [fps]
dt      = 1/fR;                 % Time resolution               [s]
time    = linspace(0,tF,tF*fR); % Time                          [s]

% Base input
L = 0.05;                       % Amplitude                     [m]
f = 1;                          % Frequency                     [Hz]
w = 2*pi*f;                     % Frequency                     [rad/s]
u_vet = L*cos(w*time');         % Displacement                  [m]

%% Simulation

%  State space model
A = [ 0               1         0       0       ;
      -(Ks+Kt)/m      -Cs/m     Ks/m    Cs/m    ;
      0               0         0       1       ;
      Ks/M            Cs/M      -Ks/M   -Cs/M   ];
B = [ 0     ;
      Kt/m  ;
      0     ;
      0     ];
C = [ 1 0 0 0 ; 
      0 0 1 0 ];
D = [0 ; 0];

sys = ss(A,B,C,D);

% Integration
[y,t,x] = lsim(sys,u_vet,time);

%% Animation

color = cool(6); % Colormap

% Unsprung mass absolute vertical position (lower center point)
y1 = y(:,1) + L0_1; 

% Sprung mass absolute vertical position (lower center point)
y2 = y(:,2) + L0_1 + L0_2; 

% Base info
baseX = [-0.3 0.3];
baseY = [u_vet u_vet];

figure
% set(gcf,'Position',[50 50 1280 720])  % YouTube: 720p
% set(gcf,'Position',[50 50 854 480])   % YouTube: 480p
set(gcf,'Position',[50 50 640 640])     % Social

hold on ; grid on ; axis equal
set(gca,'ylim',[-0.2 1.7],'xlim',[-0.8 0.8],'xtick',-0.8:0.2:0.8,'ytick',-0.2:0.2:1.8)
set(gca,'FontName','Verdana','FontSize',16)

% Create and open video writer object
v = VideoWriter('mass_spring_damper_2_dof_base.mp4','MPEG-4');
v.Quality   = 100;
v.FrameRate = fR;
open(v);

for i=1:length(time)
    cla
    
    % Unsprung mass plot
    fill([-a/2 a/2 a/2 -a/2],[y1(i) y1(i) y1(i)+h1 y1(i)+h1],color(2,:),'LineWidth',2)
    
    % Sprung mass plot
    fill([-a/2 a/2 a/2 -a/2],[y2(i) y2(i) y2(i)+h2 y2(i)+h2],color(6,:),'LineWidth',2)
    
    % Base position instant
    baseYval = baseY(i,:);
    % Base plot
    plot(baseX,baseYval,'k','linewidth',2)
    
    % Spring
    plotSpring(L0_1,L0_2,h1,u_vet,y,i)
    
    % Damper
    plotDamper(L0_1,L0_2,h1,u_vet,y,i)
  
    % Amplitude markers
    % Find steady state amplitude during second half of the simulation:
    % Unsprung mass
    y1Steady = y1(floor(end/2):end); 
    plot([-0.6 0.6],[min(y1Steady) min(y1Steady)],'k--','LineWidth',1.5)
    plot([-0.6 0.6],[max(y1Steady) max(y1Steady)],'k--','LineWidth',1.5)
    % Sprung mass
    y2Steady = y2(floor(end/2):end); 
    plot([-0.6 0.6],[min(y2Steady) min(y2Steady)],'k--','LineWidth',1.5)
    plot([-0.6 0.6],[max(y2Steady) max(y2Steady)],'k--','LineWidth',1.5)

    % Amplitude position input
    plot([-0.6 0.6],[L L],'k--','LineWidth',1.5)
    plot([-0.6 0.6],[-L -L],'k--','LineWidth',1.5)

    title('Mass spring damper - 2 DOF')
    
    xlabel('x [m]')
    ylabel('y [m]')
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);

%% Auxiliary functions

function plotSpring(L0_1,L0_2,h1,yb,y,i)

    % Unsprung mass absolute vertical position (lower center point)
    y1 = y(:,1) + L0_1; 

    % Sprung mass absolute vertical position (lower center point)
    y2 = y(:,2) + L0_1 + L0_2; 

    rodPct      = 0.11;     % Length rod percentage of L0 
    springPct   = 1/3;      % Spring pitch percentage of L0 
    spring_wid  = 3;        % Spring line width
    
    % Spring 1 (Tire) length without rods
    L1 = (y1 - yb) - 2*rodPct*L0_1;
    L2 = (y2 - (y1+h1)) - 2*rodPct*L0_2;
    
    % Spring 1 geometry 
    center1  = 0;         % Lateral position
    wid1     = 0.1;          % Width

    % Spring 2 geometry 
    center2  = -0.2;         % Lateral position
    wid2     = 0.1;          % Width

    % Spring 1
    spring1X = [ 
                center1                                     % Start
                center1                                     % rod
                center1+wid1                                % Part 1   
                center1-wid1                                % Part 2
                center1+wid1                                % Part 3
                center1-wid1                                % Part 4
                center1+wid1                                % Part 5
                center1-wid1                                % Part 6
                center1                                     % Part 7
                center1                                     % rod/End
                ];
    
	spring1Y = [ 
                yb(i)                                       % Start
                yb(i)+rodPct*L0_1                           % rod
                yb(i)+rodPct*L0_1                           % Part 1 
                yb(i)+rodPct*L0_1+springPct*L1(i)           % Part 2
                yb(i)+rodPct*L0_1+springPct*L1(i)           % Part 3
                yb(i)+rodPct*L0_1+2*springPct*L1(i)         % Part 4
                yb(i)+rodPct*L0_1+2*springPct*L1(i)         % Part 5
                yb(i)+rodPct*L0_1+3*springPct*L1(i)         % Part 6
                yb(i)+rodPct*L0_1+3*springPct*L1(i)         % Part 7
                yb(i)+2*rodPct*L0_1+3*springPct*L1(i)       % rod/End
               ]; 

    % Spring 2
    spring2X = [ 
                center2                                     % Start
                center2                                     % rod
                center2+wid2                                % Part 1   
                center2-wid2                                % Part 2
                center2+wid2                                % Part 3
                center2-wid2                                % Part 4
                center2+wid2                                % Part 5
                center2-wid2                                % Part 6
                center2                                     % Part 7
                center2                                     % rod/End
                ];
    
	spring2Y = [ 
                y1(i)+h1                                    % Start
                y1(i)+h1+rodPct*L0_2                        % rod
                y1(i)+h1+rodPct*L0_2                        % Part 1 
                y1(i)+h1+rodPct*L0_2+springPct*L2(i)        % Part 2
                y1(i)+h1+rodPct*L0_2+springPct*L2(i)        % Part 3
                y1(i)+h1+rodPct*L0_2+2*springPct*L2(i)      % Part 4
                y1(i)+h1+rodPct*L0_2+2*springPct*L2(i)      % Part 5
                y1(i)+h1+rodPct*L0_2+3*springPct*L2(i)      % Part 6
                y1(i)+h1+rodPct*L0_2+3*springPct*L2(i)      % Part 7
                y1(i)+h1+2*rodPct*L0_2+3*springPct*L2(i)    % rod/End
               ];

    % PLOT
    plot(spring1X,spring1Y,'k','LineWidth',spring_wid)
    plot(spring2X,spring2Y,'k','LineWidth',spring_wid)
        
end

function plotDamper(L0_1,L0_2,h1,~,y,i)

    % Unsprung mass absolute vertical position (lower center point)
    y1 = y(:,1) + L0_1; 

    % Sprung mass absolute vertical position (lower center point)
    y2 = y(:,2) + L0_1 + L0_2;
    
    rodLowerPct = 0.1;      % Length lower rod percentage of total gap 
    rodUpperPct = 0.4;      % Length upper rod percentage of total gap
    cylinderPct = 0.4;      % Length cylinder porcentagem of total gap
    damper_line_wid  = 3;   % Damper line width
    
    % Damper geometry
    center  = 0.2;              % Lateral position
    wid     = 0.05;             % Width

    % rod attached to unsprung mass
    rod1X = [center center];
    rod1Y = [y1+h1 y1+h1+rodLowerPct*L0_2];
    
    % Damper base cylinder - rod - base 
    cylinderX = [   
                    center-wid
                    center-wid
                    center+wid
                    center+wid
                ];
                
    cylinderY = [
                    y1(i)+h1+rodLowerPct*L0_2+cylinderPct*L0_2
                    y1(i)+h1+rodLowerPct*L0_2 
                    y1(i)+h1+rodLowerPct*L0_2 
                    y1(i)+h1+rodLowerPct*L0_2+cylinderPct*L0_2
                ];
    
    % rod attached to sprung mass
    rod2X = [center center];
    rod2Y = [y2 y2-rodUpperPct*L0_2];
    % Piston inside cylinder
    pistonX = [center-0.8*wid center+0.8*wid];
    pistonY = [y2-rodUpperPct*L0_2 y2-rodUpperPct*L0_2];
    
    % Iteration values
    rod1Yval = rod1Y(i,:);
    rod2Yval = rod2Y(i,:);
    pistonYVal = pistonY(i,:);

    % PLOT
    % rods
    plot(rod1X,rod1Yval,'k','LineWidth',damper_line_wid)
    plot(rod2X,rod2Yval,'k','LineWidth',damper_line_wid)
    % Damper parts
    plot(pistonX,pistonYVal,'k','LineWidth',damper_line_wid)
    plot(cylinderX,cylinderY,'k','LineWidth',damper_line_wid)

end
