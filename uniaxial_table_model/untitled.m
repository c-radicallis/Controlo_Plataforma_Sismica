%% State space with PL

A = [-1/tau_sv ,  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0;
k_h/A^2  , - k_h*k_pl/A^2 , 0 , 0 , 0 , 0 , -k_h/A , 0 , 0 ;
k_h/A  , - k_h*k_pl/A , 0              , 0 , 0                , 0   , -k_h , 0               , 0 ;
                0 ,                0 , 0              , 0 , 0                , 0   , 1    , 0               , 0 ;
                0 ,                0 , 0              , 0 , 0                , 0   , 0    , 1               , 0 ;
                0 ,                0 , 0              , 0 , 0                , 0   , 0    , 0               , 1 ;
                0 ,               0 , 1/mT , -k1/mT , k1/mT     ,  0              , (-cT - c1) /mT ,  c1 /mT   , 0               ;
                0 ,                0 , 0             , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1         ,(-c1-c2)/m1 , c2/m1 ;
                0 ,                0 , 0             ,         0        , k2/m2     , -k2/m2, 0                       , c2/m2     , -c2/m2];

B=[k_svk_q/tau_sv ; zeros(8,1)];

C=[zeros(1,3), 1 , zeros(1,5)];%xT
     %zeros(1,2), 1 , zeros(1,6)];%Fp

 obsv(A,C)

 %% State space without PL

 A = vpa([-1/tau_sv, 0              , 0      , 0          , 0     , 0              , 0         , 0     ;
            k_h/A , -k_h*k_pl/(A^2), 0      , 0          , 0     , -k_h           , 0         , 0     ;
                0 ,              0 , 0      , 0          , 0     , 1              , 0         , 0     ;
                0 ,              0 , 0      , 0          , 0     , 0              , 1         , 0     ;
                0 ,              0 , 0      , 0          , 0     , 0              , 0         , 1     ;
                0 ,           1/mT , -k1/mT , k1/mT      ,  0    , (-cT - c1) /mT ,  c1 /mT   , 0     ;
                0 ,            0   , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1          ,(-c1-c2)/m1, c2/m1 ;
                0 ,            0   , 0      , k2/m2      , -k2/m2, 0              , c2/m2     , -c2/m2], 500);


B=vpa([k_svk_q/tau_sv ; zeros(7,1)],500);

C=vpa([zeros(1,2), 1 , zeros(1,5)],500);%xT
     %zeros(1,2), 1 , zeros(1,6)];%Fp
%C=ones(1,9);

D=0;

ss_model = ss(double(A),double(B),double(C),D);

obs = vpa(obsv(A,C),500);

r_obsv = rank(obs);

controlability = vpa(ctrb(A,B),500);

r_controlability = rank(controlability);
