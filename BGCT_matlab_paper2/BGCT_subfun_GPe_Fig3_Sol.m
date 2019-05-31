clear;
close all;

len1 = 51;
len2 = 51;

%% ------------------------------------------------------------------
trials = 20;
v_ep2=linspace(0.0,2.0,len1)*1.0e-4;
phi_p2ext_stim=linspace(0.0,1.0,len1)*1.0e-2;
v_se = 2.2e-3;
delay = 0.05;

Pmax1_Vep2=zeros(len1,trials);
Pmax2_Vep2=zeros(len1,trials);
Pmin1_Vep2=zeros(len1,trials);
Pmin2_Vep2=zeros(len1,trials);

Phi_E = zeros(len1,trials);
Phi_P2 = zeros(len1,trials);
Phi_Xi = zeros(len1,trials);
Phi_S = zeros(len1,trials);
Phi_R = zeros(len1,trials);
V_ep2_shade = zeros(len1,trials);

Pmax1_Vep2_stim=zeros(len1,trials);
Pmax2_Vep2_stim=zeros(len1,trials);
Pmin1_Vep2_stim=zeros(len1,trials);
Pmin2_Vep2_stim=zeros(len1,trials);

Phi_E_stim = zeros(len1,trials);
Phi_P2_stim = zeros(len1,trials);
Phi_Xi_stim = zeros(len1,trials);
Phi_S_stim = zeros(len1,trials);
Phi_R_stim = zeros(len1,trials);

V_ep2_shade_stim = zeros(len1,trials);
P5_State = zeros(len1,trials);
P5_State_stim = zeros(len1,trials);
P6_Fre = zeros(len1,trials);
P6_Fre_stim = zeros(len1,trials);

for i=1:len1
    pmax1=[];
    pmax2=[];
    pmin1=[];
    pmin2=[];
    p5_state=[];
    p6_fre=[];
    VV_ep2 = [];
    
    pmax1_stim=[];
    pmax2_stim=[];
    pmin1_stim=[];
    pmin2_stim=[];
    p5_state_stim=[];
    p6_fre_stim=[];
    VV_ep2_stim = [];
    
    for j=1:trials
        tic;
        [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig2(v_se,delay,v_ep2(i));
        pmax1=[pmax1,p1];
        pmax2=[pmax2,p2];
        pmin1=[pmin1,p3];
        pmin2=[pmin2,p4];
        p5_state=[p5_state,p5];
        p6_fre=[p6_fre,p6];
        VV_ep2 = [VV_ep2,v_ep2(i) * (p5 == 1 && p6 >= 2 && p6 <= 4)];
        
        % -------------------------------------------
        [phi_e,phi_p2,phi_xi,phi_s,phi_r] = BGCT_subfun_GPe_MFR_Fig3(v_ep2(i));
        Phi_E(i,j) = phi_e;
        Phi_P2(i,j) = phi_p2;
        Phi_Xi(i,j) = phi_xi;
        Phi_S(i,j) = phi_s;
        Phi_R(i,j) = phi_r;
        
        [p1_stim,p2_stim,p3_stim,p4_stim,p5_stim,p6_stim]=BGCT_subfun_GPe_Stim_Fig3(v_se,delay,phi_p2ext_stim(i));
        pmax1_stim=[pmax1_stim,p1_stim];
        pmax2_stim=[pmax2_stim,p2_stim];
        pmin1_stim=[pmin1_stim,p3_stim];
        pmin2_stim=[pmin2_stim,p4_stim];
        p5_state_stim=[p5_state_stim,p5_stim];
        p6_fre_stim=[p6_fre_stim,p6_stim];
        VV_ep2_stim = [VV_ep2_stim,phi_p2ext_stim(i) * (p5_stim == 1 && p6_stim >= 2 && p6_stim <= 4)];
        
        % -------------------------------------------
        [phi_e_stim,phi_p2_stim,phi_xi_stim,phi_s_stim,phi_r_stim] = BGCT_subfun_GPe_MFR_Stim_Fig3(phi_p2ext_stim(i));
        Phi_E_stim(i,j) = phi_e_stim;
        Phi_P2_stim(i,j) = phi_p2_stim;
        Phi_Xi_stim(i,j) = phi_xi_stim;
        Phi_S_stim(i,j) = phi_s_stim;
        Phi_R_stim(i,j) = phi_r_stim;
        toc;
    end
    
    Pmax1_Vep2(i,:)=pmax1;
    Pmax2_Vep2(i,:)=pmax2;
    Pmin1_Vep2(i,:)=pmin1;
    Pmin2_Vep2(i,:)=pmin2;
    P5_State(i,:)=p5_state;
    P6_Fre(i,:)=p6_fre;
    V_ep2_shade(i,:) = VV_ep2;
    
    Pmax1_Vep2_stim(i,:)=pmax1_stim;
    Pmax2_Vep2_stim(i,:)=pmax2_stim;
    Pmin1_Vep2_stim(i,:)=pmin1_stim;
    Pmin2_Vep2_stim(i,:)=pmin2_stim;
    P5_State_stim(i,:)=p5_state_stim;
    P6_Fre_stim(i,:)=p6_fre_stim;
    V_ep2_shade_stim(i,:) = VV_ep2_stim;
end

save data_test_Fig3

