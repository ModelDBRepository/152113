function [pmax1,pmax2,pmin1,pmin2,state,fd]=BGCT_subfun_GPe_Stim_Fig3(v_se,delay,v_p2ext_phi_stim)
% This MATLAB code is a subfunction version of BGCT model described in the
% following paper:
%     Title:    Critical Roles of the Direct GABAergic Pallido-cortical Pathway in Controlling Absence Seizures (2015)
%     Authors: Mingming Chen, Daqing Guo*, Min Li, Tao Ma, Dunsheng Wu, Jingling Ma, Yan Cui, Yang Xia, Peng Xu, and Dezhong Yao*
%     Journal: PLoS Computational Biology
%     Emails:  twqylsf@gmail.com and dqguo@uestc.edu.cn
%
%     Parameters: "open1" and "open2" are gating parameters for the SNr-TRN
%                 pathway and SNr-SRN pathway, with 1 open and 0 close;
%                 "delay" represents GABA_B delay paramter;
%                 "v_sr" is the coupling strength of the TRN-SRN pathway
%                  v_p1xi is the coupling strength of the STN-SNr pathway

rand('state', sum(100*clock)); %random seed

dt=0.00005; % Intergation step
time=25;    % Simulation time

open1=1;    % SNr-TRN pathway, 1 open and 0 close
open2=1;    % SNr-SRN pathway, 1 open and 0 close
KK=1;       % Scale factor

% Maximum firing rate
Qmax_i=250;
Qmax_d1=65;
Qmax_d2=65;
Qmax_p1=250;
Qmax_p2=300;
Qmax_xi=500;
Qmax_s=250;
Qmax_r=250;

% Mean firing threshold
theta_i=15e-3;
theta_d1=19e-3;
theta_d2=19e-3;
theta_p1=10e-3;
theta_p2=9e-3;
theta_xi=1e-2;
theta_s=15e-3;
theta_r=15e-3;

% Coupling strength
v_ee=1.0e-3;
v_ei=1.8e-3;
v_re=5.0e-5;
v_rs=5.0e-4;
v_d1e=1.0e-3;
v_d1d1=2.0e-4;
v_d1s=1.0e-4;
v_d2e=7.0e-4;
v_d2d2=3.0e-4;
v_d2s=5e-5;
v_p1d1=1.0e-4;
v_p1p2=3.0e-5;
v_p1xi=3e-4; % STN-SNr
v_xip2=4.0e-5;
v_es=1.8e-3;
% v_se=2.2e-3;
v_sp1=open2*3.5e-5;
v_rp1=KK*open1*3.5e-5;
v_sr=8.0e-4;
%% ----------------------------------------------------------------
v_ep2 = 5.0e-5;
v_xie=0.1e-3;
v_p2xi=4.5e-4;
v_p2p2=0.75e-4;
v_p2d2=3.0e-4;

%% -----------------------------------------------------------------
% Other parameters
gamma_e=100;
alpha=50;
beta=200;
sigma=0.006;
v_sn_phi_n=2.0e-3;
% delay=0.05;

xdelay=zeros(1,ceil(delay/dt)+1);

xxs=zeros(1,time/dt);

%random initial condition
x=0.5*[0*rand;-1500+6000*rand;0.02*rand;-2*3.5*rand;0.04*rand;-0.7+1.5*rand;0.001+0.025*rand;
    -0.4+rand;0.004+0.013*rand;-0.15+0.4*rand;0.0005+0.0035*rand;-0.12+0.22*rand;-0.001+0.0055*rand;
    -0.1+0.2*rand;-0.09+0.1*rand;-4+7*rand;0.025*rand;-0.6+2*rand];

xm=zeros(18,1);


for t=1:time/dt
    
    xdelay(1:end-1)=xdelay(2:end);
    xdelay(end)=x(17);
    
    S_i=Qmax_i/(1+exp(-pi/sqrt(3)*(x(3)-theta_i)/sigma));
    S_d1=Qmax_d1/(1+exp(-pi/sqrt(3)*(x(5)-theta_d1)/sigma));
    S_d2=Qmax_d2/(1+exp(-pi/sqrt(3)*(x(7)-theta_d2)/sigma));
    S_p1=Qmax_p1/(1+exp(-pi/sqrt(3)*(x(9)-theta_p1)/sigma));
    S_p2=Qmax_p2/(1+exp(-pi/sqrt(3)*(x(11)-theta_p2)/sigma));
    S_xi=Qmax_xi/(1+exp(-pi/sqrt(3)*(x(13)-theta_xi)/sigma));
    S_s=Qmax_s/(1+exp(-pi/sqrt(3)*(x(15)-theta_s)/sigma));
    S_r=Qmax_r/(1+exp(-pi/sqrt(3)*(x(17)-theta_r)/sigma));
    S_r_lag=Qmax_r/(1+exp(-pi/sqrt(3)*(xdelay(1)-theta_r)/sigma));
    
    xm=x;
    
    k1=dt*xm(2);
    k2=dt*(xm(2)+k1/2);
    k3=dt*(xm(2)+k2/2);
    k4=dt*(xm(2)+k3);
    x(1)=xm(1)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(gamma_e^2*(-xm(1)+S_i)-2*gamma_e*xm(2));
    k2=dt*(gamma_e^2*(-xm(1)+S_i)-2*gamma_e*(xm(2)+k1/2));
    k3=dt*(gamma_e^2*(-xm(1)+S_i)-2*gamma_e*(xm(2)+k2/2));
    k4=dt*(gamma_e^2*(-xm(1)+S_i)-2*gamma_e*(xm(2)+k3));
    x(2)=xm(2)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(4);
    k2=dt*(xm(4)+k1/2);
    k3=dt*(xm(4)+k2/2);
    k4=dt*(xm(4)+k3);
    x(3)=xm(3)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(3)+v_ee*xm(1)+v_es*S_s-v_ei*S_i-v_ep2*S_p2)-(alpha+beta)*xm(4));
    k2=dt*(alpha*beta*(-xm(3)+v_ee*xm(1)+v_es*S_s-v_ei*S_i-v_ep2*S_p2)-(alpha+beta)*(xm(4)+k1/2));
    k3=dt*(alpha*beta*(-xm(3)+v_ee*xm(1)+v_es*S_s-v_ei*S_i-v_ep2*S_p2)-(alpha+beta)*(xm(4)+k2/2));
    k4=dt*(alpha*beta*(-xm(3)+v_ee*xm(1)+v_es*S_s-v_ei*S_i-v_ep2*S_p2)-(alpha+beta)*(xm(4)+k3));
    x(4)=xm(4)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(6);
    k2=dt*(xm(6)+k1/2);
    k3=dt*(xm(6)+k2/2);
    k4=dt*(xm(6)+k3);
    x(5)=xm(5)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(5)+v_d1e*xm(1)-v_d1d1*S_d1+v_d1s*S_s)-(alpha+beta)*xm(6));
    k2=dt*(alpha*beta*(-xm(5)+v_d1e*xm(1)-v_d1d1*S_d1+v_d1s*S_s)-(alpha+beta)*(xm(6)+k1/2));
    k3=dt*(alpha*beta*(-xm(5)+v_d1e*xm(1)-v_d1d1*S_d1+v_d1s*S_s)-(alpha+beta)*(xm(6)+k2/2));
    k4=dt*(alpha*beta*(-xm(5)+v_d1e*xm(1)-v_d1d1*S_d1+v_d1s*S_s)-(alpha+beta)*(xm(6)+k3));
    x(6)=xm(6)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(8);
    k2=dt*(xm(8)+k1/2);
    k3=dt*(xm(8)+k2/2);
    k4=dt*(xm(8)+k3);
    x(7)=xm(7)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(7)+v_d2e*xm(1)-v_d2d2*S_d2+v_d2s*S_s)-(alpha+beta)*xm(8));
    k2=dt*(alpha*beta*(-xm(7)+v_d2e*xm(1)-v_d2d2*S_d2+v_d2s*S_s)-(alpha+beta)*(xm(8)+k1/2));
    k3=dt*(alpha*beta*(-xm(7)+v_d2e*xm(1)-v_d2d2*S_d2+v_d2s*S_s)-(alpha+beta)*(xm(8)+k2/2));
    k4=dt*(alpha*beta*(-xm(7)+v_d2e*xm(1)-v_d2d2*S_d2+v_d2s*S_s)-(alpha+beta)*(xm(8)+k3));
    x(8)=xm(8)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(10);
    k2=dt*(xm(10)+k1/2);
    k3=dt*(xm(10)+k2/2);
    k4=dt*(xm(10)+k3);
    x(9)=xm(9)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(9)-v_p1d1*S_d1-v_p1p2*S_p2+v_p1xi*S_xi)-(alpha+beta)*xm(10));
    k2=dt*(alpha*beta*(-xm(9)-v_p1d1*S_d1-v_p1p2*S_p2+v_p1xi*S_xi)-(alpha+beta)*(xm(10)+k1/2));
    k3=dt*(alpha*beta*(-xm(9)-v_p1d1*S_d1-v_p1p2*S_p2+v_p1xi*S_xi)-(alpha+beta)*(xm(10)+k2/2));
    k4=dt*(alpha*beta*(-xm(9)-v_p1d1*S_d1-v_p1p2*S_p2+v_p1xi*S_xi)-(alpha+beta)*(xm(10)+k3));
    x(10)=xm(10)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(12);
    k2=dt*(xm(12)+k1/2);
    k3=dt*(xm(12)+k2/2);
    k4=dt*(xm(12)+k3);
    x(11)=xm(11)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(11)-v_p2d2*S_d2-v_p2p2*S_p2+v_p2xi*S_xi + v_p2ext_phi_stim)-(alpha+beta)*xm(12));
    k2=dt*(alpha*beta*(-xm(11)-v_p2d2*S_d2-v_p2p2*S_p2+v_p2xi*S_xi + v_p2ext_phi_stim)-(alpha+beta)*(xm(12)+k1/2));
    k3=dt*(alpha*beta*(-xm(11)-v_p2d2*S_d2-v_p2p2*S_p2+v_p2xi*S_xi + v_p2ext_phi_stim)-(alpha+beta)*(xm(12)+k2/2));
    k4=dt*(alpha*beta*(-xm(11)-v_p2d2*S_d2-v_p2p2*S_p2+v_p2xi*S_xi + v_p2ext_phi_stim)-(alpha+beta)*(xm(12)+k3));
    x(12)=xm(12)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(14);
    k2=dt*(xm(14)+k1/2);
    k3=dt*(xm(14)+k2/2);
    k4=dt*(xm(14)+k3);
    x(13)=xm(13)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(13)+v_xie*xm(1)-v_xip2*S_p2)-(alpha+beta)*xm(14));
    k2=dt*(alpha*beta*(-xm(13)+v_xie*xm(1)-v_xip2*S_p2)-(alpha+beta)*(xm(14)+k1/2));
    k3=dt*(alpha*beta*(-xm(13)+v_xie*xm(1)-v_xip2*S_p2)-(alpha+beta)*(xm(14)+k2/2));
    k4=dt*(alpha*beta*(-xm(13)+v_xie*xm(1)-v_xip2*S_p2)-(alpha+beta)*(xm(14)+k3));
    x(14)=xm(14)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(16);
    k2=dt*(xm(16)+k1/2);
    k3=dt*(xm(16)+k2/2);
    k4=dt*(xm(16)+k3);
    x(15)=xm(15)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(15)-v_sp1*S_p1+v_se*xm(1)-v_sr*S_r-v_sr*S_r_lag+v_sn_phi_n)-(alpha+beta)*xm(16));
    k2=dt*(alpha*beta*(-xm(15)-v_sp1*S_p1+v_se*xm(1)-v_sr*S_r-v_sr*S_r_lag+v_sn_phi_n)-(alpha+beta)*(xm(16)+k1/2));
    k3=dt*(alpha*beta*(-xm(15)-v_sp1*S_p1+v_se*xm(1)-v_sr*S_r-v_sr*S_r_lag+v_sn_phi_n)-(alpha+beta)*(xm(16)+k2/2));
    k4=dt*(alpha*beta*(-xm(15)-v_sp1*S_p1+v_se*xm(1)-v_sr*S_r-v_sr*S_r_lag+v_sn_phi_n)-(alpha+beta)*(xm(16)+k3));
    x(16)=xm(16)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*xm(18);
    k2=dt*(xm(18)+k1/2);
    k3=dt*(xm(18)+k2/2);
    k4=dt*(xm(18)+k3);
    x(17)=xm(17)+(k1+2*k2+2*k3+k4)/6;
    
    k1=dt*(alpha*beta*(-xm(17)+v_rs*S_s+v_re*xm(1)-v_rp1*S_p1)-(alpha+beta)*xm(18));
    k2=dt*(alpha*beta*(-xm(17)+v_rs*S_s+v_re*xm(1)-v_rp1*S_p1)-(alpha+beta)*(xm(18)+k1/2));
    k3=dt*(alpha*beta*(-xm(17)+v_rs*S_s+v_re*xm(1)-v_rp1*S_p1)-(alpha+beta)*(xm(18)+k2/2));
    k4=dt*(alpha*beta*(-xm(17)+v_rs*S_s+v_re*xm(1)-v_rp1*S_p1)-(alpha+beta)*(xm(18)+k3));
    x(18)=xm(18)+(k1+2*k2+2*k3+k4)/6;
    
    xxs(t)=xm(1);
end


%The following codes are used for bifurcation analysis
TD=5; % starting time used for analysis
xcut=xxs(1,TD/dt+1:end);
peaksmax=findpeaks(xcut);
peaksmax=sort(peaksmax);
peaksmin=findpeaks(-xcut);
peaksmin=sort(-peaksmin);

if (length(peaksmax)>=1)
    pmax1=peaksmax(end);
    pmax2=peaksmax(1);
    index=find(pmax1-peaksmax>=0.1);
    if (length(index)>=1)
        pmax2=peaksmax(index(end));
    end
else
    pmax1=xcut(1);
    pmax2=xcut(2);
end


if (length(peaksmin)>=1)
    pmin1=peaksmin(end);
    pmin2=peaksmin(1);
else
    pmin1=xcut(1);
    pmin2=xcut(2);
end


%The following codes are used for state analysis
state=0;
if (abs(pmin1-pmin2)>=1)
    state=1;   % SWD
elseif (abs(pmin1-pmin2)<=1 && abs(pmax1-pmin2)<=1)
    state=2;   % low firing
end

if (pmin1>248) %straution
    state=4;
end



%The following codes are used for frequency analysis
Fs=1/dt;
N_n=length(xxs(1,TD/dt+1:time/dt));

mean_xx=xxs(1,TD/dt+1:time/dt)-mean(xxs(1,TD/dt+1:time/dt));
X1=fft(mean_xx);

f1=(0:length(X1)/2-1)*Fs./(length(X1));
Power=X1.*conj(X1)/N_n;

index=find(Power(1:N_n/2)==max(Power(1:N_n/2)));

fd=f1(index); %computing the donamite frequency

