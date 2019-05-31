% This MATLAB code used to generate a figure like Fig2C, shown in the
% following paper:
%
%     Title:   Bidirectional Control of Absence Seizures by the Basal Ganglia: A
%              Computational Evidence (2014)
%     Authors: Mingming Chen, Daqing Guo*, Tiebin Wang, Wei Jing, Yang Xia, Peng Xu, 
%              Cheng Luo, Pedro A. Valdes-Sosa1, and Dezhong Yao*
%     Journal: PLoS Computational Biology
%     Emails:  twqylsf@gmail.com and dqguo@uestc.edu.cn



clear;
close all;

delay=0.05;
open1=1;
open2=1;
v_sr1=4.8e-4;
v_sr2=10.0e-4;
v_sr3=14.8e-4;
v_sr4=16.0e-4;

[tt,xx1]=BGCT_subfun(delay,v_sr1,open1,open2);
[tt,xx2]=BGCT_subfun(delay,v_sr2,open1,open2);
[tt,xx3]=BGCT_subfun(delay,v_sr3,open1,open2);
[tt,xx4]=BGCT_subfun(delay,v_sr4,open1,open2);

figure(1),
subplot(221),plot(tt,xx1),set(gca,'YDir','reverse');
subplot(222),plot(tt,xx2),set(gca,'YDir','reverse');
subplot(223),plot(tt,xx3),set(gca,'YDir','reverse');
subplot(224),plot(tt,xx4),set(gca,'YDir','reverse');
