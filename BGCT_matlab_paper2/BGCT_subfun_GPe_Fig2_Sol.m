clear;
close all;

len1 = 51;
len2 = 51;
%% ----------------------------------------------------------------
v_se = linspace(1.2,3.2,len2)*1.0e-3; %
delay = linspace(0.02,0.07,len1);
v_ep2 = 5.0e-5;

State = zeros(len1,len2);
FD = zeros(len1,len2);

for i=1:len1
    tic
    for j=1:len2
        [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig2(v_se(i),delay(j),v_ep2);
        State(i,j)=p5;
        FD(i,j)=p6;
    end
    toc
end

save data_test_Fig2

% figure(2)
% subplot(221),imshow(State,[]),set(gca,'ydir','normal'),colorbar;
% subplot(222),imshow(FD,[]),set(gca,'ydir','normal'),colorbar;
% colormap(jet);



