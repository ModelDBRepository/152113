clear;
close all;

len1 = 51;
len2 = 51;
%% ----------------------------------------------------------------
% v_p1xi=3e-4; % STN-SNr

v_ep2 = [5:0.5:7.0] * 1e-5;

v_p1xi = logspace(-1,1,len2)*3.0e-4; % 1.6
KK = linspace(0,1.1,len1);

State1 = zeros(len1,len2);
FD1 = zeros(len1,len2);

State2 = zeros(len1,len2);
FD2 = zeros(len1,len2);

State3 = zeros(len1,len2);
FD3 = zeros(len1,len2);

State4 = zeros(len1,len2);
FD4 = zeros(len1,len2);

State5 = zeros(len1,len2);
FD5 = zeros(len1,len2);


for i=1:len1
    tic
    for j=1:len2
        [p11,p21,p31,p41,p51,p61]=BGCT_subfun_GPe_Fig6_Fig7(v_p1xi(i),KK(j),v_ep2(1));  % 5.0e-5
        State1(i,j)=p51;
        FD1(i,j)=p61;
        
        [p12,p22,p32,p42,p52,p62]=BGCT_subfun_GPe_Fig6_Fig7(v_p1xi(i),KK(j),v_ep2(2));  % 5.5e-5
        State2(i,j)=p52;
        FD2(i,j)=p62;
        
        [p13,p23,p33,p43,p53,p63]=BGCT_subfun_GPe_Fig6_Fig7(v_p1xi(i),KK(j),v_ep2(3));  % 6.0e-5
        State3(i,j)=p53;
        FD3(i,j)=p63;
        
        [p14,p24,p34,p44,p54,p64]=BGCT_subfun_GPe_Fig6_Fig7(v_p1xi(i),KK(j),v_ep2(4));  % 6.5e-5
        State4(i,j)=p54;
        FD4(i,j)=p64;
        
        [p15,p25,p35,p45,p55,p65]=BGCT_subfun_GPe_Fig6_Fig7(v_p1xi(i),KK(j),v_ep2(5));  % 7.0e-5
        State5(i,j)=p55;
        FD5(i,j)=p65;
        
    end
    toc
    display(i)
end

save data_test_Fig6_Fig7

