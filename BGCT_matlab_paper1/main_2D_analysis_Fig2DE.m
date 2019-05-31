clear;
close all;

len1=31;
len2=31;

delay=linspace(0.02,0.07,len1);
v_sr=linspace(0.4,1.76,len2)*1.0e-3;
v_p1xi=3.0e-4;
open1=1;
open2=1;
KK=1;



State=zeros(len1,len2);
FD=zeros(len1,len2);

for i=1:len1
    tic
    parfor j=1:len2
        [p1,p2,p3,p4,p5,p6]=BGCT_subfun2(delay(i),v_sr(j),v_p1xi,open1,open2,KK);
        State(i,j)=p5;
        FD(i,j)=p6;
    end
    toc
end

save datatFig2DE

figure(1)
subplot(221),imshow(State,[]),set(gca,'ydir','normal'),colorbar;
subplot(222),imshow(FD,[]),set(gca,'ydir','normal'),colorbar;
colormap(jet);