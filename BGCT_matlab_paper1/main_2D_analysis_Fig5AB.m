clear;
close all;

len1=41;
len2=31;

delay=0.05;
v_sr=1.08e-3;
v_p1xi=logspace(-1,1,len1)*3.0e-4;
open1=1;
open2=1;
KK=linspace(0,1,len2);



State=zeros(len1,len2);
FD=zeros(len1,len2);

for i=1:len1
    tic
    parfor j=1:len2
        [p1,p2,p3,p4,p5,p6]=BGCT_subfun2(delay,v_sr,v_p1xi(i),open1,open2,KK(j));
        State(i,j)=p5;
        FD(i,j)=p6;
    end
    toc
end

save datatFig5AB_new

figure(1)
subplot(221),imshow(State,[0,4]),set(gca,'ydir','normal'),colorbar;
subplot(222),imshow(FD,[]),set(gca,'ydir','normal'),colorbar;
colormap(jet);