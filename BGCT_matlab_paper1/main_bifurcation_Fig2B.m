len1=31;
trails=3;

delay=linspace(0.02,0.07,len1);
v_sr=1.33*8.0e-4; %~1.06 mV s
v_p1xi=3.0e-4;
open1=1;
open2=1;
KK=1;

matlabpool local 3


Pmax1=zeros(len1,trails);
Pmax2=zeros(len1,trails);
Pmin1=zeros(len1,trails);
Pmin2=zeros(len1,trails);
FR=zeros(len1,trails);
State=zeros(len1,trails);
FD=zeros(len1,trails);

for i=1:len1
    tic
    ddelay=delay(i); 
    pmax1=[];
    pmax2=[];
    pmin1=[];
    pmin2=[];
    state=[];
    fd=[];
    parfor j=1:trails
        [p1,p2,p3,p4,p5,p6]=BGCT_subfun2(ddelay,v_sr,v_p1xi,open1,open2,KK);
        pmax1=[pmax1,p1];
        pmax2=[pmax2,p2];
        pmin1=[pmin1,p3];
        pmin2=[pmin2,p4];
        state=[state,p5];
        fd=[fd,p6];
    end
    
    Pmax1(i,:)=pmax1;
    Pmax2(i,:)=pmax2;
    Pmin1(i,:)=pmin1;
    Pmin2(i,:)=pmin2;
    State(i,:)=state;
    FD(i,:)=fd;
    toc
end

matlabpool close

delay=delay*1000;
plot(delay,Pmax1,'.r',delay,Pmax2,'.b',delay,Pmin2,'.m',delay,Pmin1,'.k')