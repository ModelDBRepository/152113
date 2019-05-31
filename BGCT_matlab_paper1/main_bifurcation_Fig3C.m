len1=31;
trails=3;

delay=0.05;
v_sr=1.2e-3;
v_p1xi=linspace(0,0.6,len1)*1.0e-3;
open1=1;
open2=0;
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
    vv_p1xi=v_p1xi(i); 
    pmax1=[];
    pmax2=[];
    pmin1=[];
    pmin2=[];
    state=[];
    fd=[];
    parfor j=1:trails
        [p1,p2,p3,p4,p5,p6]=BGCT_subfun2(delay,v_sr,vv_p1xi,open1,open2,KK);
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

v_p1xi=v_p1xi*1000;
plot(v_p1xi,Pmax1,'.r',v_p1xi,Pmax2,'.b',v_p1xi,Pmin2,'.m',v_p1xi,Pmin1,'.k')