% v_xie=0.1e-3;  v_p2xi=4.5e-4;  v_p2p2=0.75e-4;  v_p2d2=3.0e-4;  v_xip2=4.0e-5;

clear; close all;
len1 = 51;
len2 = 51;
v_ep2=linspace(0.0,2.0,len2)*1.0e-4;

for ii = 1 : 6
    if ii==1
        State_Vxie=zeros(len1,len2);   FD_Vxie=zeros(len1,len2);
        v_xie=linspace(0.0,10,len1)*0.1e-3;
        v_p2xi=4.5e-4;  v_p2p2=0.75e-4;  v_p2d2=3.0e-4;  v_xip2=4.0e-5;
        for i=1:len1
            tic
            for j=1:len2
                [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig4_Fig5(v_xie(i),v_p2xi,v_p2p2,v_p2d2,v_xip2,v_ep2(j));
                State_Vxie(i,j)=p5;
                FD_Vxie(i,j)=p6;
            end
            toc
        end
        save data_GPe_Bif_Vxie_Vep2_2015
        clear State_Vxie FD_Vxie v_xie v_p2xi v_p2p2 v_p2d2 v_xip2
        
    elseif ii==2
        State_Vp2xi=zeros(len1,len2);  FD_Vp2xi=zeros(len1,len2);
        v_p2xi=linspace(0.0,2.0,len1)*4.5e-4;
        v_xie=0.1e-3;  v_p2p2=0.75e-4;  v_p2d2=3.0e-4;  v_xip2=4.0e-5;
        for i=1:len1
            tic
            for j=1:len2
                [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig4_Fig5(v_xie,v_p2xi(i),v_p2p2,v_p2d2,v_xip2,v_ep2(j));
                State_Vp2xi(i,j)=p5;
                FD_Vp2xi(i,j)=p6;
            end
            toc
        end
        save data_GPe_Bif_Vp2xi_Vep2_2015
        clear State_Vp2xi FD_Vp2xi v_xie v_p2xi v_p2p2 v_p2d2 v_xip2
        
    elseif ii==3
        State_Vp2p2=zeros(len1,len2);  FD_Vp2p2=zeros(len1,len2);
        v_p2p2=linspace(0.0,2,len1)*0.75e-4;
        v_xie=0.1e-3;  v_p2xi=4.5e-4;  v_p2d2=3.0e-4;  v_xip2=4.0e-5;
        for i=1:len1
            tic
            for j=1:len2
                [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig4_Fig5(v_xie,v_p2xi,v_p2p2(i),v_p2d2,v_xip2,v_ep2(j));
                State_Vp2p2(i,j)=p5;
                FD_Vp2p2(i,j)=p6;
            end
            toc
        end
        save data_GPe_Bif_Vp2p2_Vep2_2015
        clear State_Vp2p2 FD_Vp2p2 v_xie v_p2xi v_p2p2 v_p2d2 v_xip2
        
    elseif ii==4
        State_Vp2d2=zeros(len1,len2);  FD_Vp2d2=zeros(len1,len2);
        v_p2d2=linspace(0.0,2,len1)*3.0e-4;
        v_xie=0.1e-3;  v_p2xi=4.5e-4;  v_p2p2=0.75e-4;  v_xip2=4.0e-5;
        for i=1:len1
            tic
            for j=1:len2
                [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig4_Fig5(v_xie,v_p2xi,v_p2p2,v_p2d2(i),v_xip2,v_ep2(j));
                State_Vp2d2(i,j)=p5;
                FD_Vp2d2(i,j)=p6;
            end
            toc
        end
        save data_GPe_Bif_Vp2d2_Vep2_2015
        clear State_Vp2d2 FD_Vp2d2 v_xie v_p2xi v_p2p2 v_p2d2 v_xip2
        
    elseif ii==5
        State_Vxip2=zeros(len1,len2);  FD_Vxip2=zeros(len1,len2);
        v_xip2=linspace(0.0,10,len1)*4.0e-5;
        v_xie=0.1e-3;  v_p2xi=4.5e-4;  v_p2p2=0.75e-4;  v_p2d2=3.0e-4;
        for i=1:len1
            tic
            for j=1:len2
                [p1,p2,p3,p4,p5,p6]=BGCT_subfun_GPe_Fig4_Fig5(v_xie,v_p2xi,v_p2p2,v_p2d2,v_xip2(i),v_ep2(j));
                State_Vxip2(i,j)=p5;
                FD_Vxip2(i,j)=p6;
            end
            toc
        end
        save data_GPe_Bif_Vxip2_Vep2_2015
        clear State_Vxip2 FD_Vxip2 v_xie v_p2xi v_p2p2 v_p2d2 v_xip2
        
    end
end




