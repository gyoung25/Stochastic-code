%Gillepsie simulations for the stochastic Lotka Volterra model

%The idea here is to record only realizations in which X wins, then average
%over only those trials to track the average total population size X+Y and
%proportion X/(X+Y)

%TRADE OFF HERE IS f=g+a=g+lambda/M, beta=alpha-b=alpha-mu/M
%In other words, added growth rate at cost of decreased defense
function stochLV()

    %Params
    g=1.1;
    d=0.1;
    alpha=0.11;
    k=5;
    M=10;
    %z=(g-d)*2*k;
    
    lambda=0;%increase in x's growth rate
    
    rho=0;
    mu=rho*lambda; %increase is y's death due to competition with x
    mu=-0.1;
    N=40000; %Number of (discrete) time steps
    
    
    %initial conditions
    x=zeros(1,N);
    y=zeros(1,N);
    tr=zeros(1,N); %"true time"
    
    x(1)=80;
    y(1)=20;
    tr(1)=0;
    
    %number of simulations over which to average
    numSims=10000;
    
    %averages to be computed over the simulations in which X wins
    xWins=0;%count the number of times X wins
    aveXx=zeros(1,N);
    aveYx=zeros(1,N);
    aveTrx=zeros(1,N);
    aveXWinTr=0;
    
    %averages to be computed over the simulations in which Y wins
    yWins=0;%count the number of times Y wins
    aveXy=zeros(1,N);
    aveYy=zeros(1,N);
    aveTry=zeros(1,N);
    aveYWinTr=0;
    
    
    for j=1:numSims
        for i=2:N

            %compute cumulative sums of rates
            %x plus
            a1=(g+lambda/M)*x(i-1);

            %x minus
            a2=a1+d*x(i-1);
            a3=a2+x(i-1)*(x(i-1)-1)/(2*k*M);
            a4=a3+alpha*x(i-1)*y(i-1)/M;

            %y plus
            a5=a4+g*y(i-1);

            %y minus
            a6=a5+d*y(i-1);
            a7=a6+y(i-1)*(y(i-1)-1)/(2*k*M);
            a8=a7+(alpha-mu/M)*x(i-1)*y(i-1)/M;

            %choose random number
            q=rand(1)*a8;
            %choose reaction to take place
            z1=(q<a1);
            z2=(q>=a1)*(q<a2);
            z3=(q>=a2)*(q<a3);
            z4=(q>=a3)*(q<a4);
            z5=(q>=a4)*(q<a5);
            z6=(q>=a5)*(q<a6);
            z7=(q>=a6)*(q<a7);
            z8=(q>=a7)*(q<=a8);

            tr(i)=tr(i-1)-log(rand(1))/a8; %time of next reaction
            x(i)=max(0,x(i-1)+z1-z2-z3-z4);
            y(i)=max(0,y(i-1)+z5-z6-z7-z8);
            
            if (x(i)==0)&&(x(i-1)~=0)
                yWinTr=tr(i);
            end
            
            if (y(i)==0)&&(y(i-1)~=0)
                xWinTr=tr(i);
            end
        end
        
        %calculate a running sum of the X and Y arrays in the cases where X
        %wins
        if y(N)==0
            aveXx=aveXx+x;
            aveYx=aveYx+y;
            aveTrx=aveTrx+tr;
            xWins=xWins+1;
            aveXWinTr=aveXWinTr+xWinTr;
            
        end
        
        %calculate a running sum of the X and Y arrays in the cases where Y
        %wins
        if x(N)==0
            aveXy=aveXy+x;
            aveYy=aveYy+y;
            aveTry=aveTry+tr;
            yWins=yWins+1;
            aveYWinTr=aveYWinTr+yWinTr;
            
        end
        
        
        if mod(j,1000)==0
            j
        end
    end
    
    
    %calculate averages when x wins
    aveXx=aveXx/xWins;
    aveYx=aveYx/xWins;
    
    aveZx=aveXx+aveYx;
    avePx=aveXx./aveZx;
    
    aveTrx=aveTrx/xWins;
    
    aveXWinTr=aveXWinTr/xWins
    xWins
    
    %calculate averages when y wins
    aveXy=aveXy/yWins;
    aveYy=aveYy/yWins;
    
    aveZy=aveXy+aveYy;
    avePy=aveXy./aveZy;

    aveTry=aveTry/yWins;
    
    aveYWinTr=aveYWinTr/yWins
    yWins
    
    %Calculate unconditional average
    aveZ=aveZx*xWins/(xWins+yWins)+aveZy*yWins/(xWins+yWins);
    aveP=avePx*xWins/(xWins+yWins)+avePy*yWins/(xWins+yWins);
    aveTr=aveXWinTr*xWins/(xWins+yWins)+aveYWinTr*yWins/(xWins+yWins);
    
    mean(aveTr)
    
%     %%%%%% plot results from x winning %%%%%%%%%%%%%
%     figure
%     hold on
%     plot(aveTrx,aveXx,'linewidth',2);
%     plot(aveTrx,aveYx,'linewidth',2);
%     xlabel('time')
%     ylabel('pop size')
%     axis([0 100 0 150])
%     legend('X','Y')
%     box on
%     set(gca,'fontsize',28)
%     
%     figure
%     hold on
%     plot(aveTrx,aveZx,'linewidth',2);
%     plot([0 100],[(g-d)*2*k*M (g-d)*2*k*M],'k--','linewidth',2)
%     xlabel('time')
%     ylabel('X+Y')
%     axis([0 100 0 150])
%     box on
%     set(gca,'fontsize',28)
%     
%     figure
%     plot(aveTrx,avePx,'linewidth',2);
%     xlabel('time')
%     ylabel('p')
%     axis([0 100 0 1]) 
%     box on
%     set(gca,'fontsize',28)
%     
%     figure
%     hold on
%     plot(aveXx,aveYx,'linewidth',2);
%     plot([0 (g-d)*2*k*M],[(g-d)*2*k*M 0],'k--','linewidth',2)
%     ylabel('Y')
%     xlabel('X')
%     axis([0 100 0 100]) 
%     box on
%     set(gca,'fontsize',28)
%     
%     %%%%%%%%%%%%%%%%%%% plot results from y winning %%%%%%%%%%%%%%%%
%     figure
%     hold on
%     plot(aveTry,aveXy,'linewidth',2);
%     plot(aveTry,aveYy,'linewidth',2);
%     xlabel('time')
%     ylabel('pop size')
%     axis([0 100 0 150])
%     legend('X','Y')
%     box on
%     set(gca,'fontsize',28)
%     
%     figure
%     hold on
%     plot(aveTry,aveZy,'linewidth',2);
%     plot([0 100],[(g-d)*2*k*M (g-d)*2*k*M],'k--','linewidth',2)
%     xlabel('time')
%     ylabel('X+Y')
%     axis([0 100 0 150])
%     box on
%     set(gca,'fontsize',28)
%     
%     figure
%     plot(aveTry,avePy,'linewidth',2);
%     xlabel('time')
%     ylabel('p')
%     axis([0 100 0 1]) 
%     box on
%     set(gca,'fontsize',28)
%     
    figure
    hold on
    plot(aveXy,aveYy,'linewidth',2);
    plot([0 (g-d)*2*k*M],[(g-d)*2*k*M 0],'k--','linewidth',2)
    ylabel('Y')
    xlabel('X')
    axis([0 100 0 100]) 
    box on
    set(gca,'fontsize',28)
    %print(fig,'-depsc','XvsY_Stoch.eps')
%     
%     %Plot unconditioned averages
%     figure
%     hold on
%     plot(aveTr,aveZ,'linewidth',2);
%     %plot([0 (g-d)*2*k*M],[(g-d)*2*k*M 0],'k--','linewidth',2)
%     ylabel('X+Y')
%     xlabel('t')
%     axis([0 100 0 200]) 
%     box on
%     set(gca,'fontsize',28)
%     %print(fig,'-depsc','ZvsT_Stoch.eps')
%     
%     figure
%     hold on
%     plot(aveTr,aveP,'linewidth',2);
%     %plot([0 (g-d)*2*k*M],[(g-d)*2*k*M 0],'k--','linewidth',2)
%     ylabel('p')
%     xlabel('t')
%     axis([0 100 0 1]) 
%     box on
%     set(gca,'fontsize',28)
%     %print(fig,'-depsc','PvsT_Stoch.eps')
%     
%     figure
%     hold on
%     plot(aveTr,aveZ,'linewidth',2);
%     %plot([0 (g-d)*2*k*M],[(g-d)*2*k*M 0],'k--','linewidth',2)
%     ylabel('X+Y')
%     xlabel('t')
%     axis([0 5 0 200]) 
%     box on
%     set(gca,'fontsize',28)
%     %print(fig,'-depsc','ZvsT_Stoch.eps')
%     
%     figure
%     hold on
%     plot(aveTr,aveP,'linewidth',2);
%     %plot([0 (g-d)*2*k*M],[(g-d)*2*k*M 0],'k--','linewidth',2)
%     ylabel('p')
%     xlabel('t')
%     axis([0 5 0 1]) 
%     box on
%     set(gca,'fontsize',28)
%     %print(fig,'-depsc','PvsT_StochClose.eps')
    

    
end