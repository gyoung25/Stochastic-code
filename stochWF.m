%Simulations of the Wright-Fisher process with periodic fitness change
function stochWF()

    %Params
    a=5;
    N=100;
    w=1*2*pi;
    phi=pi;
    b=0;
    f1=1;
    
    M=55000; %Number of (discrete) time steps
    
    for k=3:3
    %initial conditions
    x=zeros(1,M);
    
    %uncomment if varying x_0 (k should go 1 through 9)
    x(1)=k*10;
    
    %uncomment if varying phi (k should go 0 through 10)
%     x(1)=33;
%     phi=k*2*pi/10;
    
    %number of simulations over which to average
    numSims=1;
    
    %averages to be computed over the simulations in which X wins
    xWins=0;%count the number of times X wins
    aveXx=zeros(1,M);
    aveXWinT=0;
    
    %averages to be computed over the simulations in which Y wins
    yWins=0;%count the number of times Y wins
    aveXy=zeros(1,M);
    aveYWinT=0;
    
    for j=1:numSims
        for i=2:M

            %compute cumulative sums of rates
            %x plus
            %a1=(N-x(i-1))/N*f1*x(i-1)/(f1*x(i-1)+f2(i-1)*(N-x(i-1)-1));
            a1=(N-x(i-1))/N*f1*x(i-1)...
                /(f1*x(i-1)+(f1+a/N*cos(w/N^2*(i-1)+phi)+b/N)*(N-x(i-1)-1));

            %x minus
            %a2=a1+x(i-1)/N*f2(i-1)*(N-x(i-1))/(f1*(x(i-1)-1)+f2(i-1)*(N-x(i-1)));
            a2=a1+x(i-1)/N*(f1+a/N*cos(w/N^2*(i-1)+phi)+b/N)*(N-x(i-1))...
                /(f1*(x(i-1)-1)+(f1+a/N*cos(w/N^2*(i-1)+phi)+b/N)*(N-x(i-1)));

            %choose random number
            q=rand(1);
            %choose reaction to take place
            z1=(q<a1);
            z2=(q>=a1)*(q<a2);
            
            x(i)=max(0,x(i-1)+z1-z2);
            
            if (x(i)==0)&&(x(i-1)~=0)
                yWinT=i;
            end
            
            if (x(i)==N)&&(x(i-1)~=N)
                xWinT=i;
            end
        end
        
        %calculate a running sum of the X and Y arrays in the cases where X
        %wins
        if x(M)==N
            xWins=xWins+1;
            aveXx=aveXx+x;
            aveXWinT=aveXWinT+xWinT;
        end
        
        %calculate a running sum of the X and Y arrays in the cases where Y
        %wins
        if x(M)==0
            yWins=yWins+1;
            aveXy=aveXy+x;
            aveYWinT=aveYWinT+yWinT;
        end
        
        
        if mod(j,2000)==0
            j
        end
    end
    
    %calculate averages when x wins
    aveXx=aveXx/xWins;
    
    aveXWinT=aveXWinT/xWins;
    xWins
    
    %calculate averages when y wins
    aveXy=aveXy/yWins;
    
    aveYWinT=aveYWinT/yWins;
    yWins
    
    timeToFix=aveXWinT*xWins/(xWins+yWins)+aveYWinT*yWins/(xWins+yWins)
    x(1)
    phi
    k
    xWinProb=xWins/(xWins+yWins)
    end
    aveX=(aveXx*xWins+aveXy*yWins)/(xWins+yWins);
    
    %save('WF_x60_phinegPiOver2_1','xWins','yWins','aveXWinT','aveYWinT','timeToFix');
    
    w
    phi
%     figure
%     hold on
%     plot(aveXy,aveYy,'linewidth',2);
%     plot([0 90],[90 0],'k--','linewidth',2)
%     ylabel('Y')
%     xlabel('X')
%     axis([0 100 0 100]) 
%     box on
%     set(gca,'fontsize',28)
    %%%%%% plot results from x winning %%%%%%%%%%%%%
    figure
    hold on
    %plot([1:length(aveXx)],aveXx,'linewidth',2);
    %plot([1:length(aveXy)],aveXy,'linewidth',2);
    plot([1:length(aveXy)],aveX,'linewidth',2);
    xlabel('time')
    ylabel('pop size')
    axis([0 M -5 N+5])
    %legend('X wins','Y wins','Average')
    box on
    set(gca,'fontsize',28)
    
    figure
    hold on
    plot([1:M],f1+a/N*cos(w/N^2*[1:M]+phi),'linewidth',2)
    plot([0 M],[f1 f1],'k--')
    xlabel('time')
    ylabel('fitness advantage')
    axis([0 M -2*a/N+f1 2*a/N+f1])
    box on
    set(gca,'fontsize',28)
    
    function y=f2(t)
        y=f1+a/N*cos(w/N^2*t+phi)+b/N;
    end
    
end