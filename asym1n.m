%% Asymmetry of the order parameter for the symmetrical chain CN
%% computed through the auxilliary chain S3
%============================================================
% Phase differences from the main code:
%   ph14=P(1)-P(N), ph14n=P(2)-P(N-1)
%   alfa=ph14/2-ph14n/2
%============================================================
% Parameters of the auxilliary chain
%============================================================
w0=0.58;            %mean frequency
w2=0.;              % natural frequency of the middle oscillator 
sigma=0.25;         % absolute coupling disparity
s=10;               % coupling strength
num1=length(alfa);  % total length of the simulation
dt=1/365.25;        % timestep (dayly sampling in model year)
%============================================================
    %Natural frequencies of the auxilliary chain
%============================================================
    ow1=zeros(num1,2);%\omega_1 - from first oscillator
    ow3=zeros(num1,2);%\omega_3 - from last oscillator    
    odw=zeros(num1,2);%$\Delta\omega$ in the auxilliary chain
    ow2=w2;           % natural frequency of the middle oscillator
    rt=sqrt(sigma^2.*sin(alfa).*sin(alfa)+cos(alfa).*cos(alfa)-(w0-w2)*(w0-w2)/(4*s*s));
    odw(:,1)=(sigma*(w0-w2)+s*(1-sigma*sigma).*rt.*sin(2*alfa))./(2.*(sigma^2.*sin(alfa).*sin(alfa)+cos(alfa).*cos(alfa)));
    odw(:,2)=(-sigma*(w0-w2)+s*(1-sigma*sigma).*rt.*sin(2*alfa))./(2.*(sigma^2.*sin(alfa).*sin(alfa)+cos(alfa).*cos(alfa)));   
    ow1(:,1)=3/2*w0-ow2/2+odw(:,1); %positive sigma    
    ow1(:,2)=3/2*w0-ow2/2+odw(:,2); %negative sigma
    ow3(:,1)=3/2*w0-ow2/2-odw(:,1); %positive sigma
    ow3(:,2)=3/2*w0-ow2/2-odw(:,2); %negative sigma
%===========================================================
%   Coupling in the auxilliary chain: 
%   symmetry makes ph14=P(1)-P(N) to be synchronized in the neighborhood of 0
%===========================================================
        k=s; 
        dk=[sigma*s,-sigma*s];
        k12=k+dk;        %[positive sigma,negative sigma]
        k23=k-dk;        %[positive sigma,negative sigma]
%===========================================================
%   Simulation of phases in the auxilliary chain: 
%===========================================================
    odteta=zeros(num1,3,2);     %phase derivatives (t, ith oscillator, positive/negative sigma)
    oteta=zeros(num1+1,3,2);    %phases (t+1, ith oscillator, positive/negative sigma)
    oteta(1,:,:)=0;             %initialization for t=0
    or=zeros(num1,2);           %order parameter (t, positive/negative sigma)
    for n=1:num1
        for kkk_fl=1:2
        odteta(n,1,kkk_fl)=(ow1(n,kkk_fl)+k12(kkk_fl)*sin(oteta(n,2,kkk_fl)-oteta(n,1,kkk_fl)))*dt;
        odteta(n,2,kkk_fl)=(ow2+k12(kkk_fl)*sin(oteta(n,1,kkk_fl)-oteta(n,2,kkk_fl))+k23(kkk_fl)*sin(oteta(n,3,kkk_fl)-oteta(n,2,kkk_fl))).*dt;
        odteta(n,3,kkk_fl)=(ow3(n,kkk_fl)+k23(kkk_fl)*sin(oteta(n,2,kkk_fl)-oteta(n,3,kkk_fl)))*dt;  
        end
        for i=1:3
            for kkk_fl=1:2
        oteta(n+1,i,kkk_fl)=oteta(n,i,kkk_fl)+odteta(n,i,kkk_fl);
            end
        end
%===========================================================
%   Order parameter in the auxilliary chain: 
%===========================================================
        for kkk_fl=1:2
        or(n,kkk_fl)=(k12(kkk_fl)*cos(oteta(n+1,2,kkk_fl)-oteta(n+1,1,kkk_fl))+k23(kkk_fl)*cos(oteta(n+1,2,kkk_fl)-oteta(n+1,3,kkk_fl)))/(abs(k12(kkk_fl))+abs(k23(kkk_fl)));
        end
    end
%===========================================================
%   Asymmetry of the order parameter
%===========================================================
 asym=abs(or(:,1)-or(:,2));
%===========================================================
%  Phase difference to check self-consistence during synchronization
%=========================================================== 
 o1alfa=(oteta(2:num1+1,1,1)-oteta(2:num1+1,3,1))/2;%sigma>0
 o2alfa=(oteta(2:num1+1,1,2)-oteta(2:num1+1,3,2))/2;%sigma<0
%===========================================================
% Output checking self-consistence (commented after a successful check)
%===========================================================
%  figure(503)
%plot(t(num2:num),mod(2*alfa+pi,2*pi)-pi,t(num2:num),2*o1alfa,t(num2:num),2*o2alfa,'Linewidth',2),xlim([11000,13000]),ylabel('\theta_{14}')
%===========================================================
% Output asymmetry (commented when used in prediction)
%===========================================================
%  figure(5)
% semilogy(t(num2:nmax-1),asym(num2:nmax-1),'Linewidth',2),xlabel('Time (model yrs)'),ylabel('A')
 