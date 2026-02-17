%% Symmetrical equidistant chain with N Kuramoto oscillators 
%% and antisymmetrical variation of natural frequencies 
%% of the two edge oscillatiors
% May be equally used to simulate forcing applied to individual oscillators
%=======================================================================
% Constants (commented when defined in the main code)
%=======================================================================
% w0=0.58;            % mean frequency
% dt=1/365.25;        % timestep (dayly sampling in model year)
% nos=4;           % number N of oscillators in a chain
% nos2=fix(nos/2);  % M: N=2M or N=2M+1
% kap1=-0.25;       % edge coupling k(1) and k(N-1)
% kap2=0.3;         % coupling inside the chain k(i), 1<i<N-1
% p=-0.99;            % parameter of synchronization
%=======================================================================
% Variables (commented when defined in the main code)
%=======================================================================
kappa=zeros(nos-1,1);   %coupling coefficients
% w=zeros(nos,1);       %constant natural frequencies
% ======================================================================
% Initialization (commented when defined in the main code)
% ======================================================================
kappa(1)=kap1;
kappa(nos-1)=kap1;
for i=2:nos-2
    kappa(i)=kap2;
end
% if nos2*2==nos        %N=2M
%     dw=-2*min(abs(kappa))*p/(nos2-1);
%     w(1)=w0+(nos2-1)*dw/2;
%     w(nos)=w0+(nos2-1)*dw/2;
% else                  %N=2M+1
%     dw=-(2*nos2+1)*abs(kappa(1))*p/(nos2*nos2);
%     w(1)=w0+nos2*nos2*dw/nos;
%     w(nos)=w0+nos2*nos2*dw/nos;
% end
% for i=2:nos2
%     w(i)=w(i-1)-dw;
%     w(nos-i+1)=w(nos-i+2)-dw;
% end
% if 2*nos2<nos
%     w(nos2+1)=w(nos2)-dw;
% end    
disp([dw,sum(w)/nos,w0]); %Checking the mean frequency sum(w)/nos=w0
%=========================================================================
%=          Variation of natural frequencies/Forcing
%=========================================================================
fors=zeros(nos,1); % fors=[-1,1,0] indicator of variable natural frequency
fors(1)=1;                      % fors=0 - no variation, 
fors(nos)=-1;                   % fors=-1/1 - negative\positive variation
% ap=0;%max(abs(kappa))*(1-abs(p))+0.1;    % amplitude of sin variation
% nperp=10;                     % period of sin variation measured in mean periods
% lamda=w0/nperp;               % frequency of sin variation
% ======================================================================
% Parameters\variables of the evolution (commented when defined in the main code)
% ======================================================================
% per=2*pi/w0;
% num2=fix(4*pi*365.25/lamda);      % transition interval
% num=num2+fix(2.1*nperp*per/dt);   % total length of the simulation
t=zeros(num,1);                     % time in model years
ph=zeros(num,nos);                  % phases 
dph=zeros(num,nos);                 % phase derivatives
tb=2000;                            % starting time
t(1)=tb;                          
ph(1,:)=rand(nos,1)*pi/2;           % random initial phases
omega=w;                            % variable natural frequencies
% ======================================================================
% Evolution
% ======================================================================
for n=1:num-1 
        t(n+1)=t(n)+dt;
            tc=t(n)-tb;                 %current time from tb
            var=ap.*sin(lamda.*tc);     %variation of nf if present
            for i=1:nos
            omega(i)=w(i)+fors(i)*var;  %current natural frequencies 
            end
dph(n+1,1)=omega(1)-w0+kappa(1)*sin(ph(n,2)-ph(n,1));                 %dP(1)
dph(n+1,nos)=omega(nos)-w0+kappa(nos-1)*sin(ph(n,nos-1)-ph(n,nos));   %dP(N)
            for i=2:nos-1
dph(n+1,i)=omega(i)-w0+kappa(i-1)*sin(ph(n,i-1)-ph(n,i))+kappa(i)*sin(ph(n,i+1)-ph(n,i));   %dP(i)
            end
        for i=1:nos
            ph(n+1,i)=ph(n,i)+dph(n+1,i)*dt;% Phases
        end 
end
% ======================================================================
% The order parameter (commented when independently used)
% ======================================================================
rr=sum(kappa.*transp(cos(ph(num2:num,1:nos-1)-ph(num2:num,2:nos))))/sum(abs(kappa)); % the order parameter
% ======================================================================
% Output plot parameters (commented when used from the main code)
% ======================================================================
% tmin=t(num2);
% nt1=num2;
% tmax=tmin+2*per*nperp;
% nt2=find(t>tmax,1,'first');
% base=fix((ph(nt1,1)-ph(nt1,nos))/(2*pi));         %starting point mod 2pi
% ph1n=ph(nt1:nt2,1)-ph(nt1:nt2,nos);               %P(1)-P(N)
% ======================================================================
% The order parameter (commented when used from the main code)
% ======================================================================
%   rr=zeros(nt2-nt1+1,1);
%   for i=nt1:nt2                                     %
%     rr(i)=sum(transp(kappa).*cos(ph(i,1:nos-1)-ph(i,2:nos)))/sum(abs(kappa));     
%   end
% ======================================================================
% Output plot (commented when used from the main code)
% ======================================================================
% figure(9)
% subplot(2,2,1),hold on,plot(t(nt1:nt2)-t(nt1),ph1n-2*pi*base,'Linewidth',2),ylabel('\Phi'),xlabel('Time (model yrs)')
% subplot(2,2,2),hold on,plot(t(nt1:nt2)-t(nt1),rr(nt1:nt2),'Linewidth',2),ylabel('r'),xlabel('Time (model yrs)')
% subplot(2,2,3),hold on,plot(t(nt1:nt2)-t(nt1),omega(nt1:nt2,1),'Linewidth',2),ylabel('\omega'),xlabel('Time (model yrs)')
% subplot(2,2,4),hold on,plot(t(nt1:nt2)-t(nt1),ph(nt1:nt2,1)-ph(nt1:nt2,2),'Linewidth',2),ylabel('\Phi_{12}'),xlabel('Time (model yrs)')
% 