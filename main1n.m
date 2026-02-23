%% Main code: Symmetrical equidistant chain with N oscillators CN
%% and antisymmetrical sine variation of natural frequencies (nf) 
%% in the first and last oscillators (P(1) and P(N))
% uses forcing1n.m to simulate phases ph
% uses ev1n.m to find times of desynchronization events
% uses asym1n.m to compute asymmetry asym of the order parameter
% uses opt1n.n to optimize the threshold of prediction tresh
% uses predict1n.m to find the prediction error
%=================================================================
% Constants
%=================================================================
nos=4;                      % number N of oscillators in a chain
nos2=fix(nos/2);            % M: N=2M or N=2M+1
dt=1/365.25;                % timestep (dayly sampling in model year)
if 2*nos2==nos
    inos=0;                 % N=2M
else
    inos=1;                 % N=2M+1
end
w0=0.58;                        % mean frequency
per=2*pi/w0;                    % mean period
kap1=-0.25;                     % edge coupling k(1) and k(N-1)
kap2=0.25;                   % coupling inside the chain k(i), 1<i<N-1
nperp=25;                   % period of nf variation measured in mean periods
lamda=w0/nperp;             % frequency of nf variation
num2=fix(6*nperp*per/dt);       %convergence time
numint=24;                           %total cycles of simulation 
num=num2+fix(numint*nperp*per/dt);   %total time of simulation in days/steps
minev=numint;             %minimal number of extreme events, if present
%=================================================================
disp([num2,nperp*per,num,minev]) %checking times in years
%=================================================================
% Parameters
%=================================================================
pemin=-1;                                       %minimal p
pemax=1;                                        %maximal p
pp=max([abs(pemin),abs(pemax)]);
apmin=0;                                        %minimal ap
apmax=abs(fix((w0+abs(kap1*pp))*100)/100);  % max(ap) for all p with accuracy 0.01
dpe=0.01;                                       % step for p
dap=0.01;                                       % step for ap
npe=fix((pemax-pemin)/dpe);                     
nap=fix((apmax-apmin)/dap)+1;
int=fix(2*365.25);                            %D
sigma=0.25;                                   
s=10;
% tresh=-1.4;                                   %h
%=================================================================
% Variables
%=================================================================
w=zeros(nos,1);                                 %natural frequencies
upcolz=zeros(npe,nap);                          %zone of events
upcolr=zeros(npe,nap)+0.5;                          %the mean order parameter 
upcol=zeros(npe,nap);                           %the best error n+tau
upcolh=zeros(npe,nap);                          %the best threshold h
upcol4=zeros(npe,nap)-1;                        %intensity of extreme events
%=================================================================
% Simulation and prediction in cycle for (p,ap)
%=================================================================
for ine=1:npe
     p=(pemin+dpe*(ine-1));                 %p
     %=========================================
     % Natural frequencies
     %=========================================
     if inos==0 
        dw=2*kap1*p/(nos2-1);%actual dw for even N
        w(1)=w0+(nos2-1)*dw/2;
        w(nos)=w0+(nos2-1)*dw/2;
     else
      dw=(2*nos2+1)*kap1*p/(nos2*nos2);%actual dw for odd N
        w(1)=w0+nos2*nos2*dw/nos;
        w(nos)=w0+nos2*nos2*dw/nos;
     end
     for i=2:nos2
    w(i)=w(i-1)-dw;
    w(nos-i+1)=w(nos-i+2)-dw;
     end
     if inos==1
        w(nos2+1)=w(nos2)-dw;
     end
     %===========================================
     % Concervation of the sign of nf ==>   max(ap)
     %===========================================
     apma=abs(fix(w(1)*100)/100);           % max(ap) for given p with accuracy 0.01
     nap=fix((apma-apmin)/dap)+1;
     disp([p,apma,dw,w(1)])                 %checking signs and values
    for inp=1:nap                           
       ap=(apmin+dap*(inp-1));              %ap
       nhiss=0;                             %number of intermittent SD events
       forcing1n                            %simulation of phases ph
  upcolr(ine,inp)=mean(rr);                   % mean order parameter
  ph12=ph(num2+1:num,1)-ph(num2+1:num,2);       %phase difference P(1)-P(2)
  ph34=ph(num2+1:num,nos)-ph(num2+1:num,nos-1); %phase difference P(N)-P(N-1)
  ph23=ph(num2+1:num,2:nos-2)-ph(num2+1:num,3:nos-1);%phase differences P(i)-P(i+1)
  ph14=ph(num2+1:num,1)-ph(num2+1:num,nos);     %phase difference P(1)-P(N)
  ph14n=ph(num2+1:num,2)-ph(num2+1:num,nos-1); %phase difference P(2)-P(N-1)
        ev1n                                    % extreme events (EE): n12&n34 (SD), n23
  upcol4(ine,inp)=nevent*36525/(length(ph12));  % intensity of all EE per 100 years
      %========================================
% zone: 1-synch, 2 -separated SD, 3 - intermittent SD, 4 - other 
      %========================================
           if ln23>=1                        %non-SD event exists
                upcolz(ine,inp)=4;
           else                              %non-SD event does not exist
                   if ln12>=1||ln34>=1       %SD event exists
                       for nl=1:ln12-1
                           for nm=2:ln34-1
if n34(nm-1)<n12(nl)&&n34(nm)>n12(nl)&&n34(nm)<n12(nl+1)&&n34(nm+1)>n12(nl+1)                         
            nhiss=nhiss+1;               %alternating SD events exists
end
if n34(nm-1)==n12(nl) 
    nhiss=nhiss+1;               %same time SD events exists
end
                           end
                       end
    if nhiss>0&&ln12>=minev+1 %the alternation is not due to solo events
                      upcolz(ine,inp)=3; % events in 12 and 34 together
    else
                      upcolz(ine,inp)=2; % events in 12 or 34 separated
    end                      
                   else
                       upcolz(ine,inp)=1; % synchronization
                   end
           end              
%%            
disp([ap,p,upcolz(ine,inp)])                    %check zone     
disp([nevent,ln12,ln23,ln34,nhiss])             %with number of events  
      %========================================
        % Threshold optimization: precursor = asymmetry \ random
        % Comment section when use fixed threshold
      %======================================== 
%%      
            num0=fix((num-num2)/2);       %Learning interval
  ph12=ph(num2+1:num2+num0,1)-ph(num2+1:num2+num0,2);       %phase difference P(1)-P(2)
  ph34=ph(num2+1:num2+num0,nos)-ph(num2+1:num2+num0,nos-1); %phase difference P(N)-P(N-1)
  ph23=ph(num2+1:num2+num0,2:nos-2)-ph(num2+1:num2+num0,3:nos-1);%other phase differences
            ph14=ph(num2+1:num2+num0,1)-ph(num2+1:num2+num0,nos);
            ph14n=ph(num2+1:num2+num0,2)-ph(num2+1:num2+num0,nos-1);
            alfa=ph14/2-ph14n/2;    %correction for N oscillators
            asym1n                  %asymmetry function     
     fpred=log10(asym);             %asymmetry precursor, comment if random
%    fpred=rand(num0-num2+1,1);      %random precursor, comment if asymmetry
            ev1n                    %number of events in the learning interval
            opt1n                   %optimal threshold h_opt=tra
            tresh=tra;
            upcolh(ine,inp)=tresh;  %optimal threshold h_opt
        %========================================
        % Setting for the test interval: precursor = asymmetry \ random
      %========================================           
  ph12=ph(num0+1:num,1)-ph(num0+1:num,2);       %phase difference P(1)-P(2)
  ph34=ph(num0+1:num,nos)-ph(num0+1:num,nos-1); %phase difference P(N)-P(N-1)
  ph23=ph(num0+1:num,2:nos-2)-ph(num0+1:num,3:nos-1);%other phase differences
            ph14=ph(num0+1:num,1)-ph(num0+1:num,nos);
            ph14n=ph(num0+1:num,2)-ph(num0+1:num,nos-1);
            ev1n                    %number of events in the test interval
%% 
%========================================================================
%   Prediction error \epsilon=err for the test interval
%   comment section when not needed
%========================================================================
            alfa=ph14/2-ph14n/2;    %correction for N oscillators
            num0=length(alfa);      %length of test interval
            asym1n                  %asymmetry function     
            fpred=log10(asym);      %asymmetry precursor, comment if random
%           fpred=rand(num0,1);     %random precursor, comment if asymmetry  
            predict1n               %errors: errn, errtau  
            upcol(ine,inp,1)=min([errn+errtau,1.1]); %max=1.1 for the best image
            disp([tresh,int,errn,errtau]) %check prediction parameters
            disp([p,ap,upcolz(ine,inp),upcolr(ine,inp),upcol4(ine,inp),upcol(ine,inp)])
    end
end
%% Output plots 
%% Zone plot
figure(17)
subplot(2,3,1),imagesc([pemin pemax],[apmin apmax],transpose(upcolz(:,:))),xlabel('p'),ylabel('a_P'),title('N=8')
colorbar
set(gca,'YDir','normal')
subplot(2,3,2),imagesc([pemin pemax],[apmin apmax],transpose(upcolr(:,:))),xlabel('p'),ylabel('a_P'),title('r')
colorbar
set(gca,'YDir','normal')
subplot(2,3,3),imagesc([pemin pemax],[apmin apmax],transpose(upcol4(:,:))),xlabel('p'),ylabel('a_P'),title('I')
colorbar
set(gca,'YDir','normal')
colormap jet
%%
        clear('ph')
        clear('ph12')
        clear('ph23')
        clear('ph34')
        clear('ph14')
        clear('ph14n')
        clear('alfa') 
%% Prediction plot
% figure(25)
subplot(2,3,4),imagesc([pemin pemax],[apmin apmax],transpose(upcol(:,:))),xlabel('p'),ylabel('a_P'),title('\epsilon')
colorbar
colormap jet
set(gca,'YDir','normal')
subplot(2,3,5),imagesc([pemin pemax],[apmin apmax],transpose(upcolh(:,:))),xlabel('p'),ylabel('a_P'),title('h')
colorbar
colormap jet
set(gca,'YDir','normal')
%%
% subplot(1,3,3),imagesc([pemin pemax],[apmin apmax],transpose(upcol(:,:))),xlabel('p'),ylabel('a_P'),title('h=-1.4,D=2y,\sigma=0.25')
% colorbar
% set(gca,'YDir','normal')
