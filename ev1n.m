%% Times of desynchronization events   
% The phase difference ph ---> dph in [-pi,pi] with mean value in 0
% event - is a jump in dph(n+1)-dph(n) due to modulus 2pi
% synchronized solution corresponds to |dph|<pi/2
%=====================================================================
%       Phase differences from the main code: 
% ph12(t) = P(1)-P(2) - phase difference between 1st and 2nd oscillators
% ph34(t)=P(N)-P(N-1) - phase difference between Nth and (N-1)th oscillators 
% ph23(t,i)=P(i)-P(i+1) - phase difference between ith and (i+1)th oscillators 
%                         where i=2,..., N-2
%=====================================================================
    if kap1<0   % negative coupling, synchronization is near pi+2*pi*n
        dph12=mod(ph12,2*pi)-pi;
        dph34=mod(ph34,2*pi)-pi;
    else        % positive coupling, synchronization is near 2*pi*n
        dph12=mod(ph12+pi,2*pi)-pi;
        dph34=mod(ph34+pi,2*pi)-pi;
    end
    if kap2>0   % positive coupling, synchronization is near 2*pi*n     
         dph23=mod(ph23+pi,2*pi)-pi;
    else        % positive coupling, synchronization is near 2*pi*n
         dph23=mod(ph23,2*pi)-pi;
    end
 %====================================================================   
 % Number of desynchronization events 
 % in dph12 & dph34 - solitary decoherence (SD)
 % dph23 - other desynchronization
 %====================================================================
   ln23=0;
   n12=find(abs(diff(dph12))>6); % difference close to 2*pi
   n34=find(abs(diff(dph34))>6); % difference close to 2*pi
   for i=1:nos-3
      n23=find(abs(diff(dph23(:,i)))>6); % difference close to 2*pi
      ln23=length(n23)+ln23; %number of other desynchronization events
   end
      ln12=length(n12);     %number of SD events in dph12
      ln34=length(n34);     %number of SD events in dph34
      nevent=ln12+ln23+ln34; % total number of events
      clear('dph12')
      clear('dph23')
      clear('dph34')