%% Optimal threshold for prediction
%   find the threshold with minimal alarm when all SD-event are predicted
%   uses int as the alarm duration
%   uses fpred as precursor
%   uses n12 & n34 as times of SD-events
%=================================================================
ea=0.01;                        %precision of the threshold is 2*ea
tr12=max(fpred)+ea;
tr34=max(fpred)+ea;
uptr=-5;                         % initial value if nevent=0
if ln12>1
   for i=1:ln12
       if n12(i)>int 
          uptr=fix(max(fpred(n12(i)-int:n12(i)-1))./ea).*ea-ea;
          tr12=min([tr12,uptr]);
       else
          uptr=fix(fpred(1)/ea)*ea-ea; 
          tr12=min([tr12,uptr]);
       end
   end
else
        tr12=max(tr12,uptr);     % maximal value to insure errtau=0          %
end
if ln34>1
       for i=1:length(n34)
           if n34(i)>int
              uptr=fix(max(fpred(n34(i)-int:n34(i)-1)./ea)).*ea-ea;
              tr34=min([tr34,uptr]);
           else
              uptr=fix(fpred(1)/ea)*ea-ea; 
              tr34=min([tr34,uptr]);
           end   
       end
else
        tr34=max(tr34,uptr);     % maximal value to insure errtau=0          %       
end
tra=min([tr12,tr34]); %optimal threshold
