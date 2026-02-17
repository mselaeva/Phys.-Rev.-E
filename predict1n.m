%% Total prediction error of SD-events in the chain CN
%========================================================================
%Uses times of events n12 & n34 computed in ev1n.m;
%uses num0 is a total length of learning\testing intervals
%uses predictor fpred = rand\asym computed in asym1n.m
%uses int as alarm duration
%uses tresh as an alarm threshold
       alarm=zeros(num0,1); %alarm time
       nerr12=0;            %missed events in ph12=P(1)-P(2)  
       nerr34=0;            %missed events in ph34-P(N)-P(N-1)           
%===========================================================
% Alarms: an alarm is declared for [t,t+int*dt] if fpred(t)>tresh
%===========================================================
                for na=1:num0-int                                     
                    if fpred(na)>tresh% setting an alarm in random test
                     alarm(na+1:na+int)=1;
                    end
                end  
               nalarm=sum(alarm); % total alarm time
               errtau=nalarm/(num0-int); %relative alarm time in days
%===========================================================
% Missed events in ph12 & ph34: event at time t is missed if alarm(t)=0
%===========================================================                
               if ln12>=1
                    for n=1:ln12 %events in ph12 exist
                            if alarm(n12(n))==0&&n12(n)>int+1
                                nerr12=nerr12+1;%missed event 12
                            end                        
                    end
               else %no events 12
                    nerr12=0;
               end
               if ln34>=1    %events 34 exist
                for n=1:ln34
                        if alarm(n34(n))==0&&n34(n)>int+1
                       nerr34=nerr34+1;%missed event 34
                        end                   
                end
               else %no events 34
                   nerr34=0;
               end  
 %====================================================================
 % relative error of missed events
 %====================================================================
            if nevent>0
               errn=(nerr12+nerr34)/nevent; 
            else
                errn=0;
            end
 %====================================================================
 % total error 
 %====================================================================
           err=errn+errtau;
           clear('alarm')