% determines the membrane potentials of three coupled MS-oscillators
%
%   osc. 1  ---<-->---  osc. 2 ---<-->---  osc. 3
%               tau1               tau3
%
% based on PRC with equal intrinsic period, T=1.
% Between spike arrivals the membrane potential is given by
%         f = log(1+(exp(b)-1)*phi)/b
% that resets is f==1 -> 0.
% As soon as a pulse arrives the membrane potential,
% f increases instantly by epsilon .
% Time is expressed relative to the intrinsic period
% of the MS-oscillators.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Atthaphon Viriyopase                                      %
% Date:  7 July, 2009                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val=EventBased_3MS_osc(p_phi1, p_phi2, p_tau1, p_tau3, p_epsi12, p_epsi21, p_epsi32, p_epsi23, p_tEnd, p_PRC_option)
global i_e_List
global e_List_From e_List_To e_List_Int
global isDeBug

%%%%%% Parameters %%%%%%
deBugLevel=isDeBug;
b1=3;b2=3;b3=3;b=3;
tau1=p_tau1;tau3=p_tau3;
epsi12=p_epsi12;epsi21=p_epsi21;
epsi32=p_epsi32;epsi23=p_epsi23;

tend=p_tEnd;

init_phi1=p_phi1;init_phi2=p_phi2;

init_phis=[init_phi1; init_phi2; 0];

%%%%%% Initialize variables %%%%%%
% Record phases of phi1 and phi2 when phi3 spikes %
i_phisEvo=1;
phi1Evo(i_phisEvo,1)=init_phi1;
% phi2Evo(i_phisEvo,1)=init_phi2;
phi3Evo(i_phisEvo,1)=0;

curr_phis=init_phis;

% track time when an event happens %
i_t_event=1;
t_event(i_t_event,1)=0; 

i_e_List=0;
e_List_From=[];e_List_To=[];e_List_Int=[];
addEvent(1, 1, 1-init_phis(1,1));
addEvent(2, 2, 1-init_phis(2,1));
addEvent(3, 3, 1-init_phis(3,1));

% track time int. that spikes arrive at each neurons %
if isLEq(init_phis(1,1),tau1) 
    if ~(isEq(init_phis(2,1),0) && isEq(init_phis(1,1),tau1))
        addEvent(1,2,tau1-init_phis(1,1));
    end
end
if isLEq(init_phis(2,1),tau1) 
    if ~(isEq(init_phis(1,1),0) && isEq(init_phis(2,1),tau1))
        addEvent(2,1,tau1-init_phis(2,1));         
    end
end    
if isL(init_phis(2,1),tau3) addEvent(2,3,tau3-init_phis(2,1)); end

% always add this event because phi3 always begins with 0 and 0<tau3<0.5 %
addEvent(3,2,tau3-init_phis(3,1));

if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': Init: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end

while isLEq(t_event(i_t_event,1),tend)
    earliestSet=getEarliestESet(e_List_From, e_List_To, e_List_Int);
    
    % Move the data from e_List to tmp_e_List %
    tmp_e_List_From=e_List_From(earliestSet,1);tmp_e_List_To=e_List_To(earliestSet,1);tmp_e_List_Int=e_List_Int(earliestSet,1);
    
    e_List_From(earliestSet,:)=[];e_List_To(earliestSet,:)=[];e_List_Int(earliestSet,:)=[];    
    i_e_List=i_e_List-size(earliestSet,1);
    
    isN1SpNormally=0;isN2SpNormally=0;isN3SpNormally=0;
    isN2SpDriven=0;isN3SpDriven=0;

    % Record the time for this event %
    i_t_event=i_t_event+1;
    t_event(i_t_event,1)=t_event(i_t_event-1,1)+tmp_e_List_Int(1,1);
    
    % Inc. the phases of neurons up to just right before the event happen %
    curr_phis=curr_phis+tmp_e_List_Int(1,1);
    
    % Update the interval to be relatively match with the current time %
    e_List_Int=e_List_Int-tmp_e_List_Int(1,1);
    
    % Process after the events happen %
    i=1;
    while (size(tmp_e_List_From,1)>=1)
        
        % Neurons spike event %
        if isEq(tmp_e_List_From(i,1),tmp_e_List_To(i,1))
            
            % 1st neuron spikes %
            if isEq(tmp_e_List_From(i,1),1)
                isN1SpNormally=1;
                
                % Update 1st neuron spikes interval & phase %
                curr_phis(1,1)=0;
                
                % Add spike happens normally list %
                addEvent(1,1,1);                
                % Add spike arr. list %
                addEvent(1,2,tau1);                
                
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N1 sp. normally: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                
            % 2nd neuron spikes %
            elseif isEq(tmp_e_List_From(i,1),2)
                isN2SpNormally=1;
                
                % Update 2nd neuron spikes interval & phase %
                curr_phis(2,1)=0;
                
                % Add spike happens normally list %
                addEvent(2,2,1);
                % Add spike arr. list %
                addEvent(2,1,tau1);
                addEvent(2,3,tau3);

if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 sp. normally: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                
            % 3rd neuron spikes %
            else
                isN3SpNormally=1;
                
                % Update 3rd neuron spikes interval & phase %
                curr_phis(3,1)=0;
                
                % Add spike happens normally list %
                addEvent(3,3,1);                
                % Add spike arr. list %
                addEvent(3,2,tau3);
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N3 sp. normally: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                
            end
            % Delete used event in tmp_e_List %
            tmp_e_List_From(1,:)=[];tmp_e_List_To(1,:)=[];tmp_e_List_Int(1,:)=[];
            
        % Spike arriving event %
        else
            % The spike arr. at 1st neuron %
            if isEq(tmp_e_List_To(i,1),1)
                
                % The spike arr. not the same time as 1st neuron spike %
                if (isN1SpNormally==0)
                    tmpPRC=PRC(curr_phis(1,1),epsi21,b1,p_PRC_option);
                    
                    % 1st neuron suddenly spikes? %
                    if (isEq(curr_phis(1,1),1-tmpPRC)==1)
                        % Update 1st neuron spikes interval & phase %
                        curr_phis(1,1)=0;
                        
                        % Add spike happens normally list %
                        addEvent(1,1,1);                        
                        % Add spike arr. list %
                        addEvent(1,2,tau1); 
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N1 rec. sp. from ',num2str(tmp_e_List_From(i,1)),' and sp.: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                        
                    else
                        % Update 1st neuron spikes interval & phase %
                        curr_phis(1,1)=curr_phis(1,1)+tmpPRC;
                        
                        % Update the natual-spike 1st neuron in the list %
                        i_1st_neuron=findEvent(e_List_From,e_List_To,1,1);
                        cnt=size(i_1st_neuron,1);
                        
                        % There exist something to be updated %
                        if (cnt>0)
                            e_List_Int(i_1st_neuron)=1-curr_phis(1,1);
                        end        
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N1 rec. sp. from ',num2str(tmp_e_List_From(i,1)),': ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                                                
                    end                                        
                end                
                % Delete used event in tmp_e_List %
                tmp_e_List_From(1,:)=[];tmp_e_List_To(1,:)=[];tmp_e_List_Int(1,:)=[];
                
            % The spike arr. at 3rd neuron %                
            elseif isEq(tmp_e_List_To(i,1),3)

                % The spike arr. not the same time as 3rd neuron spike %
                if (isN3SpNormally==0)
                    tmpPRC=PRC(curr_phis(3,1),epsi23,b3,p_PRC_option);
                    
                    % 3rd neuron suddenly spikes? %
                    if (isEq(curr_phis(3,1),1-tmpPRC)==1)
                        isN3SpDriven=1;
                        
                        % Update 3rd neuron spikes interval & phase %
                        curr_phis(3,1)=0;
                        
                        % Add spike happens normally list %
                        addEvent(3,3,1);                        
                        % Add spike arr. list %
                        addEvent(3,2,tau3);   
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N3 rec. sp. from ',num2str(tmp_e_List_From(i,1)),' and sp.: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                        
                    else
                        % Update 3rd neuron spikes interval & phase %
                        curr_phis(3,1)=curr_phis(3,1)+tmpPRC;
                        
                        % Update the natual-spike 3rd neuron in the list %
                        i_3rd_neuron=findEvent(e_List_From,e_List_To,3,3);
                        cnt=size(i_3rd_neuron,1);
                        
                        % There exist something to be updated %
                        if (cnt>0)
                            e_List_Int(i_3rd_neuron)=1-curr_phis(3,1);
                        end 
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N3 rec. sp. from ',num2str(tmp_e_List_From(i,1)),': ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                        
                    end                    
                end                                                
                % Delete used event in tmp_e_List %
                tmp_e_List_From(1,:)=[];tmp_e_List_To(1,:)=[];tmp_e_List_Int(1,:)=[];
                
            % The spike arr. at 2nd neuron %
            else
                tmp_1_2=findEvent(tmp_e_List_From,tmp_e_List_To,1,2);
                tmp_3_2=findEvent(tmp_e_List_From,tmp_e_List_To,3,2);
                    
                % The spike arr. not the same time as 2nd neuron spike %
                if (isN2SpNormally==0)
                    cnt1=size(tmp_1_2,1);
                    cnt3=size(tmp_3_2,1);
                    
                    % spikes from 1st and 3rd neurons come? %
                    if ((cnt1==1) && (cnt3==1))
                        tmpPRC=PRC(curr_phis(2,1),epsi12+epsi32,b2,p_PRC_option);
                        
                        % 2nd neuron suddenly spikes? %
                        if (isEq(curr_phis(2,1),1-tmpPRC)==1)
                            isN2SpDriven=1;
                            % Update 2nd neuron spikes interval & phase %
                            curr_phis(2,1)=0;
                            
                            % Add spike happens normally list %
                            addEvent(2,2,1);
                            % Add spike arr. list %
                            addEvent(2,1,tau1);
                            addEvent(2,3,tau3);
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 rec. sp. from both and sp.: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                             
                        else
                            % Update 2nd neuron spikes interval & phase %
                            curr_phis(2,1)=curr_phis(2,1)+tmpPRC;
                            
                            % Update the natual-spike 2nd neuron in the list %
                            i_2nd_neuron=findEvent(e_List_From,e_List_To,2,2);
                            cnt=size(i_2nd_neuron,1);
                            
                            % There exist something to be updated %
                            if (cnt>0)
                                e_List_Int(i_2nd_neuron)=1-curr_phis(2,1);
                            end
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 rec. sp. from both: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                             
                        end
                    elseif (cnt1==1)
                        tmpPRC=PRC(curr_phis(2,1),epsi12,b2,p_PRC_option);
                        
                        % 2nd neuron suddenly spikes? %
                        if (isEq(curr_phis(2,1),1-tmpPRC)==1)
                            isN2SpDriven=1;
                            % Update 2nd neuron spikes interval & phase %
                            curr_phis(2,1)=0;
                            
                            % Add spike happens normally list %
                            addEvent(2,2,1);
                            % Add spike arr. list %
                            addEvent(2,1,tau1);
                            addEvent(2,3,tau3);
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 rec. sp. from 1 and sp.: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                            
                        else
                            % Update 2nd neuron spikes interval & phase %
                            curr_phis(2,1)=curr_phis(2,1)+tmpPRC;
                            
                            % Update the natual-spike 2nd neuron in the list %
                            i_2nd_neuron=findEvent(e_List_From,e_List_To,2,2);
                            cnt=size(i_2nd_neuron,1);
                            
                            % There exist something to be updated %
                            if (cnt>0)
                                e_List_Int(i_2nd_neuron)=1-curr_phis(2,1);
                            end
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 rec. sp. from 1: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                            
                        end
                    elseif (cnt3==1)
                        tmpPRC=PRC(curr_phis(2,1),epsi32,b2,p_PRC_option);
                        
                        % 2nd neuron suddenly spikes? %
                        if (isEq(curr_phis(2,1),1-tmpPRC)==1)
                            isN2SpDriven=1;
                            % Update 2nd neuron spikes interval & phase %
                            curr_phis(2,1)=0;
                            
                            % Add spike happens normally list %
                            addEvent(2,2,1);
                            % Add spike arr. list %
                            addEvent(2,1,tau1);
                            addEvent(2,3,tau3);
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 rec. sp. from 3 and sp.: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                            
                        else
                            % Update 2nd neuron spikes interval & phase %
                            curr_phis(2,1)=curr_phis(2,1)+tmpPRC;
                            
                            % Update the natual-spike 2nd neuron in the list %
                            i_2nd_neuron=findEvent(e_List_From,e_List_To,2,2);
                            cnt=size(i_2nd_neuron,1);
                            
                            % There exist something to be updated %
                            if (cnt>0)
                                e_List_Int(i_2nd_neuron)=1-curr_phis(2,1);
                            end
if (deBugLevel == 1) 
    disp(strcat('At ',num2str(t_event(i_t_event,1)),': N2 rec. sp. from 3: ',num2str(phi2voltage(curr_phis(1,1),b)),'_',num2str(phi2voltage(curr_phis(2,1),b)),'_',num2str(phi2voltage(curr_phis(3,1),b))));
end                            
                        end
                    else
                        disp('Strange!!!');
                    end                    
                end
                % Delete used event in tmp_e_List %
                tmp_1_2_3=[tmp_1_2;tmp_3_2];
                tmp_e_List_From(tmp_1_2_3,:)=[];tmp_e_List_To(tmp_1_2_3,:)=[];tmp_e_List_Int(tmp_1_2_3,:)=[];
            end                        
        end        
    end
    
%     if (isN3SpNormally || isN3SpDriven)
%         i_phisEvo=i_phisEvo+1;
%         phi1Evo(i_phisEvo,1)=curr_phis(1,1);
%         phi2Evo(i_phisEvo,1)=curr_phis(2,1);
%     end    
    if (isN2SpNormally || isN2SpDriven)
        i_phisEvo=i_phisEvo+1;
        phi1Evo(i_phisEvo,1)=curr_phis(1,1);
        phi3Evo(i_phisEvo,1)=curr_phis(3,1);
    end    

end

% val=[phi1Evo phi2Evo];
val=phi1Evo;

figure(11), hold on
set(gca,'FontSize', 18);
% title('phases of neurons at spiketimes of neuron 3')
ylabel('Phase')
% xlabel('spike time index of neuron 3')
xlabel('spike time index of neuron 2')
plot(phi1Evo,'b','LineWidth',4);hold on  % spikes of osc 1
% plot(phi2Evo,'c','LineWidth',2);  % spikes of osc 2
plot(phi3Evo,'c','LineWidth',2);  % spikes of osc 2
legend('neuron1','neuron3')
% maximize('all');











