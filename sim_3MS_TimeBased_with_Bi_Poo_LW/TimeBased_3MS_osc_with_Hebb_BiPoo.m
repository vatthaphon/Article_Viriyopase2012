% determines the membrane potentials of three coupled MS-oscillators
%
%   osc. 1  ---<-->---  osc. 2 ---<-->---  osc. 3
%               tau1               tau3
%
% with symmetric connectivity and equal intrinsic period, T. 
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
% Date:  17 July, 2009                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rec_SpTi1 rec_SpTi2 rec_SpTi3 rec_12_arr rec_21_arr rec_23_arr rec_32_arr o_w_12 o_w_21 o_w_23 o_w_32]...
            =TimeBased_3MS_osc_with_Hebb_BiPoo(p_phi1, p_phi2, p_phi3, p_tau1, p_tau3, p_epsi12, p_epsi21, p_epsi32, p_epsi23, p_tEnd, p_option, p_dt)
global max_cs        
%%%%%% Parameters %%%%%%
option=p_option;
b=3;            
dt=p_dt;       
tstart=0;tend=p_tEnd;        
tau1=p_tau1;tau3=p_tau3;     
epsi12=p_epsi12;epsi21=p_epsi21;
epsi32=p_epsi32;epsi23=p_epsi23;
T_0=25; % ms

PRE_EVENT=0;
POST_EVENT=1;

init_phi1=p_phi1;init_phi2=p_phi2;init_phi3=p_phi3;
%%%%%% Initialize variables %%%%%%
phi=[init_phi1; init_phi2; init_phi3];

Nbin=round(tend-tstart)/dt;
f=NaN(3,Nbin);         % init membrane potential matrix
PHI=NaN(3,2*(tend-tstart));     % init matrix with phase at spike times of MS-oscillator 3
time1=1; time2=1; time3=1;       % number of spike times
a_time1=0;a_time21=0;a_time23=0;a_time3=0;  % active index
SpTi1=NaN(1, 3*(tend-tstart));
SpTi2=NaN(1, 3*(tend-tstart));
SpTi3=NaN(1, 3*(tend-tstart));
SpTi1(1,1)=-phi(1,1);                % spike times before t=0
SpTi2(1,1)=-phi(2,1);                % spike times before t=0
SpTi3(1,1)=-phi(3,1);                % spike times before t=0
PHI(:,1)=phi';

% Variables for recording the lastest presynaptic spike arrival and 
% the lastest postsynaptic neuron spikes
is_N1_spike=NaN;        
is_N2_spike=NaN;        
is_N3_spike=NaN;        

is_E12_arr=NaN;
is_E21_arr=NaN;
is_E23_arr=NaN;
is_E32_arr=NaN;

rec_epsi12=NaN(ceil((tend-tstart)/dt)+10,1);re12_index=1;
rec_epsi21=NaN(ceil((tend-tstart)/dt)+10,1);re21_index=1;
rec_epsi23=NaN(ceil((tend-tstart)/dt)+10,1);re23_index=1;
rec_epsi32=NaN(ceil((tend-tstart)/dt)+10,1);re32_index=1;

rec_SpTi1=NaN(3*(tend-tstart), 1);rec_SpTi1_i=1;
rec_SpTi2=NaN(3*(tend-tstart), 1);rec_SpTi2_i=1;
rec_SpTi3=NaN(3*(tend-tstart), 1);rec_SpTi3_i=1;

rec_12_arr=NaN(3*(tend-tstart), 1);rec_12_arr_i=1;
rec_21_arr=NaN(3*(tend-tstart), 1);rec_21_arr_i=1;
rec_23_arr=NaN(3*(tend-tstart), 1);rec_23_arr_i=1;
rec_32_arr=NaN(3*(tend-tstart), 1);rec_32_arr_i=1;


if (isLEq(phi(1,1),tau1) && ~(isEq(phi(2,1),0) && isEq(phi(1,1),tau1)))
        a_time1=1;
end
if (isLEq(phi(2,1),tau1) &&  ~(isEq(phi(1,1),0) && isEq(phi(2,1),tau1)))
        a_time21=1;     
end    
if (isLEq(phi(2,1),tau3) &&  ~(isEq(phi(3,1),0) && isEq(phi(2,1),tau3)))
        a_time23=1;     
end    
if (isLEq(phi(3,1),tau3) && ~(isEq(phi(2,1),0) && isEq(phi(3,1),tau3)))
        a_time3=1;
end

% %% When we assume that osc 3 always starts with 0, a_time3=1;
% if (isLEq(phi(2,1),tau1) &&  ~(isEq(phi(1,1),0) && isEq(phi(2,1),tau1)))
%         a_time21=1;     
% end    
% if isL(phi(2,1),tau3) 
%         a_time23=1;     
% end
% a_time3=1;


%% main program
for t=tstart:dt:tend
    is_N1_spike=0;
    is_N2_spike=0;
    is_N3_spike=0;
    
    is_E12_arr=0;
    is_E21_arr=0;
    is_E23_arr=0;
    is_E32_arr=0;
        
    if ( (a_time1<=time1) && (a_time3<=time3) && (a_time1>0) && (a_time3>0) && ...
            (isInInt(t, t+dt, SpTi1(1,a_time1)+tau1)==1) && ...
            (isInInt(t, t+dt, SpTi3(1,a_time3)+tau3)==1))
        
        phi(2)=(exp(b*min(log(1+(exp(b)-1)*(phi(2)+dt))/b+epsi12+epsi32,1))-1)/(exp(b)-1);

        a_time1=a_time1+1;
        a_time3=a_time3+1;
        
        is_E12_arr=1;
        is_E32_arr=1;        
                
        rec_12_arr(rec_12_arr_i, 1)=t;rec_12_arr_i=rec_12_arr_i+1;
        rec_32_arr(rec_32_arr_i, 1)=t;rec_32_arr_i=rec_32_arr_i+1;
    elseif ((a_time1<=time1) && (a_time1>0) && ...
            (isInInt(t, t+dt, SpTi1(1,a_time1)+tau1)==1))
        
        phi(2)=(exp(b*min(log(1+(exp(b)-1)*(phi(2)+dt))/b+epsi12,1))-1)/(exp(b)-1); 
        
        a_time1=a_time1+1;

        is_E12_arr=1;        
        
        rec_12_arr(rec_12_arr_i, 1)=t;rec_12_arr_i=rec_12_arr_i+1;
    elseif ((a_time3<=time3) && (a_time3>0) && ...
            (isInInt(t, t+dt, SpTi3(1,a_time3)+tau3)==1))
        
        phi(2)=(exp(b*min(log(1+(exp(b)-1)*(phi(2)+dt))/b+epsi32,1))-1)/(exp(b)-1);
        
        a_time3=a_time3+1;
        
        is_E32_arr=1; 
        
        rec_32_arr(rec_32_arr_i, 1)=t;rec_32_arr_i=rec_32_arr_i+1;
    else
        phi(2)=phi(2)+dt;
    end               
    
    % if osc 1 receives pulse from osc 2
    if ((a_time21<=time2) && (a_time21>0) && ...
            (isInInt(t, t+dt, SpTi2(1,a_time21)+tau1)==1))
        
        phi(1)=(exp(b*min(log(1+(exp(b)-1)*(phi(1)+dt))/b+epsi21,1))-1)/(exp(b)-1);
                
        a_time21=a_time21+1;
        
        is_E21_arr=1;
        
        rec_21_arr(rec_21_arr_i, 1)=t;rec_21_arr_i=rec_21_arr_i+1;
    else
        phi(1)=phi(1)+dt;
    end

    % if osc 3 receives pulse from osc 2
    if ((a_time23<=time2) && (a_time23>0) && ...
            (isInInt(t, t+dt, SpTi2(1,a_time23)+tau3)==1))
        
        phi(3)=(exp(b*min(log(1+(exp(b)-1)*(phi(3)+dt))/b+epsi23,1))-1)/(exp(b)-1);
        a_time23=a_time23+1;
        
        is_E23_arr=1; 
        
        rec_23_arr(rec_23_arr_i, 1)=t;rec_23_arr_i=rec_23_arr_i+1;
    else
        phi(3)=phi(3)+dt;
    end  % if osc 1 receives pulse from osc 2

    if ~(phi(1)<1)  % neuron 1 spikes
        time1=time1+1;
        SpTi1(1,time1)=t;
        phi(1)=0;
        
        if (a_time1==0) 
            a_time1=2; 
        end
        
        is_N1_spike=1;
        
        rec_SpTi1(rec_SpTi1_i, 1)=t;rec_SpTi1_i=rec_SpTi1_i+1;
    end

    if ~(phi(2)<1)   % neuron 2 spikes
        time2=time2+1;
        SpTi2(1,time2)=t;
        phi(2)=0;
        
        if (a_time21==0) 
            a_time21=2; 
        end
        if (a_time23==0) 
            a_time23=2; 
        end        
        
        is_N2_spike=1;   
        
        rec_SpTi2(rec_SpTi2_i, 1)=t;rec_SpTi2_i=rec_SpTi2_i+1;
    end

    if ~(phi(3)<1)  % neuron 3 spikes
        phi(3)=0;
        time3=time3+1;
        SpTi3(1,time3)=t;
        PHI(1,time3)=phi(1);
        PHI(2,time3)=phi(2);
        PHI(3,time3)=phi(3);
        
        if (a_time3==0) 
            a_time3=2; 
        end        
        
        is_N3_spike=1;  
        
        rec_SpTi3(rec_SpTi3_i, 1)=t;rec_SpTi3_i=rec_SpTi3_i+1;
    end

    f(1,round(t/dt+1)) = log(1+(exp(b)-1)*phi(1))/b;
    f(2,round(t/dt+1)) = log(1+(exp(b)-1)*phi(2))/b;
    f(3,round(t/dt+1)) = log(1+(exp(b)-1)*phi(3))/b;
    
    tmpEpsi12=epsi12;
    tmpEpsi21=epsi21;
    tmpEpsi23=epsi23;
    tmpEpsi32=epsi32;
    
    deltaEpsi12Pot=0;    
    deltaEpsi21Pot=0;        
    deltaEpsi23Pot=0;
    deltaEpsi32Pot=0;    
    
    deltaEpsi12Dep=0;
    deltaEpsi21Dep=0;
    deltaEpsi23Dep=0;
    deltaEpsi32Dep=0;
    
    %% Update the coupling strengths in when the neurons spike.       
    if (is_N1_spike==1)
        deltaEpsi21Pot=adjust_W(tmpEpsi21, rec_21_arr, t, POST_EVENT, T_0); 
    end   
    if (is_N2_spike==1)        
        deltaEpsi12Pot=adjust_W(tmpEpsi12, rec_12_arr, t, POST_EVENT, T_0);
        deltaEpsi32Pot=adjust_W(tmpEpsi32, rec_32_arr, t, POST_EVENT, T_0);        
    end    
    if (is_N3_spike==1) 
        deltaEpsi23Pot=adjust_W(tmpEpsi23, rec_23_arr, t, POST_EVENT, T_0); 
    end
    
    %% Update the coupling strengths when the presynaptic input arrives.
    if (is_E12_arr==1) 
        deltaEpsi12Dep=adjust_W(tmpEpsi12, rec_SpTi2, t, PRE_EVENT, T_0); 
    end
    if (is_E21_arr==1) 
        deltaEpsi21Dep=adjust_W(tmpEpsi21, rec_SpTi1, t, PRE_EVENT, T_0); 
    end        
    if (is_E23_arr==1) 
        deltaEpsi23Dep=adjust_W(tmpEpsi23, rec_SpTi3, t, PRE_EVENT, T_0); 
    end    
    if (is_E32_arr==1) 
        deltaEpsi32Dep=adjust_W(tmpEpsi32, rec_SpTi2, t, PRE_EVENT, T_0); 
    end            
        
    epsi12=epsi12+deltaEpsi12Pot+deltaEpsi12Dep;
    epsi21=epsi21+deltaEpsi21Pot+deltaEpsi21Dep;
    epsi23=epsi23+deltaEpsi23Pot+deltaEpsi23Dep;
    epsi32=epsi32+deltaEpsi32Pot+deltaEpsi32Dep;
    
    if (epsi12<0) epsi12=0; end
    if (epsi12>max_cs) epsi12=max_cs; end
    
    if (epsi21<0) epsi21=0; end
    if (epsi21>max_cs) epsi21=max_cs; end
    
    if (epsi23<0) epsi23=0; end
    if (epsi23>max_cs) epsi23=max_cs; end
    
    if (epsi32<0) epsi32=0; end
    if (epsi32>max_cs) epsi32=max_cs; end    
    
    % Record the coupling strength
    rec_epsi12(re12_index,1)=epsi12;re12_index=re12_index+1;
    rec_epsi21(re21_index,1)=epsi21;re21_index=re21_index+1;    
    rec_epsi23(re23_index,1)=epsi23;re23_index=re23_index+1; 
    rec_epsi32(re32_index,1)=epsi32;re32_index=re32_index+1;    

    o_w_12=epsi12; 
    o_w_21=epsi21; 
    o_w_23=epsi23; 
    o_w_32=epsi32;  
end

%% Show the results
if (option==1)
    figure(51), hold on
    set(gca,'FontSize', 18);
    title('phases of neurons at spiketimes of neuron 3')
    ylabel('Phase')
    xlabel('spike time index of neuron 3')
    plot(PHI(1,:),'b','LineWidth',4);hold on  % spikes of osc 1
    plot(PHI(2,:),'c','LineWidth',2);  % spikes of osc 2
    legend('neuron1','neuron2')
elseif (option==3)
    figure(53), hold on
    set(gca,'FontSize', 18);
    ylabel('V_{membrane}')
    xlabel('time / T_0')
    
    plot(tstart:dt:tend,f(1,:),'LineWidth',3), hold on
    for i=2: time1
        plot([SpTi1(1,i) SpTi1(1,i)],[0 1],'LineWidth',3)
    end
    
    plot(tstart:dt:tend,f(2,:),'c','LineWidth',3), hold on
    for i=2: time2
        plot([SpTi2(1,i) SpTi2(1,i)],[0 1],'c','LineWidth',3)
    end
    
    plot(tstart:dt:tend,f(3,:),'r','LineWidth',3), hold on
    for i=2: time3
        plot([SpTi3(1,i) SpTi3(1,i)],[0 1],'r','LineWidth',2)
    end
    
    grid on

elseif (option==4)  
    figure(54)    
    hold on
    set(gca,'FontSize', 18);
    title('phases of neurons at spiketimes of neuron 3')
    ylabel('Phase')
    xlabel('spike time index of neuron 3')
    plot(PHI(1,:),'b','LineWidth',4);hold on  % spikes of osc 1
    plot(PHI(2,:),'c','LineWidth',2);  % spikes of osc 2
    
    lastSpTi1=SpTi1(1,time1);
    t_lin=tstart:dt:tend;
    minLastSpTi1=min(find(t_lin>=lastSpTi1));
    disp(strcat('Phi2 at last spike of Phi1: ',num2str(voltage2phi(f(2,minLastSpTi1),b)))); 
    
    if SpTi1(1,time1)>=SpTi2(1,time2)
        disp(strcat('Lag: ',num2str(SpTi1(1,time1)-SpTi2(1,time2))));        
    else
        disp(strcat('Lag: ',num2str(SpTi1(1,time1)-SpTi2(1,time2-1))));        
    end

    disp(strcat('T: ', num2str(SpTi1(1,time1)-SpTi1(1,time1-1))));
        
    figure(542) 
    hold on
    set(gca,'FontSize', 18);
    title('Dynamics of V_{Membrane}')
    ylabel('V_{membrane}')
    xlabel('time / T')
    
    t_lin=tstart:dt:tend;
    tmptstart=SpTi1(1,time1-2); 
    filteredT_lin=find(t_lin>tmptstart);
    plot(t_lin(1,filteredT_lin),f(1,filteredT_lin),'LineWidth',3), hold on
    for i=2: time1
        plot([SpTi1(1,i) SpTi1(1,i)],[0 1],'LineWidth',3)
    end
    
    plot(t_lin(1,filteredT_lin),f(2,filteredT_lin),'c','LineWidth',3), hold on
    for i=2: time2
        plot([SpTi2(1,i) SpTi2(1,i)],[0 1],'c','LineWidth',3)
    end
    
    plot(t_lin(1,filteredT_lin),f(3,filteredT_lin),'r','LineWidth',2), hold on
    for i=2: time3
        plot([SpTi3(1,i) SpTi3(1,i)],[0 1],'r','LineWidth',2)
    end
    grid on    
    axis([tmptstart tend -0.1 1.1]);  
    
elseif (option==5)
    figure(51)
    subplot(2,1,1);
    hold on
    set(gca,'FontSize', 18);        
    xlabel('spike time of neuron 3'), ylabel('Phase'), title('phases of neurons at spiketimes of neuron 3')
    plot(SpTi3(1,~isnan(SpTi3(1,:))), PHI(1,~isnan(PHI(1,:))),'b','LineWidth',4);   % spikes of osc 1
    plot(SpTi3(1,~isnan(SpTi3(1,:))), PHI(2,~isnan(PHI(2,:))),'c','LineWidth',2);   % spikes of osc 2
    xlim([0 SpTi3(1,time3)]);
    grid on

    figure(51), 
    subplot(2,1,2);
    hold on
    set(gca,'FontSize', 18);
    title('Dynamics of V_{Membrane}')
    ylabel('V_{membrane}')
    xlabel('time / T')
    
    plot(tstart:dt:tend,f(1,:),'LineWidth',3), hold on
    for i=2: time1
        plot([SpTi1(1,i) SpTi1(1,i)],[0 1],'LineWidth',3)
    end
    
    plot(tstart:dt:tend,f(2,:),'c','LineWidth',3), hold on
    for i=2: time2
        plot([SpTi2(1,i) SpTi2(1,i)],[0 1],'c','LineWidth',3)
    end
    
    plot(tstart:dt:tend,f(3,:),'r','LineWidth',3), hold on
    for i=2: time3
        plot([SpTi3(1,i) SpTi3(1,i)],[0 1],'r','LineWidth',2)
    end
    xlim([0 SpTi3(1,time3)]);
    grid on    
        
end

if (option~=0)
    figure(531)
    subplot(2,2,1)
    plot(tstart:dt:tend, rec_epsi12(~isnan(rec_epsi12(:,1)),1),'*','LineWidth',4);
    ylabel('\epsilon_{12}');xlabel('time / T_0');
    xlim([0 SpTi3(1,time3)]);
    
    subplot(2,2,2)
    plot(tstart:dt:tend, rec_epsi21(~isnan(rec_epsi21(:,1)),1),'*','LineWidth',4);
    ylabel('\epsilon_{21}');xlabel('time / T_0');
    xlim([0 SpTi3(1,time3)]);
    
    subplot(2,2,3)
    plot(tstart:dt:tend, rec_epsi23(~isnan(rec_epsi23(:,1)),1),'*','LineWidth',4);
    ylabel('\epsilon_{23}');xlabel('time / T_0');
    xlim([0 SpTi3(1,time3)]);
    
    subplot(2,2,4)
    plot(tstart:dt:tend, rec_epsi32(~isnan(rec_epsi32(:,1)),1),'*','LineWidth',4);
    ylabel('\epsilon_{32}');xlabel('time / T_0');
    xlim([0 SpTi3(1,time3)]);
end

end

function val=adjust_W(epsi, rec_Event, t, event, T_0)
    global A_pot_BiPoo A_dep_BiPoo tau_pot tau_dep
    PRE_EVENT=0;      
    
    Mul_LW=10;
        
    delta_t=(t-rec_Event(~isnan(rec_Event), 1));
    

    
    if (event==PRE_EVENT)
        gZero=0<delta_t(:,1);
        lUpper=delta_t(:,1)<Mul_LW*tau_dep/T_0;
    
%         delta_t=delta_t(0<delta_t(:,1)<Mul_LW*tau_dep/T_0, 1);
        delta_t=delta_t(gZero & lUpper, 1);
        val=Sum_F_PRE_EVENT(A_dep_BiPoo, delta_t, epsi, tau_dep/T_0);
    else
        gZero=0<delta_t(:,1);
        lUpper=delta_t(:,1)<Mul_LW*tau_pot/T_0;
        
%         delta_t=delta_t(0<delta_t(:,1)<Mul_LW*tau_pot/T_0, 1);
        delta_t=delta_t(gZero & lUpper, 1);
        val=Sum_F_POST_EVENT(A_pot_BiPoo, -delta_t, epsi, tau_pot/T_0);
    end
end

% Depression part
function val=Sum_F_PRE_EVENT(A, s, epsi, tau)
    no_Pairs=60;

	val=sum((-A.*epsi.*exp(-s./tau))./no_Pairs);
end

% Potentiation part
function val=Sum_F_POST_EVENT(A, s, epsi, tau)
	no_Pairs=60;

	val=sum((A.*epsi.*exp(s./tau))./no_Pairs);
end

function val=isInInt(t_Begin, t_End, t_Curr)
tol=0.0000001;

if ((t_Begin-tol <= t_Curr) && ( t_Curr <= t_End+tol))
    val=1;
else
    val=0;
end

end
