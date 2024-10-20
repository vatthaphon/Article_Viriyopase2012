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
function val=TimeBased_3MS_osc(p_phi1, p_phi2, p_tau1, p_tau3, p_epsi12, p_epsi21, p_epsi32, p_epsi23, p_tEnd, p_option, p_dt)
%%%%%% Parameters %%%%%%
option=p_option;
b=3;            
dt=p_dt;       
tstart=0;tend=p_tEnd;        
tau1=p_tau1;tau3=p_tau3;     
epsi12=p_epsi12;epsi21=p_epsi21;
epsi32=p_epsi32;epsi23=p_epsi23;

init_phi1=p_phi1;init_phi2=p_phi2;
%%%%%% Initialize variables %%%%%%
phi=[init_phi1; init_phi2; 0];

Nbin=round(tend-tstart)/dt;
f=NaN(3,Nbin);         % init membrane potential matrix
PHI=NaN(3,2*(tend-tstart));     % init matrix with phase at spike times of MS-oscillator 3
time1=1; time2=1; time3=1;       % number of spike times
a_time1=0;a_time21=0;a_time23=0;a_time3=0;  % active index
SpTi1(1,1)=-phi(1,1);                % spike times before t=0
SpTi2(1,1)=-phi(2,1);                % spike times before t=0
SpTi3(1,1)=-phi(3,1);                % spike times before t=0
PHI(:,1)=phi';

if (isLEq(phi(1,1),tau1) && ~(isEq(phi(2,1),0) && isEq(phi(1,1),tau1)))
        a_time1=1;
end
if (isLEq(phi(2,1),tau1) &&  ~(isEq(phi(1,1),0) && isEq(phi(2,1),tau1)))
        a_time21=1;     
end    
if isL(phi(2,1),tau3) 
        a_time23=1;     
end
a_time3=1;

% main program
for t=tstart:dt:tend

    if ( (a_time1<=time1) && (a_time3<=time3) && (a_time1>0) && (a_time3>0) && ...
            (isInInt(t, t+dt, SpTi1(1,a_time1)+tau1)==1) && ...
            (isInInt(t, t+dt, SpTi3(1,a_time3)+tau3)==1))
        
        phi(2)=(exp(b*min(log(1+(exp(b)-1)*(phi(2)+dt))/b+epsi12+epsi32,1))-1)/(exp(b)-1);

        a_time1=a_time1+1;
        a_time3=a_time3+1;
    elseif ((a_time1<=time1) && (a_time1>0) && ...
            (isInInt(t, t+dt, SpTi1(1,a_time1)+tau1)==1))
        
        phi(2)=(exp(b*min(log(1+(exp(b)-1)*(phi(2)+dt))/b+epsi12,1))-1)/(exp(b)-1); 
        
        a_time1=a_time1+1;
    elseif ((a_time3<=time3) && (a_time3>0) && ...
            (isInInt(t, t+dt, SpTi3(1,a_time3)+tau3)==1))
        
        phi(2)=(exp(b*min(log(1+(exp(b)-1)*(phi(2)+dt))/b+epsi32,1))-1)/(exp(b)-1);
        
        a_time3=a_time3+1;
    else
        phi(2)=phi(2)+dt;
    end               
    
    % if osc 1 receives pulse from osc 2
    if ((a_time21<=time2) && (a_time21>0) && ...
            (isInInt(t, t+dt, SpTi2(1,a_time21)+tau1)==1))
        
        phi(1)=(exp(b*min(log(1+(exp(b)-1)*(phi(1)+dt))/b+epsi21,1))-1)/(exp(b)-1);
        a_time21=a_time21+1;
    else
        phi(1)=phi(1)+dt;
    end  % if osc 1 receives pulse from osc 2

    % if osc 3 receives pulse from osc 2
    if ((a_time23<=time2) && (a_time23>0) && ...
            (isInInt(t, t+dt, SpTi2(1,a_time23)+tau3)==1))
        
        phi(3)=(exp(b*min(log(1+(exp(b)-1)*(phi(3)+dt))/b+epsi23,1))-1)/(exp(b)-1);
        a_time23=a_time23+1;
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
    end % neuron 1 spikes

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
    end % neuron 2 spikes

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
    end % neuron 3 spikes

    f(1,round(t/dt+1)) = log(1+(exp(b)-1)*phi(1))/b;
    f(2,round(t/dt+1)) = log(1+(exp(b)-1)*phi(2))/b;
    f(3,round(t/dt+1)) = log(1+(exp(b)-1)*phi(3))/b;
    
end % for t-loop

val=PHI(1,:)';
i_NaN=isnan(val);
val(i_NaN,:)=[];

if (option==1)
    figure(11), hold on
    set(gca,'FontSize', 18);
    title('phases of neurons at spiketimes of neuron 3')
    ylabel('Phase')
    xlabel('spike time index of neuron 3')
    plot(PHI(1,:),'b','LineWidth',4);hold on  % spikes of osc 1
    plot(PHI(2,:),'c','LineWidth',2);  % spikes of osc 2
    legend('neuron1','neuron2')
elseif (option==2)
    figure(12), hold on
    set(gca,'ytick',[1 2 3],'FontSize', 18);
    title('Spike times')
    ylabel('neuron index')
    xlabel('time / T')
    axis([tstart tend 0 4])
    for i=2:length(SpTi)
        plot(SpTi(1,i),1,'.','MarkerSize',30);  % spikes of osc 1
        plot(SpTi(2,i),2,'c.','MarkerSize',30);  % spikes of osc 2
        plot(SpTi(3,i),3,'r.','MarkerSize',30);  % spikes of osc 3
    end;
    legend('neuron1','neuron2','neuron3')
    
elseif (option==3)
    figure(13), hold on
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
%     legend('Neuron 1','Neuron 2','Neuron 3');
    grid on
elseif (option==4)  
    figure(14)    
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
    
    
    figure(24) 
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
end

end

function val=isInInt(t_Begin, t_End, t_Curr)
tol=0.0000001;
% if ((t_Begin < t_Curr) && ( t_Curr <= t_End+tol))
%     val=1;
% else
%     val=0;
% end

if ((t_Begin-tol <= t_Curr) && ( t_Curr <= t_End+tol))
    val=1;
else
    val=0;
end   
end
