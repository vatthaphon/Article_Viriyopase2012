clc
clear all
close all

tend=50;

phi1_min=0;phi1_max=phi1_min;phi1len=1;
phi2_min=0.948;phi2_max=phi2_min;phi2len=1;

tau1_min=0.05;tau1_max=tau1_min;tau1len=1;
tau3_min=0.4;tau3_max=tau3_min;tau3len=1;
epsi12_min=0.1;epsi12_max=epsi12_min;epsi12len=1;
epsi21_min=0.1;epsi21_max=epsi21_min;epsi21len=1;
epsi23_min=0.1;epsi23_max=epsi23_min;epsi23len=1;
epsi32_min=0.1;epsi32_max=epsi32_min;epsi32len=1;

isSymTau=0; % 0 No, 1 Yes
isSymE=0;

dt=0.0001;      % 0.000001<0.0001
figOption=3;    % 0=no figure, 1=phase info., 2=raster plot, 3=voltage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi1_lin=linspace(phi1_min,phi1_max,phi1len);
phi2_lin=linspace(phi2_min,phi2_max,phi2len);
tau1_lin=linspace(tau1_min,tau1_max,tau1len);
tau3_lin=linspace(tau3_min,tau3_max,tau3len);
epsi12_lin=linspace(epsi12_min,epsi12_max,epsi12len);
epsi21_lin=linspace(epsi21_min,epsi21_max,epsi21len);
epsi32_lin=linspace(epsi32_min,epsi32_max,epsi32len);
epsi23_lin=linspace(epsi23_min,epsi23_max,epsi23len);

i_end_phi1=1;

% i<->phi1, j<->phi2, k<->tau1, l<->tau3, m<->epsi12, n<->epsi21, o<->epsi32, p<->epsi23 %
for i=1:1:phi1len  % Phi1
    init_phi1=phi1_lin(1,i);
    for j=1:1:phi2len  % Phi2
        init_phi2=phi2_lin(1,j);
        if isSymTau
            for k=1:1:tau1len
                tau1=tau1_lin(1,k);
                tau3=tau1;
                if isSymE
                    for m=1:1:epsi12len
                        epsi12=epsi12_lin(1,m);
                        epsi21=epsi12;
                        epsi32=epsi12;
                        epsi23=epsi12;
retStr=TimeBased_3MS_osc(init_phi1, init_phi2, tau1, tau3, epsi12, epsi21, epsi32, epsi23, tend, figOption, dt);                        

rows=size(retStr,1);
if isEq(retStr(rows-1,1),retStr(rows,1))
    if isLEq(0.00001,retStr(rows,1))
        end_phi1(i_end_phi1,1)=retStr(rows,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
    else
        end_phi1(i_end_phi1,1)=0;end_phi1(i_end_phi1,2)=0;
    end
else
    end_phi1(i_end_phi1,1)=retStr(rows-1,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
end
i_end_phi1=i_end_phi1+1;

                    end
                else
                    for m=1:1:epsi12len
                        epsi12=epsi12_lin(1,m);
                        epsi21=epsi21_lin(1,m);
                        epsi32=epsi32_lin(1,m);
                        epsi23=epsi23_lin(1,m);
retStr=TimeBased_3MS_osc(init_phi1, init_phi2, tau1, tau3, epsi12, epsi21, epsi32, epsi23, tend, figOption, dt);                        

rows=size(retStr,1);
if isEq(retStr(rows-1,1),retStr(rows,1))
    if isLEq(0.00001,retStr(rows,1))
        end_phi1(i_end_phi1,1)=retStr(rows,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
    else
        end_phi1(i_end_phi1,1)=0;end_phi1(i_end_phi1,2)=0;
    end
else
    end_phi1(i_end_phi1,1)=retStr(rows-1,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
end
i_end_phi1=i_end_phi1+1;
                    
                    end
                    
                end
            end
        else
            for k=1:1:tau1len
                tau1=tau1_lin(1,k);
                for l=1:1:tau3len
                    tau3=tau3_lin(1,l);                
                    if isSymE
                        for m=1:1:epsi12len
                            epsi12=epsi12_lin(1,m);
                            epsi21=epsi12;
                            epsi32=epsi12;
                            epsi23=epsi12;                        
retStr=TimeBased_3MS_osc(init_phi1, init_phi2, tau1, tau3, epsi12, epsi21, epsi32, epsi23, tend, figOption, dt);                        

rows=size(retStr,1);
if isEq(retStr(rows-1,1),retStr(rows,1))
    if isLEq(0.00001,retStr(rows,1))
        end_phi1(i_end_phi1,1)=retStr(rows,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
    else
        end_phi1(i_end_phi1,1)=0;end_phi1(i_end_phi1,2)=0;
    end
else
    end_phi1(i_end_phi1,1)=retStr(rows-1,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
end
i_end_phi1=i_end_phi1+1;                      
                        end
                    else
                        for m=1:1:epsi12len
                            epsi12=epsi12_lin(1,m);
                            epsi21=epsi21_lin(1,m);
                            for o=1:1:epsi32len
                            epsi32=epsi32_lin(1,o);
                            epsi23=epsi23_lin(1,o);                        
retStr=TimeBased_3MS_osc(init_phi1, init_phi2, tau1, tau3, epsi12, epsi21, epsi32, epsi23, tend, figOption, dt);                        

rows=size(retStr,1);
if isEq(retStr(rows-1,1),retStr(rows,1))
    if isLEq(0.00001,retStr(rows,1))
        end_phi1(i_end_phi1,1)=retStr(rows,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
    else
        end_phi1(i_end_phi1,1)=0;end_phi1(i_end_phi1,2)=0;
    end
else
    end_phi1(i_end_phi1,1)=retStr(rows-1,1);end_phi1(i_end_phi1,2)=retStr(rows,1);
end
i_end_phi1=i_end_phi1+1;                      
                            end
                        end                        
                        
                    end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strcat= ['[1,3   driven + 2   driven    ]',' 1,3 vs 2:',num2str(tau1_min),' new period:', num2str(2*tau1_min)];
strcat1=['[1,3 n-driven + 2   driven (U)]',' 1,3 vs 2:',num2str(0),'<lag<',num2str(tau1_min),' new period:', num2str(1-alpha(epsi12_min,3)-beta(epsi12_min,3)*tau1_min+beta(epsi12_min,3)*phi2_min)];
strcat2=['[1,3 n-driven + 2   driven (L)]',' 1,3 vs 2:',num2str(tau1_min),'<lag<',num2str(1-tau1_min),' new period:', num2str(1-alpha(epsi12_min,3)-2*beta(epsi12_min,3)*tau1_min)];
% strcat1=['[1,3 n-driven + 2   driven]',' 1,3 vs 2:',num2str(phi2_min),' new period:', num2str(1-alpha(epsi12_min,3)-beta(epsi12_min,3)*tau1_min+beta(epsi12_min,3)*phi2_min)];
% strcat2=['[1,3 n-driven + 2   driven]',' 1,3 vs 2:',num2str(phi2_min),' new period:', num2str(1-alpha(epsi12_min,3)-2*beta(epsi12_min,3)*tau1_min)];
% strcat1=['[1,3 n-driven + 2   driven]',' 1,3 vs 2:phi2 new period:', num2str(1-alpha(epsi12_min,3)-beta(epsi12_min,3)*tau1_min+beta(epsi12_min,3)*phi2_min)];
% strcat2=['[1,3 n-driven + 2   driven]',' 1,3 vs 2:phi2 new period:', num2str(1-alpha(epsi12_min,3)-2*beta(epsi12_min,3)*tau1_min)];
strcat3=['[1,3   driven + 2 n-driven    ]',' 1,3 vs 2:',num2str(tau1_min),' new period:', num2str(1-alpha(2*epsi12_min,3)-2*beta(2*epsi12_min,3)*tau1_min)];
disp(strcat);
disp(strcat1);
disp(strcat2);
disp(strcat3);


