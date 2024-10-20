clc
clear all
format long

phi1_min=0;phi1_max=0.99;phi1len=50;
phi2_min=0;phi2_max=0.99;phi2len=50;

tau1_min=0.01;tau1_max=0.49;tau1len=60;
tau3_min=0.01;tau3_max=0.49;tau3len=60;
epsi12_min=0.01;epsi12_max=0.01;epsi12len=1;
epsi21_min=0.01;epsi21_max=0.01;epsi21len=1;
epsi32_min=0.01;epsi32_max=0.01;epsi32len=1;
epsi23_min=0.01;epsi23_max=0.01;epsi23len=1;

isSymTau=0; % 0 No, 1 Yes
isSymE=1;
i_end_phi1=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi1_lin=linspace(phi1_min,phi1_max,phi1len);
phi2_lin=linspace(phi2_min,phi2_max,phi2len);
tau1_lin=linspace(tau1_min,tau1_max,tau1len);
tau3_lin=linspace(tau3_min,tau3_max,tau3len);
epsi12_lin=linspace(epsi12_min,epsi12_max,epsi12len);
epsi21_lin=linspace(epsi21_min,epsi21_max,epsi21len);
epsi32_lin=linspace(epsi32_min,epsi32_max,epsi32len);
epsi23_lin=linspace(epsi23_min,epsi23_max,epsi23len);

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
                        
                        if (i_end_phi1==16)
                            init_phi1
                            init_phi2
                            tau1
                            epsi12
                        end                        
                        i_end_phi1=i_end_phi1+1;

                    end
                else
                end
            end
        else
        end
    end
end

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
if (i_end_phi1==16)
    init_phi1
    init_phi2
    tau1
	epsi12
end                        
i_end_phi1=i_end_phi1+1;

                    end
                else
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
if (i_end_phi1==16)
    init_phi1
    init_phi2
    tau1
    epsi12
end                        
    i_end_phi1=i_end_phi1+1;                        
                        end
                    else
                    end
                end
            end
        end
    end
end
 