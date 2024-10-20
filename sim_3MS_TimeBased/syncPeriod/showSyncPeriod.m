clc
clear all
close all

filename='tend_50_P1_0.00_0.99_10_P2_0.00_0.99_10_P3_0.00_0.00_1_T1_0.25_0.25_1_T3_0.25_0.25_1_E1_0.01_0.99_20_E3_0.01_0.99_20_isSCS_1_isST_1_tDt_0.0010_sp_.txt';

fid = fopen(filename,'rt');          
if (fid < 0)
    error('could not open file');   
end

descp=regexp(filename, '_', 'split');
tEnd=str2double(descp(1,2));
phi1Begin=str2double(descp(1,4));phi1End=str2double(descp(1,5));phi1Len=str2double(descp(1,6));
phi2Begin=str2double(descp(1,8));phi2End=str2double(descp(1,9));phi2Len=str2double(descp(1,10));
phi3Begin=str2double(descp(1,12));phi3End=str2double(descp(1,13));phi3Len=str2double(descp(1,14));
tau1Begin=str2double(descp(1,16));tau1End=str2double(descp(1,17));tau1Len=str2double(descp(1,18));
tau3Begin=str2double(descp(1,20));tau3End=str2double(descp(1,21));tau3Len=str2double(descp(1,22));
epsi1Begin=str2double(descp(1,24));epsi1End=str2double(descp(1,25));epsi1Len=str2double(descp(1,26));

if (phi1Len==1) 
    phi1_lin=phi1Begin; 
else
    phi1_lin=linspace(phi1Begin,phi1End,phi1Len);
end

if (phi2Len==1) 
    phi2_lin=phi2Begin; 
else
    phi2_lin=linspace(phi2Begin,phi2End,phi2Len);
end

if (phi3Len==1) 
    phi3_lin=phi1Begin; 
else
    phi3_lin=linspace(phi3Begin,phi3End,phi3Len);
end

if (tau1Len==1) 
    tau1_lin=tau1Begin; 
else
    tau1_lin=linspace(tau1Begin,tau1End,tau1Len);
end

if (epsi1Len==1) 
    epsi1_lin=epsi1Begin; 
else
    epsi1_lin=linspace(epsi1Begin,epsi1End,epsi1Len);
end

whichNeuron=0;
i1=1;i2=1;i3=1;i4=1;
sp1=zeros(phi1Len*phi2Len*phi3Len*tau1Len*epsi1Len,2);
sp2=zeros(phi1Len*phi2Len*phi3Len*tau1Len*epsi1Len,2);
sp3=zeros(phi1Len*phi2Len*phi3Len*tau1Len*epsi1Len,2);
phi=zeros(phi1Len*phi2Len*phi3Len*tau1Len*epsi1Len,2);

tline = fgetl(fid);
while ischar(tline)
    if strcmp(tline,'Spike1:')
        step=1;
    elseif strcmp(tline,'Spike2:')
        step=2;
    elseif strcmp(tline,'Spike3:')
        step=3;
    elseif strcmp(tline,'EvoPhi:')
        step=4;
    end

    tokens=regexp(tline, ';', 'split');    
    l=length(tokens);
    
    if (length(tokens)>1)
        if (step==1)
            lastI=str2double(tokens(1,l-1));
            prevlI=str2double(tokens(1,l-2));
            sp1(i1,1)=prevlI;sp1(i1,2)=lastI;
            i1=i1+1;
        elseif (step==2)
            lastI=str2double(tokens(1,l-1));
            prevlI=str2double(tokens(1,l-2));
            sp2(i2,1)=prevlI;sp2(i2,2)=lastI;
            i2=i2+1;
        elseif (step==3)
            lastI=str2double(tokens(1,l-1));
            prevlI=str2double(tokens(1,l-2));
            sp3(i3,1)=prevlI;sp3(i3,2)=lastI;
            i3=i3+1;
        elseif (step==4)
            lastI=str2double(tokens(1,l-1));
            prevlI=str2double(tokens(1,l-2));
            phi(i4,1)=prevlI;phi(i4,2)=lastI;
            i4=i4+1;           
        end
    end
        
    tline = fgetl(fid);
end

index=1;i_NP=1;s_index=1;
for i=1:1:phi1Len
    for j=1:1:phi2Len
        for k=1:1:phi3Len
            for l=1:1:tau1Len
                for m=1:1:epsi1Len
                    if (isEq(phi(index,1),0) && isEq(phi(index,2),0))
                        newPeriod(i_NP)=sp1(index,2)-sp1(index,1);
                        
                        if (newPeriod(i_NP)<0.6) && (newPeriod(i_NP)>0.4)
                            s_init(s_index,1)=phi1_lin(i);
                            s_init(s_index,2)=phi2_lin(j);
                            s_init(s_index,3)=phi3_lin(k);
                            s_init(s_index,4)=tau1_lin(l);
                            s_init(s_index,5)=epsi1_lin(m);

                            s_index=s_index+1;
                        end                        
                        i_NP=i_NP+1;
                    end
                    index=index+1;
                end
            end
        end
    end
end

% plot(newPeriod,'+');

fclose(fid);
                       
