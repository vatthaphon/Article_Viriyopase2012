clear all
close all
clc

filename='tend_50_P1_0.00_0.99_100_P2_0.50_0.50_1_P3_0.00_0.00_1_T1_0.13_0.13_1_T3_0.13_0.13_1_E1_0.04_0.04_1_E3_0.04_0.04_1_isSymCS_1_isSymTau_1_tDt_0.0010_.txt';
% filename='tend_50_P1_0.00_0.99_100_P2_0.25_0.25_1_P3_0.00_0.00_1_T1_0.13_0.13_1_T3_0.13_0.13_1_E1_0.04_0.04_1_E3_0.04_0.04_1_isSymCS_1_isSymTau_1_tDt_0.0010_.txt';
% filename='tend_50_P1_0.00_0.99_100_P2_0.75_0.75_1_P3_0.00_0.00_1_T1_0.13_0.13_1_T3_0.13_0.13_1_E1_0.04_0.04_1_E3_0.04_0.04_1_isSymCS_1_isSymTau_1_tDt_0.0010_.txt';
% filename='tend_50_P1_0.00_0.99_100_P2_0.50_0.50_1_P3_0.00_0.00_1_T1_0.31_0.31_1_T3_0.31_0.31_1_E1_0.04_0.04_1_E3_0.04_0.04_1_isSymCS_1_isSymTau_1_tDt_0.0010_.txt';
% filename='tend_50_P1_0.00_0.99_100_P2_0.50_0.50_1_P3_0.00_0.00_1_T1_0.46_0.46_1_T3_0.46_0.46_1_E1_0.04_0.04_1_E3_0.04_0.04_1_isSymCS_1_isSymTau_1_tDt_0.0010_.txt';

ip = fopen(filename,'rt');          
if (ip < 0) 
    error('could not open file');   
end;

descp=regexp(filename, '_', 'split');
tEnd=str2double(descp(1,2));
phi1Begin=str2double(descp(1,4));phi1End=str2double(descp(1,5));phi1Len=str2double(descp(1,6));
phi2Begin=str2double(descp(1,8));phi2End=str2double(descp(1,9));phi2Len=str2double(descp(1,10));
phi3Begin=str2double(descp(1,12));phi3End=str2double(descp(1,13));phi3Len=str2double(descp(1,14));
tau1Begin=str2double(descp(1,16));tau1End=str2double(descp(1,17));tau1Len=str2double(descp(1,18));
tau3Begin=str2double(descp(1,20));tau3End=str2double(descp(1,21));tau3Len=str2double(descp(1,22));
epsi1Begin=str2double(descp(1,24));epsi1End=str2double(descp(1,25));epsi1Len=str2double(descp(1,26));
epsi3Begin=str2double(descp(1,28));epsi3End=str2double(descp(1,29));epsi3Len=str2double(descp(1,30));
isSymCS=str2double(descp(1,32));isSymTau=str2double(descp(1,34));

lines=textscan(ip,'%s',2*phi1Len+3,'delimiter','\n');
fclose(ip);

if (isSymCS==1) && (isSymTau==1)
    phi1_neuron_lin=linspace(phi1Begin,phi1End,phi1Len);
    
    zeroValueCnt=zeros(tau1Len-1,epsi1Len-1);    
    
    index=2;
    for i=1:1:phi1Len  % Phi1
        tokensSpTime=regexp(lines{1,1}{index,1}, ';', 'split');
        tokensEvoPhi1=regexp(lines{1,1}{index+phi1Len+2,1}, ';', 'split');
        index=index+1;

        col=size(tokensSpTime,2)-1;
        tmpSPT=NaN(1,col);
        tmpEP=NaN(1,col);
        
        for i2=1:1:col
            tmpSPT(1,i2)=str2double(tokensSpTime{i2});
            tmpEP(1,i2)=str2double(tokensEvoPhi1{i2});
        end
        
        figure(1)
        hold on
        set(gca,'FontSize', 18);
        if (tmpEP(1,i2)==0)          
            plot(tmpSPT,tmpEP);
        else 
            plot(tmpSPT,tmpEP,'r');
        end
        axis([0 50 0 1]);
        xlabel('time','fontsize',20),ylabel('phi1','fontsize',20),title('Evolution of phi1','fontsize',20);        
    end
else
    'Not yet implement'
end;


return;