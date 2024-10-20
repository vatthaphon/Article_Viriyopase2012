function draw_all_line_separation
% clc
% clear all
% close all

b=3;


lw=10;
%% 1,2,3 non driven
%%%%%% 1-tau < phi2 < 1 %%%%%%
A=exp(b)-1;

epsi_l=linspace(0.01,0.21,1000);

alpha_epsi=alpha(epsi_l, b);
alpha_2epsi=alpha(2.*epsi_l, b);
A_epsi=A.*alpha_epsi;
A_2epsi=A.*alpha_2epsi;
theta_c1=theta_c(epsi_l, b);

G=A_2epsi+A_epsi;
H=A_epsi-A_2epsi;
F=alpha_epsi+A_epsi+A_2epsi-alpha_2epsi;

tau2=(theta_c1-1+F./G)./(1-H./G);

hold on
z=ones(1,size(epsi_l,2));
plot3(epsi_l,tau2,z,'r','LineWidth',lw);

%% 1,3 driven + 2 driven
%%%%%% tau=phi2 %%%%%%%
epsi_lin=linspace(0.01,0.21,1000);

tau=theta_c(epsi_lin,b)/2;
hold on
z=ones(1,size(epsi_lin,2));
plot3(epsi_lin,tau,z,'k','LineWidth',lw);
plot(epsi_lin,tau,'k','LineWidth',lw);

%% 1,3 non driven + 2 driven
%%%%%% 0<phi2<tau %%%%%%%
epsi_lin=linspace(0.01,0.21,1000);
lb_tau=(1-alpha(epsi_lin,b))./2;

hold on
z=ones(1,size(epsi_lin,2));
plot3(epsi_lin,lb_tau,z,'k','LineWidth',lw);
plot(epsi_lin,lb_tau,'k','LineWidth',lw);

%% 1,3 driven + 2 non driven
%%%%%% phi2=tau %%%%%%%
epsi_lin=linspace(0.01,0.21,1000);

tau1=(1-alpha(2*epsi_lin,b)-theta_c(epsi_lin,b))./(2*beta(2*epsi_lin,b));
tau2=theta_c(2*epsi_lin,b)./2;

cols=size(tau1,2);
for i=1:1:cols
    if (tau1(1,i)<=tau2(1,i))
        less_tau(1,i)=tau1(1,i);
    else
        less_tau(1,i)=tau2(1,i);
    end
end

hold on
z=ones(1,size(epsi_lin,2));
plot3(epsi_lin,less_tau,z,'k','LineWidth',lw);
plot(epsi_lin,less_tau,'k','LineWidth',lw);

%%%%%% 1-tau<phi2<1 %%%%%%%
epsiMTZero_lin=linspace(0.01,(1/(2*b))*log(2),1000);
epsiLTZero_lin=linspace((1/(2*b))*log(2),0.21,1000);
A=exp(b)-1;
D_MTZero=exp(2.*epsiMTZero_lin.*b);
D_LTZero=exp(2.*epsiLTZero_lin.*b);

tau2_MTZero=(A-D_MTZero+1)./(2.*A);

tau2_LTZero=(A-D_LTZero+1)./(2.*A);

hold on
z=ones(1,size(epsiMTZero_lin,2));
plot3(epsiMTZero_lin,tau2_MTZero,z,'y','LineWidth',lw);
plot3(epsiLTZero_lin,tau2_LTZero,z,'k','LineWidth',lw);

% axis([0.01 0.21 0.01 0.49]);

end

