clc
clear all
close all

format long
global max_cs A_pot_BiPoo A_dep_BiPoo tau_pot tau_dep

%% Initialize
tend=50;

phi1=0.3;
phi2=0.0;
phi3=0.8;

% tau1=0.3;
% tau3=tau1;

tau1=0.05;
tau3=0.375;

% epsi12=0.15;
% epsi21=epsi12;
% epsi23=epsi12;
% epsi32=epsi12;

epsi12=0.0;
epsi21=0.21;
epsi23=0.21;
epsi32=0.21;

max_cs=0.21;
% A_pot_BiPoo=0.777;
% A_dep_BiPoo=0.273;
A_pot_BiPoo=0;
A_dep_BiPoo=0;
tau_pot=16.8; % ms
tau_dep=33.7; % ms

dt=0.0001;      % 0.000001<0.0001
figOption=3;    % 0=no figure, 1=phase info., 2=raster plot, 3=voltage

%% Execution
[rec_SpTi1 rec_SpTi2 rec_SpTi3 rec_12_arr rec_21_arr rec_23_arr rec_32_arr o_w_12 o_w_21 o_w_23 o_w_32]...
=TimeBased_3MS_osc_with_Hebb_BiPoo(phi1, phi2, phi3, tau1, tau3, epsi12, epsi21, epsi32, epsi23, tend, figOption, dt);                        

maximize('all');

%% Cross checking
% w_12=check(rec_12_arr, rec_SpTi2, epsi12);
% disp(strcat(num2str(o_w_12), '==?', num2str(w_12)));

% w_21=check(rec_21_arr, rec_SpTi1, epsi21);
% disp(strcat(num2str(o_w_21), '==?', num2str(w_21)));
% 
% w_23=check(rec_23_arr, rec_SpTi3, epsi23);
% disp(strcat(num2str(o_w_23), '==?', num2str(w_23)));
% 
% w_32=check(rec_32_arr, rec_SpTi2, epsi32);
% disp(strcat(num2str(o_w_32), '==?', num2str(w_32)));


