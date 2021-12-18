% Call this function to plot the figures:
% Figure 1, COVID-19 transmission dynamics in relation to the introduction of NPIs,
% derived from the full model.
% Figure 2, Comparisons between actual values and model projections for two indicators of
% reduced efficacy of tracing and testing.

function [] = plot_HK_dynamics_soc11_1_sliding(ps, sys_par, par)
%% model parameters
names = ps.Properties.VariableNames
%generate 100 random number 
no = height(ps);%the number of row in table ps
a = sys_par.burnIn;
%a = 600000;
b = no;
%b = 600000;
sample = 50;
r = round((b-a).*rand(sample,1)) + a;%size: 100*1

%[beta, Lambda, tau, Tc, Tqr, Tqr1]; 

beta = ps.beta(r);
tau = ps.tau(r);
Tqr = ps.Tqr(r);
Tqr1 = ps.Tqr(r);
q = ps.q(r);
sigma = ps.sigma(r);
Tc = ps.Tc(r);
ct_b = ps.ct_b(r);
thld = ps.thld(r);
cr = ps.cr(r);
dur = ps.dur(r);
rp_d = ps.rp_d(r);
thld_d = ps.thld_d(r);
p2 = ps.p2(r);
sd_t2 = ps.sd_t2(r);
rp = ps.rp(r);
alpha = ps.alpha(r);
sd_t3 = ps.sd_t3(r);
delratio = ps.delratio(r);
deltrans = ps.deltrans(r);
delratio2 = ps.delratio2(r);
deltrans2 = ps.deltrans2(r);
sd_fm = ps.sd_fm(r);
maxdelay = ps.maxdelay(r);
sd_t1f = ps.sd_t1f(r);


%par.p1 = 0.8;
S2 = par.N2;
E21 = 2; I21 = 0; R21 = 0; QE21 = 0; QI21 = 0; H21 = 0;
E22 = 0; I22 = 0; R22 = 0; QE22 = 0; QI22 = 0; HE22 = 0; HI22 = 0; H22 = 0; 
CHE22 = 0; CHI22 = 0; CHR22 = 0; CH22 = 0; CI = 0; 
INC22 = 0; QINC22 = 0; HINC22 = 0; CHINC22 = 0;


totalt = 90;
%% current situation
par.sc=0;
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % add the variable, check ln61
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp(i,:) = detected_imp(1:end-1);
    local(i,:) = detected_local(1:end-1);
    localnl(i,:) = detected_local_nolink(1:end-1); 
    udet_imp(i,:) = undetected_imp(1:end-1);
    local_delay(i,:) = detected_delay(1:end-1);
    true_local(i,:) = total_local(1:end-1);
end

par.sc=1; % without any NPIs
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % check ln74
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp2(i,:) = detected_imp(1:end-1);
    local2(i,:) = detected_local(1:end-1);
    localnl2(i,:) = detected_local_nolink(1:end-1); 
    udet_imp2(i,:) = undetected_imp(1:end-1);
    true_local2(i,:) = total_local(1:end-1);
end

par.sc=2; % without 2nd 3rd and 4th social distancing tightening
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % go change the ode45 
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp3(i,:) = detected_imp(1:end-1);
    local3(i,:) = detected_local(1:end-1);
    localnl3(i,:) = detected_local_nolink(1:end-1); 
    udet_imp3(i,:) = undetected_imp(1:end-1);
    true_local3(i,:) = total_local(1:end-1);
end

par.sc=3; % without 3rd and 4th social distancing tightening
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % go change the ode45 
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp4(i,:) = detected_imp(1:end-1);
    local4(i,:) = detected_local(1:end-1);
    localnl4(i,:) = detected_local_nolink(1:end-1); 
    udet_imp4(i,:) = undetected_imp(1:end-1);
    true_local4(i,:) = total_local(1:end-1);
end

par.sc=4; % without 4th social distancing tightening and without isolation boosting
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % go change the ode45 
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp5(i,:) = detected_imp(1:end-1);
    local5(i,:) = detected_local(1:end-1);
    localnl5(i,:) = detected_local_nolink(1:end-1); 
    udet_imp5(i,:) = undetected_imp(1:end-1);
    true_local5(i,:) = total_local(1:end-1);
end

par.sc=5; % without 4th social distancing tightening but with isolation boosting
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % go change the ode45 
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp6(i,:) = detected_imp(1:end-1);
    local6(i,:) = detected_local(1:end-1);
    localnl6(i,:) = detected_local_nolink(1:end-1); 
    udet_imp6(i,:) = undetected_imp(1:end-1);
    true_local6(i,:) = total_local(1:end-1);
end

par.sc=6; % all measures without contact tracing boosting
for  i=1:length(r)
    % single strain
    %[beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
    % model 7 will include sd_t2(i)
    par.p = [beta(i) tau(i) Tqr(i) Tqr1(i) q(i) sigma(i) Tc(i) ct_b(i) thld(i) cr(i) dur(i) rp_d(i) thld_d(i) p2(i) sd_t2(i) rp(i) alpha(i) sd_t3(i) delratio(i) deltrans(i) delratio2(i) deltrans2(i) sd_fm(i) maxdelay(i) sd_t1f(i)]; % go change the ode45 
    [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(i);
    imp7(i,:) = detected_imp(1:end-1);
    local7(i,:) = detected_local(1:end-1);
    localnl7(i,:) = detected_local_nolink(1:end-1); 
    udet_imp7(i,:) = undetected_imp(1:end-1);
    true_local7(i,:) = total_local(1:end-1);
end



%% plot the dynamics of cases
%load('virus_hk07.mat');
virus_hk = load('virus_hk07.mat');
virus = virus_hk.virus_hk07;
case_dd = load('delay.mat');
X = case_dd.X;
figure();clf
subplot(3,1,1);
hold on

totaltime = totalt;
t = [1:totaltime+1];
tt = t(1:end-2);
n = length(tt);
current_t = 62; % 62-> 15 Aug

%% with all restrictions
local_total = local+localnl;
local_delay_mean = mean(local_delay);
local_total_mean = mean(local_total);
% for revision
var(sum(local_total(1:50,18:24),2))^0.5

local_total_bound = prctile(local_total,[2.5 97.5]);
XX = [tt(1),  tt,  tt(n),  fliplr(tt)];
YY = [local_total_bound(1,1), local_total_bound(2,:), local_total_bound(1,n), fliplr(local_total_bound(1,:))];
h1=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
set(h1,'facealpha',.5);
b1 = plot(tt(1:current_t),local_total_mean(1:current_t),'r-','linewidth',2.5) % day 1 to day 43 (29 July)
b5 = plot(tt(1:current_t),virus(tt(1:current_t),3),'ko','MarkerFaceColor','r');


current_t = current_t;
%local_total_lk = local;
%local_total_lk_mean = mean(local_total_lk);
%local_total_lk_bound = prctile(local_total_lk,[2.5 97.5]);
%b1 = plot(tt(1:current_t),local_total_lk_mean(1:current_t),'g-','linewidth',2.5) % day 1 to day 43 (29 July)
% t=25 => 11 July
rp = round(rp);
t0 = 25; %No NPIs
local_total2 = local2+localnl2;
local_total_mean2 = mean(local_total2);
local_total_bound2 = prctile(local_total2,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_total_bound2(1,t0), local_total_bound2(2,t0:end), local_total_bound2(1,n), fliplr(local_total_bound2(1,t0:end))];
h2=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
b2 = plot(tt(t0:end),local_total_mean2(t0:end),'r:','linewidth',2.5);
set(h2,'facealpha',.5);

t0 = 25; %T1=25
local_total3 = local3+localnl3;
local_total_mean3 = mean(local_total3);
local_total_bound3 = prctile(local_total3,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_total_bound3(1,t0), local_total_bound3(2,t0:end), local_total_bound3(1,n), fliplr(local_total_bound3(1,t0:end))];
h4=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
b4 = plot(tt(t0:end),local_total_mean3(t0:end),'r:','linewidth',2.5);
set(h4,'facealpha',.5)

t0 = 33; %T1,2
local_total4 = local4+localnl4;
local_total_mean4 = mean(local_total4);
local_total_bound4 = prctile(local_total4,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_total_bound4(1,t0), local_total_bound4(2,t0:end), local_total_bound4(1,n), fliplr(local_total_bound4(1,t0:end))];
h6=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
b6 = plot(tt(t0:end),local_total_mean4(t0:end),'r:','linewidth',2.5)
set(h6,'facealpha',.5);

%t0 = 43 + rp;
t0 = 35;  %T1,2,3
local_total5 = local5+localnl5;
local_total_mean5 = mean(local_total5);
local_total_bound5 = prctile(local_total5,[2.5 97.5]);%b7 = plot(tt(43:end),local_total_mean5(43:end),'g:','linewidth',3)

XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_total_bound5(1,t0), local_total_bound5(2,t0:end), local_total_bound5(1,n), fliplr(local_total_bound5(1,t0:end))];
%h7=fill(XX,YY,[228/255 237/255 218/255],'Linestyle','none');
h7=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
b7 = plot(tt(t0:end),local_total_mean5(t0:end),'r:','linewidth',2.5)
set(h7,'facealpha',.5);

t0 = 37;
local_total7 = local7+localnl7;
local_total_mean7 = mean(local_total7);
local_total_bound7 = prctile(local_total7,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_total_bound7(1,t0), local_total_bound7(2,t0:end), local_total_bound7(1,n), fliplr(local_total_bound7(1,t0:end))];
%h8=fill(XX,YY,[228/255 237/255 218/255],'Linestyle','none');
%b8 = plot(tt(t0:end),local_total_mean7(t0:end),'color',[101/255 171/255 9/255],'linestyle',':','linewidth',2.5)
h8=fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
b8 = plot(tt(t0:end),local_total_mean7(t0:end),'r:','linewidth',2.5);
set(h8,'facealpha',.5);

%% 
set(gcf,'position',[400,350,900,600]);
ymax = 500;
ymaxl = 300;
ll = 0:0.1:ymaxl;
%plot(3*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
%plot(25*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
%plot(29*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
%plot(31*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
%plot(37*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
%plot(43*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
%plot(45*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);

% title('Local cases in Hong Kong');
xlim([1,60]);
ylim([0,ymaxl]);
%xlabel('Date','fontsize',15);
ylabel('Daily reported local cases','fontsize',15);
%%lgd = legend('CI:95%','mean','mean','Observed','Location','best');
lgd = legend([b5 b1 b2],{'Observed','Predicted, with all NPIs','With T1; T1,2; T1,2,3; T1-4'});
%lgd.FontSize = 15;
set(gca,'FontSize',15);
date0 = ([1 6 11 16 21 26 31 36 41 46 51 56 60]);
% date = ({'06/17','06/22','06/27','07/02','07/07','07/12','07/17','07/22','07/27','08/01','08/06','08/11','08/15'});
date = {};
set(gca,'xtick',date0,'XTickLabel',date)
xtickangle(90);
%ax = gca;
%ax.XGrid = 'on';
%ax.YGrid = 'on';
box on;

%% Figure2 Recovered true epidemic dynamics
subplot(3,1,2);
hold;
% What are those 6 i s?
for i = 1:6
if i == 1
    true_local_mean = mean(true_local);
    true_local_bound = prctile(true_local,[2.5 97.5]);
    t0 = 1;
end
if i == 2 % without NPIs
    true_local_mean = mean(true_local2);
    true_local_bound = prctile(true_local2,[2.5 97.5]);
    t0 = 25-2;
end
if i == 3 % without T2-T4
    true_local_mean = mean(true_local3);
    true_local_bound = prctile(true_local3,[2.5 97.5]);
    t0 = 29;
end
if i == 4  % without T3-T4
    true_local_mean = mean(true_local4);
    true_local_bound = prctile(true_local4,[2.5 97.5]);
    t0 = 29;
end
if i == 5  % without T4 without boosting
    true_local_mean = mean(true_local5);
    true_local_bound = prctile(true_local5,[2.5 97.5]);
    %t0 = 43 + rp;
    t0 = 37;
end
%if i == 6
%    true_local_mean = mean(true_local6);  % without T4 with boosting
%    true_local_bound = prctile(true_local6,[2.5 97.5]);
%    %t0 = 43 + rp;
%    t0 = 39;
%end
if i == 6
    true_local_mean = mean(true_local7);  % all measures without contact tracing boosting
    true_local_bound = prctile(true_local7,[2.5 97.5]);
    %t0 = 43 + rp;
    t0 = 39;
end

%b8 = plot(tt(t0:end),local_total_mean6(t0:end),'g:','linewidth',2.2)
%XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
%YY = [local_total_bound6(1,t0), local_total_bound6(2,t0:end), local_total_bound6(1,n), fliplr(local_total_bound6(1,t0:end))];
%h8=fill(XX,YY,[228/255 237/255 218/255],'Linestyle','none');

    XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
    YY = [true_local_bound(1,t0), true_local_bound(2,t0:end), true_local_bound(1,n), fliplr(true_local_bound(1,t0:end))];
if i == 1
    ht = fill(XX,YY,[208/255 227/255 252/255],'Linestyle','none');
    b1 = plot(tt(t0:current_t),true_local_mean(t0:current_t),'b-','linewidth',2.5) % day 1 to day 43 (29 July)
    %b1 = plot(tt(t0:end),true_local_mean(t0:end),'b-','linewidth',2.5) % day 1 to day 43 (29 July)

elseif i < 6 
    ht = fill(XX,YY,[208/255 227/255 252/255],'Linestyle','none');
    b2 = plot(tt(t0:current_t),true_local_mean(t0:current_t),'b:','linewidth',2.5) % day 1 to day 43 (29 July)
    %b2 = plot(tt(t0:end),true_local_mean(t0:end),'b:','linewidth',2.5) % day 1 to day 43 (29 July)

else
    %ht = fill(XX,YY,[228/255 237/255 218/255],'Linestyle','none');
    %b3 = plot(tt(t0:current_t),true_local_mean(t0:current_t),'color',[101/255 171/255 9/255],'linestyle',':','linewidth',2.5) % day 1 to day 43 (29 July)
    ht = fill(XX,YY,[208/255 227/255 252/255],'Linestyle','none');
    b3 = plot(tt(t0:current_t),true_local_mean(t0:current_t),'b:','linewidth',2.5) % day 1 to day 43 (29 July)
    %b3 = plot(tt(t0:end),true_local_mean(t0:end),'b:','linewidth',2.5) % day 1 to day 43 (29 July)

end
set(ht,'facealpha',.5);
end
plot(3*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(25*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(29*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(31*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(37*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(43*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(46*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
xlim([1,60]);
ylim([0,300]);
lgd = legend([b1 b2],{'Predicted, with all NPIs','With T1; T1,2; T,1,2,3; T1-4'});
%xlabel('Date','fontsize',15);
ylabel('Daily local infections','fontsize',15);
date0 = ([1 6 11 16 21 26 31 36 41 46 51 56 60]);
%date = ({'06/17','06/22','06/27','07/02','07/07','07/12','07/17','07/22','07/27','08/01','08/06','08/11','08/15'});
date = {};
set(gca,'FontSize',15);
set(gca,'xtick',date0,'XTickLabel',date)
%xtickangle(90);
box on;

%% Figure3
subplot(3,1,3);
hold;
%Should pass the same variables#######20210203
%plot_Rt2curves_import_HK9_10(ps, sys_par, par, r);
%plot_Rt2curves_import_HK9_10_fm(ps, sys_par, par, r);
plot_Rt2curves_import_HK11_1_fm(ps, sys_par, par, r);

%% Figure
% Proportion of EpiLink vs Without EpiLink
figure;
subplot(2,1,1);
hold;

t0 = 17;
t00 = 20;
link_prop = 100 - local./(local+localnl).*100;
link_prop_mean = mean(link_prop);
true_link_bound = prctile(link_prop,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [true_link_bound(1,t0), true_link_bound(2,t0:end), true_link_bound(1,n), fliplr(true_link_bound(1,t0:end))];
hl = fill(XX,YY,[208/255 227/255 252/255],'Linestyle','none');
set(hl,'facealpha',.5);
b10 = plot(tt,link_prop_mean,'b-','linewidth',2);
% col4: cases with epi-link; col5: cases without epi-link
% actual_prop: Percentage of cases without epi-link (used to represent contract tracing inefficiency)
actual_prop = 100 - virus(19:current_t,4)./(virus(19:current_t,4)+virus(19:current_t,5)).*100;
window = 5;
for i=1:length(virus)
%  delay(i) = sum(local_delay_mean(i+15-6:i+15))/7*100;
  if i >= round(window/2) && i <= length(virus) - round(window/2) + 1
      actual_prop_avg(i) = sum(virus(i-(round(window/2)-1):i+(round(window/2)-1),4))./sum(virus(i-(round(window/2)-1):i+(round(window/2)-1),4)+virus(i-(round(window/2)-1):i+(round(window/2)-1),5));
  end 
  if i < round(window/2)
      actual_prop_avg(i) = sum(virus(1:i+(round(window/2)-1),4))./sum(virus(1:i+(round(window/2)-1),4)+virus(1:i+(round(window/2)-1),5));
  end
  if i > length(virus) - round(window/2) + 1
      actual_prop_avg(i) = sum(virus(i-(round(window/2)-1):end,4))./sum(virus(i-(round(window/2)-1):end,4)+virus(i-(round(window/2)-1):end,5));
  end
end
actual_prop_avg = (1-actual_prop_avg).*100;
%actual_prop_avg = movmean(actual_prop,7);
b11_dot = plot(tt(21:current_t-1),actual_prop(3:end-1),'r.','MarkerSize',12);
b11 = plot(tt(21:current_t-1),actual_prop_avg(21:61),'-','LineWidth',2,'Color',[1 0 0]);
%b11 = plot(tt(21:current_t-1),actual_prop_avg(3:end-1),'-','LineWidth',2,'Color',[1 0 0]);
%b11 = plot(tt(19:current_t),actual_prop_avg,'-','LineWidth',1.5,'MarkerEdgeColor',[120/255 171/255 48/255]);

ymaxl = 100;
ll = 0:0.1:ymaxl;
plot(3*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(25*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(29*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(31*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(37*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(43*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(46*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);

lgd = legend([b10 b11_dot],{'Model output','Observerd'});
xlim([t00-1,61-1]);
ylim([0,60]);
%xlabel('Date','fontsize',15);
ylabel('Percentage without epi-link','fontsize',15);
date0 = t0+2:2:60;
%date0 = ([1 6 11 16 21 26 31 36 41 46 51 56 61]);
%date = ({'06/17','06/22','06/27','07/02','07/07','07/12','07/17','07/22','07/27','08/01','08/06','08/11','08/16'});
date = {};
set(gca,'FontSize',15);
set(gca,'xtick',date0,'XTickLabel',date)
xtickangle(90);
box on;


%% Figure
%% plot the percentage of cases confirmed after symptom onset 
%figure();clf
%hold;
subplot(2,1,2);
hold;
t00 = 20;
%DelayP = 1 - local_delay;
%Asym = 0.17;
%DelayAsym = (1-DelayP).*0.17;
%local_delay = local_delay + DelayAsym;
local_delay_mean = mean(local_delay);
local_delay_bound = prctile(local_delay,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_delay_bound(1,t0), local_delay_bound(2,t0:end), local_delay_bound(1,n), fliplr(local_delay_bound(1,t0:end))]*100;
%h=fill(XX,YY,[208/255 227/255 252/255],'Linestyle','none');
%set(h,'facealpha',.5);
loc_p = local_delay_mean*100;
loc_p_list = [0 0 0];
window = 5;
for i = window:length(loc_p)-(window-1)
  loc_p_list(i) = mean(loc_p(i-(window-1):i+(window-1)));
end
loc_p_list = [loc_p_list 0 0 0];
%b10 = plot(tt(t00-3:end),loc_p(t00-3:end),'b-','linewidth',2);
b10 = plot(tt(t00-3:end),loc_p_list(t00-(window-1):end),'b-','linewidth',2) %sliding winder

% make an error for moving average
loc_p = local_delay*100;
%loc_p_list = [0 0 0];
loc_p_list = zeros(sample,length(loc_p(1,:)));
window = 5;
for i = window:length(loc_p(1,:))-(window-1)
  loc_p_list(:,i) = mean(loc_p(:,i-(window-1):i+(window-1)),2);
end

local_delay_bound = prctile(loc_p_list,[2.5 97.5]);
XX = [tt(t0),  tt(t0:end),  tt(n),  fliplr(tt(t0:end))];
YY = [local_delay_bound(1,t0), local_delay_bound(2,t0:end), local_delay_bound(1,n), fliplr(local_delay_bound(1,t0:end))];
h=fill(XX,YY,[208/255 227/255 252/255],'Linestyle','none');
set(h,'facealpha',.5);
%loc_p_list = [loc_p_list 0 0 0];
%b10 = plot(tt(t00-3:end),loc_p_list(t00-(window-1):end),'b-','linewidth',2) %sliding winder

window = 3;
window_2 = floor(window/2);
for i=1:51
%  delay(i) = sum(local_delay_mean(i+15-6:i+15))/7*100;
 %delay(i) = sum(local_delay_mean(i+t00-(window-2):i+t00+(window-2)).*(local(i+t00-(window-2):i+t00+(window-2))+localnl(i+t00-(window-2):i+t00+(window-2)))/(sum((local(i+t00-(window-2):i+t00+(window-2))+localnl(i+t00-(window-2):i+t00+(window-2))))))*100;
  delay(i) = sum(local_delay_mean(i+t00-(window_2):i+t00+(window_2)).*(local(i+t00-(window_2):i+t00+(window_2))+localnl(i+t00-(window_2):i+t00+(window_2)))/(sum((local(i+t00-(window_2):i+t00+(window_2))+localnl(i+t00-(window_2):i+t00+(window_2))))))*100;
 if i<2
      %delay(i) = sum(local_delay_mean(t00:i+t00+(window-2)).*(local(t00:i+t00+(window-2))+localnl(t00:i+t00+(window-2)))/(sum((local(t00:i+t00+(window-2))+localnl(t00:i+t00+(window-2))))))*100;
      delay(i) = sum(local_delay_mean(t00:i+t00+floor(window/2)).*(local(t00:i+t00+floor(window/2))+localnl(t00:i+t00+floor(window/2)))/(sum((local(t00:i+t00+floor(window/2))+localnl(t00:i+t00+floor(window/2))))))*100;
 end
end
%plot(t0:66,delay)
actual_delay_avg = movmean(X,3)*100;
b11_dot = plot(tt(21:current_t-2),X(3:end)*100,'r.','MarkerSize',12);
b11 = plot(tt(21:current_t-2),virus(21:current_t-2,6),'-','LineWidth',2,'Color',[1 0 0]);

ymaxl = 100;
ll = 0:0.1:ymaxl;
plot(3*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(25*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(29*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(31*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(37*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(43*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);
plot(46*ones(length(ll),1), 0:0.1:ymaxl,'color',[0.3 0.3 0.3],'linestyle','-','linewidth',1);

lgd = legend([b10 b11_dot],{'Model output','Observerd'});
xlim([t00-1,61-1]);
ylim([30 100]);
xlabel('Date','fontsize',15);
ylabel('Percentage with confirmation delay','fontsize',15);
date0 = t0+2:2:60;
date = ({'05/07','07/07','09/07','11/07','13/07','15/07','17/07','19/07','21/07','23/07','25/07','27/07','29/07','31/07','02/08','04/08','06/08','08/08','10/08','12/08','14/08'});
%date0 = ([1 6 11 16 21 26 31 36 41 46 51 56 61]);
%date = ({'06/17','06/22','06/27','07/02','07/07','07/12','07/17','07/22','07/27','08/01','08/06','08/11','08/16'});
set(gca,'FontSize',15);
set(gca,'xtick',date0,'XTickLabel',date)
xtickangle(90);
box on;



function [detected_imp detected_local detected_local_nolink undetected_imp detected_delay total_local] = get_sir_import(sid)
x0 = [S2; E21; I21; R21; QE21; QI21; H21; E22; I22; R22; QE22; QI22; HE22; HI22; H22; CHE22; CHI22; CHR22; CH22; CI; INC22; QINC22; HINC22; CHINC22];
starttime = 0;
totaltime = 90;
[Sol] = ode45(@odef_3pop_imp_soc11_1,[starttime totaltime] ,x0,[],par);

tx = linspace( starttime+1, starttime+totaltime, totaltime-starttime);
x = deval(Sol, tx)';

%group2 for Other provinces
s2 = x(:,1);

%Imported
e21 = x(:,2);
i21 = x(:,3);
r21 = x(:,4);
qe21 = x(:,5);
qi21 = x(:,6);
h21 = x(:,7);

e22 = x(:,8);    %exposed 
i22 = x(:,9);    %infectious
r22 = x(:,10);   %recovered
qe22 = x(:,11);  %quarantined exposed
qi22 = x(:,12);  %quarantined infectious
he22 = x(:,13);  %confirmed exposed --- Local with Epi-Link (before symptom onset)
hi22 = x(:,14);  %confirmed infectious --- Local with Epi-Link (after symptom onset)
h22 = x(:,15);   %confirmed --- Local without Epi-Link
che22 = x(:,16);  %confirmed exposed --- Local with Epi-Link (before symptom onset)
chi22 = x(:,17);  %confirmed infectious --- Local with Epi-Link (before and after symptom onset)
chr22 = x(:,18);  %confirmed recovered --- Local with Epi-Link (after symptom onset; back traced)
ch22 = x(:,19);   %confirmed --- Local without Epi-Link
ci = x(:,20);    %cumulative cases

inc22 = x(:,21);
qinc22 = x(:,22);
hinc22 = x(:,23); 
chinc22 = x(:,24);


%% Get Tqr

%get parameters
lambda1 = 1/par.p(4); % tqre 
lambda3 = 1/par.p(5); % tqri 
lambda2 = 1/par.p(2); % tau 
infperiod = par.p(7) - par.p(3);
presymp = par.p(6);

%detected_imp_ratio
i1_list_avg = [5.0    4.0    3.8    8.2   10.6   10.4   13.0   13.0    7.2    4.4    4.8    2.4    7.4    9.0    9.6   11.0   12.2    9.8    9.0    9.0    8.4    8.0    7.2    7.8    9.0    9.0    8.8    7.2    7.2    5.8    9.2    9.6    9.4   10.0   10.8    7.4    7.4   11.8   12.0   13.3];

factor = 0.1*(1+((1:totaltime)/30)*0.18);
exempt_prop = 1;
for tid = 1:totaltime
if tid > length(i1_list_avg)
   %if round((t - floor(t))*1000) < 10
   e1(tid) = mean(i1_list_avg)*exempt_prop;
   i1(tid) = mean(i1_list_avg)*factor(tid)*exempt_prop;
   %end 
else
   %if round((t - floor(t))*1000) < 10
   e1(tid) = i1_list_avg(tid);  
   i1(tid) = i1_list_avg(tid)*factor(tid);   
   %end
end
end

%detected_imp = [0; diff(h21)];
detected_imp = e1';
undetected_imp = i1';
%detected_imp_ratio = qi21./(i21+qi21);
%detected_imp_shift = [detected_imp_ratio(3:end); detected_imp_ratio(end); detected_imp_ratio(end)];
%detected_imp = ((1/par.p(3)).*(e21+qe21)).*detected_imp_shift;

%detected_local_ratio
%step1. calculate daily newly qe2 and qi2
%newQe2_delay = lambda1*e22_d;
%newqi2_delay = lambda3.*i22_d; % ??P6

%step2. for Qe group, return how many of them will be reported d days later
%[ra qe2_before qi2_before]= getAsymProp( tqre, tqri, tau, Tc, sigma, dlay)
%Qe2_delay_Asym = (1-qe2_before).*newQe2_delay;

%reported_local link(+)
detected_local = [0; diff(che22 + chi22 + chr22)];
%reported_local link(-)
detected_local_nolink = [0; diff(ch22)];
%total quarantined
TotalCases = qe22+qi22+he22+hi22+h22;
TotalCases_M = qe22+qi22; %Total cases in mobility restriction


boost_day = 31;
del_durat = 5;
k = 1; %a scaling factor

%% t = 1-50;
%% time shift with stepwise change
for dd = 1:totaltime
maxdelay_dd = maxdelay(sid);
if dd == 29
  xxx = 1;
end
if dd >= boost_day + deltrans(sid)
 delrat = delratio(sid);
 maxdelay1 = delrat*maxdelay(sid);
 if dd-(boost_day+deltrans(sid)) <= del_durat  
     %maxdelay_dd = mean(maxdelay) - ((mean(maxdelay)-maxdelay1).*1./(1+exp(-0.5.*((dd-((boost_day+mean(deltrans))))))));
     delta_t = dd-(boost_day+deltrans(sid));
     maxdelay_dd = maxdelay(sid)-((maxdelay(sid)-maxdelay1).*(1 - 1./(exp(0.3*delta_t))));
 else
     maxdelay_dd = maxdelay1;
 end
end

if (dd == 49)
   xxx = 1;
end
if dd >= 46 + deltrans2(sid) % Consider the effect of the second case isolation improvement on 1 August 
 delrat = delratio2(sid);
 maxdelay2 = delrat*maxdelay_dd;
 if dd-(46+mean(deltrans)) <= del_durat
     %maxdelay_dd = maxdelay_dd - ((maxdelay_dd-maxdelay2).*1./(1+exp(-0.5.*((dd-((46+mean(deltrans2))))))));
     delta_t2 = dd-(46+deltrans(sid));
     maxdelay_dd = maxdelay_dd - ((maxdelay_dd-maxdelay2).*(1 - 1./(exp(0.3*delta_t2))));
 else
     maxdelay_dd = maxdelay2;
 end
end

%% Case isolation delay
if (TotalCases(dd))<thld_d(sid);
 b = 0;
else
 % exponential    
 b = maxdelay_dd*(1-exp(-((TotalCases(dd))-thld_d(sid))*rp_d(sid))); 
 %b = maxdelay*(1-exp(-((TotalCases)-thld_d)*rp_d)); 
end
b = 0.5 + b; % 0.5 is the minimum testing time


%% Contact tracing delay
if dd>=16
if (TotalCases_M(dd))<thld(sid);
 a = 0;
else
 % exponential
 a = maxdelay_dd*(1-exp(-((TotalCases_M(dd))-thld(sid))*ct_b(sid)));
end
else
     a = 0;
end
a = a + b*0.5; % Case isolation delay will cause contact tracing delay 
tqr = Tqr(sid)+a;
%dlay(dd) = tqr + b - tau(sid); % Confirmation time minus tau
dlay(dd) = tqr + b + tau(sid); % Confirmation time
end

%pro = presymp/infperiod;
pro2 = exp(-((tau(sid)+sigma(sid))^(-1)).*dlay)';


%Asymptomatic proportion
%pro = presymp/infperiod;
%maxdlay = dlay;
%pro = presymp./(presymp+maxdlay');
%detected_delay = [0; diff([0; 0; 0; chi22(1:end-3)].*(1-pro)+chr22)./diff(chi22+che22+chr22)];
detected_delay = [0; diff(chi22.*(1-pro2)+chr22)./diff(chi22+che22+chr22)];

total_local = [0; diff(ci)];
%detected_delay = [0; diff((chi22-chinc22)+chr22)./diff(chi22+che22+chr22)];


end


end

