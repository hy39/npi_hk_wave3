%% SEIR model with 2 populations
%% Correct the time for second tightening event
%% Timing of the social distancing tightening
%% Threshod for reporting should be the total hospitalized cases

function [xdot t1]= odef_3pop_imp11_1(t,x,par)
%group2 for Other provinces
s2 = x(1);
e21 = x(2); %Imported
i21 = x(3);
r21 = x(4);
qe21 = x(5);
qi21 = x(6);
h21 = x(7);

e22 = x(8); %Local clusters
i22 = x(9);
r22 = x(10);
qe22 = x(11);
qi22 = x(12);
he22 = x(13);
hi22 = x(14);
h22 =  x(15);
che22 = x(16);
chi22 = x(17);
chr22 = x(18);
ch22 =  x(19);
ci = x(20);

% Epidemiological parameters
% pa.p = [beta, tau, Tqr, Tqr1, q, sigma, Tc, ct_b];
beta = par.p(1);
tau = par.p(2);
Tqr = par.p(3);
Tqr1 = 0.05;
%Tqr1 = 0.02;
q = par.p(5);
sigma = par.p(6);
Tc = par.p(7);
ct_b = par.p(8);
thld = par.p(9);
cr = par.p(10);
dur = par.p(11);
rp_d = par.p(12);
thld_d = par.p(13);
p2 = par.p(14);
sd_t2 = par.p(15);
rp = par.p(16); %risk perception
alpha = par.p(17);
sd_t3 = par.p(18);
delratio = par.p(19);
deltrans = par.p(20);
delratio2 = par.p(21);
deltrans2 = par.p(22);
sd_fm = par.p(23);
maxdelay = par.p(24);
sd_t1f = par.p(25);
%ct_b = 0.05;
%cr = 1.5;
gamma = 1/(Tc-tau);

% Other parameters
N2 = par.N2;
p1 = par.p1;

tid = round(t)+1;
t1 = tid;

%% A simple SIR model without birth rates
e1 = 0;
i1 = 0;
i1_list = [8, 4, 3, 1, 3, 30, 16, 2, 14, 3, 1, 2, 4, 2, 28, 9, 5, 11, 8, 16, 5, 5, 8, 6, 12, 8, 11, 8, 5, 4, 8, 4, 25, 7, 3, 11, 8, 8, 7, 25, 4, 8, 5, 4, 3, 1, 1, 0, 5, 3, 4, 8, 2, 9, 2, 1, 1, 4, 2, 7, 4, 13, 1];
i1_list_avg = movmean(i1_list,5);
%i1_list_avg = [5.0    4.0    3.8    8.2   10.6   10.4   13.0   13.0    7.2    4.4    4.8    2.4    7.4    9.0    9.6   11.0   12.2    9.8    9.0    9.0    8.4    8.0    7.2    7.8    9.0    9.0    8.8    7.2    7.2    5.8    9.2    9.6    9.4   10.0   10.8    7.4    7.4   11.8   12.0   13.3];

%% Imported cases
exempt_prop = 1;
factor = 0.1*(1+(t/30)*0.18); % 10% of visitors are exempted; monthly increase rate 18%
if tid > length(i1_list_avg)
   %if round((t - floor(t))*1000) < 10
   e1 = mean(i1_list_avg)*exempt_prop;
   i1 = mean(i1_list_avg)*factor*exempt_prop;
   %end 
else
   %if round((t - floor(t))*1000) < 10
   e1 = i1_list_avg(tid);  
   i1 = i1_list_avg(tid)*factor;   
   %end
end

tqr1 = Tqr1;

%% Reporting delay
%maxdelay = 12.4; % 10 or 12.4 are better than 4
%% linear improvement
%maxdelay = 10;
%if t >= 37  % Consider the effect of the first case isolation improvement on 23 July
 %%maxdelay = 3;
% delay_a = delratio*maxdelay;
% delay_b = deltrans;
% newdelay1 = maxdelay + ((delay_a - maxdelay)/delay_b)*(t-37); 
% maxdelay = newdelay1;
% if newdelay1 < delay_a
%    maxdelay = delay_a;
% end
%end
%if t >= 46
 %delay2_b = deltrans2;
 %newdelay2 = maxdelay + ((delay2_a - maxdelay)/delay2_b)*(t-37); 
 %maxdelay = newdelay2;
 %if newdelay2 < delay2_a
 %   maxdelay = delay2_a;
 %end
%end

%% time shift with stepwise change
%maxdelay = 12.4;
del_durat = 5;
boost_day = 31;
if t >= boost_day + deltrans
 maxdelay1 = delratio*maxdelay;
 if t-(boost_day+deltrans) <= del_durat
     delta_t = t-(boost_day+deltrans);
     maxdelay = maxdelay-((maxdelay-maxdelay1).*(1 - 1./(exp(0.3*delta_t))));
 else
     maxdelay = maxdelay1;
 end
end

if t >= 46 + deltrans2 % Consider the effect of the second case isolation improvement on 1 August 
 maxdelay2 = delratio2*maxdelay;
 if t-(46+deltrans) <= del_durat
     delta_t2 = t-(46+deltrans);
     maxdelay = maxdelay-((maxdelay-maxdelay2).*(1 - 1./(exp(0.3*delta_t2))));
 else
     maxdelay = maxdelay2;
 end
end
 


%% Case isolation delay
TotalCases = qe22+qi22+he22+hi22+h22; 
TotalCases_M = qe22+qi22; %Total cases in mobility restriction
if t>=16
if (TotalCases)<thld_d;
    b = 0;
else
 % exponential    
% b = maxdelay*(1-exp(-((qe22+qi22)-thld_d)*rp_d)); %rp_d determines the quantile
b = maxdelay*(1-exp(-((TotalCases)-thld_d)*rp_d)); 
% b = maxdelay*(1-exp(-((he22+hi22+h22)-thld_d)*ct_b*rp_d)); 
% linear
 %b = 0 + ((qe22+qi22)-thld_d)*rp_d;
end
else
    b = 0;
end
dlay = 0.5 + b;

gamma2 = gamma;
%% Case isolation delay for non-traced; no contact tracing delay
%if  dlay > (Tc-tau) || t > boost_day
if  dlay > (Tc-tau)
   gamma2 = 1/(dlay+sigma);
end


%% Contact tracing delay

threshold = thld; % Test different threshold
%TotalCases = he22+hi22+h22;
%TotalCases = qe22 + qi22;
if t>=16
if (TotalCases_M)<threshold;
 a = 0;
else
 % linear    
 %a = 0 + ((qe22+qi22)-threshold)*ct_b;
 % exponential
 
 a = maxdelay*(1-exp(-((TotalCases_M)-threshold)*ct_b));
end
else
 a = 0;
end
a = a + dlay*0.5; % Case isolation delay will cause contact tracing delay 
tqr = Tqr+a;

    
if (t<3) 
    cont_ratio = 1;
end

if (t>=3 && t<25 + rp) % Relaxation  
    %logistic curve
    relax = 1+(cr-1)*1/(1+exp(alpha*((t-3)-dur)));
    cont_ratio = relax;
end
if (t>=25 + rp) % Scenario1 = Back to Strong contact tracing; Tightening 1
    T1 = 25 + rp;
    relax = 1+(cr-1)*1/(1+exp(alpha*((T1-3)-dur)));
    if (rp<0)
        if t <= 25
            tight1 = relax + ((1*sd_t1f-relax)/abs(rp))*(t-(25+rp));
        end
        if t > 25
            tight1 = 1*sd_t1f;
        end
    end
    if (rp>=0)
      tight1 = 1*sd_t1f;
    end
    cont_ratio = tight1;
end
if (t>=29 + rp) % Scenario2 = Tightening 2
    T2 = 29 + rp;
    if (rp<0)
        if t < 29
            tight2 = tight1 + ((sd_t2-tight1)/abs(rp))*(t-(29+rp));
        end
        if t >= 29
            tight2 = sd_t2;
        end
    end
    if (rp>=0)
      tight2 = sd_t2;
    end
    cont_ratio = tight2;
end

if (t >= 37) % 23 July, mandatory facial mask 
  tight_fm = sd_fm;
  cont_ratio = cont_ratio - tight_fm;
end

% when t=43, 29 July, social distancing begins
if (t>=43 + rp) % Scenario3 = Tightening 3
    % 
    if (rp<0)
        if t<43
            tight3 = tight2 + ((sd_t3-tight2)/abs(rp))*(t-(43+rp));
        end
        if t>=43
            tight3 = sd_t3;
        end
    end
    if (rp>=0)
        tight3 = sd_t3;
    end
    cont_ratio = tight3 - tight_fm; % add the effect of wearing mask
end


hosp_day = 7;% 14->9->7
%% Quarantined only infect local with link
beta = beta * cont_ratio;
s2_dot =  - beta*s2/N2*(i21 + i22) - q*beta*s2/N2*(qi21+qi22);
% Imported
e21_dot = e1 - (1/tau)*e21 - (1/tqr1)*e21;   
i21_dot = i1 + (1/tau)*e21 - gamma*i21 -  0*i21;
r21_dot = gamma*i21;
qe21_dot = (1/tqr1)*(e21) - (1/tau)*qe21 - (1/dlay)*qe21;
qi21_dot = 0*(i21) + (1/tau)*qe21 - (1/dlay)*qi21 - gamma*qi21 ;
h21_dot = (1/dlay)*(qe21 + qi21) + gamma*qi21;
%Local
%High R, Strong Q   
e22_dot = beta*s2/N2*(i21+i22+q*qi21+q*qi22) - (1/tau)*e22 - (1/tqr)*e22;
i22_dot = (1/tau)*e22 - (1-p1)*gamma*i22 - p1*gamma2*i22 -(1/tqr)*i22;
r22_dot = (1-p1)*gamma*i22;
qe22_dot = (1/tqr)*e22 - (1/tau)*qe22 - (1/dlay)*qe22;
qi22_dot = (1/tqr)*i22 + (1/tau)*qe22 - gamma*qi22 - (1/dlay)*qi22;
he22_dot = (1/dlay)*qe22 - 1/hosp_day*he22;
hi22_dot = (1/dlay)*qi22 + p1*p2*gamma2*i22 - 1/hosp_day*hi22;
h22_dot =  (p1)*gamma2*i22 - 1/hosp_day*h22; 
che22_dot = (1/dlay)*qe22;   % With link: exposed
chi22_dot = (1/dlay)*qi22;   % With link: infectious
chr22_dot = p1*p2*gamma2*i22; % With link: recovered
ch22_dot =  (p1*(1-p2))*gamma2*i22;  % Without link
ci_dot =  beta*s2/N2*i21 + beta*s2/N2*i22 + q*beta*s2/N2*(qi21+qi22); % Cumulative cases
xdot = [s2_dot; e21_dot; i21_dot; r21_dot; qe21_dot; qi21_dot; h21_dot;
        e22_dot; i22_dot; r22_dot; qe22_dot; qi22_dot; he22_dot; hi22_dot; 
        h22_dot; che22_dot; chi22_dot; chr22_dot; ch22_dot; ci_dot;];
    
end




