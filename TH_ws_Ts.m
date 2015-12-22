M=8.8393;
C=8.8962e3;
Cpa=1006.43;
Cpw=1800;
G=0.1;
winf=0.003;
Tinf=291.15;
K=100;

w0=0.004;
T0=293.15;
Td=273.15+24;
wd=0.009302;

mdotH=0.5;mdotL=0.0;
wsH=0.05; wsL=0.002;
TsH=310;TsL=280;

dt=0.2;

mdot=0.2;ws=0.012;Ts=301.15;

w=w0;T=T0;

w_traj=[];T_traj=[];
mdot_traj=[];Ts_traj=[];ws_traj=[];
Td_traj=[];wd_traj=[];
C_traj=[];

intw=0;intT=0;
epsilon = 1e-5;

set_flag=0;

steps=250;

for i=1:steps
    w_traj=[w_traj w];
    T_traj=[T_traj T];

    if Ts<TsH && Ts>TsL && ws<wsH && ws>wsL  % Anti-windup
        intw=intw+dt*(w-wd);
        intT=intT+dt*(T-Td);
    end
    
    % Controllers
    Ts=295-8*(T-Td)-0.3*intT;
    ws=0.02-30*(w-wd)-2*intw;
%     mdot = 0.2;
    
    
    Ts=satu(Ts,TsH,TsL);
    ws=satu(ws,wsH,wsL);
    mdot=satu(mdot, mdotH, mdotL);
    
    
    
    Ts_traj=[Ts_traj Ts];ws_traj=[ws_traj ws];mdot_traj=[mdot_traj mdot];
    Td_traj=[Td_traj Td];wd_traj=[wd_traj wd];
    
%     if i==steps/2
%         Td=290;
%         wd=0.01;
%         set_flag=1;
%     end
    
%     if i>1000 && abs(w_traj(i)-w_traj(i-1))+abs(T_traj(i)-T_traj(i-1))<=epsilon && set_flag>0
%         break;
%     end

    C = Cpa+w*Cpw;
    wdot = 1/M*(mdot*(ws-w)-G*(w-winf));
    Tdot = 1/(M*C)*(mdot*Cpa*(Ts-T)+mdot*Cpw*(ws*Ts-w*T)-K*(T-Tinf)-M*wdot);
    w = w+dt*wdot;
    T = T+dt*Tdot; 
    
end


time = 0:dt:(i-1)*dt;

figure(1);
subplot(2,1,1);
plot(time, w_traj);hold on;
plot(time, wd_traj,'--');
xlabel('Time(s)');ylabel('w');

subplot(2,1,2);
plot(time, T_traj);hold on;
plot(time, Td_traj,'--');
xlabel('Time(s)');ylabel('T(K)');

    
figure(2);
subplot(2,1,1);
plot(time, ws_traj);
xlabel('Time(s)'); ylabel('ws');

subplot(2,1,2);
plot(time, Ts_traj);
xlabel('Time(s)'); ylabel('Ts(K)');


