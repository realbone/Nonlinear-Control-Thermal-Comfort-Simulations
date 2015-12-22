clear all;


syms mdot;
syms ms;
syms ws;
syms w;
syms Ts;
syms T;


M=8.8393;
C=8.8962e3;
Cpa=1006.43;
Cpw=1800;
G=0.1;
winf=0.003;
Tinf=291.15;
K=100;

% Define the setpoint
w_d=0.009302;
T_d=297.15;

% Define and solve for steady state input
% mdot_d=2;
% ws_d=((mdot_d+G)*w_d-G*winf)/mdot_d;
% Ts_d=(mdot_d*Cpa*T_d+mdot_d*Cpw*w_d*T_d+ K*T_d-K*Tinf)/(mdot_d*Cpa+mdot_d*Cpw*ws_d);


ws_d= 0.02;
mdot_d= G*(w_d-winf)/(ws_d-w_d);
Ts_d=(mdot_d*Cpa*T_d+mdot_d*Cpw*w_d*T_d+ K*T_d-K*Tinf)/(mdot_d*Cpa+mdot_d*Cpw*ws_d);

% Derivatives
wdot = 1/M*(mdot*(ws-w)-G*(w-winf));
Tdot = 1/(M*Cpa+M*w*Cpw)*(mdot*Cpa*(Ts-T)+mdot*Cpw*(ws*Ts-w*T)-K*(T-Tinf)-M*wdot);


j11=subs(diff(wdot, w),[w,ws,mdot],[w_d, ws_d, mdot_d]);
j12=0;
j21=subs(diff(Tdot, w),[w,ws,mdot,Ts,T],[w_d,ws_d,mdot_d,Ts_d,T_d]);
j22=subs(diff(Tdot, T),[w,ws,mdot,Ts,T],[w_d,ws_d,mdot_d,Ts_d,T_d]);

A=[j11 j12;j21 j22]


b11=subs(diff(wdot, mdot), [w,ws,mdot],[w_d, ws_d, mdot_d]);
b12=subs(diff(wdot, ws), [w,ws,mdot],[w_d, ws_d, mdot_d]);
b13=0;
b21=subs(diff(Tdot, mdot),[w,ws,mdot,Ts,T],[w_d,ws_d,mdot_d,Ts_d,T_d]);
b22=subs(diff(Tdot, ws),[w,ws,mdot,Ts,T],[w_d,ws_d,mdot_d,Ts_d,T_d]);
b23=subs(diff(Tdot, Ts),[w,ws,mdot,Ts,T],[w_d,ws_d,mdot_d,Ts_d,T_d]);

B=[b11 b13;b21 b23]


Cc=diag([1 1]);
D=zeros(2,2);



% Passivity
TFF=tf(ss(A,B,Cc,D));
omega = logspace(-4,10,1000);

mineig=[];
for i=1:1
    tff=evalfr(TFF+TFF',j*omega(i));
    eigvalue=eig(tff);
    mineig=[mineig min(eigvalue)];
end

loglog(omega,mineig)
xlabel('w');
ylabel('Smallest Singular Value of G(jw)+G^T(-jw)')
    
    
    
    

