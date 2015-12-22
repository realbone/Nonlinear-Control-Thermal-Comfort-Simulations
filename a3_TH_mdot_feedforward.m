clear all;
M=1.19*8;
C=8.8962e3;
Cpa=1006.43;
Cpw=1800;
G=0.1;
winf=0.003;
Tinf0=273.15+23;
K=100;

w0=0.004;
T0=Tinf0;
Td=273.15+24;
wd=0.009302;

mdotH=0.25;mdotL=0.0;
wsH=0.05; wsL=0.002;
TsH=301;TsL=280;

dt=2;

ws_d= 0.013; 
mdot_d= G*(wd-winf)/(ws_d-wd);
Ts_d=(mdot_d*Cpa*Td+mdot_d*Cpw*wd*Td+ K*Td-K*Tinf0)/(mdot_d*Cpa+mdot_d*Cpw*ws_d);

mdot=mdot_d;ws=ws_d;Ts=Ts_d;

w=w0;T=T0;

w_traj=[];T_traj=[];
mdot_traj=[];Ts_traj=[];ws_traj=[];
Td_traj=[];wd_traj=[];
C_traj=[];Tinf_traj=[];

intw=0;intT=0;
epsilon = 1e-5;
steps=200;

for i=1:steps
    w_traj=[w_traj w];
    T_traj=[T_traj T];

    if Ts<TsH || Ts>TsL || mdot<mdotH || mdot>mdotL  % Anti-windup
        intw=intw+dt*(w-wd);
        intT=intT+dt*(T-Td);
    end
    
    % Controllers

     mdot=mdot_d - 200*((w-wd)+0.003*(T-Td));%-10*(intw+0.01*intT);
    
    
    Ts=satu(Ts,TsH,TsL);
    ws=satu(ws,wsH,wsL);
    mdot=satu(mdot, mdotH, mdotL);
    Tinf=Tinf0+0*sin(0.2*i);
    
    Ts_traj=[Ts_traj Ts];ws_traj=[ws_traj ws];mdot_traj=[mdot_traj mdot];
    Td_traj=[Td_traj Td];wd_traj=[wd_traj wd];Tinf_traj=[Tinf_traj Tinf];
    
    C = Cpa+w*Cpw;
    wdot = 1/M*(mdot*(ws-w)-G*(w-winf));
    Tdot = 1/(M*C)*(mdot*Cpa*(Ts-T)+mdot*Cpw*(ws*Ts-w*T)-K*(T-Tinf)-M*wdot);
    w = w+dt*wdot;
    T = T+dt*Tdot;  
end


time = 0:dt:(i-1)*dt;

fontsize=15;

% Plots that's nor moving!
figure(1);
subplot(2,2,1)
pbaspect([2,1,1]);
plot(time,w_traj,'LineWidth',1.5);hold on;
plot(time, wd_traj,'--');
set(gca,'fontsize',fontsize)
xlabel('Time(s)','FontSize',fontsize);ylabel('w','FontSize',fontsize);
legend('w','w_d','location','southeast');
% hgexport(gcf,'1_w.eps')

subplot(2,2,2)
pbaspect([2,1,1]);
plot(time,T_traj,'LineWidth',1.5);hold on;
plot(time, Td_traj,'--');
plot(time(1:i), Tinf_traj(1:i),'LineWidth',1.5,'Color','r');hold on;
legend('T','T_d','T_{\infty}','location','southeast');
set(gca,'fontsize',fontsize)
xlabel('Time(s)','FontSize',fontsize);ylabel('T (K)','FontSize',fontsize);
% hgexport(gcf,'1_T.eps')


subplot(2,2,3)
pbaspect([2,1,1]);
plot(time, mdot_traj,'LineWidth',1.5);hold on;
plot(time, repmat(mdotH,1,size(time,2)),'r','LineWidth',1.5);
xlabel('Time(s)','FontSize',fontsize);yl=ylabel('$\dot{m}\ (kg/s)$','FontSize',fontsize);
lg=legend('$\dot{m}$', '$\dot{m}_{H}$','location','southeast');
set(gca,'fontsize',fontsize)
set(lg, 'Interpreter', 'latex');
set(yl, 'Interpreter', 'latex');
axis([0 400 0 0.3])





% 
% 
% red(:,:,1)=ones(100,100);red(:,:,2)=zeros(100,100);red(:,:,3)=zeros(100,100);
% blue(:,:,1)=zeros(100,100);blue(:,:,2)=zeros(100,100);blue(:,:,3)=ones(100,100);
% green(:,:,1)=zeros(100,100);green(:,:,2)=ones(100,100);green(:,:,3)=zeros(100,100);
% 
% 
% 
% for i=1:steps
% figure(1);
% subplot(3,2,1);
% pbaspect([2,1,1]);
% plot(time(1:i), w_traj(1:i),'LineWidth',1.5);hold on;
% plot(time, wd_traj,'--');
% set(gca,'fontsize',15)
% xlabel('Time(s)','FontSize',fontsize);ylabel('w','FontSize',fontsize);
% 
% subplot(3,2,2);
% im=(w_traj(i)-w0)/(wd-w0)*(green-blue)+blue;
% 
% imshow(im);
% title(strcat('Humidity Level::',num2str(w_traj(i))),'FontSize',fontsize);
% 
% 
% subplot(3,2,3);
% pbaspect([2,1,1]);
% plot(time(1:i), T_traj(1:i),'LineWidth',1.5);hold on;
% plot(time(1:i), Tinf_traj(1:i),'LineWidth',1.5,'Color','r');hold on;
% legend('T','T_{\infty}','Location','southeast');
% plot(time, Td_traj,'--');
% set(gca,'fontsize',15)
% 
% xlabel('Time(s)','FontSize',fontsize);ylabel('T(K)','FontSize',fontsize);
% 
% subplot(3,2,4);
% im=(T_traj(i)-T0)/(Td-T0)*(red-blue)+blue;
% imshow(im);
% title(strcat('Temperature Level::',num2str(T_traj(i))),'FontSize',fontsize);
% 
% subplot(3,2,5);
% pbaspect([2,1,1])
% plot(time(1:i), mdot_traj(1:i),'LineWidth',1.5);hold on;
% plot(time, repmat(mdotH,1,size(time,2)),'r');
% xlabel('Time(s)','FontSize',fontsize);ylabel('mdot','FontSize',fontsize);
% 
% 
% 
% set(gca,'fontsize',15)
% if i==1
%     pause();
% end
% 
% if i>1
%     MM(i-1)=getframe(gcf);
% end
% 
% end
% 
% movie2avi(MM,'mdot_feedforward.avi');
% 
% 
% 
% 
