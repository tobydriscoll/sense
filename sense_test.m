clear all
f_ = @(t,u,v) [v(1)*(u(2)-u(1));
             u(1)*(v(2)-u(3))-u(2);
             u(1)*u(2)-v(3)*u(3)];
tspan = [0 10];
u0 = [1;1;1];
p = [10;28;8/3];
[t,soln,u,Ju,Jp] = sense(f_,tspan,u0,p);
subplot 221
plot3(soln(:,1),soln(:,2),soln(:,3))
subplot 222
plot3(u(:,1),u(:,2),u(:,3))
subplot 223
plot3(u(:,4),u(:,5),u(:,6))
subplot 224
plot3(u(:,7),u(:,8),u(:,9))
%%
clear all
f = @(t,u,v) [u(1)*(v(1)-v(2)*u(2));
               u(2)*(v(3)*u(1)-v(4))];
tspan = [0 100];
u0 = [10;10];
p = [0.1;0.02;0.02;0.4];
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 221
plot(u(:,1),u(:,2))
xlabel('prey population')
ylabel('predator population')
title('sensitivity to growth rate of prey')
subplot 222
plot(u(:,3),u(:,4))
xlabel('prey population')
ylabel('predator population')
title('sensitivitiy to kill rate of predators')
subplot 223
plot(u(:,5),u(:,6))
xlabel('prey population')
ylabel('predator population')
title('sensitivity to growth rate of predators')
subplot 224
plot(u(:,7),u(:,8))
xlabel('prey population')
ylabel('predator population')
title('sensitivity to death rate of predators')
%%
clear all
f = @(t,u,v) [u(1)*(v(1)-0.02*u(2));
                u(2)*(v(2)*u(1)-0.4)];
tspan = [0 140];
u0 = [10;10];
p = [0.1;0.02];
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 121
plot(t,[soln(:,1),soln(:,2)])
subplot 122
plot(u(:,1),u(:,2)), hold on
plot(u(:,3),u(:,4))
%%
clear all
f = @(t,u,v) [u(1)*(v-0.02*u(2));
                u(2)*(0.02*u(1)-0.4)];
tspan = [0 140];
u0 = [10;10];
p = 0.1;
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 121
plot(t,[soln(:,1),soln(:,2)])
subplot 122
plot(u(:,1),u(:,2))
%% 
clear all
f = @(t,u,v) u*v(1)*v(2);
tspan = [0 10];
u0 = 1;
p = [0.1,3];
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 121
plot(t,soln);
subplot 122
plot(t,u(:,1)), hold on
plot(t,u(:,2))
%%
clear all
f = @(t,u,v) u*v;
tspan = [0 10];
u0 = 1;
p = 3;
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 121
plot(t,soln);
subplot 122
plot(t,u);
%%
clear all
f_ = @(t,u,v) [v(1)*(u(2)-u(1));
             u(1)*(v(2)-u(3))-u(2);
             u(1)*u(2)-v(3)*u(3)];
tspan = [0 10];
u0 = [1;1;1];
p = [10;28;8/3];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,soln,u,Ju,Jp] = sense(f_,tspan,u0,p,@ode23,options);
subplot 221
plot3(soln(:,1),soln(:,2),soln(:,3))
subplot 222
plot3(u(:,1),u(:,2),u(:,3))
subplot 223
plot3(u(:,4),u(:,5),u(:,6))
subplot 224
plot3(u(:,7),u(:,8),u(:,9))
%%
clear all
f = @(t,u,v) [u(2);
             -(v(1)/v(2)) *sin(u(1))];
tspan = [0 10];
u0 = [0;6];
p = [98.1;1];
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 211
plot(t,soln(:,2))
subplot 212
plot(t,[u(:,2),u(:,4)])
%%
clear all
f = @(t,u,v) [v(1)-v(2)*u(1)-v(3)*u(2)*u(1);
             v(3)*u(2)*u(1) - (v(4)-v(2))*u(2);
             v(4)*u(2)-v(2)*u(3)];
tspan = [0 20];
u0 =[99;1;0];
p = [0;0;.01;0.1];
[t,soln,u,Ju,Jp] = sense(f,tspan,u0,p);
subplot 221
plot(t,soln)
title('S-I-R Model')
subplot 222
plot(t,[u(:,1),u(:,2),u(:,3),u(:,4)])
title('Sensitivity of Susceptible Population')
subplot 223
plot(t,[u(:,5),u(:,6),u(:,7),u(:,8)])
title('Sensitivity of Infected Population')
subplot 224
plot(t,[u(:,9),u(:,10),u(:,11),u(:,12)])
title('Sensitivity of Recovered Population')