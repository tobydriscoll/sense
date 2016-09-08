clear all
f = @(t,u,v) [u(1)*(v(1)-v(2)*u(2));
               u(2)*(v(3)*u(1)-v(4))];
tspan = [0 100];
u0 = [10;10];
p = [0.1;0.02;0.02;0.4];
[t,soln,par_sense,var_sense,Ju,Jp] = sense(f,tspan,u0,p);
subplot 221
plot(par_sense(:,1),par_sense(:,2))
xlabel('prey population')
ylabel('predator population')
title('sensitivity to growth rate of prey')
subplot 222
plot(par_sense(:,3),par_sense(:,4))
xlabel('prey population')
ylabel('predator population')
title('sensitivitiy to kill rate of predators')
subplot 223
plot(par_sense(:,5),par_sense(:,6))
xlabel('prey population')
ylabel('predator population')
title('sensitivity to growth rate of predators')
subplot 224
plot(par_sense(:,7),par_sense(:,8))
xlabel('prey population')
ylabel('predator population')
title('sensitivity to death rate of predators')
%%
clear all
f = @(t,u,v) [-1*v(1)*u(2)*u(1);
              v(1)*u(2)*u(1)-v(2)*u(2);
              v(2)*u(2)];
tspan = [0 30];
u0 =[99;1;0];
p = [.05;0.1];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,soln,par_sense,Ju,Jp] = sense(f,tspan,u0,p);
subplot 221
plot(t,soln)
title('S-I-R Model')
subplot 222
plot(t,[par_sense(:,1),par_sense(:,2)])
title('Sensitivity of Susceptible Population')
subplot 223
plot(t,[par_sense(:,3),par_sense(:,4)])
title('Sensitivity of Infected Population')
subplot 224
plot(t,[par_sense(:,5),par_sense(:,6)])
title('Sensitivity of Recovered Population')