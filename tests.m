% Tests of ODE sensitivity

%% Linear, constant-coefficients
A = [3 0 -1;1 1 0;2 -2 2];
f = @(t,u,p) p*A*u;
t = linspace(0,2,50)';
p = -6;
[t,u,z] = sense(f,t,[-1;0;1],p);

% exact answer
B = [ p*A, 0*A; A, p*A ];  % system for w = [u;z] 
[t,w] = ode45(@(t,w) B*w,t,[-1;0;1;0;0;0]);

assert( norm(squeeze(z)-w(:,4:6),inf) < 1e-8 )

%% Lorenz system
f = @(t,u,v) [v(1)*(u(2)-u(1));
             u(1)*(v(2)-u(3))-u(2);
             u(1)*u(2)-v(3)*u(3)];
tspan = [0 2];
u0 = [1;1;1];
p = [10;28;8/3];
opt = odeset('reltol',1e-7,'abstol',1e-7);
[t,u,Z] = sense(f,tspan,u0,p,@ode113,opt);
finalZ = [
   0.054540202470694  -0.895125686886151  -5.021772002532813
   0.310395767098313  -0.822574382021570  -9.757372302834186
   0.160477095899202   2.145918785312579  -2.200955088434970
   ];
assert( norm(squeeze(Z(end,:,:))-finalZ) < 1e-5 )

%% Lotka-Volterra
f = @(t,u,v) [u(1)*(v(1)-v(2)*u(2));
               u(2)*(v(3)*u(1)-v(4))];
tspan = [0 100];
u0 = [10;10];
p = [0.1;0.02;0.02;0.4];
opt = odeset('reltol',1e-7,'abstol',1e-7);

[t,u,Z] = sense(f,tspan,u0,p,@ode45,opt);

%% SIR model
f = @(t,u,v) [-1*v(1)*u(2)*u(1);
              v(1)*u(2)*u(1)-v(2)*u(2);
              v(2)*u(2)];
tspan = [0 30];
u0 = [99;1;0];
p = [.05;0.1];
opt = odeset('reltol',1e-7,'abstol',1e-7);

[t,u,Z] = sense(f,tspan,u0,p,@ode45,opt);

%% Van der Pol
f = @(t,y,mu) [ y(2); mu*(1-y(1)^2)*y(2)-y(1) ];
tspan = [0 100];
mu = 100;
u0 = [2;0];
[t,u,Z,Ju] = sense(f,tspan,u0,mu,@ode45);
for i = 1:size(t,1)
    lambda(i,:) = eig( Ju(t(i),u(i,:),mu) )';
end
plot(t,lambda)
