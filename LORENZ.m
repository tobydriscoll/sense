f_ = @(t,u,v) [v(1)*(u(2)-u(1));
             u(1)*(v(2)-u(3))-u(2);
             u(1)*u(2)-v(3)*u(3)];
tspan = [0 13];
u0 = [1;1;1];
p = [10;28;8/3];
[t,soln,p_sense,v_sense,Ju,Jp] = sense(f_,tspan,u0,p);
x = soln(:,1);
y = soln(:,2);
z = soln(:,3);
%%
clf
subplot 311
hold on
plot(t,x,'b')
plot(t,[x+p_sense(:,1),x-p_sense(:,1)],'--b')
subplot 312
hold on
plot(t,y,'r')
plot(t,[y+p_sense(:,2),x-p_sense(:,2)],'--r')
subplot 313
hold on
plot(t,z,'g')
plot(t,[z+p_sense(:,3),z-p_sense(:,3)],'--g')
%%
clf
subplot 311
hold on
plot(t,x,'b')
plot(t,[x+p_sense(:,4),x-p_sense(:,4)],'--b')
subplot 312
hold on
plot(t,y,'r')
plot(t,[y+p_sense(:,5),x-p_sense(:,5)],'--r')
subplot 313
hold on
plot(t,z,'g')
plot(t,[z+p_sense(:,6),z-p_sense(:,6)],'--g')