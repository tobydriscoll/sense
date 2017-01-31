function [t,solution,par_sensitivities,var_sensitivities,Ju,Jp] = sense(fun,tspan,u0,param,solver,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SENSE Utilizes automatic differentiation to find the sensitivity of a
%system of first-order differential equations to the parameters of the
%system.
%INPUTS:
%   fun     dudt = fun(t,u,v) (n x 1)
%   tspan   time span of interest
%   u0      initial values of the variables for fun (n x 1)
%   param   values of parameters for fun (m x 1)
%   solver  solver to be used by the function. if empty, uses ode45 as
%           default
%   opts    MATLAB ODE solver options. If empty, uses default options as
%           defined by MATLAB
%OUTPUTS:
%   t       time range given by ode solver
%   soln    solution to dudt = fun(t,u,v) for initial value u0 and timespan
%           tspan (n x length(t))
%   par_sensitivities   sensitivity to fun with respect to each parameter,
%                       formatted as (1:n,:) = sensitivity of fun to first
%                       parameter, (n+1:2*n,:) = sensitivity of fun to 2nd param,
%                       etc. (m*n x length(t))
%   Ju      Jacobian with respect to the variables. i.e., first row is
%           derivative of each component of f with respect to x, second is
%           each component wrt y, etc.
%   Jp      Jacobian wrt parameters. i.e., first row is derivative of each
%           component of f wrt first parameter, second is each component wrt
%           2nd, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6
    opts = odeset;  % use built-in defaults
    if nargin < 5
        solver = @ode45;
    end
    if ischar(solver), solver = str2func(solver); end    
end

nvar = length(u0);
npar = length(param);
initsense = reshape(eye(nvar),[],1);
if nvar>1 && npar>1
    v0 = zeros(nvar*(2+2*npar+2*nvar),1);
    v0(1:nvar) = u0;
    v0(end-length(initsense)+1:end) = initsense;
    f = @(t,u) odefun1(t,u,param,fun,nvar,npar);
    [t,u] = solver(f,tspan,v0,opts);
    solution = u(:,1:nvar);
    par_sensitivities = u(:,nvar+1:nvar+nvar*npar);
    Ju_ = u(:,nvar+nvar*npar+1:nvar+nvar*npar+nvar*nvar);
    Jp_ = u(:,nvar+nvar*npar+nvar*nvar+1:nvar+nvar*npar+nvar*nvar+nvar*npar);
    var_sensitivities = u(:,end-nvar^2+1:end);
    Ju = zeros(nvar,nvar,length(t));
    Jp = zeros(nvar,npar,length(t));
    for i = 1:length(t)
        r = 1;
        for j = 1:nvar
            for k = 1:nvar
                Ju(j,k,i) = Ju_(i,r);
                r = r+1;
            end
        end
        
    end
    for i = 1:length(t)
        r = 1;
        for j = 1:nvar
            for k = 1:npar
                Jp(j,k,i) = Jp_(i,r);
                r = r+1;
            end
        end
    end
elseif nvar == 1 && npar >1
    v0 = zeros(3+2*npar,1);
    v0(1) = u0;
    f = @(t,u) odefun2(t,u,param,fun,1,npar);
    [t,u] = solver(f,tspan,v0,opts);
    solution = u(:,1);
    par_sensitivities = u(:,2:1+npar);
    Ju = u(:,2+npar)';
    Jp = u(3+npar:2+npar+npar);
    var_sensitivities = u(end-nvar^2+1:end);
elseif nvar > 1 && npar == 1
    v0 = zeros(nvar*(4+nvar),1);
    v0(1:nvar) = u0;
    f = @(t,u) odefun3(t,u,param,fun,nvar,1);
    [t,u] = solver(f,tspan,v0,opts);
    solution = u(:,1:nvar);
    par_sensitivities = u(:,nvar+1:2*nvar);
    Ju_ = u(:,2*nvar+1:2*nvar+nvar*nvar);
    Ju = zeros(nvar,nvar,length(t));
    for i = 1:t
        r = 1;
        for j = 1:nvar
            for k = 1:nvar
                Ju(k,j,i) = Ju_(i,r);
                r = r+1;
            end
        end
    end
    Jp = u(:,2*nvar+nvar*nvar+1:2*nvar+nvar*nvar+nvar)';
    var_sensitivities = u(:,end-nvar^2+1:end);
elseif nvar == 1 && npar == 1
    v0 = zeros(6,1);
    v0(1) = u0;
    f = @(t,u) odefun4(t,u,param,fun,nvar,npar);
    [t,u] = solver(f,tspan,v0,opts);
    solution = u(:,1);
    par_sensitivities = u(:,2);
    Ju = u(:,3);
    Jp = u(:,4);
    var_sensitivities = u(:,5);
end
end
function dvdt = odefun1(t,v,p,f,n,m)
% >1 variable, >1 parameter
derivs = eye(n+m);
dvdt = zeros(length(v),1);
vars = valder(n);
pars = valder(m);
for j = 1:n
    vars(j) = valder(v(j),derivs(:,j));
end
for j = 1:m
    pars(j) = valder(p(j),derivs(:,n+j));
end
fun = f(t,vars,pars);

jacvar = zeros(n,n);
for i = 1:n
    for j = 1:n
        jacvar(i,j) = fun(i).der(j);
    end
end
jacpar = zeros(n,m);
for i = 1:n
    for j = 1:m
        jacpar(i,j) = fun(i).der(n+j);
    end
end
paramder = reshape(v(n+1:n+n*m),n,m);
for i = 1:n
    dvdt(i) = fun(i).val;
end
V = reshape(jacvar*paramder+jacpar,[],1);
dvdt(n+1:n+n*m) = V;
dvdt(n+n*m+1:n+n*m+n*n) = reshape(jacvar,[],1);
dvdt(n+n*m+n*n+1:end-n^2-n) = reshape(jacpar,[],1);
dvdt(end-n^2+1:end) = reshape(jacvar*reshape(v(end+1-n^2:end),n,n),[],1);
end
function dvdt = odefun2(t,v,p,f,n,m)
% 1 variable, >1 paramter
derivs = eye(n+m);
vars = valder(v(1),derivs(:,1));
pars = valder(m);
for i = 1:m
    pars(i) = valder(p(i),derivs(:,1+i));
end
dvdt = zeros(length(v),1);
fun = f(t,vars,pars);
dvdt(1) = fun.val;
jacvar = fun.der(1);
jacpar = zeros(1,m);
for i = 1:m
    jacpar(i) = fun.der(1+i);
end
paramder = reshape(v(n+1:n+n*m),n,m);
V = reshape(jacvar*paramder+jacpar,[],1);
dvdt(n+1:n+n*m) = V;
dvdt(n+n*m+1) = jacvar;
dvdt(n+n*m+2:end-n^2) = reshape(jacpar,[],1);
dvdt(end-n^2+1:end) = reshape(jacvar*reshape(v(end-n^2+1:end),n,n),[],1);
end
function dvdt = odefun3(t,v,p,f,n,m)
% >1 variable, 1 parameter
derivs = eye(n+m);
vars = valder(n);
pars = valder(p,derivs(:,end));
for i = 1:n
    vars(i) = valder(v(i),derivs(:,i));
end
fun  = f(t,vars,pars);
dvdt = zeros(length(v),1);
for i = 1:n
    dvdt(i) = fun(i).val;
end
jacvar = zeros(n,n);
jacpar = zeros(n,1);
for i = 1:n
    for j = 1:n
        jacvar(i,j) = fun(i).der(j);
    end
end
for i = 1:n
    jacpar(i) = fun(i).der(end);
end
paramder = v(n+1:n+n*m);
V = reshape(jacvar*paramder+jacpar,[],1);
dvdt(n+1:n+n*m) = V;
dvdt(n+n*m+1:n+n*m+n*n) = reshape(jacvar,[],1);
dvdt(n+n*m+n*n+1:n+n*m+n*n+n*m) = reshape(jacpar,[],1);
dvdt(end-n^2+1:end) = reshape(jacvar*(reshape(v(end-n^2+1:end),n,n)),[],1);
end
function dvdt = odefun4(t,v,p,f,n,m)
dvdt = zeros(4,1);
var = valder(v(1),[1 0]);
par = valder(p,[0 1]);
fun = f(t,var,par);
jacvar = fun.der(1);
jacpar = fun.der(2);
dvdt(1) = fun.val;
dvdt(2) = jacvar*v(2)+jacpar;
dvdt(3) = jacvar;
dvdt(4) = jacpar;
dvdt(5) = jacvar*v(5);
end