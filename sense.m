function [t,solution,J,jacvar,jacpar] = sense(fun,tspan,u0,param,solver,opts)
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
    if nargin < 4 || nargin(fun) == 2
        param = [];
    end
end

npar = length(param);
nvar = length(u0);

v0 = [ u0; zeros(nvar*npar,1) ];
[t,v] = solver(@odefun,tspan,v0,opts);
solution = v(:,1:nvar).';
nt = length(t);
J = reshape( v(:,nvar+(1:nvar*npar)).',[nvar npar nt] ); 
jacvar = @dfun_dvar;
jacpar = @dfun_dpar;

    function [Au,dudt] = dfun_dvar(t,u,p)
        u_diff = valder(u,eye(nvar));
        f = fun(t,u_diff,p);
        Au = f.der;
        dudt = f.val;
    end

    function Ap = dfun_dpar(t,u,p)
        p_diff = valder(p,eye(npar));
        f = fun(t,u,p_diff);
        Ap = cat(1,f.der);
    end

function dvdt = odefun(t,v)
    inputs = mat2cell(v,[nvar nvar*npar]);
    [u,z] = deal(inputs{:}); 
    Z = reshape(z,nvar,npar);
    
    [Au,dudt] = dfun_dvar(t,u,param);
    Ap = dfun_dpar(t,u,param);
    
    dZdt = Au*Z + Ap;
    
    dvdt = [ dudt; dZdt(:) ];

end

end
