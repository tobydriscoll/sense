classdef odesol
    
    properties
        dudt      % function defining ODE
        init      % initial value
        p         % parameter value
        tspan     % independent variable range
        dfdu      % Jacobian of dudt WRT u
        dfdp      % Jacobian of dudt WRT parameters  
        opts      % solver options
        solver    % selected solver
    end
    
    properties (Access=private)
        solution  % DEVAL-type solution struct
    end
    
    properties (Dependent=true)
        t         % times selected by the integrator
        u         % solution
        stats     % solver statistics
        nvar      % number of scalar variables
        npar      % number of scalar parameters
    end
    
    methods
        function s = odesol(dudt,tspan,u0,p,solver,opts)
            if nargin < 6
                opts = odeset;  % use built-in defaults
                if nargin < 5
                    solver = @ode45;
                end
                if nargin < 4 || nargin(fun) == 2
                    p = [];
                end
            end
            
            if ischar(solver), solver = str2func(solver); end
            
            npar = length(p);
            nvar = length(u0);

            odefun = @(t,u) dudt(t,u,p);
            s.solution = solver(odefun,tspan,u0,opts);

            function Au = df_du(t,u,p)
                u_diff = valder(u,eye(nvar));
                y = s.dudt(t,u_diff,p);
                Au = y.der;
            end
                
            function Ap = df_dp(t,u,p)
                p_diff = valder(p,eye(npar));
                y = s.dudt(t,u,p_diff);
                Ap = cat(1,y.der);
            end
            
            s.dudt = dudt;
            s.init = u0;
            s.p = p;
            s.tspan = tspan([1 end]);
            s.dfdu = @df_du;
            s.dfdp = @df_dp;
            s.opts = opts;
            s.solver = solver;          
        end
        
        function t = get.t(s)
            t = s.solution.x;
        end
        
        function u = get.u(s)
            u = s.solution.y;
        end
        
        function stat = get.stats(s)
            stat = s.solution.stats;
        end
        
        function n = get.nvar(s)
            n = length(s.init);
        end
        
        function m = get.npar(s)
            m = length(s.p);
        end
        
        function varargout = eval(s,t,varargin)
            [varargout{1:nargout}] = deval(s.solution,t,varargin{:}).';
        end
        
        function varargout = plot(s,varargin)
            f = @(t) s.eval(t,varargin{:});
            [varargout{1:nargout}] = fplot(f,s.tspan);
        end
        
        function varargout = phaseplot(s,k)
            if nargin < 2
                k = 1:s.nvar;
            end
            if length(k)<2 || length(k)>3
                error('Must have 2 or 3 variables for a phase plot.')
            end
            x = @(t) s.eval(t,k(1));
            y = @(t) s.eval(t,k(2));
            if length(k)==2
                [varargout{1:nargout}] = fplot(x,y,s.tspan);
            elseif length(k)==3
                z = @(t) s.eval(t,k(3));
                [varargout{1:nargout}] = fplot3(x,y,z,s.tspan);
            end
        end

        function Z = sense(s)
            np = s.npar;
            nv = s.nvar;

            z0 = zeros(nv*np,1);
            [~,z] = s.solver(@odefun,s.t,z0,s.opts);
            nt = length(s.t);
            Z = reshape( z,[nt nv np] );
            
            function dzdt = odefun(t,z)
                Z = reshape(z,nv,np);
                u = s.eval(t);
                
                Au = s.dfdu(t,u,s.p);
                Ap = s.dfdp(t,u,s.p);
    
                dZdt = Au*Z + Ap;
                dzdt = dZdt(:);
            end
        end
        
        function ev = jaceig(s)
            ev = NaN(length(s.t),s.nvar);
            for i = 1:length(s.t)
                ev(i,:) = eig( s.dfdu(s.t(i),s.u(:,i),s.p) ).';
            end
        end
        
        function lambda = adjoint(s,phi)
            % Solve the adjoint equation.           
            v0 = zeros(size(s.init));
            lambda = odesol(@odefun,s.tspan([2 1]),v0,s.p,s.solver,s.opts);
            function dvdt = odefun(t,v,p)
                Au = s.dfdu(t,s.eval(t),p);
                dvdt = phi(t) - Au'*v;
            end
            lambda.tspan = s.tspan;
        end
        
        function v = adjsense(s,g)          
            
            function g_u = dgdu(t)
                u = s.eval(t);
                u = valder(u,eye(s.nvar));
                G = g(t,u,s.p);
                g_u = G.der.';
            end
            lambda = s.adjoint(@dgdu);
            
            function q = integrand(t)
                u = s.eval(t);
                f_p = s.dfdp(t,u,s.p);
                p = valder(s.p,eye(s.npar));
                G = g(t,u,p);
                g_p = G.der;
                q = g_p - lambda.eval(t)*f_p;
            end
            v = integral(@integrand,s.tspan(1),s.tspan(2),...
                'arrayvalued',true);

        end
            
    end
        
    
end
