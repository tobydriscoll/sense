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
        
        function varargout = eval(s,t,varargin)
            [varargout{1:nargout}] = deval(s.solution,t,varargin{:}).';
        end

        function Z = sense(s)
            npar = length(s.p);
            nvar = length(s.init);

            z0 = zeros(nvar*npar,1);
            [~,z] = s.solver(@odefun,s.t,z0,s.opts);
            nt = length(s.t);
            Z = reshape( z,[nt nvar npar] );
            
            function dzdt = odefun(t,z)
                Z = reshape(z,nvar,npar);
                u = s.eval(t);
                
                Au = s.dfdu(t,u,s.p);
                Ap = s.dfdp(t,u,s.p);
    
                dZdt = Au*Z + Ap;
                dzdt = dZdt(:);
            end
        end
        
        function lambda = jaceig(s)
            lambda = NaN(length(s.t),length(s.u));
            for i = 1:length(s.t)
                lambda(i,:) = eig( s.dfdu(s.t(i),s.u(i,:),s.p) ).';
            end
        end
        
        
            
        end
    end
        
    
end
