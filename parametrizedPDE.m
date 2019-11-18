classdef parametrizedPDE
    %
    %    - grad.(k grad u) + (a . grad u) = f  in Omega
    %                                  u  = 0  on dOmega
    %
    %                       +-------+-------+
    %                       |       |       |
    %                       |  x_1  |  x_4  |
    %  Omega=[0,1]^2 ,      |       |       |
    %  with 4 subdomains:   +-------+-------+
    %                       |       |       |
    %                       |  x_2  |  x_3  |
    %                       |       |       |
    %                       +-------+-------+
    %
    %  Diffusion field k = k(x) takes the value x_i on the subdomain
    %  Omega_i. The advection field a = a(x) is a rotating vector field
    %  proportional to x_5. The 5D parameter x=(x_1,...,x_5) is such 
    %  that
    %               0.05 <= x_i <= 1     i=1..4
    %               100  <= x_5 <= 100
    %  When it is random, the parameter X is uniformly distributed.
    %
    %  Solution : u(x) in R^n
    %  Quantity of interest: y = y(x) = int_Omega3 u(x) dx
    %
    %
    % Author: Olivier Zahm (olivier.zahm@inria.fr)
    % Date: November 2019
    
    properties
        
        n  % Size of the solution vector u(x)
        
        K1 % Diffusion matrix on Omega1
        K2 % Diffusion matrix on Omega2
        K3 % Diffusion matrix on Omega3
        K4 % Diffusion matrix on Omega4
        Ad % Advection matrix
        
        b  % Right-hand side
        q  % Extractor of the quantity of interest
        
    end
    
    
    methods
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function self = parametrizedPDE()
            
            % Number of element per spatial dimension (must be even!)
            nelem = 80;
            
            % Create the algebraic problem with the Finite Element Method
            self = initFEM(self,nelem);
            
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function X = randX(self)
            % Random draw of the parameter uniformly on its interval
            
            % Draw the four diffusion coefficients
            Kmin = 0.05;
            Kmax = 1;
            K = Kmin + rand(4,1)*(Kmax-Kmin);
            
            % Draw the advection coefficient
            amin = -100;
            amax = +100;
            a = amin + rand(1,1)*(amax-amin);
            
            % Assemble the parameter
            X = [K(1);K(2);K(3);K(4);a];
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function u = u(self,x)
            % Compute the PDE solution associated with x
            
            % If x is not specified, draw it randomly
            if nargin<2
                x = self.randX();
            end
            
            % Assemble the parameter-dependent operator
            A = self.A(x);
            
            % Solve the system
            u = A\self.b;
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function A = A(self,x)
            % Compute the parameter dependent matrix A(x)
            
            % If x is not specified, draw it randomly
            if nargin<2
                x = self.randX();
            end
            A = self.K1*x(1) + self.K2*x(2) + self.K3*x(3) + self.K4*x(4)  + x(5)*self.Ad;
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function plotsol(self,u)
            % Plot the solution u
            
            % If u is not specified, draw it randomly
            if nargin<2
                u = self.u();
            end
            
            % Reshape the solution for post-precessing
            nelem = sqrt(size(u,1));
            u = reshape(u,nelem,nelem);
            
            % Add the contour of 0 (homogeneous boundary conditions)
            Z = zeros(nelem,1);
            ZZ= zeros(1,nelem+2);
            u = [ZZ; [Z,u,Z] ; ZZ];
            
            % Plot the solution
            h = 1/(nelem+1);
            x = 0:h:1;
            surface(x,x,u,'EdgeColor','none')
            colorbar
            axis([0 1 0 1])
            axis square
            
            view([0,90])
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function plotPDEsetup(self,X)
            % Didactic plot the PDE setup: aiffusion coefficient, advection
            %field, source term and one realization of the solution.
            
            % Draw a random parameter
            if nargin<2
                X = self.randX();
                a = X(5);
            end
            
            clf
            set(groot,'defaulttextinterpreter','latex');  
            
            %%% Plot the diffusion coefficient:
            subplot(2,2,1) 
            Kmin = 0.05;
            Kmax = 1;
            
            colors = parula(200);
            colors = colors( ceil((X(1:4)-Kmin)/(Kmax-Kmin)*200),: );
            rectangle('Position',[0 0.5 0.5 0.5],'FaceColor',colors(1,:));
            text(0.25,0.75,'$\kappa_1=x_1$','FontSize',15,'HorizontalAlignment','center')
            rectangle('Position',[0 0 0.5 0.5],'FaceColor',colors(2,:));
            text(0.25,0.25,'$\kappa_2=x_2$','FontSize',15,'HorizontalAlignment','center')
            rectangle('Position',[0.5 0 0.5 0.5],'FaceColor',colors(3,:));
            text(0.75,0.25,'$\kappa_3=x_3$','FontSize',15,'HorizontalAlignment','center')
            rectangle('Position',[0.5 0.5 0.5 0.5],'FaceColor',colors(4,:));
            text(0.75,0.75,'$\kappa_4=x_4$','FontSize',15,'HorizontalAlignment','center')
            
            axis([0 1 0 1])
            axis square
            title('Diffusion coefficient $\kappa$')
            
            %%% Plot the advection field:
            subplot(2,2,2) 
            h=1/10;
            [x,y] = meshgrid(0:h:1,0:h:1);
            coef=0.007;
            u = (-y+1/2)*a*coef;
            v = (-1/2+x)*a*coef;
            
            quiver(x,y,u,v,0)
            axis([0 1 0 1])
            axis square
            title('Advection field $a$')
            
            %%% Plot the source term:
            subplot(2,2,3) 
            rectangle('Position',[0 0 1 1])
            rectangle('Position',[0.05 0.4 0.15 0.2])
            rectangle('Position',[0.2 0.4 0.15 0.2])
            text(0.125,0.5,'$+1$','HorizontalAlignment','center')
            text(0.275,0.5,'$-1$','HorizontalAlignment','center')
            text(0.5,0.5,'$0$','HorizontalAlignment','center')
            
            axis([0 1 0 1])
            axis square
            title('Source term $f$')
            
            %%% Plot the solution:
            subplot(2,2,4) 
            self.plotsol(self.u(X))
            colorbar('off') 
            axis([0 1 0 1])
            axis square
            
            Xchar = arrayfun(@(X)[num2str(X,2) ', '],X','UniformOutput',0);
            Xchar = [Xchar{:}];
            title(['Solution $u(x)$ for $x=(' Xchar(1:(end-2)) ')$'])
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function self = initFEM(self,nelem)
            % Initialize the model.
            
            %%% Number of elements, mesh size, nb of degrees of freedom
            nelem = ceil(nelem/2)*2; % check nelem is even
            h = 1/nelem; % mesh size
            ndof = nelem-1; % number of degrees of freedom per dimension
            self.n = ndof^2; % total number of degrees of freedom
            
            %%% Build the diffusion matrices K1,K2,K3,K4
            %
            D_left = [2*ones(nelem/2-1,1) ; 1 ; zeros(nelem/2-1,1)];
            Dul_left = [-ones(nelem/2-1 ,1 ) ; zeros(nelem/2-1,1)];
            %
            K_left = spdiags(D_left,0,ndof,ndof);
            K_left = K_left + spdiags(Dul_left,-1,ndof,ndof);
            K_left = K_left + spdiags(Dul_left,-1,ndof,ndof)';
            K_left = K_left/h;
            %
            K_right = flip(fliplr(K_left));
            
            D_left = [2/3*ones(nelem/2-1,1) ; 1/3 ; zeros(nelem/2-1,1)];
            Dul_left = [1/6*ones(nelem/2-1 ,1 ) ; zeros(nelem/2-1,1)];
            %
            M_left = spdiags(D_left,0,ndof,ndof);
            M_left = M_left + spdiags(Dul_left,-1,ndof,ndof);
            M_left = M_left + spdiags(Dul_left,-1,ndof,ndof)';
            M_left = M_left*h;
            %
            M_right = flip(fliplr(M_left));
            
            self.K1 = kron(K_left,M_right) + kron(M_left,K_right);
            self.K2 = kron(K_left,M_left) + kron(M_left,K_left);
            self.K3 = kron(K_right,M_left) + kron(M_right,K_left);
            self.K4 = kron(K_right,M_right) + kron(M_right,K_right);
            
            
            %%% Build the advection matrix:
            D = 2/3*[1:ndof] *h^2;
            Dud = 1/12*(2*[1:(ndof-1)]+1) *h^2;
            
            M_lin = spdiags(D',0,ndof,ndof);
            M_lin = M_lin + spdiags(Dud',-1,ndof,ndof);
            M_lin = M_lin + spdiags(Dud',-1,ndof,ndof)';
            
            M_aff1 = M_lin - 1/2*(M_left+M_right);
            M_aff2 =-M_lin + 1/2*(M_left+M_right);
            
            Ad_aff = 1/2*spdiags(ones(ndof-1,1),-1,ndof,ndof);
            Ad_aff = Ad_aff - 1/2*spdiags(ones(ndof-1,1),-1,ndof,ndof)';
            
            self.Ad = kron(Ad_aff,M_aff1) + kron(M_aff2,Ad_aff);
            
            %%% Build the right-hand side b
            funX = @(x) (x>0.05).*(x<0.2)*h - (x>0.2).*(x<=0.35)*h ;
            funY = @(x) (x>=0.4).*(x<=0.6)*h;
            Fx = funX(h*[1:ndof]');
            Fy = funY(h*[1:ndof]');
            self.b = kron(Fx,Fy);
            
            %%% Build the observation operator
            self.q = 4*kron(sum(M_right) , sum(M_left))';
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
    end %endMethods
end %endClass
