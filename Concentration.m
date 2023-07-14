classdef Concentration
    % free actin and polymerized actin concentration
    %

    properties
        N = 32
        field
        initial_concentration = 3.3*10^-8
        concentration_randomness = 0
        depolymerization_rate
        polymerization_rate
        diffusivity

    end

    methods
        function obj = Concentration(beadSimulation)
            % G- and F-actin concentration fields and their properties 
            %   
            obj.field = zeros(beadSimulation.N,beadSimulation.N,2);
            obj.field(:,:,1)=beadSimulation.initial_concentration*(1+obj.concentration_randomness...
                *randn(beadSimulation.N,beadSimulation.N));
            switch beadSimulation.depolymerization_constant
                case 0
                    obj.depolymerization_rate = 0;
                otherwise
                    obj.depolymerization_rate = 1/beadSimulation.depolymerization_constant;
            end
            obj.polymerization_rate = beadSimulation.polymerization_constant;
            obj.diffusivity = beadSimulation.diffusivity;
        end



        
        function [obj,midpoint_C]=updateConcentration(obj,delta_function,beadSimulation)
            % Time step the concentration field
            a = beadSimulation.a;
            dt = beadSimulation.dt;
            c = obj.field(:,:,1); % G-actin concentration
            cp = obj.field(:,:,2); %F-actin concentration

            % the following section uses Fourier Transform to solve
            % the pde

            % the equation for G-actin
            w=c+obj.polymerization_rate.*(dt/2)*delta_function.*c...
                +obj.depolymerization_rate*(dt/2)*cp;
            w=fft(w,[],1);
            w=fft(w,[],2);
            uu=a.*w;
            uu=ifft(uu,[],2);
            uu=real(ifft(uu,[],1));

            % the equation for F-actin
            uup=cp-obj.polymerization_rate.*(dt/2)*delta_function.*c-obj.depolymerization_rate*(dt/2)*cp;

            midpoint_C(:,:,1) = uu; midpoint_C(:,:,2) = uup;

            % solve for the next full timestep
            % G-actin equation
            w=c+obj.polymerization_rate.*dt*delta_function.*uu...
                +(dt/2)*(obj.diffusivity)*laplacian(beadSimulation)...
                +(dt/(beadSimulation.depolymerization_constant))*uup;
            w=fft(w,[],1);
            w=fft(w,[],2);
            uuu=a.*w;
            uuu=ifft(uuu,[],2);
            uuu=real(ifft(uuu,[],1)); 
            % F-actin equation
            uuup=cp-obj.polymerization_rate.*dt*delta_function.*uu...
                -obj.depolymerization_rate*dt*uup;

            % update the concentration field
            obj.field(:,:,1) = uuu; obj.field(:,:,2) = uuup;

            function Delta_u=laplacian(beadSimulation)
                % Second order centered Laplacian for periodic domain
                % N  = number of points in each space direction
                %             im = [beadSimulation.N,1:(beadSimulation.N-1)]; %:= circular version of i-1
                %             ip = [2:beadSimulation.N,1];     %:= circular version of i+1
                u = beadSimulation.concentration.field(:,:,1);
                Delta_u=(u(beadSimulation.ip,:,:)+u(beadSimulation.im,:,:)+u(:,beadSimulation.ip,:)+u(:,beadSimulation.im,:)-4*u)/(beadSimulation.h*beadSimulation.h);

            end
        end

        
        
        % this function is currently not being used
        function f=sk(u,g,N,h)
            % Help skew.m get a disretized (u . nabla)u with nice symmetries
            % N  = number of points in each space direction

            
            im = [N,1:(N-1)]; %:= circular version of i-1
            ip = [2:N,1];     %:= circular version of i+1
            f=((u(ip,:,1)+u(:,:,1)).*g(ip,:)...
                -(u(im,:,1)+u(:,:,1)).*g(im,:)...
                +(u(:,ip,2)+u(:,:,2)).*g(:,ip)...
                -(u(:,im,2)+u(:,:,2)).*g(:,im))/(4*h);
        end

    end
end

