classdef BeadSimulation
    %BEADSIMULATIN Summary of this class goes here
    %   Detailed explanation goes here

    properties
       
        diffusivity = 1.3e-7; % diffusion coeff (13 um^2/sec);
        delta_function_mode = 6; %only 6 works right now
        bead_distance=0.025

        simulation_ID = string(datetime('now','format','yyyyMMdd'));
        simulation_mode = '1';

        L = 0.04 % Box size
        N = 32; % Number of grid cells
        h = 0.04/32; % Grid spacing
        
        bead_size = 0.04/32;
        bead_count = 1; % Number of particles
        polymerization_constant = 2e-5;
        repulsion_coefficient = 5e-2;%Interaction/repulsion coefficient
        interaction_range=0.04/32*2; %Interaction/repulsion range
        alpha = 1e5; %get rid of alpha
        depolymerization_constant = 100; % 100s? Depolymerization timescale
        initial_speed = 5e-5;
        mass = 5e-5/(3e4); %?% effective mass coefficient  (um^3/s)
        inertia = (0.04/32)/(7e3)*5e-5;% effective rotational inertia coefficient   (um^2/s)
        initial_concentration = 3.3*10^-8
        time_max = 6; % Run until time
        dt = 0.01; % Time step instability in speed when dt=>0.04
        clockmax = ceil(6/0.01);
        sample_count = 30 % how often save data/how many video frames
        rng_state
        
        beadList
        
        a=ones(32,32);
        concentration 
        k_c = 1;
        force_type = 'Lennard-Jones';
        model = 'stuck beads'
        ip = [(2:32),1]; % Grid index shifted left
        im = [32,(1:(32-1))]; % Grid index shifted right
        kp = [(2:1),1]; % IB index shifted left
        km = [1,(1:(1-1))]; % IB index shifted right
        coupling_constant=-2e-5*ones(1,1);



    end

    methods
        function obj = BeadSimulation(varargin)
            % either go through and assign prameters from the set to this
            % or just assing the parameter set to something that we'll
            % reference
            % Set default values
            defaults = {'diffusivity', 1.3*10^-7;...
                'delta_function_mode', 6;...
                'simulation_mode', '1';...
                'L', 0.04;...
                'N', 32;...
                'h',0.04/32;...
                'dt',0.01;...
                'bead_size', 0.04/32;...
                'bead_count', 1;...
                'clockmax',1000;
                'polymerization_constant', 2e-5;...
                'repulsion_coefficient', 1;
                'interaction_range', 0.04/32;
                'depolymerization_constant', 100;
                'initial_speed', 25;
                'mass', 5e-5/(3e4);
                'inertia', (0.04/32)/(7e3)*(5e-5);
                'time_max', 6;
                'sample_count',30;
                'force_type', 'Lennard-Jones';
                'model', 'stuck beads';
                'initial_concentration', 3.3*10^-8;
                'beadList',StuckBead();
                'bead_distance',0.025
                };

            % Parse inputs
            p = inputParser;
            p.CaseSensitive = false;
            for i = 1:1:size(defaults,1)
                addParameter(p, defaults{i,1}, defaults{i,2})
            end
            parse(p, varargin{:})
            
            % Assign parsed values to object properties
            fields = fieldnames(p.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = p.Results.(fields{i});
            end

            rng shuffle; rng_struct = rng();
            obj.rng_state = rng_struct.Seed; % setting run
%             obj.simulation_ID = str2double(strcat(num2str(obj.bead_count), num2str(obj.N), obj.simulation_mode, num2str(obj.rng_state))); % Simulation number

            obj.h = obj.L/obj.N; % Grid spacing
            obj.a = create_a(obj);
            obj.ip = [(2:obj.N),1]; % Grid index shifted left
            obj.im = [obj.N,(1:(obj.N-1))]; % Grid index shifted right
            obj.kp = [(2:obj.bead_count),1]; % IB index shifted left
            obj.km = [obj.bead_count,(1:(obj.bead_count-1))]; % IB index shifted right
            obj.coupling_constant=-obj.polymerization_constant*ones(obj.bead_count,1);
            obj.dt = (0.01/(obj.L/32)^2)*obj.h^2; % Time step instability in speed when dt=>0.04
            obj.clockmax = ceil(obj.time_max/obj.dt);
            obj.k_c = obj.bead_size/obj.h;
%             obj.interaction_range =  obj.bead_size*2;
            obj.concentration = Concentration(obj);

        end

        function a = create_a(obj)
            a = ones(obj.N);
            for m1=0:(obj.N-1)
                for m2=0:(obj.N-1)
                    t=(pi/obj.N)*[m1;m2];
                    s=sin(t);
                    a(m1+1,m2+1)=a(m1+1,m2+1)...
                        /(1+(obj.dt/2)*(obj.diffusivity)*(4/(obj.h^2))*(s'*s));
                end
            end

        end


        function run(obj,my_logger,my_movie_maker)
            
            my_logger.logInfo(obj);

            for time_step=1:obj.clockmax
                
                % advance the beads to time+(dt/2)
                obj = halfTimeStepBeads(obj,obj.concentration.field);
                
                % update the concentration field to time+dt
                [obj.concentration,half_step_concentration]=obj.concentration.updateConcentration(obj.spread6,obj);

                 % advance the beads to time+(dt)
                obj = halfTimeStepBeads(obj,half_step_concentration);

                if (mod(time_step,my_logger.sampling_step) == 0)

                    my_movie_maker = my_movie_maker.makeMovieFrame(obj, time_step);
                    
                    my_logger = my_logger.logData(obj);

                    my_logger = my_logger.nextEntry();

                    my_logger.progressMessage()
                    

                end
            end

            my_movie_maker.closeMovie;

            my_logger.saveData()
            
            
        end

        function obj = halfTimeStepBeads(obj,u)
            % get the position of the beads
            for i=1:obj.bead_count
                obj.beadList(i).position = mod(obj.beadList(i).position,obj.L);
            end
                % update velocity
                obj=obj.repulsion;
                
                direction = -unitGrad(obj,u(:,:,1));
                direction = interp6(obj,direction);

                % calculate the concentration at bead position
                concentration = interp6(obj,u);
                concentration = concentration(:,1);
                % calculate the gradient at the position of beads
                concentration_gradient = grad(obj,u(:,:,1));
                concentration_gradient = interp6(obj,concentration_gradient);

                for i=1:obj.bead_count
                    % update velocity
                    obj.beadList(i)=obj.beadList(i).updateVelocity('concentration',concentration(i,:),...
                        'dt',obj.dt/2,'direction',direction(i,:));
                    % update position
                    obj.beadList(i)=obj.beadList(i).updatePosition(obj.dt/2);

                    % update the anglar velocity
                    obj.beadList(i)=obj.beadList(i).updateAngularVelocity(concentration_gradient(i,:),obj.dt/2);
                    % update the orientation
                    obj.beadList(i)=obj.beadList(i).updateOrientation(obj.dt/2);
                end


                
                

        end

        function value_at_X = vecInterpolate(obj,u,X)
            % Interpolated the grid point value to position of the bead
            % 
            value_at_X=zeros(obj.bead_count,2);

            % Get body position relative to grid % now scaling by 'k_c'
            switch nargin
                case 3
                    [central_index,w] = makeSpreadMatrix(obj,X);
                case 2
                    [central_index,w] = makeSpreadMatrix(obj);
            end
            
            for k=1:obj.bead_count

                i1 = indexCalculation(obj,central_index(k,1));
                i2 = indexCalculation(obj,central_index(k,2));

                ww = w(:,:,k);
                UU = [sum(sum((obj.k_c).*ww.*u(i1,i2,1))), sum(sum((obj.k_c).*ww.*u(i1,i2,2)))]; %Interpolate
                value_at_X(k,:)= UU;

            end
        end

        function delta_function = vecSpread(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            c=1/(obj.bead_size^2);
            delta_function=zeros(obj.N);

            [central_index,w] = makeSpreadMatrix(obj);
            
            
            for k=1:obj.bead_count

                i1 = indexCalculation(obj,central_index(k,1));
                i2 = indexCalculation(obj,central_index(k,2));

                ww = w(:,:,k);
                delta_function(i1,i2)=delta_function(i1,i2)+(c*obj.coupling_constant(k))*ww;

            end
        end

        function C =grad(obj,u)
            % calculate the vector gradient and return it as one array
            % [ux, uy] = gradient(u(:,:,1));
            % u = obj.concentration.field(:,:,1);

            ux = (u(obj.ip,:,1)-u(obj.im,:,1))/(2*obj.h);
            uy = (u(:,obj.ip,1)-u(:,obj.im,1))/(2*obj.h);
            C(:,:,1) = ux;C(:,:,2) = uy;
        end

        function [central_index,w] = makeSpreadMatrix(obj,X,varargin)

            % Get body position relative to grid % now scaling by 'k_c'
            switch nargin
                case 2
                    s = X/obj.h;
                case 1
                    s=obj.beadList.position/obj.h; 
            end

            central_index=floor(s);
            r=(s-central_index)/obj.k_c;

            w=vecPhiX(r(:,1)).*vecPhiY(r(:,2),obj.k_c);%Evaluate delta function % returns...
            w = permute(w, [1,3,2]); %Reogranize, this is quite fast


        end

        function U=interp6(obj,u)
%             u = obj.concentration.field;
            X = ones(obj.bead_count,2);
            for i = 1:obj.bead_count
                X(i,:) = obj.beadList(i).position;
            end

            U=zeros(size(X,1),2);
            s=X/obj.h; % Get body position relative to grid
            i=floor(s);
            r=s-i;
            w=phi6X(r(:,1)).*phi6Y(r(:,2));

            w = permute(w, [1,3,2]); %Reogranize, this is quite fast

            for k=1:size(X,1)
                i1=mod((i(k,1)-3):(i(k,1)+2),obj.N)+1; % Find adjacent fluid cells
                i2=mod((i(k,2)-3):(i(k,2)+2),obj.N)+1;

                ww = w(:,:,k);
                UU = [sum(sum(ww.*u(i1,i2,1))), sum(sum(ww.*u(i1,i2,2)))]; %Interpolate
                U(k,:)= UU;
            end
        end

            function f=spread6(obj)
%                 global N h;
                F = obj.coupling_constant;
                X = ones(obj.bead_count,2);
                for i = 1:obj.bead_count
                    X(i,:) = obj.beadList(i).position;
                end

                c=1/(obj.h*obj.h);
                f=zeros(obj.N,obj.N);
                s=X/obj.h; % Get body position relative to grid
                i=floor(s);
                r=s-i;
                w=phi6X(r(:,1)).*phi6Y(r(:,2));

                w = permute(w, [1,3,2]); %Reogranize, this is quite fast

                for k=1:obj.bead_count
                    i1=mod((i(k,1)-3):(i(k,1)+2),obj.N)+1; % Find adjacent fluid cells
                    i2=mod((i(k,2)-3):(i(k,2)+2),obj.N)+1;

                    ww = w(:,:,k);
                    f(i1,i2)=f(i1,i2)+(c*F(k))*ww;
                end
            end

            function n=unitGrad(obj,u)

                n1=(u(obj.ip,:)-u(obj.im,:))./(2*obj.h);
                n2=(u(:,obj.ip)-u(:,obj.im))./(2*obj.h);
                norm = sqrt(n1.*n1+n2.*n2)+eps;

                n1 = n1./norm; n2 = n2./norm;
                n(:,:,1)=n1; n(:,:,2)=n2;
            end




        function w=vecPhiX(r,k_c)
            % 6 point delta function in y

            w=zeros(6,size(r,1),6);
            % calculate the 6 cells from 6pt delta function
            [w0,~] = deltaCK(r);
            % make it 6x6xNb
            ii=0;
            for i = 1:k_c:6*k_c
                ii = ii +1;
                w(i:i+k_c-1,:,:)=repmat(w0(:,ii),[k_c,1,6*k_c]);
            end
        end



        function w=vecPhiY(r,k_c)
            % 6 point delta function in x
            w=zeros(6*k_c,size(r,1),6*k_c);
            % calculate the 6 cells from 6pt delta function
            [w0,~] = deltaCK(r);
            ii = 0;
            % make it 6x6xNb
            for i = 1:k_c:6*k_c
                ii = ii +1;
                w(:,:,i:i+k_c-1)=repmat(w0(:,ii)',[6*k_c,1,k_c]);
            end
        end


        function [w,flag] = deltaCK(r)
            %deltaCK:  6pt delta function with K=constant second moment

            w=zeros(length(r),6);

            K_delta=(59/60)*(1-sqrt(1-(3220/3481)));

            Alpha=28;
            beta=(9/4)-(3/2)*(K_delta+r.^2)+((22/3)-7*K_delta)*r-(7/3)*r.^3;
            GAMMA=(1/4)*(((161/36)-(59/6)*K_delta+5*K_delta^2)*(1/2)*r.^2 + (-(109/24)+5*K_delta)*(1/3)*r.^4 + (5/18)*r.^6 );

            discr=beta.^2-4*Alpha*GAMMA;
            flag=0;
            if(any(discr<0))
                flag=1
            end

            pm3=(-beta+sign((3/2)-K_delta)*sqrt(discr))/(2*Alpha);
            w(:,6)=pm3;
            w(:,5)= -3*pm3 - (1/16) + (1/8)*(K_delta+r.^2) + (1/12)*(3*K_delta-1)*r + (1/12)*r.^3;
            w(:,4)=  2*pm3 + (1/4)                   +  (1/6)*(4-3*K_delta)*r -  (1/6)*r.^3;
            w(:,3)=  2*pm3 + (5/8)  - (1/4)*(K_delta+r.^2);
            w(:,2)= -3*pm3 + (1/4)                   -  (1/6)*(4-3*K_delta)*r +  (1/6)*r.^3;
            w(:,1)=    pm3 - (1/16) + (1/8)*(K_delta+r.^2) - (1/12)*(3*K_delta-1)*r - (1/12)*r.^3;

        end

        

        function n=unit_grad(obj)
            
            u = obj.concentration.field(:,:,1);

            n1=(u(obj.ip,:)-u(obj.im,:))./(2*obj.h);
            n2=(u(:,obj.ip)-u(:,obj.im))./(2*obj.h);
            norm = sqrt(n1.*n1+n2.*n2)+eps;

            n1 = n1./norm; n2 = n2./norm;
            n(:,:,1)=n1; n(:,:,2)=n2;
        end


        

        function indices = indexCalculation(obj,middle_index)
           
            switch obj.delta_function_mode
                case 6
                    indices=mod((middle_index-(3*obj.k_c)+1):(middle_index+(3*obj.k_c)),obj.N); %Find affected cells
                case 4
                    indices = mod((middle_index-2):(middle_index+1),obj.N)+1;
            end
        
        end

        function obj = repulsion(obj)
            % update the positions based on periodic boundaries
            X = ones(obj.bead_count,2);
            for i = 1:obj.bead_count
                X(i,:) = obj.beadList(i).position;
            end
            X = mod(X,obj.L);
            acceleration_vec = 0.*X;

            % forceType 1: inspired by a paper on modeling attraction
            % and repulsion in (social) collective behavior
            % forceType 2: (chemical) Lennard-Jones force
            switch obj.force_type
                case '1/r'
                    for j = 1:size(X,1)
                        [Idx,D] = rangesearch(X(:,:),X(j,:),obj.interaction_range);
                        Idx = Idx{1};
                        if nnz(Idx ~= j) ~=0
                            D = D{1};
                            D = D(Idx ~= j);Idx = Idx(Idx ~= j);
                            acceleration_vec(j,:)= - obj.repulsion_coefficient*sum((X(Idx,:)-X(j,:))./(D'.^2),1);
                        end
                    end
                case 'Lennard-Jones'
                    for j = 1:size(X,1)
                        % Lennard-Jones Potential
                        [Idx,D] = periodicDist(X(j,:),X(:,:),obj.L,obj.interaction_range);
                        if nnz(Idx ~= j) ~=0
                            D = D(Idx ~= j);Idx = Idx(Idx ~= j);
                            D(D <= obj.h/2 ) = obj.h/2;
                            acceleration_vec(j,:)=(- obj.repulsion_coefficient*(2.4e-41)*sum((X(Idx,:)-X(j,:)).*...
                                ((2./D.^(14))-(1./D.^8)),1));
                        end
                    end
            end
            for i=1:obj.bead_count
                obj.beadList(i).acceleration = acceleration_vec(i,:);
            end    
            function [Idx2,D2] = periodicDist(ZI,ZJ,L,interaction_range,varargin)
            % calculate the periodic distance
            switch nargin
                case 4
                    dx1 = abs(ZI(:,1)-ZJ(:,1));
                    dx1(dx1 > L/2) = L - dx1(dx1 > L/2);
                    dy1 = abs(ZI(:,2)-ZJ(:,2));
                    dy1(dy1 > L/2) = L - dy1(dy1 > L/2);

                    D2 = sqrt(dy1.^2+dx1.^2);
                    Idx2 = 1:length(ZJ);
                    Idx2 = Idx2(D2 <= interaction_range);
                    D2 = D2(D2 <= interaction_range);
                case 3
                    dx1 = abs(ZI(:,1)-ZJ(:,1));
                    dx1(dx1 > L/2) = L - dx1(dx1 > L/2);
                    dy1 = abs(ZI(:,2)-ZJ(:,2));
                    dy1(dy1 > L/2) = L - dy1(dy1 > L/2);

                    D2 = sqrt(dy1.^2+dx1.^2);
                    Idx2 = 1:length(ZJ);
            end
            end

        end

        
%         
%         function delta_C_vec = calculateDeltaC(obj,query_distance)
%             % calculate the concentration difference across the beads
%             X = obj.beadList.position;
%             u = obj.concentration.field(:,:,1);
%             delta_C_vec = ones(1,size(X,2));
%             
%             for indx = 1:size(X,2)
% 
%                 x_left =  X(indx,1)-query_distance;
%                 x_right =  X(indx,1)+query_distance;
%                 y = X(indx,2);
%                 c_left = vecInterpolate(u,[x_left,y]);
%                 c_right = vecInterpolate(u,[x_right,y]);
%                 delta_C_vec(1,indx) = c_left - c_right;
% 
%             end
% 
%         end


    end
end


