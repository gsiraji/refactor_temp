classdef MovieMaker
    %MOVIEMAKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frame_rate = 5 % Video frame rate
        frame_count = 30% number of video frames
        step4video = 10000 % Video frame frequency^-1
        movie = 0 %0: no movie, 3: 2 movies G&F
        v
        movie_frame_number = 1
        grid % ygrid
        video_ID
    end
    
    methods
        function obj = MovieMaker(movie_parameters,beadSimulation)
            %MOVIEMAKER Construct an instance of this class
            %   Detailed explanation goes here clockmax,movie,video_ID,frame_rate
            obj.grid = repmat((0:1:(beadSimulation.N-1)).*beadSimulation.h,[beadSimulation.N,1]); % ygrid
            obj.movie = movie_parameters.movie;
            obj.frame_count = beadSimulation.sample_count;
            obj.step4video = ceil(beadSimulation.clockmax/beadSimulation.sample_count);% Video frame frequency^-1
            if obj.movie ~=0
                obj.video_ID = strjoin({movie_parameters.name,...
                    num2str(movie_parameters.trial_number),...
                    movie_parameters.variable_of_interest},'_');
            end

        end
        
        function obj = makeMovie(obj)
            movieTitle = obj.video_ID;
            switch obj.movie
                case 1
                    % only make the G-actin movie
                    obj.v = VideoWriter(strcat('G',movieTitle),'MPEG-4');
                    obj.v.Quality = 100;
                    obj.v.FrameRate = obj.frame_rate;
                    open(obj.v)
                    % only make the F-actin movie
                case 2
                    obj.v = VideoWriter(strcat('F',movieTitle),'MPEG-4');
                    obj.v.Quality = 100;
                    obj.v.FrameRate = obj.frame_rate;
                    open(obj.v)
                case 3
                    % make both the G- and F-actin movies
                    v1 = VideoWriter(strcat('G',movieTitle),'MPEG-4');
                    v1.Quality = 100;
                    v1.FrameRate = obj.frame_rate;
                    open(v1)
                    v2 = VideoWriter(strcat('F',movieTitle),'MPEG-4');
                    v2.Quality = 100;
                    v2.FrameRate = obj.frame_rate;
                    open(v2)
                    obj.v = {v1, v2};
                case 0
                    obj.v = 0;
            end
        end

        function obj = makeMovieFrame(obj,beadSimulation, Step)

            u = beadSimulation.concentration.field;
            for i=1:beadSimulation.bead_count
                X(i,:) = beadSimulation.beadList(i).position;
            end
            L = beadSimulation.L;
            dt = beadSimulation.dt;

            if obj.movie > 0
                switch obj.movie
                    case 3

                        figure(1)
                        plotFrame(obj,u(:,:,1),X,L,dt,Step)
                        F1(obj.movie_frame_number) = getframe(gcf);
                        writeVideo(obj.v{1},F1(obj.movie_frame_number));
                        set(gca,  'CLim', [0 40]);
                        close all;

                        figure(2)
                        plotFrame(obj,u(:,:,2),X,L,dt,Step)
                        F2(obj.movie_frame_number) = getframe(gcf);
                        writeVideo(obj.v{2},F2(obj.movie_frame_number));
                        set(gca,  'CLim', [0 40]);
                        close all;
                    otherwise
                        figure(1)
                        plotFrame(obj,u(:,:,obj.movie),X,L,dt,Step)
                        set(gca,  'CLim', [32.98 33]);
                        F1(obj.movie_frame_number) = getframe(gcf);
                        writeVideo(obj.v,F1(obj.movie_frame_number));
                        
                        close all;
                end
                
            end

        end
        
        function plotFrame(obj,u,X,L,dt,Step)

            
            pcolor(10^4*obj.grid',10^4*obj.grid,10^9*u)
            hold on
            Xplt = mod(X,L);
            plot(10^4*Xplt(:,1),10^4*Xplt(:,2),'ko','MarkerSize',7,'MarkerFaceColor','k')
%             plot(10^4*Xplt(1,1),10^4*Xplt(1,2),'ro','MarkerFaceColor','r','MarkerSize',5)

            axis([0,10^4*L,0,10^4*L])
            hcb=colorbar;
            H = ylabel(hcb, 'uM');
            H.Rotation = 270;
            H.Position = [4.0813   20.0000         0];
            xlabel('um'); ylabel('um');
            title(['t is ' num2str(dt*Step/60,'%.1f') ' min']);
            axis equal
            axis manual
            set(gca, 'fontsize', 18);
            drawnow
        end
        function obj = closeMovie(obj)

            if (obj.movie > 0)
                switch obj.movie
                    case 3
                        for movie_indx = 1:length(obj.v)
                            close(obj.v{movie_indx})
                        end
                    otherwise
                        close(obj.v)
                end
            end

        end

    end
end








