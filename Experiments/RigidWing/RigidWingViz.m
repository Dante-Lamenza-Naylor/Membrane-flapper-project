function RigidWingViz(basket,params,show_vid,name)

set(0,'DefaultFigureVisible','off');

N = params{1,2};
sigma = params{4,2};
amp = params{6,2};
U   = params{7,2};

% Completely unpack basket
j = 1;
t_sol           = basket.LTM(1,:);
y_sol           = basket.LTM(2,:);
dydt_sol        = basket.LTM(3,:);
dydt2_sol       = basket.LTM(4,:);
tension_sol     = basket.LTM(5,:);
fluid_load_sol  = basket.LTM(6,:);


% Avoid issue of empty cells
nonEmptyCells = ~cellfun('isempty', y_sol); 
y_sol = y_sol(nonEmptyCells);
xs = cos(pi*(0:N)/N)';

% Video stuff
f = figure('Position',[10,10,900,400]);
tiledlayout(2,2,"TileSpacing","tight");
ax1=nexttile(1,[2,1]);
ax2=nexttile(2);
ax3=nexttile(4);

v = VideoWriter(name,'MPEG-4');
v.FrameRate = 30;
v.Quality = 100;
% v.LosslessCompression = true;
open(v);
skip = 20;
skipcount = 0;

% Run video
for j = 2:length(y_sol)
    skipcount = skipcount + 1;
    figure(f);

    if show_vid && mod(skipcount,skip) == 0
        
        t = t_sol{j};
        memys = y_sol{j};
        fluid_load  = fluid_load_sol{j};
        analpresh = AnalyticWingPressure(N,t,2*amp,0,sigma,U);
        analpresh(end) = 0; % Dealing with singularity...

        ghostys = zeros(N+1,1) + amp*sin(sigma*t);
         
        %Membrane
        plot(ax1,xs,memys,'o-black','LineWidth',2,'MarkerSize',5);
        hold(ax1,'on');
        plot(ax1,xs,ghostys,'.-','LineWidth',2)
        hold(ax1,'off');

        % plot(ax2,xs,real(fluid_load./max(abs(fluid_load))),'.-');
        plot(ax2,xs,real(fluid_load),'.-');
        hold(ax2,'on')
        % plot(ax2,xs,analpresh./max(abs(analpresh)),'.-');
        plot(ax2,xs,analpresh,'.-');
        hold(ax2,'off')
        legend(ax2,'numerical','analytic')

       
        plot(ax3,xs,...
            real(fluid_load)-(analpresh),...
            '.-')

        xlabel(ax1,'x/H_m');
        ylabel(ax1,'y/H_m');

        axis(ax1,[-1.1 1.1 -1.1 1.1]);
        axis(ax2,[-1.1 1.1 -100 100]);
        axis(ax3,[-1.1 1.1 -10 10]);

        title(ax1,"wing");
        title(ax2,"fluid loading");
        title(ax3,"loading difference");

        % txt = text(0.1, -0.025, 'Initial Value');
        % set(txt, 'String', num2str(basket.LTM{1,j}))
        % 
        % txt = text(0.1, -0.07, 'Initial Value');
        % set(txt, 'String', num2str(basket.LTM{5,j}))

        % Write to video
        frame = getframe(f);
        writeVideo(v,frame);

    end
end

% Close video
close(v);
set(0,'DefaultFigureVisible','on');
close all;

