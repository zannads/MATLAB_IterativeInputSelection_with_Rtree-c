function draw_R2(R2, ax)
    %UNTITLED4 Summary of this function goes here
    %   Detailed explanation goes here
    if nargin == 1
        figure; 
        ax = gca;
    end
    
    bar(ax, [R2(1, :); diff(R2, 1, 1)] );
    hold(ax, 'on');
    plot(ax,  R2, '-o', 'LineWidth', 2 );
    ylim(ax, [0,1]);
    yticks(ax, 0:0.2:1 );
    ylabel(ax, 'R2 [-]', 'FontSize', 14 );
    grid(ax, 'on');
    xlabel(ax, 'Position of the selected variable', 'FontSize', 14 )
    xticks(ax, 1:length(R2) );
end

