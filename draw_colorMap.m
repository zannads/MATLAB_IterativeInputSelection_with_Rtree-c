function draw_colorMap(X, R2)
    figure; 
    tiledlayout(2,1);
    nexttile;
    imagesc(X); colormap('Jet')
    ylabel('Input ranking', 'FontSize', 16);
    xticks( 1:size(X,2) )
    nexttile; 
    plot(max(R2), 'k', 'LineWidth', 2); grid on;
    xlim([0.5 size(X,2)+0.5] );
    xticks( 1:size(X,2) )
    ylim( [fix(min(max(R2))*100)/100 fix(max(max(R2))*100+1)/100]);
    ylabel('Model performance R2', 'FontSize', 16);
    xlabel('IIS run');
end
