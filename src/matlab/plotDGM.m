function [] = plotDGM( I, color, sz, plotDiag )
    hold on;
    if size(I, 1) == 0
        return;
    end
    %Plot a persistence diagram
    if nargin < 2
        color = 'b';
    end
    if nargin < 3
        sz = 20;
    end
    if nargin < 4
        plotDiag = 1;
    end
    
    if plotDiag
        axMin = min(I);
        axMax = max(I);
        axRange = axMax - axMin;
        a = max(axMin - axRange/5, 0);
        b = axMax + axRange/5;
        plot([a, b], [a, b], 'k');
    end
    scatter(I(:, 1), I(:, 2), sz, color, 'fill');
    xlabel('Birth Time');
    ylabel('Death Time');
end

