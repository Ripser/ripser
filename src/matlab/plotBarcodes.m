function [] = plotBarcodes( I )
    hold on;
    s = (max(I(:)) - min(I(:)))/size(I, 1);
    for ii = 1:size(I, 1)
        b = I(ii, 1);
        d = I(ii, 2);
        plot([b, d], [s*ii, s*ii], 'color', 'k');
        scatter([b, d], [s*ii, s*ii], 20, 'k', 'fill');
    end
end

