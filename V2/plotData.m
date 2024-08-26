function plotData(matrix, colnames, ylabel_text)
    numrows = size(matrix, 1);
    submatrices = cell(1, numrows / 3);
    for i = 1:length(submatrices)
        start_row = (i - 1) * 3 + 1;
        submatrices{i} = matrix(start_row:start_row + 2, :);
    end
    num_submatrices = numel(submatrices);
    figure;
    for i = 1:num_submatrices
        for col = 1:size(matrix, 2)
            subplot(2, 4, col);
            hold on;
            plot(35:10:55, submatrices{i}(:, col), 'o-', 'DisplayName', sprintf('Doppler Offset: %.2f', 2.5 + (i-1)*2.5));
            xlabel('CN0-dBHz');
            ylabel(ylabel_text, 'Interpreter','latex');
            title(sprintf('%s', colnames{col}));
            grid on;
            legend('Location', 'northoutside');
        end
    end
end