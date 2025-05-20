function plotSens(sens)
    % Ensure the input is the correct size
    % if ndims(sens) ~= 4 || size(sens, 4) ~= 10
    %     error('Input must be a 1x1x1x1x1x10 array');
    % end

    % Create a new figure
    figure;
    % Define the number of rows and columns
    numRows = 2;
    numCols = 5;
    hMargin = 0.01; % Horizontal margin between images
    vMargin = 0.02; % Vertical margin between rows     % Calculate width and height for each subplot
    % Calculate width and height for each subplot
    subplotWidth = (1 - (numCols + 1) * hMargin) / numCols;
    subplotHeight = (1 - (numRows + 1) * vMargin) / numRows;

    % Loop through each image and plot it in a 2x5 grid
    for i = 1:10
        % Calculate the position for the current subplot
        row = floor((i-1) / numCols);
        col = mod((i-1), numCols);
       position = [hMargin + col * (subplotWidth + hMargin), ...
                    1 - vMargin - (row + 1) * (subplotHeight + vMargin), ...
                    subplotWidth, subplotHeight];
        % Create a subplot with specified position
        subplot('Position', position);
        imshow(squeeze(abs(sens(:,:,:,i))), []); % Display the image
        % title(['Image ', num2str(i)]); % Add a title
    end
end