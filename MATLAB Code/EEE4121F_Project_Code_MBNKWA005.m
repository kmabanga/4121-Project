clear; clc; close all;

% Normalized decision matrix
% Criteria: Cost per Byte, Available Bandwidth, Delay
normalized_matrix = [
    15/60 2/60 50/50; % 3G (UMTS)
    15/15 15/60 50/140% WLAN (802.11n)
    15/50 60/60 50/100; % 4G (LTE)
    ];

% Lambda value: Weighting factor for combining WSM and WPM
lambda = 0.5;

% Sample size
N = 1000;

%  Priority Cases
%   Low (1-3)
%   Moderate (4-6)
%   High (9-10)
%   No priority (Full Range) (1-10)
cases = {[1, 3];  [4, 7];  [9, 10]; [1, 10]};

% Array to store selection count for each RAT
selection_counts = zeros(4, 3);

% Criterion labels
% Change the criterion labels
criteria_labels = {'Cost', 'Bandwidth', 'Delay'};

% RAT labels
% change RAT-labels to available RATs
RAT_labels = {'3G (UMTS)', 'WLAN (802.11n)', '4G (LTE)'};

% Colors for the cases
colors = ['b', 'r', 'g', 'y'];

% 7. Effect of the Users' Weight Assigned to the Price on RAT-Selection Decisions
% 8. Effect of the Users' Weight Assigned to the Available Data Rate on RAT-Selection Decisions
% 9. Effect of the Users' Weight Assigned to the Delay on RAT-Selection Decisions

% Loop over each criterion
for criterion = 1:3
    % Loop over each case
    for c = 1:4
        % Initialize counts to zero for the current case
        count_1 = 0; count_2 = 0; count_3 = 0;

        % Loop over the sample size (N/4 iterations)
        for i = 1:N / 4
            % Generate random user weights based on the current criterion and case
            if criterion == 1  % Criterion 1 = Cost per Byte
                user_weights = [randi(cases{c}), randi([1, 10]), randi([1, 10])];
            elseif criterion == 2  % Criterion 2 = Available Data Rate
                user_weights = [randi([1, 10]), randi(cases{c}), randi([1, 10])];
            else  % Criterion 3 = Delay
                user_weights = [randi([1, 10]), randi([1, 10]), randi(cases{c})];
            end

            % Normalize user weights
            normalized_weights = user_weights / sum(user_weights);

            % Compute WASPAS scores
            % RAT-1 = 3G (UMTS)
            row_1 = normalized_matrix(1, :);
            WSM_1 = dot(row_1, normalized_weights);
            WPM_1 = prod(row_1 .^ normalized_weights);
            WASPAS_1 = lambda * WSM_1 + (1 - lambda) * WPM_1;

            % RAT-2 = WLAN (802.11n)
            row_2 = normalized_matrix(2, :);
            WSM_2 = dot(row_2, normalized_weights);
            WPM_2 = prod(row_2 .^ normalized_weights);
            WASPAS_2 = lambda * WSM_2 + (1 - lambda) * WPM_2;

            % RAT-3 = 4G (LTE)
            row_3 = normalized_matrix(3, :);
            WSM_3 = dot(row_3, normalized_weights);
            WPM_3 = prod(row_3 .^ normalized_weights);
            WASPAS_3 = lambda * WSM_3 + (1 - lambda) * WPM_3;

            % Select the RAT with the highest WASPAS score
            [~, idx] = max([WASPAS_1, WASPAS_2, WASPAS_3]);

            % Update selection count for the selected RAT
            if idx == 1
                count_1 = count_1 + 1; % increment count for 3G (UMTS)
            elseif idx == 2
                count_2 = count_2 + 1; % increment count for WLAN (802.11n)
            else
                count_3 = count_3 + 1; % increment count for 4G (LTE)
            end
        end

        % Store counts for the current priority level
        selection_counts(c, :) = [count_1, count_2, count_3];
    end

    % Plot the results for the current criterion
    figure;
    hold on;

    % Create grouped bar chart
    bar_handles = bar(categorical(RAT_labels), selection_counts', 'grouped');

    % Set color for each priority level
    for c = 1:4
        bar_handles(c).FaceColor = colors(c);
    end

    % Add labels on top of bars
    for i = 1:length(bar_handles)
        x = bar_handles(i).XEndPoints;
        y = bar_handles(i).YEndPoints;
        labels = string(bar_handles(i).YData);
        text(x, y, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    xlabel('RATs');
    ylabel('No. of calls');

    title(['WASPAS Selection Frequency for Different ', criteria_labels{criterion}, ' Weight Ranges']);

    legend({'Low Priority', 'Moderate Priority', 'High Priority', 'No Priority'});

    grid on;

    hold off;
end

% 10. Effect of the Proportion of the WSM and the WPM on RAT-Selection Decisions

% Lambda values to test
%   0 (WPM only)
%   0.25 (25% WSM, 75% WPM)
%   0.5 (Balanced WSM and WPM)
%   0.75 (75% WSM, 25% WPM)
%   1 (WSM only)
lambda_values = [0, 0.25, 0.5, 0.75, 1];

% Array to store counts for each RAT for each lambda case
selection_counts = zeros(5, 3);

% Iterate over each lambda case
for l = 1:5
    % Initialize the counts to zero
    count_1 = 0; count_2 = 0; count_3 = 0;

    % Current lambda value
    lambda = lambda_values(l);

    for i = 1:N / 5
        % Generate random user weights
        user_weights = [randi([1, 10]), randi([1, 10]), randi([1, 10])];

        % Normalize user weights
        normalized_weights = user_weights / sum(user_weights);

        % Compute WASPAS scores
        % RAT-1 = 3G (UMTS)
        row_1 = normalized_matrix(1, :);
        WSM_1 = row_1(1) * normalized_weights(1) + row_1(2) * normalized_weights(2) + row_1(3) * normalized_weights(3);
        %   -1 * normalized_weights(i) - cost criteria
        %    1 * normalized_weights(i) - benefit criteria
        WPM_1 = (row_1(1)^(1 * normalized_weights(1))) * (row_1(2)^(normalized_weights(2))) * (row_1(3) ^ (1 * normalized_weights(3)));
        WASPAS_1 = lambda * WSM_1 + (1 - lambda) * WPM_1;

        % RAT-2 = WLAN (802.11n)
        row_2 = normalized_matrix(2, :);
        WSM_2 = row_2(1) * normalized_weights(1) + row_2(2) * normalized_weights(2) + row_2(3) * normalized_weights(3);
        %   -1 * normalized_weights(i) - cost criteria
        %    1 * normalized_weights(i) - benefit criteria
        WPM_2 = (row_2(1)^(1 * normalized_weights(1))) * (row_2(2)^(normalized_weights(2))) * (row_2(3) ^ (1 * normalized_weights(3)));
        WASPAS_2 = lambda * WSM_2 + (1 - lambda) * WPM_2;

        % RAT-3 = 4G (LTE)
        row_3 = normalized_matrix(3, :);
        WSM_3 = row_3(1) * normalized_weights(1) + row_3(2) * normalized_weights(2) + row_3(3) * normalized_weights(3);
        % change according to your criteria
        %   -1 * normalized_weights(i) - cost criteria
        %    1 * normalized_weights(i) - benefit criteria
        WPM_3 = (row_3(1)^(1 * normalized_weights(1))) * (row_3(2)^(normalized_weights(2))) * (row_3(3) ^ (1 * normalized_weights(3)));
        WASPAS_3 = lambda * WSM_3 + (1 - lambda) * WPM_3;
        % Select the best RAT
        [~, idx] = max([WASPAS_1, WASPAS_2, WASPAS_3]);

        % Update count for the selected RAT
        if idx == 1
            count_1 = count_1 + 1; % increment count for 3G
        elseif idx == 2
            count_2 = count_2 + 1; % increment count for 4G
        else
            count_3 = count_3 + 1; % increment count for WLAN
        end
    end

    % Store counts for the current lambda case
    selection_counts(l, :) = [count_1, count_2, count_3];
end

% Plot the results
colors = ['b', 'c', 'r', 'g', 'y'];

figure;

hold on;

% Create grouped bar chart
bar_handles = bar(categorical(RAT_labels), selection_counts', 'grouped');

% Set color for each lambda case
for l = 1:5
    bar_handles(l).FaceColor = colors(l);
end
% Add labels on top of bars
for i = 1:length(bar_handles)
    x = bar_handles(i).XEndPoints;
    y = bar_handles(i).YEndPoints;
    labels = string(bar_handles(i).YData);
    text(x, y, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end


xlabel('RATs');
ylabel('No. of calls');

title('WASPAS Selection Frequency for Different Lambda Values');

legend('${\lambda_1} = 0$','${\lambda_2} = 0.25$','${\lambda_3} = 0.5$','${\lambda_4} = 0.75$','${\lambda_5} = 1$', 'Interpreter','latex');

grid on;

hold off;

