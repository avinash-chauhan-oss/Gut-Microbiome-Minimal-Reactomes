% PLOT 1: Interaction Types Across All Communities (6 Diets)
%
% This script calculates the relative growth change (alpha) for 1326 
% pairwise communities, classifies them into 6 interaction types based 
% on a +/- 0.1 threshold, and plots the WT vs MRM stacked bar charts.


% --- Loading Data ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
% Monoculture File
min_dir = 'C:\Users\Avinash Chauhan\agora_2_models\minReactModels\';
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 
fprintf('Loading datasets...\n');
T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietTitles = {'High Fiber Diet', 'Mediterranean Diet', 'Unhealthy Diet', 'Vegetarian Diet', 'Western Diet', 'No Diet'};

nDiets = numel(dietNames);
nPairs = height(T_comm);

% Creating a dictionary/map for fast monoculture lookups
if any(strcmp(T_mono.Properties.VariableNames, 'ModelNames'))
    modelList = T_mono.ModelNames;
else
    modelList = T_mono.Model;
end
modelMap = containers.Map(modelList, 1:numel(modelList));

% --- Classification Setup ---
% Colors matched exactly to the provided visualization guide
customColors = [
    0.17 0.63 0.17;  % 1: Mutualism (Green)
    0.09 0.75 0.81;  % 2: Commensalism (Teal)
    0.60 0.60 0.60;  % 3: Neutralism (Grey)
    0.95 0.76 0.20;  % 4: Amensalism (Yellow)
    0.90 0.57 0.22;  % 5: Parasitism/Exploitation (Orange)
    0.84 0.15 0.16   % 6: Competition (Red)
];

interactionLabels = {'Mutualism', 'Commensalism', 'Neutralism', 'Amensalism', 'Parasitism', 'Competition'};
threshold = 0.1; 

fprintf('Calculating alpha values and classifying %d pairs across 6 diets...\n', nPairs);

% --- Process and Plots for Each Diet ---
for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    
    % Initialize counts: Rows: 1=WT, 2=Reduced. Cols: 1 to 6 (Interaction Types)
    dietCounts = zeros(2, 6); 
    
    for p = 1:nPairs
        % Safely splitting the community string
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        idx1 = modelMap(m1);
        idx2 = modelMap(m2);
        
        % --- A. Wild-Type (WT) Classification ---
        WT_m1_mono = T_mono{idx1, [prefix '_WT']};
        WT_m2_mono = T_mono{idx2, [prefix '_WT']};
        WT_m1_comm = T_comm{p, ['WT_M1_Growth_' diet]};
        WT_m2_comm = T_comm{p, ['WT_M2_Growth_' diet]};
        
        if ~isnan(WT_m1_comm) && WT_m1_comm >= 0 && ~isnan(WT_m2_comm) && WT_m2_comm >= 0
            WT_a1 = calcAlpha(WT_m1_comm, WT_m1_mono);
            WT_a2 = calcAlpha(WT_m2_comm, WT_m2_mono);
            typeWT = classifyInteraction(WT_a1, WT_a2, threshold);
            dietCounts(1, typeWT) = dietCounts(1, typeWT) + 1;
        end
        
        % --- B. Reduced (MRM) Classification ---
        MRM_m1_mono = T_mono{idx1, [prefix '_MRM']};
        MRM_m2_mono = T_mono{idx2, [prefix '_MRM']};
        MRM_m1_comm = T_comm{p, ['MRM_M1_Growth_' diet]};
        MRM_m2_comm = T_comm{p, ['MRM_M2_Growth_' diet]};
        
        if ~isnan(MRM_m1_comm) && MRM_m1_comm >= 0 && ~isnan(MRM_m2_comm) && MRM_m2_comm >= 0
            MRM_a1 = calcAlpha(MRM_m1_comm, MRM_m1_mono);
            MRM_a2 = calcAlpha(MRM_m2_comm, MRM_m2_mono);
            typeMRM = classifyInteraction(MRM_a1, MRM_a2, threshold);
            dietCounts(2, typeMRM) = dietCounts(2, typeMRM) + 1;
        end
    end
    
    % --- Generating the Stacked Bar Plot for this Diet ---
    fig = figure('Name', ['Interaction Types: ' diet], 'Color', 'w', 'Position', [100+(d*20), 100+(d*20), 650, 500]);
    
    % Creating stacked bar chart
    b = bar([1, 2], dietCounts, 0.45, 'stacked');
    
    % Applying custom colors
    for i = 1:6
        b(i).FaceColor = customColors(i, :);
        b(i).EdgeColor = 'none'; 
    end
    
    % Formatting to match the provided layout template
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'WT', 'Reduced'}, 'FontSize', 12, 'TickDir', 'out');
    ylabel('no. of pairwise communities', 'FontSize', 13);
    
    title(sprintf('Interaction types across all %d communities\n(%s)', nPairs, dietTitles{d}), 'FontSize', 14, 'FontWeight', 'normal');
    
    % Adding Legend (Reversed so Mutualism is at the top of the box)
    leg = legend(flip(b), flip(interactionLabels), 'Location', 'northeastoutside');
    legend boxoff;
    
    % Cleaning up axes
    box off;
    set(gca, 'YGrid', true, 'GridColor', [0.9 0.9 0.9]);
    
    % Saving image safely using 'print'
    export_filename = fullfile(base_dir, sprintf('Plot1_Interactions_%s.png', diet));
    print(fig, export_filename, '-dpng', '-r300');
end

fprintf('All 6 plots generated and saved successfully.\n');


% HELPER FUNCTIONS


function alpha = calcAlpha(gc, gm)
    % Robustly calculates the relative growth change.
    % Handles edge cases where monoculture growth is practically zero.
    if gm < 1e-6
        if gc > 1e-6
            alpha = 1.0; % Positive growth gained from absolute zero 
        else
            alpha = 0.0; % Zero growth to zero growth remains Neutral
        end
    else
        alpha = (gc - gm) / gm;
    end
end

function type = classifyInteraction(a1, a2, th)
    % Categorizes alpha into +, 0, or -
    function c = getCat(a, t)
        if a > t
            c = 1;
        elseif a < -t
            c = -1;
        else
            c = 0;
        end
    end

    c1 = getCat(a1, th);
    c2 = getCat(a2, th);

    % Maps (c1, c2) combinations to interaction types as per the provided table
    if c1 == 1 && c2 == 1
        type = 1; % Mutualism (+, +)
    elseif (c1 == 1 && c2 == 0) || (c1 == 0 && c2 == 1)
        type = 2; % Commensalism (+, 0)
    elseif c1 == 0 && c2 == 0
        type = 3; % Neutralism (0, 0)
    elseif (c1 == -1 && c2 == 0) || (c1 == 0 && c2 == -1)
        type = 4; % Amensalism (-, 0)
    elseif (c1 == 1 && c2 == -1) || (c1 == -1 && c2 == 1)
        type = 5; % Parasitism / Exploitation (+, -)
    elseif c1 == -1 && c2 == -1
        type = 6; % Competition (-, -)
    else
        type = 3; % Safety Fallback
    end
end