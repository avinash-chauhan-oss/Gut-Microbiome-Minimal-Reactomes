% PLOT 3: Interaction Composition Across Diets
%
% Generates a single figure with side-by-side stacked bars (WT left, 
% Reduced right) for each diet to compare the overall shift in interaction 
% types across different nutritional environments.


% --- Loading Data ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
min_dir  = fullfile(base_dir, 'minReactModels');
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 

fprintf('Loading datasets for multi-diet plot...\n');
T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);


dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietLabels = {'Diet 1\n(HF)', 'Diet 2\n(MED)', 'Diet 3\n(UNH)', 'Diet 4\n(VEG)', 'Diet 5\n(WES)', 'Diet 6\n(NOD)'};

nDiets = numel(dietNames);
nPairs = height(T_comm);

if any(strcmp(T_mono.Properties.VariableNames, 'ModelNames'))
    modelList = T_mono.ModelNames;
else
    modelList = T_mono.Model;
end
modelMap = containers.Map(modelList, 1:numel(modelList));

% --- Classification Setup ---
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

% Master matrix: [WT_Diet1, MRM_Diet1, WT_Diet2, MRM_Diet2, ...]
% Rows: 12 bars total (2 per diet). Cols: 6 interaction types.
masterCounts = zeros(nDiets * 2, 6);

fprintf('Calculating alpha values and classifying %d pairs across %d diets...\n', nPairs, nDiets);

% --- Processing Data ---
for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    
    wt_idx  = (d - 1) * 2 + 1;
    mrm_idx = (d - 1) * 2 + 2;
    
    for p = 1:nPairs
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        i1 = modelMap(m1);
        i2 = modelMap(m2);
        
        % --- WT Classification ---
        gm1_wt = T_mono{i1, [prefix '_WT']};
        gm2_wt = T_mono{i2, [prefix '_WT']};
        gc1_wt = T_comm{p, ['WT_M1_Growth_' diet]};
        gc2_wt = T_comm{p, ['WT_M2_Growth_' diet]};
        
        if ~isnan(gc1_wt) && gc1_wt >= 0 && ~isnan(gc2_wt) && gc2_wt >= 0
            a1_wt = calcAlpha(gc1_wt, gm1_wt);
            a2_wt = calcAlpha(gc2_wt, gm2_wt);
            typeWT = classifyInteraction(a1_wt, a2_wt, threshold);
            masterCounts(wt_idx, typeWT) = masterCounts(wt_idx, typeWT) + 1;
        end
        
        % --- MRM Classification ---
        gm1_mrm = T_mono{i1, [prefix '_MRM']};
        gm2_mrm = T_mono{i2, [prefix '_MRM']};
        gc1_mrm = T_comm{p, ['MRM_M1_Growth_' diet]};
        gc2_mrm = T_comm{p, ['MRM_M2_Growth_' diet]};
        
        if ~isnan(gc1_mrm) && gc1_mrm >= 0 && ~isnan(gc2_mrm) && gc2_mrm >= 0
            a1_mrm = calcAlpha(gc1_mrm, gm1_mrm);
            a2_mrm = calcAlpha(gc2_mrm, gm2_mrm);
            typeMRM = classifyInteraction(a1_mrm, a2_mrm, threshold);
            masterCounts(mrm_idx, typeMRM) = masterCounts(mrm_idx, typeMRM) + 1;
        end
    end
end

% --- Plotting ---
fig = figure('Name', 'Plot 3: Multi-Diet Composition', 'Color', 'w', 'Position', [100, 100, 1000, 500]);

% To get grouped AND stacked bars in MATLAB, we format the X coordinates
% Manually define X positions to group the 2 bars per diet tightly together,
% while leaving a larger gap between different diets.
x_pos = zeros(1, nDiets * 2);
current_x = 1;
for d = 1:nDiets
    x_pos((d-1)*2 + 1) = current_x;       % WT bar position
    x_pos((d-1)*2 + 2) = current_x + 0.8; % MRM bar position (close to WT)
    current_x = current_x + 2.5;          % Jump to next diet block
end

% Draw the stacked bars
b = bar(x_pos, masterCounts, 0.8, 'stacked');

% Applying custom colors
for i = 1:6
    b(i).FaceColor = customColors(i, :);
    b(i).EdgeColor = 'none'; 
end

% Formatting X-Axis
% Placing the tick marks exactly in the center of the two bars for each diet
tick_positions = zeros(1, nDiets);
for d = 1:nDiets
    tick_positions(d) = (x_pos((d-1)*2 + 1) + x_pos((d-1)*2 + 2)) / 2;
end

set(gca, 'XTick', tick_positions, 'XTickLabel', dietLabels, 'FontSize', 11, 'TickDir', 'out');
xlabel('WT | reduced', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('no. of pairwise communities', 'FontSize', 13);
title('Interaction composition across diets (WT left, reduced right)', 'FontSize', 14, 'FontWeight', 'normal');

% Adding Legend (Reversed so Mutualism is at the top)
leg = legend(flip(b), flip(interactionLabels), 'Location', 'northeastoutside');
legend boxoff;

% Cleaning up axes
box off;
set(gca, 'YGrid', true, 'GridColor', [0.9 0.9 0.9]);

% Saving image
export_filename = fullfile(base_dir, 'Plot3_Interaction_Composition_AllDiets.png');
print(fig, export_filename, '-dpng', '-r300');

fprintf('\nPlot 3 generated and saved successfully to:\n%s\n', export_filename);


% HELPER FUNCTIONS

function alpha = calcAlpha(gc, gm)
    if gm < 1e-6
        if gc > 1e-6
            alpha = 1.0; 
        else
            alpha = 0.0; 
        end
    else
        alpha = (gc - gm) / gm;
    end
end

function type = classifyInteraction(a1, a2, th)
    function c = getCat(a, t)
        if a > t, c = 1; elseif a < -t, c = -1; else, c = 0; end
    end
    c1 = getCat(a1, th); c2 = getCat(a2, th);

    if c1 == 1 && c2 == 1, type = 1;
    elseif (c1 == 1 && c2 == 0) || (c1 == 0 && c2 == 1), type = 2;
    elseif c1 == 0 && c2 == 0, type = 3;
    elseif (c1 == -1 && c2 == 0) || (c1 == 0 && c2 == -1), type = 4;
    elseif (c1 == 1 && c2 == -1) || (c1 == -1 && c2 == 1), type = 5;
    elseif c1 == -1 && c2 == -1, type = 6;
    else, type = 3; 
    end
end