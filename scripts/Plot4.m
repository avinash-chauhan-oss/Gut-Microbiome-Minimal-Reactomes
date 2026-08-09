% PLOT 4: Community-wise interaction heatmap (52x52)
%
% Generates a 52x52 heatmap for each diet.
% Upper triangle: Wild-Type (WT) interaction classifications
% Lower triangle: Minimal Reduced Model (MRM) classifications
% Diagonal: Left blank (white)


% --- Loading Data ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
min_dir  = 'C:\Users\Avinash Chauhan\agora_2_models\minReactModels\';
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 

fprintf('Loading datasets for heatmap generation...\n');
T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietTitles = {'High Fiber Diet', 'Mediterranean Diet', 'Unhealthy Diet', 'Vegetarian Diet', 'Western Diet', 'No Diet'};

nDiets = numel(dietNames);
nPairs = height(T_comm);

% Extracting and sort the 52 models alphabetically to guarantee consistent axes
if any(strcmp(T_mono.Properties.VariableNames, 'ModelNames'))
    modelList = sort(T_mono.ModelNames);
else
    modelList = sort(T_mono.Model);
end
nModels = numel(modelList);

% Creating a dictionary/map for fast 1-52 index lookups
modelMap = containers.Map(modelList, 1:nModels);

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

fprintf('Generating 52x52 Interaction Heatmaps across 6 diets...\n');

% --- Process and Plotting for Each Diet ---
for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    
    % Initializing 52x52 matrix with NaNs (NaNs will render as white)
    H = nan(nModels, nModels);
    
    for p = 1:nPairs
        % Safely split the community string
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        idx1 = modelMap(m1);
        idx2 = modelMap(m2);
        
        % Ensuring consistent upper/lower mapping (i < j for upper triangle)
        i = min(idx1, idx2);
        j = max(idx1, idx2);
        
        % --- A. Wild-Type (WT) Classification -> Upper Triangle ---
        gm1_wt = T_mono{idx1, [prefix '_WT']};
        gm2_wt = T_mono{idx2, [prefix '_WT']};
        gc1_wt = T_comm{p, ['WT_M1_Growth_' diet]};
        gc2_wt = T_comm{p, ['WT_M2_Growth_' diet]};
        
        if ~isnan(gc1_wt) && gc1_wt >= 0 && ~isnan(gc2_wt) && gc2_wt >= 0
            a1_wt = calcAlpha(gc1_wt, gm1_wt);
            a2_wt = calcAlpha(gc2_wt, gm2_wt);
            typeWT = classifyInteraction(a1_wt, a2_wt, threshold);
            H(i, j) = typeWT; % Assign to Upper Triangle
        end
        
        % --- B. Reduced (MRM) Classification -> Lower Triangle ---
        gm1_mrm = T_mono{idx1, [prefix '_MRM']};
        gm2_mrm = T_mono{idx2, [prefix '_MRM']};
        gc1_mrm = T_comm{p, ['MRM_M1_Growth_' diet]};
        gc2_mrm = T_comm{p, ['MRM_M2_Growth_' diet]};
        
        if ~isnan(gc1_mrm) && gc1_mrm >= 0 && ~isnan(gc2_mrm) && gc2_mrm >= 0
            a1_mrm = calcAlpha(gc1_mrm, gm1_mrm);
            a2_mrm = calcAlpha(gc2_mrm, gm2_mrm);
            typeMRM = classifyInteraction(a1_mrm, a2_mrm, threshold);
            H(j, i) = typeMRM; % Assign to Lower Triangle
        end
    end
    
    % --- Plotting the Heatmap ---
    fig = figure('Name', ['Plot 4 Heatmap: ' diet], 'Color', 'w', 'Position', [100, 100, 800, 650]);
    
    % Create the heatmap image
    h_img = imagesc(H);
    
    % Apply the exact custom colormap
    colormap(customColors);
    
    % Set discrete color limits so 1 maps to color 1, 2 to color 2, etc.
    caxis([0.5 6.5]); 
    
    % Make the NaNs (the diagonal) completely transparent/white
    set(h_img, 'AlphaData', ~isnan(H));
    
    % Formatting the axes to match the visual guide
    axis square;
    set(gca, 'XTick', [], 'YTick', []); % Remove tick labels for blocky look
    
    xlabel('partner organism', 'FontSize', 12);
    ylabel('focal organism', 'FontSize', 12);
    
    title(sprintf('Community-wise interaction type\n(upper = WT, lower = reduced)\n%s', dietTitles{d}), ...
        'FontSize', 13, 'FontWeight', 'normal');
    
    % --- Custom Legend ---
    % Create dummy patches for the legend to match the exact colors
    hold on;
    h_legend_patches = zeros(6, 1);
    for leg_idx = 1:6
        h_legend_patches(leg_idx) = patch(NaN, NaN, customColors(leg_idx, :), 'EdgeColor', 'none');
    end
    hold off;
    
    leg = legend(h_legend_patches, interactionLabels, 'Location', 'northeastoutside');
    legend boxoff;
    
    % Save image automatically
    export_filename = fullfile(base_dir, sprintf('Plot4_Heatmap_%s.png', diet));
    print(fig, export_filename, '-dpng', '-r300');
end

fprintf('\nAll 6 heatmaps generated and saved successfully.\n');


% HELPER FUNCTIONS


function alpha = calcAlpha(gc, gm)
    % Robustly calculates the relative growth change.
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
    % Interaction is symmetric, categorization maps directly to 1-6
    function c = getCat(a, t)
        if a > t, c = 1; elseif a < -t, c = -1; else, c = 0; end
    end
    c1 = getCat(a1, th); 
    c2 = getCat(a2, th);

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