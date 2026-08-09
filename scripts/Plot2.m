% PLOT 2: Monoculture vs. Pair Growth, per member, by diet

% Generates a 2x6 grid of scatter plots.
% X-axis: Monoculture Growth
% Y-axis: Growth inside the pairwise community
% Top Row: Wild-Type (WT) models (Blue)
% Bottom Row: Minimal Reduced Models (MRM) (Orange/Red)


% --- Loading Data ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
% Monoculture summary CSV
min_dir = 'C:\Users\Avinash Chauhan\agora_2_models\minReactModels\';
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 

fprintf('Loading datasets for scatter plots...\n');
T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietTitles = {'Diet 1 (HF)', 'Diet 2 (MED)', 'Diet 3 (UNH)', 'Diet 4 (VEG)', 'Diet 5 (WES)', 'Diet 6 (NOD)'};

nDiets = numel(dietNames);
nPairs = height(T_comm);

% Creating a dictionary/map for fast monoculture lookups
if any(strcmp(T_mono.Properties.VariableNames, 'ModelNames'))
    modelList = T_mono.ModelNames;
else
    modelList = T_mono.Model;
end
modelMap = containers.Map(modelList, 1:numel(modelList));

% --- Setup Figure ---
% Creates a wide figure to fit all 6 columns cleanly
fig = figure('Name', 'Plot 2: Monoculture vs Pair Growth', 'Color', 'w', 'Position', [50, 100, 1600, 600]);

fprintf('Extracting growth coordinates and plotting...\n');

% --- Process Data and Plotting ---
for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    
    % Preallocate arrays for X (monoculture) and Y (community) coordinates
    % 1326 pairs = 2652 individual data points per diet
    x_wt = zeros(nPairs * 2, 1);  y_wt = zeros(nPairs * 2, 1);
    x_mrm = zeros(nPairs * 2, 1); y_mrm = zeros(nPairs * 2, 1);
    
    idx_wt = 1; idx_mrm = 1;
    
    for p = 1:nPairs
        % Safely split the community string
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        i1 = modelMap(m1);
        i2 = modelMap(m2);
        
        % --- Extracting WT Data ---
        gm1_wt = T_mono{i1, [prefix '_WT']};
        gm2_wt = T_mono{i2, [prefix '_WT']};
        gc1_wt = T_comm{p, ['WT_M1_Growth_' diet]};
        gc2_wt = T_comm{p, ['WT_M2_Growth_' diet]};
        
        if ~isnan(gm1_wt) && ~isnan(gc1_wt) && gc1_wt >= 0
            x_wt(idx_wt) = gm1_wt; y_wt(idx_wt) = gc1_wt; idx_wt = idx_wt + 1;
        end
        if ~isnan(gm2_wt) && ~isnan(gc2_wt) && gc2_wt >= 0
            x_wt(idx_wt) = gm2_wt; y_wt(idx_wt) = gc2_wt; idx_wt = idx_wt + 1;
        end
        
        % --- Extracting MRM Data ---
        gm1_mrm = T_mono{i1, [prefix '_MRM']};
        gm2_mrm = T_mono{i2, [prefix '_MRM']};
        gc1_mrm = T_comm{p, ['MRM_M1_Growth_' diet]};
        gc2_mrm = T_comm{p, ['MRM_M2_Growth_' diet]};
        
        if ~isnan(gm1_mrm) && ~isnan(gc1_mrm) && gc1_mrm >= 0
            x_mrm(idx_mrm) = gm1_mrm; y_mrm(idx_mrm) = gc1_mrm; idx_mrm = idx_mrm + 1;
        end
        if ~isnan(gm2_mrm) && ~isnan(gc2_mrm) && gc2_mrm >= 0
            x_mrm(idx_mrm) = gm2_mrm; y_mrm(idx_mrm) = gc2_mrm; idx_mrm = idx_mrm + 1;
        end
    end
    
    % Trim preallocated arrays to actual data length
    x_wt = x_wt(1:idx_wt-1);   y_wt = y_wt(1:idx_wt-1);
    x_mrm = x_mrm(1:idx_mrm-1); y_mrm = y_mrm(1:idx_mrm-1);
    
    % ----------------------------------------------------
    % Plot WT (Top Row)
    % ----------------------------------------------------
    subplot(2, nDiets, d);
    % Scatter points (small, light blue, slightly transparent)
    scatter(x_wt, y_wt, 8, [0.35 0.55 0.85], 'filled', 'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Diagonal reference line (y = x)
    max_val_wt = max([x_wt; y_wt; 0.1]);
    plot([0 max_val_wt*1.05], [0 max_val_wt*1.05], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    hold off;
    
    % Formatting
    title(dietTitles{d}, 'FontSize', 11, 'FontWeight', 'normal');
    if d == 1, ylabel('WT pair growth', 'FontSize', 11, 'FontWeight', 'bold'); end
    
    axis equal;
    xlim([0 max_val_wt*1.05]); 
    ylim([0 max_val_wt*1.05]);
    set(gca, 'TickDir', 'out', 'box', 'off');
    
    % ----------------------------------------------------
    % Plot MRM (Bottom Row)
    % ----------------------------------------------------
    subplot(2, nDiets, d + nDiets);
    % Scatter points (small, orange/red, slightly transparent)
    scatter(x_mrm, y_mrm, 8, [0.85 0.40 0.30], 'filled', 'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Diagonal reference line (y = x)
    max_val_mrm = max([x_mrm; y_mrm; 0.1]);
    plot([0 max_val_mrm*1.05], [0 max_val_mrm*1.05], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    hold off;
    
    % Formatting
    if d == 1, ylabel('reduced pair growth', 'FontSize', 11, 'FontWeight', 'bold'); end
    xlabel('mono growth', 'FontSize', 11);
    
    axis equal;
    xlim([0 max_val_mrm*1.05]); 
    ylim([0 max_val_mrm*1.05]);
    set(gca, 'TickDir', 'out', 'box', 'off');
end

% --- Main Title & Output ---
sgtitle('Monoculture vs pair growth, per member — above line = grows faster in the community', 'FontSize', 14, 'FontWeight', 'bold');

% Automatically save the generated figure
export_filename = fullfile(base_dir, 'Plot2_Monoculture_vs_Pair_Growth.png');
print(fig, export_filename, '-dpng', '-r300');

fprintf('\nScatter plot grid generated and saved successfully to:\n%s\n', export_filename);