% PLOT 7: Interaction composition by metabolic overlap (ALL 6 DIETS)
%
% Computes the WT Jaccard similarity (shared reactions / union) for all 
% 1326 pairs ONCE, bins them into low / mid / high overlap, and loops 
% through all 6 diets to plot a 3x3 grid of scaled pie charts.


% --- Directories & Data Loading ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
min_dir  = fullfile(base_dir, 'minReactModels\');
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 

fprintf('Loading datasets for Jaccard analysis...\n');
T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietTitles = {'Diet 1 (High Fiber)', 'Diet 2 (Mediterranean)', 'Diet 3 (Unhealthy)', ...
              'Diet 4 (Vegetarian)', 'Diet 5 (Western)', 'Diet 6 (No Diet)'};

nDiets = numel(dietNames);
nPairs = height(T_comm);

model_files = dir(fullfile(base_dir, '*.mat'));
model_files = model_files(~[model_files.isdir]);
nModels = numel(model_files);

% --- Pre-compute Jaccard Similarities ---
fprintf('Loading %d base .mat models to compute exact reaction Jaccard similarities...\n', nModels);
modelRxns = cell(nModels, 1);
modelNames = cell(nModels, 1);

for m = 1:nModels
    [~, name] = fileparts(model_files(m).name);
    modelNames{m} = name;
    
    S = load(fullfile(base_dir, model_files(m).name));
    fn = fieldnames(S);
    mod = S.(fn{1});
    if iscell(mod), mod = mod{1}; end
    
    modelRxns{m} = strtrim(cellstr(mod.rxns));
end
modelMap = containers.Map(modelNames, 1:nModels);

if any(strcmp(T_mono.Properties.VariableNames, 'ModelNames'))
    monoList = T_mono.ModelNames;
else
    monoList = T_mono.Model;
end
monoMap = containers.Map(monoList, 1:numel(monoList));

jaccardValues = nan(nPairs, 1);
fprintf('Calculating Jaccard overlaps for all %d pairs...\n', nPairs);

for p = 1:nPairs
    names = strsplit(T_comm.Community{p}, '_AND_');
    m1 = strrep(names{1}, '.mat', ''); 
    m2 = strrep(names{2}, '.mat', '');
    
    if isKey(modelMap, m1) && isKey(modelMap, m2)
        idx1_base = modelMap(m1);
        idx2_base = modelMap(m2);

        % Compute Exact Jaccard: |Intersection| / |Union|
        rxns1 = modelRxns{idx1_base};
        rxns2 = modelRxns{idx2_base};
        interCount = numel(intersect(rxns1, rxns2));
        unionCount = numel(union(rxns1, rxns2));
        jaccardValues(p) = interCount / unionCount;
    end
end

% Compute global tertile thresholds for Low / Mid / High bins (ignoring any NaNs)
validJaccards = jaccardValues(~isnan(jaccardValues));
jaccardBinsEdges = prctile(validJaccards, [33.3, 66.7]);

% --- Classification Setup ---
customColors = [
    0.17 0.63 0.17;  % 1: Mutualism (Green)
    0.09 0.75 0.81;  % 2: Commensalism (Teal)
    0.60 0.60 0.60;  % 3: Neutralism (Grey)
    0.95 0.76 0.20;  % 4: Amensalism (Yellow)
    0.90 0.57 0.22;  % 5: Exploitation (Orange)
    0.84 0.15 0.16   % 6: Competition (Red)
];
interactionLabels = {'Mutualism', 'Commensalism', 'Neutralism', 'Amensalism', 'Exploitation', 'Competition'};
threshold = 0.1; 

% Flat text labels to prevent MATLAB from falling back to 1, 2, 3
rowLabels = {'Synergistic Interaction', 'Neutral', 'Antagonistic Interaction'};
colLabels = {'Low', 'Mid', 'High'};

fprintf('Processing interactions across all 6 diets...\n');

% --- Main Loop for All Diets ---
for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    dietTitle = dietTitles{d};
    
    % Initialize 3x3 Grid
    gridData = num2cell(zeros(3, 3, 6), 3);
    
    for p = 1:nPairs
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        if ~isKey(monoMap, m1) || ~isKey(monoMap, m2), continue; end
        
        idx1_mono = monoMap(m1);
        idx2_mono = monoMap(m2);
        
        gm1 = T_mono{idx1_mono, [prefix '_WT']};
        gm2 = T_mono{idx2_mono, [prefix '_WT']};
        gc1 = T_comm{p, ['WT_M1_Growth_' diet]};
        gc2 = T_comm{p, ['WT_M2_Growth_' diet]};
        
        if isnan(gc1) || gc1 < 0 || isnan(gc2) || gc2 < 0, continue; end
        
        a1 = calcAlpha(gc1, gm1);
        a2 = calcAlpha(gc2, gm2);
        itype = classifyInteraction(a1, a2, threshold);
        
        % Map Jaccard to Bin
        jval = jaccardValues(p);
        if isnan(jval), continue; end
        
        if jval <= jaccardBinsEdges(1)
            jBin = 1;
        elseif jval <= jaccardBinsEdges(2)
            jBin = 2;
        else
            jBin = 3;
        end
        
        % Map Interaction Type to Meta-class Row
        if itype == 1 || itype == 2
            metaRow = 1; % Synergistic
        elseif itype == 3
            metaRow = 2; % Neutral
        else
            metaRow = 3; % Antagonistic
        end
        
        gridData{metaRow, jBin}(itype) = gridData{metaRow, jBin}(itype) + 1;
    end
    
    % --- Plotting ---
    fig = figure('Name', ['Plot 7: ' dietTitle], 'Color', 'w', 'Position', [100, 100, 1000, 800]);

    ax = axes('Position', [0.22, 0.15, 0.55, 0.75]);
    hold on;

    xlim([0.5, 3.5]); ylim([0.5, 3.5]);
    set(ax, 'YDir', 'reverse');
    set(ax, 'XTick', 1:3, 'XTickLabel', colLabels, 'FontSize', 12, 'FontWeight', 'bold');
    set(ax, 'YTick', 1:3, 'YTickLabel', rowLabels, 'FontSize', 12, 'FontWeight', 'bold');
    set(ax, 'TickLength', [0 0], 'box', 'off');

    xlabel('Jaccard similarity bin', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Metabolic overlap (Jaccard) vs Interaction type — %s', dietTitle), ...
        'FontSize', 15, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [0.5, 1.05, 0]);

    % Find Max Count for Scaling
    maxCount = 0;
    for r = 1:3
        for c = 1:3
            % Flatten the array to a 1D vector before summing
            counts = gridData{r,c};
            maxCount = max(maxCount, sum(counts(:)));
        end
    end

    maxRadius = 0.35;
    for r = 1:3
        for c = 1:3
            counts = gridData{r,c};
            counts = counts(:); % Flatten array here too
            total = sum(counts);
            
            % Now total is guaranteed to be a single scalar integer
            if total > 0
                radius = maxRadius * sqrt(total / maxCount);
                drawScaledPie(c, r, radius, counts, customColors);
            end
        end
    end
    hold off;

    % --- Custom Legend ---
    ax_leg = axes('Position', [0.82, 0.3, 0.15, 0.4]);
    axis(ax_leg, [0, 1, 0, 7.5]);
    hold on;
    for i = 1:6
        patch([0 0.15 0.15 0], [i i i+0.4 i+0.4], customColors(7-i,:), 'EdgeColor', 'none');
        text(0.20, i+0.2, interactionLabels{7-i}, 'FontSize', 11, 'VerticalAlignment', 'middle');
    end
    text(0, 6.8, 'Interaction type', 'FontSize', 12, 'FontWeight', 'bold');
    axis off;
    hold off;

    export_filename = fullfile(base_dir, sprintf('Plot7_Jaccard_%s.png', prefix));
    print(fig, export_filename, '-dpng', '-r300');
    
    close(fig);
    fprintf('Completed and saved: Plot7_Jaccard_%s.png\n', prefix);
end

fprintf('\nAll 6 Jaccard plots generated successfully!\n');


% HELPER FUNCTIONS


function drawScaledPie(xCenter, yCenter, radius, counts, colors)
    total = sum(counts);
    if total == 0, return; end
    angles = 2 * pi * (counts / total);
    startAngle = pi/2; 
    
    for i = 1:length(counts)
        if counts(i) > 0
            endAngle = startAngle - angles(i);
            theta = linspace(startAngle, endAngle, 30);
            x = [xCenter, xCenter + radius * cos(theta), xCenter];
            y = [yCenter, yCenter - radius * sin(theta), yCenter]; 
            patch(x, y, colors(i,:), 'EdgeColor', 'w', 'LineWidth', 0.5);
            startAngle = endAngle;
        end
    end
    
    text(xCenter, yCenter, num2str(total), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');
end

function alpha = calcAlpha(gc, gm)
    if gm < 1e-6
        if gc > 1e-6, alpha = 1.0; else, alpha = 0.0; end
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
    else, type = 3; end
end