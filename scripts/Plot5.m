
%
% This script builds a 6x6 transition matrix for each diet tracking how
% the classification of every pairwise community shifts when minimized.
% It uses a custom patch-based algorithmic renderer to draw smooth 
% Sankey ribbons since standard MATLAB lacks a native Sankey function.

% --- Professional Aesthetic Overrides ---
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultFigureColor', 'w');

% --- Loading Data ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
min_dir  = fullfile(base_dir, 'minReactModels');
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 

T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietTitles = {'High Fiber Diet', 'Mediterranean Diet', 'Unhealthy Diet', 'Vegetarian Diet', 'Western Diet', 'No Diet'};

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
    0.17 0.63 0.17;  
    0.10 0.74 0.61;  
    0.50 0.50 0.50;  
    0.89 0.65 0.00;  
    0.96 0.43 0.26;  
    0.84 0.15 0.16   
];

interactionLabels = {'Mutualism', 'Commensalism', 'Neutralism', 'Amensalism', 'Parasitism', 'Competition'};
threshold = 0.1; 

for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    
    T_matrix = zeros(6, 6);
    
    for p = 1:nPairs
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        idx1 = modelMap(m1);
        idx2 = modelMap(m2);
        
        gm1_wt = T_mono{idx1, [prefix '_WT']};
        gm2_wt = T_mono{idx2, [prefix '_WT']};
        gc1_wt = T_comm{p, ['WT_M1_Growth_' diet]};
        gc2_wt = T_comm{p, ['WT_M2_Growth_' diet]};
        
        gm1_mrm = T_mono{idx1, [prefix '_MRM']};
        gm2_mrm = T_mono{idx2, [prefix '_MRM']};
        gc1_mrm = T_comm{p, ['MRM_M1_Growth_' diet]};
        gc2_mrm = T_comm{p, ['MRM_M2_Growth_' diet]};
        
        if ~isnan(gc1_wt) && gc1_wt >= 0 && ~isnan(gc2_wt) && gc2_wt >= 0 && ...
           ~isnan(gc1_mrm) && gc1_mrm >= 0 && ~isnan(gc2_mrm) && gc2_mrm >= 0
            
            a1_wt = calcAlpha(gc1_wt, gm1_wt);
            a2_wt = calcAlpha(gc2_wt, gm2_wt);
            typeWT = classifyInteraction(a1_wt, a2_wt, threshold);
            
            a1_mrm = calcAlpha(gc1_mrm, gm1_mrm);
            a2_mrm = calcAlpha(gc2_mrm, gm2_mrm);
            typeMRM = classifyInteraction(a1_mrm, a2_mrm, threshold);
            
            T_matrix(typeWT, typeMRM) = T_matrix(typeWT, typeMRM) + 1;
        end
    end
    
    % --- Plotting Sankey Diagram ---
    fig = figure('Name', ['Plot 5 Sankey: ' diet], 'Color', 'w', 'Position', [100, 100, 1000, 600]);
    
    drawSankey(T_matrix, interactionLabels, customColors);
    
    % Cleaned up title grammar
    title(sprintf('Interaction transitions from WT to Reduced\n%s', dietTitles{d}), ...
        'FontSize', 15, 'FontWeight', 'normal');
    
    export_filename = fullfile(base_dir, sprintf('Plot5_%s.png', prefix));
    set(fig, 'InvertHardcopy', 'off');
    print(fig, export_filename, '-dpng', '-r300');
    close(fig);
end

% --- HELPER FUNCTIONS ---

function drawSankey(T, labels, colors)
    totalFlow = sum(T(:));
    if totalFlow == 0
        text(0.5, 0.5, 'No Data Available', 'HorizontalAlignment', 'center', 'FontSize', 14);
        axis off; return;
    end
    
    figHeight = 100;
    gapSize = 4; 
    nNodes = 6;
    totalGapSpace = gapSize * (nNodes - 1);
    scale = (figHeight - totalGapSpace) / totalFlow;
    
    xLeft = 0.2; nodeWidth = 0.05; xLeftEnd = xLeft + nodeWidth;
    xRight = 0.8; xRightStart = xRight - nodeWidth;
    
    leftNodeHeights = sum(T, 2) * scale;
    rightNodeHeights = sum(T, 1) * scale;
    
    leftNodeY = zeros(nNodes, 1);
    rightNodeY = zeros(nNodes, 1);
    
    currY = figHeight;
    for i = 1:nNodes
        leftNodeY(i) = currY;
        currY = currY - leftNodeHeights(i) - gapSize;
    end
    
    currY = figHeight;
    for j = 1:nNodes
        rightNodeY(j) = currY;
        currY = currY - rightNodeHeights(j) - gapSize;
    end
    
    hold on;
    
    curFlowLeft = leftNodeY;
    curFlowRight = rightNodeY;
    
    xNorm = linspace(0, 1, 50);
    xCurve = xLeftEnd + xNorm * (xRightStart - xLeftEnd);
    cosineInterp = (1 - cos(pi * xNorm)) / 2;
    
    for i = 1:nNodes
        for j = 1:nNodes
            if T(i, j) > 0
                flowH = T(i, j) * scale;
                
                yTopL = curFlowLeft(i);
                yBotL = yTopL - flowH;
                yTopR = curFlowRight(j);
                yBotR = yTopR - flowH;
                
                yTopCurve = yTopL + (yTopR - yTopL) * cosineInterp;
                yBotCurve = yBotL + (yBotR - yBotL) * cosineInterp;
                
                patch([xCurve, fliplr(xCurve)], [yTopCurve, fliplr(yBotCurve)], ...
                      colors(i, :), 'FaceAlpha', 0.45, 'EdgeColor', 'none');
                  
                curFlowLeft(i) = yBotL;
                curFlowRight(j) = yBotR;
            end
        end
    end
    
    for i = 1:nNodes
        if leftNodeHeights(i) > 0
            rectangle('Position', [xLeft, leftNodeY(i) - leftNodeHeights(i), nodeWidth, leftNodeHeights(i)], ...
                      'FaceColor', colors(i, :), 'EdgeColor', 'none');
            text(xLeft - 0.02, leftNodeY(i) - leftNodeHeights(i)/2, labels{i}, ...
                 'HorizontalAlignment', 'right', 'FontSize', 12);
        end
        
        if rightNodeHeights(i) > 0
            rectangle('Position', [xRightStart, rightNodeY(i) - rightNodeHeights(i), nodeWidth, rightNodeHeights(i)], ...
                      'FaceColor', colors(i, :), 'EdgeColor', 'none');
            text(xRight + 0.02, rightNodeY(i) - rightNodeHeights(i)/2, labels{i}, ...
                 'HorizontalAlignment', 'left', 'FontSize', 12);
        end
    end
    
    text(xLeft + nodeWidth/2, figHeight + 5, 'WT', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
    text(xRight - nodeWidth/2, figHeight + 5, 'Reduced', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
    
    axis off;
    ylim([-5, figHeight + 10]);
    xlim([0, 1]);
    hold off;
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
    else, type = 3; 
    end
end