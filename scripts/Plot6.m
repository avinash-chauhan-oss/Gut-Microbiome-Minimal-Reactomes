% PLOT 6: Phylum-pair interaction composition (ALL 6 DIETS)
%
% Creates a grid of pie charts for WT (top) and MRM (bottom).
% Loops through all 6 diets and saves 6 separate high-res images.
% - X/Y Axes: Bacterial Phyla (dynamically generated from data)
% - Pie Size: Proportional to the number of communities in that pairing.


% --- Loading Data ---
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
comm_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');
min_dir  = 'C:\Users\Avinash Chauhan\agora_2_models\minReactModels\';
mono_csv = fullfile(min_dir, 'minReactModels_WideSummary.csv'); 

fprintf('Loading datasets for Phylum-Pair plots...\n');
T_comm = readtable(comm_csv, 'PreserveVariableNames', true);
T_mono = readtable(mono_csv, 'PreserveVariableNames', true);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietPrefixes = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
dietTitles = {'Diet 1 (High Fiber)', 'Diet 2 (Mediterranean)', 'Diet 3 (Unhealthy)', ...
              'Diet 4 (Vegetarian)', 'Diet 5 (Western)', 'Diet 6 (No Diet)'};

nDiets = numel(dietNames);
nPairs = height(T_comm);

if any(strcmp(T_mono.Properties.VariableNames, 'ModelNames'))
    modelList = T_mono.ModelNames;
else
    modelList = T_mono.Model;
end
modelMap = containers.Map(modelList, 1:numel(modelList));

% --- Taxonomy Setup ---
% Getting the exact biological mapping for the 52 models
taxonomyMap = getPhylumMapping(); 

% Dynamically build the Phyla list from the mapped models
phylaList = unique(values(taxonomyMap));
phylaList = sort(phylaList); % Sort alphabetically for consistent axes
nPhyla = length(phylaList);
phylaMap = containers.Map(phylaList, 1:nPhyla);

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

fprintf('Processing all %d diets. This will generate 6 figures...\n', nDiets);

% --- Main Loop for All Diets ---
for d = 1:nDiets
    diet = dietNames{d};
    prefix = dietPrefixes{d};
    dietTitle = dietTitles{d};
    
    % Grids to hold counts for THIS diet
    gridWT  = num2cell(zeros(nPhyla, nPhyla, 6), 3);
    gridMRM = num2cell(zeros(nPhyla, nPhyla, 6), 3);
    
    for p = 1:nPairs
        names = strsplit(T_comm.Community{p}, '_AND_');
        m1 = strrep(names{1}, '.mat', ''); 
        m2 = strrep(names{2}, '.mat', '');
        
        i1 = modelMap(m1);
        i2 = modelMap(m2);
        
        if isKey(taxonomyMap, m1) && isKey(taxonomyMap, m2)
            p1 = phylaMap(taxonomyMap(m1));
            p2 = phylaMap(taxonomyMap(m2));
        else
            continue; 
        end
        
        % Enforcing lower-triangle sorting (Row >= Col)
        row_idx = max(p1, p2);
        col_idx = min(p1, p2);
        
        % --- WT Classification ---
        gm1_wt = T_mono{i1, [prefix '_WT']};
        gm2_wt = T_mono{i2, [prefix '_WT']};
        gc1_wt = T_comm{p, ['WT_M1_Growth_' diet]};
        gc2_wt = T_comm{p, ['WT_M2_Growth_' diet]};
        
        if ~isnan(gc1_wt) && gc1_wt >= 0 && ~isnan(gc2_wt) && gc2_wt >= 0
            a1 = calcAlpha(gc1_wt, gm1_wt);
            a2 = calcAlpha(gc2_wt, gm2_wt);
            typeWT = classifyInteraction(a1, a2, threshold);
            gridWT{row_idx, col_idx}(typeWT) = gridWT{row_idx, col_idx}(typeWT) + 1;
        end
        
        % --- MRM Classification ---
        gm1_mrm = T_mono{i1, [prefix '_MRM']};
        gm2_mrm = T_mono{i2, [prefix '_MRM']};
        gc1_mrm = T_comm{p, ['MRM_M1_Growth_' diet]};
        gc2_mrm = T_comm{p, ['MRM_M2_Growth_' diet]};
        
        if ~isnan(gc1_mrm) && gc1_mrm >= 0 && ~isnan(gc2_mrm) && gc2_mrm >= 0
            a1 = calcAlpha(gc1_mrm, gm1_mrm);
            a2 = calcAlpha(gc2_mrm, gm2_mrm);
            typeMRM = classifyInteraction(a1, a2, threshold);
            gridMRM{row_idx, col_idx}(typeMRM) = gridMRM{row_idx, col_idx}(typeMRM) + 1;
        end
    end
    
    % Find Maximum Count for Global Scaling (Constant across all diets, but calculated safely here)
    maxCount = 0;
    for r = 1:nPhyla
        for c = 1:nPhyla
            maxCount = max([maxCount, sum(gridWT{r,c}), sum(gridMRM{r,c})]);
        end
    end
    
    % --- Plotting ---
    fig = figure('Name', ['Plot 6: ' dietTitle], 'Color', 'w', 'Position', [100, 50, 1000, 1200]);
    
    % Draw WT Grid (Top Half)
    axes('Position', [0.1, 0.55, 0.7, 0.4]); 
    drawPhylumGrid(gridWT, phylaList, maxCount, customColors, sprintf('Phyla-wise distribution of interaction types — %s (WT)', dietTitle));
    
    % Draw MRM Grid (Bottom Half)
    axes('Position', [0.1, 0.05, 0.7, 0.4]); 
    drawPhylumGrid(gridMRM, phylaList, maxCount, customColors, sprintf('Phyla-wise distribution of interaction types — %s (reduced)', dietTitle));
    
    % Adding Legend 
    axes('Position', [0.82, 0.4, 0.15, 0.2]);
    hold on;
    xlim([0, 1]); % Explicit limit
    ylim([0, 7.5]); % Explicit limit
    for i = 1:6
        patch([0 0.2 0.2 0], [i i i+0.5 i+0.5], customColors(7-i,:), 'EdgeColor', 'none');
        text(0.25, i+0.25, interactionLabels{7-i}, 'FontSize', 11, 'VerticalAlignment', 'middle');
    end
    text(0, 6.8, 'Interaction type', 'FontSize', 12, 'FontWeight', 'bold');
    axis off;
    hold off;
    
    export_filename = fullfile(base_dir, sprintf('Plot6_PhylumPair_%s.png', prefix));
    print(fig, export_filename, '-dpng', '-r300');
    
    % Close figure to prevent memory buildup in loop
    close(fig); 
    
    fprintf('Completed and saved: Plot6_PhylumPair_%s.png\n', prefix);
end

fprintf('\nAll 6 Phylum-Pair plots generated successfully!\n');


% HELPER FUNCTIONS


function drawPhylumGrid(gridData, phylaList, maxCount, colors, titleStr)
    n = length(phylaList);
    hold on;
    
    xlim([0.5, n+0.5]); ylim([0.5, n+0.5]);
    set(gca, 'YDir', 'reverse');
    
    set(gca, 'XTick', 1:n, 'XTickLabel', phylaList, 'XTickLabelRotation', 30);
    set(gca, 'YTick', 1:n, 'YTickLabel', phylaList);
    set(gca, 'TickLength', [0 0], 'FontSize', 11, 'box', 'off');
    
    title(titleStr, 'FontSize', 13, 'FontWeight', 'bold');
    
    maxRadius = 0.45; % Max radius to fit inside a 1x1 grid cell
    
    for r = 1:n
        for c = 1:r % Lower triangle only
            counts = gridData{r,c};
            total = sum(counts);
            if total > 0
                radius = maxRadius * sqrt(total / maxCount);
                drawScaledPie(c, r, radius, counts, colors);
            end
        end
    end
    axis equal;
    xlim([0.5, n+0.5]); ylim([0.5, n+0.5]);
end

function drawScaledPie(xCenter, yCenter, radius, counts, colors)
    total = sum(counts);
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
    
    % --- Text color changed to 'k' (black) here ---
    text(xCenter, yCenter, num2str(total), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 9, 'FontWeight', 'bold');
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

function taxMap = getPhylumMapping()
    % Exact biological classification for the 52 AGORA models
    taxMap = containers.Map();
    taxMap('Acidaminococcus_fermentans_DSM_20731') = 'Firmicutes';
    taxMap('Akkermansia_muciniphila_ATCC_BAA_835') = 'Verrucomicrobia';
    taxMap('Alistipes_finegoldii_DSM_17242') = 'Bacteroidetes';
    taxMap('Alistipes_putredinis_DSM_17216') = 'Bacteroidetes';
    taxMap('Alistipes_shahii_WAL_8301') = 'Bacteroidetes';
    taxMap('Bacteroides_caccae_ATCC_43185') = 'Bacteroidetes';
    taxMap('Bacteroides_coprocola_M16_DSM_17136') = 'Bacteroidetes';
    taxMap('Bacteroides_eggerthii_DSM_20697') = 'Bacteroidetes';
    taxMap('Bacteroides_fragilis_NCTC_9343') = 'Bacteroidetes';
    taxMap('Bacteroides_intestinalis_341_DSM_17393') = 'Bacteroidetes';
    taxMap('Bacteroides_stercoris_ATCC_43183') = 'Bacteroidetes';
    taxMap('Bacteroides_thetaiotaomicron_VPI_5482') = 'Bacteroidetes';
    taxMap('Bacteroides_uniformis_ATCC_8492') = 'Bacteroidetes';
    taxMap('Bifidobacterium_adolescentis_ATCC_15703') = 'Actinobacteria';
    taxMap('Bifidobacterium_bifidum_PRL2010') = 'Actinobacteria';
    taxMap('Bifidobacterium_longum_NCC2705') = 'Actinobacteria';
    taxMap('Bilophila_wadsworthia_3_1_6') = 'Proteobacteria';
    taxMap('Blautia_wexlerae_DSM_19850') = 'Firmicutes';
    taxMap('Citrobacter_amalonaticus_Y19') = 'Proteobacteria';
    taxMap('Clostridium_clostridioforme_CM201') = 'Firmicutes';
    taxMap('Collinsella_tanakaei_YIT_12063') = 'Actinobacteria';
    taxMap('Coprococcus_catus_GD_7') = 'Firmicutes';
    taxMap('Desulfovibrio_desulfuricans_subsp_desulfuricans_DSM_642') = 'Proteobacteria';
    taxMap('Desulfovibrio_piger_ATCC_29098') = 'Proteobacteria';
    taxMap('Dialister_invisus_DSM_15470') = 'Firmicutes';
    taxMap('Dorea_longicatena_DSM_13814') = 'Firmicutes';
    taxMap('Enterobacter_cloacae_EcWSU1') = 'Proteobacteria';
    taxMap('Enterococcus_faecium_TX1330') = 'Firmicutes';
    taxMap('Escherichia_coli_str_K_12_substr_MG1655') = 'Proteobacteria';
    taxMap('Eubacterium_rectale_M104_1') = 'Firmicutes';
    taxMap('Faecalibacterium_prausnitzii_L2_6') = 'Firmicutes';
    taxMap('Fusobacterium_varium_ATCC_27725') = 'Fusobacteria';
    taxMap('Haemophilus_parainfluenzae_T3T1') = 'Proteobacteria';
    taxMap('Helicobacter_pylori_26695') = 'Proteobacteria';
    taxMap('Klebsiella_pneumoniae_pneumoniae_MGH78578') = 'Proteobacteria';
    taxMap('Lactobacillus_mucosae_LM1') = 'Firmicutes';
    taxMap('Lactobacillus_reuteri_SD2112_ATCC_55730') = 'Firmicutes';
    taxMap('Megasphaera_elsdenii_DSM_20460') = 'Firmicutes';
    taxMap('Parabacteroides_distasonis_ATCC_8503') = 'Bacteroidetes';
    taxMap('Parabacteroides_johnsonii_DSM_18315') = 'Bacteroidetes';
    taxMap('Phascolarctobacterium_succinatutens_YIT_12067') = 'Firmicutes';
    taxMap('Prevotella_ruminicola_23') = 'Bacteroidetes';
    taxMap('Pseudoflavonifractor_capillosus_strain_ATCC_29799') = 'Firmicutes';
    taxMap('Roseburia_intestinalis_L1_82') = 'Firmicutes';
    taxMap('Roseburia_inulinivorans_DSM_16841') = 'Firmicutes';
    taxMap('Ruminococcus_bromii_L2_63') = 'Firmicutes';
    taxMap('Ruminococcus_callidus_ATCC_2776001') = 'Firmicutes';
    taxMap('Ruminococcus_torques_ATCC_27756') = 'Firmicutes';
    taxMap('Shigella_flexneri_2002017') = 'Proteobacteria';
    taxMap('Streptococcus_salivarius_JIM8777') = 'Firmicutes';
    taxMap('Subdoligranulum_variabile_DSM_15176') = 'Firmicutes';
    taxMap('Veillonella_atypica_ACS_049_V_Sch6') = 'Firmicutes';
end