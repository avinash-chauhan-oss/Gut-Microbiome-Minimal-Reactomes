% Monoculture52.m
%
% Loads the 52 base WT models and 52 MT (minimal) models across all 6 diets.
% Applies the specific diet constraints, runs FBA to calculate monoculture 
% growth, and compiles the results into 'monoculture52.csv'.


fprintf('\n----------------------------------------------------\n');
fprintf('SIMULATING & COMPILING MONOCULTURE 52 CSV\n');
fprintf('----------------------------------------------------\n\n');

%% Directories & Setup
base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
mt_dir   = fullfile(base_dir, 'minReactModels');
diet_dir = fullfile(base_dir, 'Diets');
output_csv = fullfile(base_dir, 'monoculture52.csv');

changeCobraSolver('gurobi', 'all', 0);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
dietShort = {'HF', 'MED', 'UNH', 'VEG', 'WES', 'NOD'};
nDiets = numel(dietNames);

% Getting the list of all 52 original WT model names
model_files = dir(fullfile(base_dir, '*.mat'));
model_files = model_files(~[model_files.isdir]);
nModels = numel(model_files);

%% Loading Diet Definitions
fprintf('Loading diet constraint files...\n');
DietTables = cell(nDiets, 1);
for d = 1:nDiets
    file = fullfile(diet_dir, [dietNames{d} '.txt']);
    if exist(file, 'file')
        T = readtable(file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
        rxns = strtrim(string(T{:,1}));
        rxns = strrep(rxns, '[u]', '(e)'); % Standardize extracellular compartment
        DietTables{d}.rxns = rxns;
        DietTables{d}.lbs  = T{:,2};
    end
end

%% Initializing Results Table
ModelNames = cell(nModels, 1);
for m = 1:nModels
    [~, name] = fileparts(model_files(m).name);
    ModelNames{m} = name;
end

WideTable = table(ModelNames, 'VariableNames', {'Model'});

% Preallocate columns for the 6 diets
for d = 1:nDiets
    WideTable.([dietShort{d} '_WT_Growth']) = nan(nModels, 1);
    WideTable.([dietShort{d} '_MT_Growth']) = nan(nModels, 1);
end

%% Main Simulation Loop
fprintf('Simulating WT and MT models across %d diets...\n', nDiets);

for m = 1:nModels
    modelName = ModelNames{m};
    
    % --- Load WT Model ---
    try
        S_wt = load(fullfile(base_dir, model_files(m).name));
        fn = fieldnames(S_wt);
        wt_model = S_wt.(fn{1});
        if iscell(wt_model), wt_model = wt_model{1}; end
        wt_model.rxns = strtrim(wt_model.rxns);
    catch
        fprintf('Failed to load WT model: %s\n', modelName);
        continue;
    end
    
    % --- Loop Through Diets ---
    for d = 1:nDiets
        dietName = dietNames{d};
        prefix = dietShort{d};
        
        % Evaluating WT Model
        wt_diet_model = applyDietConstraints(wt_model, DietTables{d});
        wt_sol = optimizeCbModel(wt_diet_model, 'max', 'one');
        if ~isempty(wt_sol) && wt_sol.stat == 1 && wt_sol.f > 0
            WideTable{m, [prefix '_WT_Growth']} = wt_sol.f;
        else
            WideTable{m, [prefix '_WT_Growth']} = 0;
        end
        
        % Evaluating MT (Minimal) Model
        mt_file = fullfile(mt_dir, dietName, [modelName '.mat']);
        if exist(mt_file, 'file')
            try
                S_mt = load(mt_file, 'minimalModel');
                if isfield(S_mt, 'minimalModel')
                    mt_model = S_mt.minimalModel;
                    mt_diet_model = applyDietConstraints(mt_model, DietTables{d});
                    mt_sol = optimizeCbModel(mt_diet_model, 'max', 'one');
                    
                    if ~isempty(mt_sol) && mt_sol.stat == 1 && mt_sol.f > 0
                        WideTable{m, [prefix '_MT_Growth']} = mt_sol.f;
                    else
                        WideTable{m, [prefix '_MT_Growth']} = 0;
                    end
                end
            catch
                % Silently handle load errors for MT models
            end
        end
    end
    fprintf('Completed %d/%d: %s\n', m, nModels, modelName);
end

%% Save to CSV
writetable(WideTable, output_csv);
fprintf('\n----------------------------------------------------\n');
fprintf('Success! Monoculture FBA data written to:\n%s\n', output_csv);
fprintf('----------------------------------------------------\n');


% HELPER FUNCTION: Apply Diet Constraints

function model = applyDietConstraints(model, dietTable)
    if isempty(dietTable)
        return;
    end
    
    % Scale lower bounds to per-hour format
    dietLB = dietTable.lbs ./ 24; 
    modelRxnsLower = strrep(lower(string(model.rxns)), '[u]', '(e)');
    
    [tf, loc] = ismember(lower(dietTable.rxns), modelRxnsLower);
    if any(tf)
        valid = loc > 0;
        matchedIdx = loc(valid);
        matchedRxns = model.rxns(matchedIdx);
        matchedLB = dietLB(valid);
        
        % Apply the constraints
        model = changeRxnBounds(model, matchedRxns, matchedLB, 'l');
        model = changeRxnBounds(model, matchedRxns, 1000, 'u');
    end
end