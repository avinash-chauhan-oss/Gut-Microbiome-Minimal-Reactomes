% Pairwise Community Analysis (WT & MT)
%
% UPDATES & FIXES:
% Unit Scaling Fix: Diet bounds are now correctly divided by 24 
% (DietTables{d}.lbs(k) / 24) to match the monoculture simulations, 
% fixing the 24x growth bias in the community scatter plots.


OverallTimer = tic;

fprintf('\n------------------------------------------------------\n');
fprintf('STEP 2 : PAIRWISE COMMUNITY ANALYSIS (WT & MT)\n');
fprintf('------------------------------------------------------\n\n');


% DIRECTORIES

base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
diet_dir = fullfile(base_dir,'Diets');
mt_dir   = fullfile(base_dir, 'minReactModels');
output_csv = fullfile(base_dir,'Community_Results_Wide_52Models.csv');


% SOLVER

changeCobraSolver('gurobi','LP',0);


% LOADING MODEL FILES

model_files = dir(fullfile(base_dir,'*.mat'));
model_files = model_files(~[model_files.isdir]);
nModels = numel(model_files);

fprintf('Total models found : %d\n', nModels);


% CREATING ALL PAIRWISE COMBINATIONS

Pairs = nchoosek(1:nModels, 2);
nPairs = size(Pairs, 1);
fprintf('Total pairwise communities : %d\n\n', nPairs);


% DIETS SETUP

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
nDiets = numel(dietNames);
DietTables = cell(nDiets, 1);

for d = 1:nDiets
    if strcmp(dietNames{d}, 'No_Diet')
        DietTables{d}.rxns = [];
        DietTables{d}.lbs  = [];
        continue;
    end
    
    file = fullfile(diet_dir,[dietNames{d} '.txt']);
    T = readtable(file, 'FileType','text', 'Delimiter','\t', 'ReadVariableNames',false);

    rxns = string(T{:,1});
    lbs  = T{:,2};
    rxns = strtrim(rxns);
    rxns = strrep(rxns,'[u]','(e)');

    DietTables{d}.rxns = rxns;
    DietTables{d}.lbs  = lbs;
end
fprintf('All diet files loaded successfully.\n\n');


% RESULT TABLE HEADERS

headers = {'Community'};
for d = 1:nDiets
    dn = dietNames{d};
    headers{end+1} = sprintf('WT_Rxns_%s', dn);
    headers{end+1} = sprintf('WT_Growth_%s', dn);
    headers{end+1} = sprintf('WT_M1_Growth_%s', dn);
    headers{end+1} = sprintf('WT_M2_Growth_%s', dn);
    headers{end+1} = sprintf('MRM_Rxns_%s', dn);
    headers{end+1} = sprintf('MRM_Growth_%s', dn);
    headers{end+1} = sprintf('MRM_M1_Growth_%s', dn);
    headers{end+1} = sprintf('MRM_M2_Growth_%s', dn);
end


% PARALLEL POOL

if isempty(gcp('nocreate'))
    parpool('local');
end


% MAIN PARFOR LOOP

WorkerResults = cell(nPairs, 1);
fprintf('\nStarting Parallel Processing...\n');

parfor p = 1:nPairs
    
    model1_MT = [];
    model2_MT = [];
    changeCobraSolver('gurobi','LP',0);

    idx1 = Pairs(p,1);
    idx2 = Pairs(p,2);

    file1 = model_files(idx1).name;
    file2 = model_files(idx2).name;

    [~,name1] = fileparts(file1);
    [~,name2] = fileparts(file2);

    communityName = [name1 '_AND_' name2];
    
    LocalRow = cell(1, length(headers));
    LocalRow{1} = communityName;

    %% ---------------------------------------------------------
    % Loading WT Models
    %% ---------------------------------------------------------
    try
        S1 = load(fullfile(base_dir, file1)); fn1 = fieldnames(S1); model1_WT = S1.(fn1{1});
        S2 = load(fullfile(base_dir, file2)); fn2 = fieldnames(S2); model2_WT = S2.(fn2{1});
        if iscell(model1_WT), model1_WT = model1_WT{1}; end
        if iscell(model2_WT), model2_WT = model2_WT{1}; end
        
        model1_WT.rxns = strtrim(model1_WT.rxns);
        model2_WT.rxns = strtrim(model2_WT.rxns);
        
        bio1_WT = model1_WT.rxns(model1_WT.c ~= 0);
        bio2_WT = model2_WT.rxns(model2_WT.c ~= 0);

        if isempty(bio1_WT) || isempty(bio2_WT)
            WorkerResults{p}=LocalRow;
            continue;
        end
        
        removeFields = { ...
        'C','ctrs','d','dsense','osense','osenseStr',...
        'metCharge','metCharges','metFormulas','metNames','metCHEBIID',...
        'metHMDBID','metKEGGID','metPubChemID','metSmile','metSmiles',...
        'metInChIString','metInchiString','metSEEDID','rules','grRules',...
        'rxnGeneMat','genes','geneNames','subSystems','proteinClasses',...
        'comments','rxnConfidenceScores','citations','ecNumbers'};

        for f = 1:numel(removeFields)
            if isfield(model1_WT,removeFields{f})
                model1_WT = rmfield(model1_WT,removeFields{f});
            end
            if isfield(model2_WT,removeFields{f})
                model2_WT = rmfield(model2_WT,removeFields{f});
            end
        end
        
        WT_Community = createMultipleSpeciesModel({model1_WT; model2_WT}, {bio1_WT{1}; bio2_WT{1}}, 'mergeGenesFlag', false);
        
        if ~isfield(WT_Community,'d')
                WT_Community.d = zeros(size(WT_Community.S,1),1);
        end
        if ~isfield(WT_Community,'dsense')
                WT_Community.dsense = repmat('E',size(WT_Community.S,1),1);
        end
        
        WT_Community.c(:) = 0;
        bioIdx1_WT = find(strcmp(WT_Community.rxns,['model1_' bio1_WT{1}]));
        bioIdx2_WT = find(strcmp(WT_Community.rxns,['model2_' bio2_WT{1}]));

        if isempty(bioIdx1_WT)
            bioIdx1_WT = find(endsWith(WT_Community.rxns,['model1_' bio1_WT{1}]));
        end
        if isempty(bioIdx2_WT)
            bioIdx2_WT = find(endsWith(WT_Community.rxns,['model2_' bio2_WT{1}]));
        end

        if isempty(bioIdx1_WT) || isempty(bioIdx2_WT)
            WorkerResults{p}=LocalRow;
            continue;
        end

        bioIdx1_WT = bioIdx1_WT(1);
        bioIdx2_WT = bioIdx2_WT(1);

        WT_Community.c(:)=0;
        WT_Community.c(bioIdx1_WT)=1;
        WT_Community.c(bioIdx2_WT)=1;
        WT_Community.rxns = strtrim(WT_Community.rxns);
       
    catch ME
        WorkerResults{p}=LocalRow;
        continue;
    end

    col = 2; 
    
    %% ---------------------------------------------------------
    % Loop Over Diets
    %% ---------------------------------------------------------
    for d = 1:nDiets
        dietName = dietNames{d};
        
        WT_Rxns = NaN; WT_G = NaN; WT_M1 = NaN; WT_M2 = NaN;
        MT_Rxns = NaN; MT_G = NaN; MT_M1 = NaN; MT_M2 = NaN;
        
        %% WT Optimization
        try
            WT_Current = WT_Community;
            
            if ~strcmp(dietName, 'No_Diet')
                if ~isempty(DietTables{d}.rxns)
                    for k = 1:length(DietTables{d}.rxns)
                        target = char(DietTables{d}.rxns(k));
                        target = strrep(target,'(e)','[u]');
                        idx = strcmp(strtrim(WT_Current.rxns),target);

                        if any(idx)
                            rxnList = WT_Current.rxns(idx);
                            for r = 1:numel(rxnList)
                                rxnID = rxnList(r);
                                % Dividing by 24 to match monoculture scale
                                % and to convert to per-hour basis
                                WT_Current = changeRxnBounds(WT_Current, rxnID, DietTables{d}.lbs(k) / 24, 'l');
                                WT_Current = changeRxnBounds(WT_Current, rxnID, 1000, 'u');
                            end
                        end
                    end
                end
            end
            
            if ~isfield(WT_Current,'d') || isempty(WT_Current.d)
                WT_Current.d = zeros(size(WT_Current.S,1),1);
            end
            if ~isfield(WT_Current,'dsense') || isempty(WT_Current.dsense)
                WT_Current.dsense = repmat('E',size(WT_Current.S,1),1);
            end
            
            solWT = optimizeCbModel(WT_Current);

            if ~isempty(solWT) && solWT.stat == 1
                WT_G = solWT.f;
                WT_Rxns = nnz(abs(solWT.x)>1e-6);
                WT_M1 = solWT.x(bioIdx1_WT);
                WT_M2 = solWT.x(bioIdx2_WT);
            end
        catch ME
        end
        
        %% Loading Pre-minimized MT Models & Optimize
        try
            mtFile1 = fullfile(mt_dir, dietName, [name1 '.mat']);
            mtFile2 = fullfile(mt_dir, dietName, [name2 '.mat']);
            
            if ~exist(mtFile1,'file') || ~exist(mtFile2,'file')
              col = col + 8;
              continue;
            end
            
            S1 = load(mtFile1);
            if isfield(S1,'minimalModel'), model1_MT = S1.minimalModel; else, continue; end
            while iscell(model1_MT), model1_MT = model1_MT{1}; end
            
            S2 = load(mtFile2);
            if isfield(S2,'minimalModel'), model2_MT = S2.minimalModel; else, continue; end
            while iscell(model2_MT), model2_MT = model2_MT{1}; end
            
            model1_MT.rxns = strtrim(model1_MT.rxns);
            model2_MT.rxns = strtrim(model2_MT.rxns);
            
            bio1_MT = find(model1_MT.c~=0);
            bio2_MT = find(model2_MT.c~=0);

            if isempty(bio1_MT) || isempty(bio2_MT), continue; end

            bioRxn1 = model1_MT.rxns{bio1_MT(1)};
            bioRxn2 = model2_MT.rxns{bio2_MT(1)};
            
            coreFields = {'rxns','mets','S','lb','ub','c','b','csense','description'};
            allF = fieldnames(model1_MT);
            for k = 1:numel(allF)
                if ~ismember(allF{k},coreFields)
                    if isfield(model1_MT,allF{k}), model1_MT = rmfield(model1_MT,allF{k}); end
                    if isfield(model2_MT,allF{k}), model2_MT = rmfield(model2_MT,allF{k}); end
                end
            end
            
            MT_Community = createMultipleSpeciesModel({model1_MT; model2_MT}, {bioRxn1; bioRxn2}, 'mergeGenesFlag', false);
            
            if ~isfield(MT_Community,'d')
                MT_Community.d = zeros(size(MT_Community.S,1),1);
            end
            if ~isfield(MT_Community,'dsense')
                MT_Community.dsense = repmat('E',size(MT_Community.S,1),1);
            end
            
            MT_Community.c(:) = 0;
            bioIdx1_MT = find(strcmp(MT_Community.rxns,['model1_' bioRxn1]));
            bioIdx2_MT = find(strcmp(MT_Community.rxns,['model2_' bioRxn2]));

            if isempty(bioIdx1_MT) || isempty(bioIdx2_MT)
                bioIdx1_MT = find(endsWith(MT_Community.rxns,['model1_' bioRxn1]));
                bioIdx2_MT = find(endsWith(MT_Community.rxns,['model2_' bioRxn2]));
                if isempty(bioIdx1_MT), bioIdx1_MT = find(strcmp(MT_Community.rxns,['model1_' bioRxn1])); end
                if isempty(bioIdx2_MT), bioIdx2_MT = find(strcmp(MT_Community.rxns,['model2_' bioRxn2])); end
            end

            if isempty(bioIdx1_MT) || isempty(bioIdx2_MT), continue; end
            
            bioIdx1_MT = bioIdx1_MT(1);
            bioIdx2_MT = bioIdx2_MT(1);
            
            MT_Community.c(:) = 0;
            MT_Community.c(bioIdx1_MT) = 1;
            MT_Community.c(bioIdx2_MT) = 1;
            MT_Community.rxns = strtrim(MT_Community.rxns);
            
            MT_Current = MT_Community;
            
            % Applying Diet to MT
            if ~strcmp(dietName, 'No_Diet')
                if ~isempty(DietTables{d}.rxns)
                    for k = 1:length(DietTables{d}.rxns)
                        target = char(DietTables{d}.rxns(k));
                        target = strrep(target,'(e)','[u]');
                        idx = strcmp(strtrim(MT_Current.rxns),target);
                        
                        if any(idx)
                            rxnList = MT_Current.rxns(idx);
                            for r = 1:numel(rxnList)
                                rxnID = rxnList(r);
                                % Dividing by 24 to match monoculture scale
                                MT_Current = changeRxnBounds(MT_Current, rxnID, DietTables{d}.lbs(k) / 24, 'l');
                                MT_Current = changeRxnBounds(MT_Current, rxnID, 1000, 'u');
                            end
                        end
                    end
                end
            end
            
            solMT = optimizeCbModel(MT_Current);

            if ~isempty(solMT) && solMT.stat == 1
                MT_G = solMT.f;
                MT_Rxns = nnz(abs(solMT.x)>1e-6);
                MT_M1 = solMT.x(bioIdx1_MT);
                MT_M2 = solMT.x(bioIdx2_MT);
            end
            
        catch ME
        end
        
        %% Storing row data
        LocalRow{col}   = WT_Rxns;
        LocalRow{col+1} = WT_G;
        LocalRow{col+2} = WT_M1;
        LocalRow{col+3} = WT_M2;
        LocalRow{col+4} = MT_Rxns;
        LocalRow{col+5} = MT_G;
        LocalRow{col+6} = MT_M1;
        LocalRow{col+7} = MT_M2;
        
        col = col + 8;
    end
    
    WorkerResults{p} = LocalRow;
    fprintf('Processed Pair %d/%d : %s (%.2f%% Complete)\n', p, nPairs, communityName, 100*p/nPairs);
end


% COLLECTING RESULTS & SAVE

fprintf('\n---------------------------------------------------\n');
fprintf('Collecting Results & Saving CSV...\n');
fprintf('---------------------------------------------------\n');

Results = vertcat(WorkerResults{:});
ResultTable = cell2table(Results, 'VariableNames', headers);

writetable(ResultTable, output_csv);
fprintf('\nCSV written successfully to:\n%s\n', output_csv);

Elapsed = toc(OverallTimer);
fprintf('\n---------------------------------------------------\n');
fprintf('Total execution time : %.2f minutes\n', Elapsed/60);
fprintf('---------------------------------------------------\n');