function magic_gene_set_reproduction(output_path)
% This function aims to reproduce Figure3: 
% Hsiao and Chiu et al. Differential network analysis reveals the genome-wide landscape of estrogen receptor modulation in hormonal cancers. Scientific Reports. 2016;6:23035. doi:10.1038/srep23035.  

% In order to extend the use of the tool designed for gene level
% interaction analysis to gene set level interaction analysis, the methods
% outlined in the Methods section of the paper and the Supplementary
% document was implemented. 
% In order to obtain data for the gene sets (CGP, TFT, OS, GO, CB),
% a gene set content matrix, which quantifies the presence of a gene within
% the gene set, is combined with the z-transformed gene expression data by
% means of calculating the inner product. The resulting matrix contains the
% degree of activation of a gene set in a sample. This matrix is then used
% to perform MAGIC analysis and obtain the resulting DGRN. 

    %% Load CGP(Chemical and genetic perturbations) gene sets
    % " Gene sets represent expression signatures of genetic and 
    % chemical perturbations. A number of these gene sets come
    % in pairs: xxx_UP (and xxx_DN) gene set representing genes 
    % induced (and repressed) by the perturbation.": Molecular Signatures Database (MSigDB)
    
    json_path = 'Datasets\c2.cgp.v2023.2.Hs.json';
    jsondata = jsondecode(fileread(json_path));
    
    cgp_sets = cell(1, numel(fieldnames(jsondata)));
    geneSets = fieldnames(jsondata);
    
    for i = 1:numel(geneSets)
        geneSetStruct = struct((geneSets{i}), jsondata.(geneSets{i}).geneSymbols);
        cgp_sets{i} = geneSetStruct;
    end
    
    clear i geneSetStruct geneSets json_path jsondata
    
    %% Transcription Factor Targets 
    % "All transcription factor target prediction gene sets. 
    % Combined superset of both GTRD prediction methods and legacy sets.": Molecular Signatures Database (MSigDB) v3.1
    
    json_path = 'Datasets\c3.tft.v2023.2.Hs.json';
    jsondata = jsondecode(fileread(json_path));
    
    tft_sets = cell(1, numel(fieldnames(jsondata)));
    geneSets = fieldnames(jsondata);
    
    for i = 1:numel(geneSets)
        geneSetStruct = struct((geneSets{i}), jsondata.(geneSets{i}).geneSymbols);
        tft_sets{i} = geneSetStruct;
    end
    
    clear i geneSetStruct geneSets json_path jsondata
    
    %% Gene Ontology gene sets
    % "Gene sets that contain genes annotated by the same ontology term. 
    % The C5 collection is divided into two subcollections, the first 
    % derived from the Gene Ontology resource (GO) which contains BP, CC, and MF 
    % components and a second derived from the Human Phenotype Ontology (HPO). 
    % details" Molecular Signatures Database (MSigDB) v3.1
    json_path = 'Datasets\c5.go.v2023.2.Hs.json';
    jsondata = jsondecode(fileread(json_path));
    
    go_sets = cell(1, numel(fieldnames(jsondata)));
    geneSets = fieldnames(jsondata);
    
    for i = 1:numel(geneSets)
        geneSetStruct = struct((geneSets{i}), jsondata.(geneSets{i}).geneSymbols);
        go_sets{i} = geneSetStruct;
    end
    
    clear i geneSetStruct geneSets json_path jsondata
    
    %% (OS) Oncogenic Singature gene sets
    % "Gene sets that represent signatures of cellular pathways 
    % which are often dis-regulated in cancer. The majority of 
    % signatures were generated directly from microarray data 
    % from NCBI GEO or from internal unpublished profiling experiments
    % involving perturbation of known cancer genes": Molecular Signatures Database (MSigDB) v3.1
    
    json_path = 'Datasets\c6.all.v2023.2.Hs.json';
    jsondata = jsondecode(fileread(json_path));
    
    os_sets = cell(1, numel(fieldnames(jsondata)));
    geneSets = fieldnames(jsondata);
    
    for i = 1:numel(geneSets)
        geneSetStruct = struct((geneSets{i}), jsondata.(geneSets{i}).geneSymbols);
        os_sets{i} = geneSetStruct;
    end
    
    clear i geneSetStruct geneSets json_path jsondata
    
    %% Load Cytogenetic Bands (CB)
    % "Gene sets corresponding to human chromosome 
    % cytogenetic bands.": Molecular Signatures Database (MSigDB) v3.1
    
    json_path = 'Datasets\c1.all.v2023.2.Hs.json';
    jsondata = jsondecode(fileread(json_path));
    
    cb_sets = cell(1, numel(fieldnames(jsondata)));
    geneSets = fieldnames(jsondata);
    
    for i = 1:numel(geneSets)
        geneSetStruct = struct((geneSets{i}), jsondata.(geneSets{i}).geneSymbols);
        cb_sets{i} = geneSetStruct;
    end
    
    clear i geneSetStruct geneSets json_path jsondata
    %% Preprocessing [as defined in Methods]
    % small or oversized gene sets (<20 or >500 genes) removed [except CB]
    cgp_sets = cgp_sets(cellfun(@(x) numel(x) >= 20 && numel(x) <= 500, cgp_sets));
    go_sets = go_sets(cellfun(@(x) numel(x) >= 20 && numel(x) <= 500, go_sets));
    os_sets = os_sets(cellfun(@(x) numel(x) >= 20 && numel(x) <= 500, os_sets));
    tft_sets = tft_sets(cellfun(@(x) numel(x) >= 20 && numel(x) <= 500, tft_sets));
    %% Load Gene Expression data for Breast Cancer
    fid = fopen("Datasets\GSE2034-22071.txt", 'r');
    
    % read the first row to get the variable names
    header = fgetl(fid);
    variableNames = strsplit(header, '\t');
    fclose(fid);
    
    data = readtable("Datasets\GSE2034-22071.txt", 'ReadVariableNames', false, 'Delimiter', '\t', 'HeaderLines', 1);
    %% Map Affymetrix probe ID to gene symbol
       
    third_rows = height(data)/3;
    data = data(1:third_rows, : );
    % get gene names after data split
    geneNames = data{:, 1};
    % load affymetrix annotations
    gene_annotation = readtable("Datasets\GPL96-57554.txt",'Delimiter', '\t');
    desired_cols = {'ID','GeneSymbol'};
    filt_gene_annotation = gene_annotation(:,desired_cols);
    % create map for easy lookup
    probeIDs = cellfun(@(x) strrep(x, '''', ''), filt_gene_annotation.ID, 'UniformOutput', false);
    geneSymbols = cellfun(@(x) strrep(x, '''', ''), filt_gene_annotation.GeneSymbol, 'UniformOutput', false);
    id2symbol_map = containers.Map(probeIDs, geneSymbols);
    % replace probe IDs with gene symbols
    for i = 1:numel(geneNames)
       geneNames{i} = id2symbol_map(geneNames{i});
    end
    
    probeIDs_filt = filt_gene_annotation.ID;
    
    % replace probe IDs with gene names 
    for i = 1:numel(probeIDs_filt)
       if isKey(id2symbol_map, probeIDs_filt{i})
           filt_gene_annotation.ID{i} = id2symbol_map(probeIDs_filt{i});
       end
    end
    
    
    data.Properties.VariableNames = variableNames;
    data_matrix = double(data{:, 2:end});
    
    clear fid header variableNames ans
    
    
    %% Read Modulator Data
    % contains information about ER status of each sample
    mod_data = readtable("Datasets\modulator_data.txt");
    
    columns_to_keep = {'GEOAsscessionNumber', 'ERStatus'};
    all_vars = mod_data.Properties.VariableNames;
    vars_to_remove = setdiff(all_vars, columns_to_keep);
    
    mod_data_clean = removevars(mod_data,vars_to_remove);
    clear mod_data vars_to_remove all_vars columns_to_keep
    
    %% Define GROUP
    % binary vector that defines which samples are ER+ and which are ER-
    group = zeros(1, numel(data.Properties.VariableNames)-1);
    modulator_sign = contains(mod_data_clean.ERStatus, '+');
    
    [~, idx] = ismember(mod_data_clean.GEOAsscessionNumber, data.Properties.VariableNames(2:end));
    
    group(idx(modulator_sign)) = 1;
    
    clear modulator_sign idx
    
    %% Sample Size
    count_plus = sum(contains(mod_data_clean{:,2}, "+"));
    count_minus = sum(contains(mod_data_clean{:,2}, "-"));
    
    equ_sam_size = (count_plus + count_minus)/2;
    disp(equ_sam_size);
    
    clear count_minus count_plus
    %% Z-transformation as outlined in Methods
    
    % Mean and Strandard Deviation for each gene
    mu = mean(data_matrix, 2); % col-wise
    s = std(data_matrix, 0, 2); % col-wise
    
    % z-value matrix
    Z = zeros(size(data_matrix));
    
    % z-transformation
    for k = 1:size(data_matrix,1)
        Z(k, group==0) = (data_matrix(k, group==0) - mu(k)) ./ s(k); % ER-
        Z(k, group==1) = (data_matrix(k, group==1) - mu(k)) ./ s(k); % ER+
    end
    
    %% Gene Set Content Matrix
    categories = {cb_sets, cgp_sets, go_sets, os_sets, tft_sets};
    % all gene sets
    gene_sets = cat(2, categories{:});
    gene_sets = cellfun(@struct2cell, gene_sets, 'UniformOutput', false);
    
    
    % gene set content matrix
    S = length(gene_sets); % Number of gene sets
    K = length(geneNames); % Number of genes
    G = zeros(S, K); 
    
    % for each gene set
    for s = 1:S
        n_s = length(gene_sets{s}); % number of genes in gene set 
    
        % for each each gene in gene set
        for k = 1:K
            gene = geneNames{k}; 
            % index of gene in gene set s
            gene_idx = find(strcmp(gene, gene_sets{s}));
            if ~isempty(gene_idx)
                % If gene is in gene set s, set G(s,k) = 1/n_s, else 0
                G(s,k) = 1 / n_s;
            end
        end
    end
    %% Gene Set Enrichment Score
    % not using dot(A,B) as MATLAB treats matrix dot product like the
    % collection of vectors, so identical dimensions for A and B needed. 
    
    % The * operator achieves the desired matrix
    ES_M = G*Z;
    
    % Filter based on L_0.05 criterion
    ES_std = std(ES_M,0,2);
    boundary = 1.96 * ES_std;
    
    l_boundary = ES_M - boundary;
    u_boundary = ES_M + boundary;
    
    within = (ES_M>= l_boundary) & (ES_M <= u_boundary);
    samples = size(Z,2);
    % Filter out if a gene set falls within boundary in more than 80% of
    % samples
    percentage_within = sum(within,2)/samples;
    filtered_gene_sets_indices = find(percentage_within <= 0.8);
    
    % remove the corresponding rows
    ES_M(filtered_gene_sets_indices, :) = [];
    
    % new matrix with filtered gene sets
    nonzero_rows = any(ES_M, 2);
    ES_M_filt = ES_M(nonzero_rows, :);
    
    
    %% Group Vector Approximation
    data = ES_M_filt;
    group = group;
    bonf = 1;
    equ_sam_size = equ_sam_size;
    p_cutoff = 0.05;
    mod_score_cutoff = 0.6;
    output_filename = fullfile(output_path,'ER-MGSIN');
    [p1, p0, mod_score, adj_mat] = MAGIC(ES_M_filt, group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff, output_filename);
    %%
    % for CGP, GO, OS, use Kappa statistic to cluster similar genes and assign
    % to one "functionally representing gene set"
    network = graph(adj_mat); % network(i,j) = [1 OR 2 OR 0 OR -2 OR -1]
    
    % extract genes which are ER-modulated
    gene1_nodes = network.Edges.EndNodes(:,1);
    gene2_nodes = network.Edges.EndNodes(:,2);
    % replace indices with gene names
    
    
    % extract modulation
    % 2: M=1 specific positive correlation; 1: M=1 specific negative correlation
    % -1: M=0 specific positive correlation; -2: M=0 specific negative correlation
    weights = network.Edges.Weight;
    
    
    %%
    % unique gene numbers (since we will create a new adjacency matrix)
    gene_names = unique([gene1_nodes; gene2_nodes]);
    
    % new adjacency matrix
    % init symmetric adjacency matrix for only the modulated pairs
    num_nodes = length(gene_names);
    new_adj_mat = zeros(num_nodes);
    
    
    % weights 
    for k = 1:length(gene1_nodes)
        i = find(gene_names == gene1_nodes(k));
        j = find(gene_names == gene2_nodes(k));
        new_adj_mat(i, j) = weights(k);
        new_adj_mat(j, i) = weights(k); % symmetry
    end
    
    % handle problem with multiple probes mapping onto a single gene
    % leads to self loops in network
    new_adj_mat(logical(eye(size(new_adj_mat)))) = 0;
    % new graph w/o self loops
    network2 = graph(new_adj_mat);
    % removing the self loops leads to isolated nodes
    % remove isolated nodes from the graph manually
    zero_degree_nodes = find(degree(network2) == 0);
    network2 = rmnode(network2, zero_degree_nodes);
    
    
    % plot
    figure;
    plot(network2, 'NodeLabel', {}, 'NodeFontSize', 5, 'LineWidth',0.2, 'Layout','circle', 'NodeColor',[0.5 0.7 1]);
    axis off;
    title("ER-MGSIN");
    %% Output File
    % network.tsv
        % Col 1 -> target
        % Col 2 -> regulator
        % Col 3 -> condition ???
        % Col 4 -> Weight
    
    
    output_file = fullfile(output_path,'network.tsv');
    fid = fopen(output_file, 'w');
    
    % header
    fprintf(fid, 'Target\tRegulator\tCondition\tWeight\n');
    target = gene1_nodes;
    regulator = gene2_nodes;
    w = weights; 
    
    % data
    for k = 1:length(gene1_nodes)
        % condition column
        if w(k)>0
            condition ='Modulated';
        else 
            condition = 'Not Modulated';
        end
        fprintf(fid, '%d\t%d\t%s\t%d\n', target(k), regulator(k), condition, w(k));
    end
    
    
    fclose(fid);
    
    disp(['File ', output_file, ' created successfully!']);

end


