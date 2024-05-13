function reproduction_gse2034(output_path)
    %% Load Data
    % gene expression table with each row representing a gene and each column
    % representing the expression in a sample
    fid = fopen("Datasets\GSE2034-22071.txt", 'r');
    
    % read the first row to get the variable names
    header = fgetl(fid);
    variableNames = strsplit(header, '\t');
    
    fclose(fid);
    
    data = readtable("Datasets\GSE2034-22071.txt", 'ReadVariableNames', false, 'Delimiter', '\t', 'HeaderLines', 1);
    
    %% sample of data [only run this section if memory runs out for full dataset]
    
    third_rows = height(data)/3;
    data = data(1:third_rows, : );
    
    %% Map Affymetrix probe ID to gene symbol
    % NOTE: multiple probes can refer to the same gene. This is something that
    % needs to be handled later
    
    
    % get gene names after data split
    geneNames = data{:, 1};
    % load affymetrix annotations
    gene_annotation = readtable("Datasets\GPL96-57554.txt",'Delimiter', '\t');
    desired_cols = {'ID','GeneOntologyBiologicalProcess','GeneOntologyMolecularFunction','GeneSymbol'};
    filt_gene_annotation = gene_annotation(:,desired_cols);
    % create map for easy lookup
    probeIDs = cellfun(@(x) strrep(x, '''', ''), filt_gene_annotation.ID, 'UniformOutput', false);
    geneSymbols = cellfun(@(x) strrep(x, '''', ''), filt_gene_annotation.GeneSymbol, 'UniformOutput', false);
    id2symbol_map = containers.Map(probeIDs, geneSymbols);
    % replace probe IDs with gene symbols
    for n = 1:numel(geneNames)
       geneNames{n} = id2symbol_map(geneNames{n});
    end
    
    probeIDs_filt = filt_gene_annotation.ID;
    
    % r probe IDs with gene names 
    for j = 1:numel(probeIDs_filt)
       if isKey(id2symbol_map, probeIDs_filt{j})
           filt_gene_annotation.ID{j} = id2symbol_map(probeIDs_filt{j});
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
    
    %% Sample size average
    count_plus = sum(contains(mod_data_clean{:,2}, "+"));
    
    count_minus = sum(contains(mod_data_clean{:,2}, "-"));
    
    equ_sam_size = (count_plus + count_minus)/2;
    disp(equ_sam_size);
    
    clear count_minus count_plus
    
    %% Input for Algorithm
    data = data_matrix;
    group = group;
    bonf = 1;
    equ_sam_size = equ_sam_size;
    p_cutoff = 0.05;
    mod_score_cutoff = 0.6;
    output_filename = fullfile(output_path,'modulation_gse2034');
    
    clear data_matrix mod_data_clean
    %% MAGIC
    [p1, p0, mod_score, adj_mat] = MAGIC(data, group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff, output_filename);
    
    %% Network graph [Not TF Filtered]
    % merge gene pairs
    % nodes: genes
    % edges: modulated interaction
    
    network = graph(adj_mat); % network(i,j) = [1 OR 2 OR 0 OR -2 OR -1]
    
    % extract genes which are ER-modulated
    gene1_nodes = network.Edges.EndNodes(:,1);
    gene2_nodes = network.Edges.EndNodes(:,2);
    % replace indices with gene names
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    
    % extract modulation
    % 2: M=1 specific positive correlation; 1: M=1 specific negative correlation
    % -1: M=0 specific positive correlation; -2: M=0 specific negative correlation
    weights = network.Edges.Weight;
    
    % unique gene numbers (since we will create a new adjacency matrix)
    gene_names = unique([gene1_names; gene2_names]);
    
    % new adjacency matrix
    % init symmetric adjacency matrix for only the modulated pairs
    num_nodes = length(gene_names);
    new_adj_mat = zeros(num_nodes);
    
    
    % weights 
    for k = 1:length(gene1_names)
        l = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_adj_mat(l, j) = weights(k);
        new_adj_mat(j, l) = weights(k); % symmetry
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
    gene_names(zero_degree_nodes)=[];
    
    % properties
    edge_colors  = assignEdgeColors(network2.Edges.Weight);
    % adjust sizes of nodes according to number of neighbors
    node_degree = degree(network2);
    scaling_factor = 1;  % customizable to preference
    node_sizes = scaling_factor * node_degree;
    
    % plot
    figure;
    plot(network2, 'NodeLabel', gene_names, 'MarkerSize',node_sizes, 'NodeFontSize', 5, 'LineWidth',1.3, 'EdgeColor', edge_colors);
    axis off;
    title("FIG2A - TF Unfiltered");
    
    %% Subnetwork for node with highest degree [Not TF Filtered]
    % for the important subnetworks, the hub genes are identified using the
    % 3-sigma rule. this section is demonstrative, and takes a simplified
    % approach, only identifying the most connected node
    
    % identify node with highest degree and its connections
    [max_degree, max_idx] = max(node_degree);
    connected_nodes = neighbors(network2, max_idx);
    
    subnetwork_nodes = unique([max_idx; connected_nodes]);
    % select submatrix
    subnetwork_adj_mat = adjacency(network2);
    subnetwork_adj_mat = subnetwork_adj_mat(subnetwork_nodes,subnetwork_nodes);
    
    % construct subnetwork
    subnetwork = graph(subnetwork_adj_mat);
    
    subweights = subnetwork.Edges.Weight;
    subedge_colors = assignEdgeColors(subweights);
    
    subnode_degree = degree(subnetwork);
    subnode_sizes = scaling_factor * subnode_degree;
    
    % plot
    figure;
    plot(subnetwork, 'NodeLabel', gene_names(subnetwork_nodes), 'MarkerSize', subnode_sizes, 'EdgeColor', subedge_colors, 'LineWidth', 1.3);
    axis off;
    
    title("Subnetwork Unfiltered");
    
    %% Output File - Unfiltered
    % network.tsv
        % Col 1 -> target
        % Col 2 -> regulator
        % Col 3 -> condition ???
        % Col 4 -> Weight
    
    
    output_file = fullfile(output_path,'unfiltered_network.tsv');
    fid = fopen(output_file, 'w');
    
    % header
    fprintf(fid, 'Target\tRegulator\tCondition\tWeight\n');
    target = gene1_names;
    regulator = gene2_names;
    w = weights; 
    
    % data
    for k = 1:length(gene1_names)
        % condition column
        if w(k)>0
            condition ='Modulated';
        else 
            condition = 'Not Modulated';
        end
        fprintf(fid, '%s\t%s\t%s\t%d\n', target{k}, regulator{k}, condition, w(k));
    end
    
    
    fclose(fid);
    
    disp(['File ', output_file, ' created successfully!']);
    
    
    %% Network - [TF Filtered]
    % merge gene pairs
    % nodes: genes
    % edges: modulated interaction
    
    % import list of tf genes
    tf_list = importdata("Datasets\allTFs_hg38.txt");
    
    gene_table=table(gene1_names,gene2_names,weights);
    % check if gene1_names or gene2_names are TF genes
    tf_gene_indices = ismember(gene_table.gene1_names, tf_list) | ismember(gene_table.gene2_names, tf_list);
    
    % prune rows that don't contain TF genes in at least one column
    gene_table_pruned = gene_table(tf_gene_indices, :);
    gene_table_pruned = unique(gene_table_pruned); % remove duplicates
    pruned_gene_names = unique([gene_table_pruned.gene1_names;gene_table_pruned.gene2_names]);
    
    
    % create new graph
    pruned_network = graph(gene_table_pruned.gene1_names, gene_table_pruned.gene2_names, gene_table_pruned.weights);
    
    % properties
    pruned_node_degree = degree(pruned_network);
    scaling_factor = 3;  % customizable to preference
    pruned_node_sizes = scaling_factor * pruned_node_degree;
    pruned_edge_colors = assignEdgeColors(gene_table_pruned.weights);
    
    figure;
    
    plot(pruned_network, 'EdgeColor', pruned_edge_colors, 'NodeLabel', pruned_gene_names, 'MarkerSize',pruned_node_sizes, 'LineWidth',1.3, 'NodeFontSize', 5);
    axis off;
    title("FIG2A - TF Filtered Genes");
    
    %% Finding hub genes with 3-sigma rule
    % the 3 sigma rule was not defined in the paper. Instead, only genes with 20+ connections were considered.
    % Here, the 3 sigma rule is being a good practice
    mean_degree = mean(pruned_node_degree);
    std_degree = std(pruned_node_degree);
    
    hub_gene_idx  = find(pruned_node_degree>mean_degree + 3*std_degree);
    hub_genes=pruned_gene_names(hub_gene_idx);
    
    %% Plot subnetworks for hub genes
    annotation_id = filt_gene_annotation.ID;
    
    
    for j = 1:numel(hub_genes)
        hub_gene_idx = find(strcmp(pruned_gene_names, hub_genes{j}));
        
        % identify neighbors of hub gene
        connected_nodes = neighbors(pruned_network, hub_gene_idx);
        subnetwork_nodes = unique([hub_gene_idx; connected_nodes]);
        
        % create adjacency matrix from netwrok and select revelant submatrix
        subnetwork_adj_mat = adjacency(pruned_network);
        subnetwork_adj_mat = subnetwork_adj_mat(subnetwork_nodes,subnetwork_nodes);
        
      
        subnetwork = graph(subnetwork_adj_mat);
    
        subnetwork_gene_names = pruned_gene_names(subnetwork_nodes);
        
        % properties
        subweights = subnetwork.Edges.Weight;
        subedge_colors = assignEdgeColors(subweights);
        subnode_degree = degree(subnetwork);
        subnode_sizes = scaling_factor * subnode_degree;
        node_labels = cell(numel(subnetwork_gene_names), 1);
        
        %[~, gene_idx_table] = ismember(subnetwork_gene_names, annotation_id);
        %biological_process_subnetwork = biological_process(gene_idx_table);
        %for j = 1:numel(subnetwork_gene_names)
        %    node_labels{j} = [subnetwork_gene_names{j}, ': ', biological_process_subnetwork{j}];
        %end
    
        % plot subnetwork
        figure;
        plot(subnetwork, 'NodeLabel', subnetwork_gene_names, 'MarkerSize', subnode_sizes, 'EdgeColor', subedge_colors, 'LineWidth', 1.3);
        axis off;
        title(['Subnetwork for Hub Gene: ', hub_genes{j}]);
    end
    
    
    %% Output TF Filtered
    
    % network.tsv
        % Col 1 -> target
        % Col 2 -> regulator
        % Col 3 -> condition ???
        % Col 4 -> Weight
    
    
    output_file = fullfile(output_path,'network.tsv');
    fid = fopen(output_file, 'w');
    
    % header
    fprintf(fid, 'Target\tRegulator\tCondition\tWeight\n');
    target = pruned_network.Edges.EndNodes(:,1);
    regulator = pruned_network.Edges.EndNodes(:,2);
    w = weights; 
    
    % data
    for k = 1:length(pruned_network.Edges.EndNodes)
        % condition column
        if w(k)>0
            condition ='Modulated';
        else 
            condition = 'Not Modulated';
        end
        fprintf(fid, '%s\t%s\t%s\t%d\n', target{k}, regulator{k}, condition, w(k));
    end
    
    
    fclose(fid);
    
    disp(['File ', output_file, ' created successfully!']);
    
    %% Scatter Plots of Gene Pair with highest Modulation Score
    % we find the gene pair with the highest mod_score value (delta_I_adj) 
    % we costruct scatter plots for the gene expression under each condition [modulator expressed or not]
    % find raw correlation under each condition
    
    % split original data into two separate tables, each containing samples
    % under only one of the conditions
    sam1 = find(group==1);
    sam0 = find(group==0);
    
    input1 = data(:,sam1);
    input1 = table(input1);
    input1 = [geneNames, input1]; % concat. geneNames for easy matching
    input2 = data(:,sam0);
    input2 = table(input2);
    input2 = [geneNames, input2];
    
    
    % to extract raw correlation under each condition separately
    % read output of algorithm, which contains:
    % modulation score for each identified gene pair
    % raw correlation under each condition for each identified gene pair
    out_network = readtable([output_filename, '.txt']);
    
    % find gene pair with highest modulation score (6th column in every generated file)
    % indexing columns by name is problematic since the name endings are
    % adjusted for every database (sample size)
    max_mod_score = max(out_network{:,6});
    max_mod_idx = find(out_network{:,6} == max_mod_score);
    %%
    mod_pair = [out_network.Gene1(max_mod_idx);out_network.Gene2(max_mod_idx)];
    max_mod_pair = geneNames(mod_pair);
    
    % raw corr. under each condition for the pair
    raw_corrm1 = out_network{:,3}(max_mod_idx);
    raw_corrm0 = out_network{:,4}(max_mod_idx);
    
    
    %% Scatter under M=0
    
    gene1_expression_m0_idx = strcmp(input1{:,1}, max_mod_pair{1});
    gene1_expression_m0 = input1(gene1_expression_m0_idx,2:end);
    gene1_expression_m0 = gene1_expression_m0{:,:};
    
    gene2_expression_m0_idx = strcmp(input1{:,1}, max_mod_pair{2});
    gene2_expression_m0 = input1(gene2_expression_m0_idx,2:end);
    gene2_expression_m0=gene2_expression_m0{:,:};
    
    figure;
    % massive variability in data here, outlier?
    scatter(gene1_expression_m0,gene2_expression_m0, 'Marker','.');
    
    xlabel(max_mod_pair{1});
    ylabel(max_mod_pair{2});
    title("Scatter Plot of Raw correlation under M=0");
    txt = {'Raw correlation = ', raw_corrm0, 'Modulation score = ', max_mod_score};
    annotation('textbox',[.9 .5 .1 .2], ...
        'String',txt,'EdgeColor','none')
    %% Scatter under M=1
    gene1_expression_m1_idx = strcmp(input2{:,1}, max_mod_pair{1});
    gene1_expression_m1 = input2(gene1_expression_m1_idx,2:end);
    gene1_expression_m1 = gene1_expression_m1{:,:};
    
    gene2_expression_m1_idx = strcmp(input2{:,1}, max_mod_pair{2});
    gene2_expression_m1 = input2(gene2_expression_m1_idx,2:end);
    gene2_expression_m1=gene2_expression_m1{:,:};
    
    figure;
    scatter(gene1_expression_m1,gene2_expression_m1, 'Marker','.');
    
    xlabel(max_mod_pair{1});
    ylabel(max_mod_pair{2});
    title("Scatter Plot of Raw correlation under M=1");
    txt = {'Raw Correlation = ', raw_corrm1, 'Modulation score = ', max_mod_score};
    annotation('textbox',[.9 .5 .1 .2], ...
        'String',txt,'EdgeColor','none')
    
    %% [Fig2B] Scatter plots of AKR1C1-LPL Gene pair
    % This was done separately, as pair AKR1C1-LPL did not result as the pair
    % with the highest modulation score (expected: 0.8..)
    % find indices of the pair, as genes are referred to by their indices in
    % the output file
    
    idx_geneA = find(strcmp(geneNames, 'AKR1C1'));
    idx_geneB = find(strcmp(geneNames, 'LPL'));
    % find rows in the output file where both genes appear
    pair_idx = find(out_network.Gene1 == idx_geneA & out_network.Gene2 == idx_geneB(1) | ...
        out_network.Gene1 == idx_geneA & out_network.Gene2 == idx_geneB(2) | ...
        out_network.Gene1 == idx_geneB(1) & out_network.Gene2 == idx_geneA | ...
        out_network.Gene1 == idx_geneB(2) & out_network.Gene2 == idx_geneA);
    
    % obtain maximum modulation score (there can be multiple due to the probe-gene symbol map)
    pair_mod_score = out_network{:,6}(pair_idx);
    [max_pair_mod_score, max_pair_mod_idx] = max(pair_mod_score); % find which index gives maximum
    % use index of pair that gives maximum to also obtain raw correlation under
    % each condition
    pair_raw_corrm0 = out_network{:,3}(pair_idx(max_pair_mod_idx));
    pair_raw_corrm1 = out_network{:,4}(pair_idx(max_pair_mod_idx));
    
    % scatter plot for raw correlation under no modulation
    AKR1C1_expression_m0 = input1(idx_geneA,2:end);
    AKR1C1_expression_m0 = AKR1C1_expression_m0{:,:};
    
    LPL_expression_m0 = input1(idx_geneB(max_pair_mod_idx),2:end);
    LPL_expression_m0=LPL_expression_m0{:,:};
    
    figure;
    scatter(AKR1C1_expression_m0,LPL_expression_m0,'Marker', '.');
    
    xlabel('AKR1C1');
    ylabel('LPL');
    title("Scatter Plot of Raw correlation between AKR1C1-LPL under M=0");
    txt = {'Raw Correlation = ', pair_raw_corrm0, 'Modulation score = ', max_pair_mod_score};
    annotation('textbox',[.9 .5 .1 .2], ...
        'String',txt,'EdgeColor','none')
    
    
    % scatter plot for raw correlation under modulation
    AKR1C1_expression_m1 = input2(idx_geneA,2:end);
    AKR1C1_expression_m1 = AKR1C1_expression_m1{:,:};
    
    LPL_expression_m1 = input2(idx_geneB(max_pair_mod_idx),2:end);
    LPL_expression_m1=LPL_expression_m1{:,:};
    
    figure;
    scatter(AKR1C1_expression_m1,LPL_expression_m1,'Marker', '.');
    
    xlabel('AKR1C1');
    ylabel('LPL');
    title("Scatter Plot of Raw correlation between AKR1C1-LPL under M=1");
    txt = {'Raw Correlation = ', pair_raw_corrm1, 'Modulation score = ', max_pair_mod_score};
    annotation('textbox',[.9 .5 .1 .2], ...
        'String',txt,'EdgeColor','none')
    
    %Conclusion:
    % We see different values than the study, starting from the modulation
    % score
    
    %% [Fig2C] Subnetwork and functional annotations of NRN1 + modulated pairs
    idx_NRN1 = find(strcmp(geneNames, 'NRN1'));
    % could not find
    
    %% [Fig2D] Subnetwork and functional annotations of SFPR1 + modulated pairs
    idx_SFPDR1 = find(strcmp(geneNames, 'SFPDR1'));
    % could not find
    %% Functions
    
    function edge_colors = assignEdgeColors(weights)
        edge_colors = zeros(length(weights), 3);
        for i = 1:length(weights)
            if weights(i) == 2 || weights(i) == -1
                edge_colors(i, :) = [0, 0, 1];  % Blue for M=1 specific positive correlation or M=0 specific positive correlation
            elseif weights(i) == 1 || weights(i) == -2
                edge_colors(i, :) = [1, 0, 0];  % Red for M=1 specific negative correlation or M=0 specific negative correlation
            end
        end
    end
    
    %% MAGIC funct
    function [p1,p0,mod_score,adj_mat] = ...
    MAGIC(data,group,bonf,equ_sam_size,p_cutoff,mod_score_cutoff,output_filename)
    
    % MATLAB tool for modulated gene/gene set interaction (MAGIC) analysis
    % 
    % 
    % MAGIC(DATA,GROUP,BONF,EQU_SAM_SIZE,P_CUTOFF,MOD_SCORE_CUTOFF,OUTPUT_FILENAME)
    % identifies differentially correlated gene (or gene set) pairs modulated
    % by states of a modulator; i.e., pair of genes that is correlated
    % specifically in one state of the modulator (M). All possible combinations
    % of genes deposited in DATA are tested. Take pair of gene i and j for
    % example, correlation coefficients of gene i and j are separately
    % calculated in samples with M=1 and samples with M=0. The correlation
    % coefficients are Fisher transformed to a sample-size-free domain and
    % tested for significance of their difference in the absolute manner (the
    % modulation test). To ensure biologically meaningful change between the
    % correlation coefficients, inverse Fisher transformation is utilized to
    % convert the Fisher transformed coefficients back to the domain with a
    % user-defined equivalent sample size. The modulation score measures the
    % difference of transformed correlation coefficients. Gene (or gene set)
    % pairs that meet the criteria on p-value from modulation test and
    % modulation score are defined as modulated interaction pairs. The MAGIC
    % tool outputs three matrices: modulation p-values, modulation scores, and
    % adjacency matrix of the modulated interaction network, as well as a .txt
    % file that can be used to generate the modulated interaction network by
    % the Cytoscape software.
    % 
    % 
    % Description of the input parameters:
    % 
    % DATA is a K-by-N numeric matrix (in double precision), which contains the
    % expression profiles of K genes (or enrichment scores of K gene sets) in N
    % samples. DATA should not contain NaNs.
    % 
    % GROUP is an N-length numeric vector that defines binary states of the
    % modulator in N samples. GROUP can contain only 0s and 1s.
    % 
    % BONF is set as 1 to perform Bonferroni correction to the number of
    % testing (nchoosek(K,2)). To analyze raw p-values, BONF should be set to
    % 0. Suggested setting: 1
    % 
    % EQU_SAM_SIZE is a numeric value denoting user-assigned sample size at
    % which correlation coefficients from two sample sizes (i.e., number of
    % samples with M=1 and M=0) are compared; that is, the modulation scores
    % are calculated at the sample size of EQU_SAM_SIZE. Suggested setting:
    % average of number of samples with M=1 and M=0
    % 
    % P_CUTOFF is the threshold on raw (or Bonferroni corrected, when BONF = 1)
    % p-value to define "statistically" significant modulated interaction.
    % Suggested setting: 0.05
    % 
    % MOD_SCORE_CUTOFF is the threshold on modulation score to define
    % "biologically" significant modulated interaction. MOD_SCORE_CUTOFF must
    % be a positive numeric value. Suggested setting: 0.6
    % 
    % OUTPUT_FILENAME is a string specifying a filename for the output .txt
    % file. If set as 'NA', no output txt file will be generated.
    % 
    % 
    % Description of the outputs:
    % 
    % P1 is a K-by-K symmetric matrix, with elements of p-value from the
    % modulation test. Significant P1(i,j) (typically < 0.05) means that genes
    % i and j are strongly (either positively or negatively) correlated
    % specifically in M=1 samples.
    % 
    % P0 is a K-by-K symmetric matrix, denoting the significance of strong
    % correlation specifically in M=0 samples.
    % 
    % MOD_SCORE is a K-by-K symmetric matrix of modulation scores. Larger
    % positive elements have stronger correlation in M=1 samples compared to
    % M=0.
    % 
    % ADJ_MAT is a K-by-K symmetric adjacency matrix, of which a non-zero
    % element ADJ_MAT(i,j) denotes a modulated interaction pair of i and j
    % (i.e., the i-j edge in the modulated interaction network.
    % 
    % When output_filename is specified with any string except for 'NA', a
    % Cytoscape compatible 'output_filename.txt' will be generated. The .txt
    % file can be imported to Cytoscape for construction, visualization, and
    % analyses of the modulated interaction network.
    % 
    % 
    % Reference: The MAGIC tool is for academic purposes only and all rights
    % are reserved. To reference the MAGIC algorithm or the tool, please cite
    % the paper: Hsiao and Chiu et al. Differential network analysis reveals
    % the genome-wide landscape of estrogen receptor modulation in hormonal
    % cancers. Scientific Reports. 2016;6:23035. doi:10.1038/srep23035.
    % Mathematical details and biological applications of MAGIC can be found
    % in this paper.
    % 
    % Enjoy!
    
    tic
    
    num_gene = size(data,1);
    sam1 = find(group==1);
    sam0 = find(group==0);
    num_sam1 = length(sam1);
    num_sam0 = length(sam0);
    
    disp(sam1);
    
    
    % calculation of raw correlation coefficients
    [Corr1 Corr1_p] = corrcoef(data(:,sam1)');
    z1 = 0.5*log((1+Corr1)./(1-Corr1));
    CS1 = sqrt(num_sam1-3)*z1; % Fisher-transformed correlation
    
    [Corr0 Corr0_p] = corrcoef(data(:,sam0)');
    z0 = 0.5*log((1+Corr0)./(1-Corr0));
    CS0 = sqrt(num_sam0-3)*z0; % Fisher-transformed correlation
    
    corr_diff_fisher = abs(CS1) - abs(CS0);
    
    % p-value from the modulation test p-value
    p1 = 1-(0.5+erf(corr_diff_fisher/2)-0.5*sign(corr_diff_fisher).*erf(corr_diff_fisher/2).*erf(corr_diff_fisher/2)); % right-tail
    p0 = (0.5+erf(corr_diff_fisher/2)-0.5*sign(corr_diff_fisher).*erf(corr_diff_fisher/2).*erf(corr_diff_fisher/2)); % left-tail
    
    % Bonferroni correction
    if bonf==1
        p1 = p1*nchoosek(num_gene,2);
        p0 = p0*nchoosek(num_gene,2);
    end
    
    p1(p1>1) = 1;
    p0(p0>1) = 1;
    
    % inverse Fisher transformation to N = equ_sam_size %
    z1_b = 1/sqrt(equ_sam_size-3)*CS1;
    r1_b = (exp(2*z1_b)-1)./(exp(2*z1_b)+1);
    z0_b = 1/sqrt(equ_sam_size-3)*CS0;
    r0_b = (exp(2*z0_b)-1)./(exp(2*z0_b)+1);
    mod_score = abs(r1_b)-abs(r0_b);
    
    % identification of M=1 specific interaction pairs
    [row1 col1] = find((p1<=p_cutoff).*(mod_score>=mod_score_cutoff).*(triu(ones(num_gene),1)));
    id1 = find((p1<=p_cutoff).*(mod_score>=mod_score_cutoff).*(triu(ones(num_gene),1)));
    
    % identification of M=0 specific interaction pairs
    [row0 col0] = find((p0<=p_cutoff).*(mod_score<=-mod_score_cutoff).*(triu(ones(num_gene),1)));
    id0 = find((p0<=p_cutoff).*(mod_score<=-mod_score_cutoff).*(triu(ones(num_gene),1)));
    
    % adjacency matrix
    % 2: M=1 specific positive correlation; 1: M=1 specific negative correlation
    % -1: M=0 specific positive correlation; -2: M=0 specific negative correlation
    adj_mat = zeros(num_gene);
    adj_mat(id1(Corr1(id1)>0)) = 2;
    adj_mat(id1(Corr1(id1)<0)) = 1;
    adj_mat(id0(Corr0(id0)>0)) = -1;
    adj_mat(id0(Corr0(id0)<0)) = -2;
    for i=1:(num_gene-1)
        adj_mat((i+1):num_gene,i) = adj_mat(i,(i+1):num_gene);
    end
    
    time_used = toc;
    
    % export Cytoscape .txt file
    if ~strcmp(output_filename,'NA')
        fid = fopen(sprintf('%s.txt',output_filename),'w');
        fprintf(fid, ['Gene1' '\t' 'Gene2' '\t' sprintf('Raw corr in M=1 (N=%s)',num2str(num_sam1)) ...
            '\t' sprintf('Raw corr in M=0 (N=%s)',num2str(num_sam0))  '\t' 'P-value' '\t' ...
            sprintf('Modulation score (N=%s)',num2str(equ_sam_size)) '\n']);
        fprintf(fid, '%d \t %d \t %.3f \t %.3f \t %.2e \t %.3f \n', [col1 row1 Corr1(id1) Corr0(id1) p1(id1) mod_score(id1)]');
        fprintf(fid, '%d \t %d \t %.3f \t %.3f \t %.2e \t %.3f \n', [col0 row0 Corr1(id0) Corr0(id0) p0(id0) mod_score(id0)]');
        fclose(fid);
        disp(sprintf('\n\nSuccess! MAGIC analysis comes true!\n\n%s.txt has been generated.\n\nComputation time: %.2f seconds.\n',output_filename,time_used));
    else
        disp(sprintf('\n\nSuccess! MAGIC analysis comes true!\n\nComputation time: %.2f seconds.\n',time_used));
    end
    
    end
end