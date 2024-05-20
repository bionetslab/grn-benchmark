function downstream_analysis_rand(output_path)   
%% Analysis of Algorithm
    % Comments explaining the original procedure have been left out for
    % legibility of method used for downstream analysis
    % Original construction with comments can be found:
        % reproduction_reference_data.m
    %% Data
    input1 = readtable("Datasets\out_CD8_exhausted.tsv", 'Delimiter', '\t', 'FileType','text');
    input2 = readtable("Datasets\out_Macrophages.tsv",'Delimiter', '\t', 'FileType','text');
    input2_Excl = input2(:, 2:end);
    geneNames = input1{:,1};
    input2_Excl.Properties.VariableNames = strrep(input2_Excl.Properties.VariableNames, "Hackathon_April2024", '');
    total_data = [input1, input2_Excl];
    data_matrix = double(total_data{:, 2:end});
    
    %% Modulation Vector
    num_vars_input1 = size(input1,2)-1;
    num_vars_input2 = size(input2_Excl, 2);
    
    % original group vector
    group = [zeros(1,num_vars_input1), ones(1,num_vars_input2)];
    
    % randomized group vector
    % using random permutations of original vector to ensure same number of
    % occurrences of 1 and 0 in both vectors
    % a change in the number of occurrences conflicts with the equ_sam_size
    % parameter which would skew the algorithm along an entirely different
    % dimension. Here, an attempt at isolating the influence of the group
    % vector is crucial.
    randomized_group = group(randperm(length(group)));
    % Hamming distance as similarity measure, to ensure that the orignial
    % vector and the randomized vector are not too similar
        % Calculates number of positions at which corresponding elements btw
        % the two vectors are dissimilar
        % Low value indicates greater similarity;
    hamming_dist = sum(group ~= randomized_group);
    disp(['Hamming distance (group and randomized group): ', num2str(hamming_dist)]);
    
    %% Sample size average
    equ_sam_size = (num_vars_input1 + num_vars_input1)/2;
    
    %% Algorithm Input
    data = data_matrix;
    group = group;
    bonf = 1;
    equ_sam_size = equ_sam_size;
    p_cutoff = 0.05;
    mod_score_cutoff = 0.6;
    output_filename_original = fullfile(output_path, 'modulation');
    %% MAGIC algorithm
    [p1, p0, mod_score, adj_mat] = MAGIC(data, group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff, output_filename_original);
    
    %% Differential Gene Regulatory Network with Original Group vector 
    network = graph(adj_mat); 
    
    gene1_nodes = network.Edges.EndNodes(:,1);
    gene2_nodes = network.Edges.EndNodes(:,2);
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    weights = network.Edges.Weight;
    gene_names = unique([gene1_names; gene2_names]);
    
    num_nodes = length(gene_names);
    new_adj_mat = zeros(num_nodes);
    
    for k = 1:length(gene1_nodes)
        i = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_adj_mat(i, j) = weights(k);
        new_adj_mat(j, i) = weights(k); 
    end
    
    edge_colors = assignEdgeColors(weights);
    new_network = graph(new_adj_mat);
    
    node_degree = degree(new_network);
    scaling_factor = 2; 
    node_sizes = scaling_factor * node_degree;
    
    figure;
    plot(new_network, 'NodeLabel', gene_names, 'MarkerSize',node_sizes, 'EdgeColor',edge_colors, 'LineWidth',1.3, 'NodeFontSize', 5);
    axis off;
    title("Differential Gene Regulatory Network using normal grouping");
    %% Modules for normal group
    bins = conncomp(new_network);
    bin_nodes = cell(max(bins), 1);
    
    % group nodes according to bins
    % each cell in bin_nodes contains indices of all nodes that share a bin
    for l = 1:max(bins)
        bin_nodes{l} = find(bins == l);
    end
    
    % do not consider pairs
    % pairs are not strictly uninformative, but this decision was made
    % due to the fact given that denser communities were detected
    valid_bins = find(cellfun(@numel, bin_nodes) > 2);
    bin_nodes = bin_nodes(valid_bins);
    
    num_clusters = length(bin_nodes);
    num_cols = ceil(sqrt(num_clusters));
    num_rows = ceil(num_clusters/num_cols);
    cluster_gene_names_all = cell(num_clusters,1);
    
    figure;
    % build subnetworks for all the disconnected subnetworks
    for k = 1:num_clusters
        % subnetwork for current bin
        subnetwork_nodes = bin_nodes{k};
        cluster_gene_names = gene_names(subnetwork_nodes);
        cluster_gene_names_all{k} = cluster_gene_names;
        %cluster_gene_bio = bio_processes{subnetwork_nodes};
    
        subnetwork_adj_mat = new_adj_mat(subnetwork_nodes, subnetwork_nodes);
    
        subnetwork = graph(subnetwork_adj_mat);
    
        subplot(num_rows, num_cols, k);
        plot(subnetwork, 'Layout', 'force', 'NodeLabel',cluster_gene_names,'NodeFontSize',6);
        axis off;
        title(['Module ', num2str(k)]);
    
    end
    sgtitle('Modules of Differential Regulatory Gene Network with normal grouping');
    
    
    
    %% MAGIC with randomized group vector
    randomized_group = randomized_group;
    output_filename_original = fullfile(output_path, 'modulation_rand');;
    
    [p1r, p0r, mod_score_r, adj_mat_r] = MAGIC(data, randomized_group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff, output_filename_rand);
    
    %% Differential Gene Regulatory Network with Randomized Group vector 
    network_rand = graph(adj_mat_r); 
    
    gene1_nodes = network_rand.Edges.EndNodes(:,1);
    gene2_nodes = network_rand.Edges.EndNodes(:,2);
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    weights = network_rand.Edges.Weight;
    gene_names = unique([gene1_names; gene2_names]);
    
    num_nodes = length(gene_names);
    new_adj_mat_r = zeros(num_nodes);
    
    for k = 1:length(gene1_nodes)
        m = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_adj_mat_r(m, j) = weights(k);
        new_adj_mat_r(j, m) = weights(k); 
    end
    
    edge_colors = assignEdgeColors(weights);
    new_network_rand = graph(new_adj_mat_r);
    
    node_degree = degree(new_network_rand);
    scaling_factor = 2; 
    node_sizes = scaling_factor * node_degree;
    
    figure;
    plot(new_network_rand, 'NodeLabel', gene_names, 'MarkerSize',node_sizes, 'EdgeColor',edge_colors, 'LineWidth',1.3, 'NodeFontSize', 5);
    axis off;
    %% DISCUSSION: Randomized Group Vector
    % Network is not generated when group vector is randomized (Sample not
    % corresponding with sample). Contains no modulated pairs. Checking the new
    % generated modulation score matrix:
    
    disp(max(mod_score_r, [], 'all')); % 0.5390 or some variation that is under 0.6
    % for reference, the maximum value in the original modulation score matrix:
    disp(max(mod_score, [], 'all')); % 0.7284 or some variation that is above 0.6
    % Our mod_score_cutoff is the recommended value of 0.6. That is to say,
    % that when running the MAGIC algorithm with the randomized group vector,
    % no biologically significant modulated pairs are found, and therefore a
    % network can not be constructed. 
    
    %%
    % If we adjust the modulation score cutoff to be a lower value:
    mod_score_cutoff_4 = 0.4;
    output_filename_original = fullfile(output_path, 'modulation_rand_mod4');
    [p1r, p0r, mod_score_r, adj_mat_r] = MAGIC(data, randomized_group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff_4, output_filename_rand_mod4);
    
    network_rand = graph(adj_mat_r); 
    
    gene1_nodes = network_rand.Edges.EndNodes(:,1);
    gene2_nodes = network_rand.Edges.EndNodes(:,2);
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    weights = network_rand.Edges.Weight;
    gene_names = unique([gene1_names; gene2_names]);
    
    num_nodes = length(gene_names);
    new_adj_mat_r = zeros(num_nodes);
    
    for k = 1:length(gene1_nodes)
        i = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_adj_mat_r(i, j) = weights(k);
        new_adj_mat_r(j, i) = weights(k); 
    end
    
    edge_colors = assignEdgeColors(weights);
    new_network_rand = graph(new_adj_mat_r);
    
    node_degree = degree(new_network_rand);
    scaling_factor = 2; 
    node_sizes = scaling_factor * node_degree;
    
    figure;
    plot(new_network_rand, 'NodeLabel', gene_names, 'MarkerSize',node_sizes, 'EdgeColor',edge_colors, 'LineWidth',1.3, 'NodeFontSize', 5);
    axis off;
    %% 
    %% Modules
    bins = conncomp(new_network_rand);
    bin_nodes = cell(max(bins), 1);
    
    % group nodes according to bins
    % each cell in bin_nodes contains indices of all nodes that share a bin
    for l = 1:max(bins)
        bin_nodes{l} = find(bins == l);
    end
    
    % do not consider pairs
    % pairs are not strictly uninformative, but this decision was made
    % due to the fact given that denser communities were detected
    valid_bins = find(cellfun(@numel, bin_nodes) > 2);
    bin_nodes = bin_nodes(valid_bins);
    
    num_clusters = length(bin_nodes);
    num_cols = ceil(sqrt(num_clusters));
    num_rows = ceil(num_clusters/num_cols);
    cluster_gene_names_all_r = cell(num_clusters,1);
    
    figure;
    % build subnetworks for all the disconnected subnetworks
    for k = 1:num_clusters
        % subnetwork for current bin
        subnetwork_nodes = bin_nodes{k};
        cluster_gene_names = gene_names(subnetwork_nodes);
        cluster_gene_names_all_r{k} = cluster_gene_names;
        %cluster_gene_bio = bio_processes{subnetwork_nodes};
    
        subnetwork_adj_mat = new_adj_mat(subnetwork_nodes, subnetwork_nodes);
    
        subnetwork = graph(subnetwork_adj_mat);
    
        subplot(num_rows, num_cols, k);
        plot(subnetwork, 'Layout', 'force', 'NodeLabel',cluster_gene_names,'NodeFontSize',6);
        axis off;
        title(['Module ', num2str(k)]);
    
    end
    sgtitle('Modules of Differential Regulatory Gene Network with randomized grouping');
    
    %% Similarity: Genes in all clusters
    unique_n = vertcat(cluster_gene_names_all{:});
    unique_rand = vertcat(cluster_gene_names_all_r{:});
    intersection = intersect(unique_n, unique_rand);
    union_set = union(unique_n, unique_rand);
    
    % Jaccard Similarity: similarity by comparing intersection to union
    % 0: no common elements
    % 1: identical
    jaccard_similarity = numel(intersection) / numel(union_set);
    
    fprintf('Jaccard Similarity: %.4f\n', jaccard_similarity);
    %% Similartiy: Genes per cluster
    
    %empty matrix to store Jaccard similarities
    num_clusters_norm = length(cluster_gene_names_all);
    num_clusters_rand = length(cluster_gene_names_all_r);
    jaccard_matrix = zeros(num_clusters_norm, num_clusters_rand);
    
    % Jaccard similarity for each pair of clusters
    for n = 1:num_clusters_norm
        for j = 1:num_clusters_rand
            % intersection and union
            intersection = intersect(cluster_gene_names_all{n}, cluster_gene_names_all_r{j});
            union_set = union(cluster_gene_names_all{n}, cluster_gene_names_all_r{j});
            
            % Jaccard similarity
            jaccard_similarity = numel(intersection) / numel(union_set);
            jaccard_matrix(n, j) = jaccard_similarity;
        end
    end
    
    
    %% Functions
    
    function edge_colors = assignEdgeColors(weights)
        edge_colors = zeros(length(weights), 3);
        for w = 1:length(weights)
            if weights(w) == 2 || weights(w) == -1
                edge_colors(w, :) = [0, 0, 1];  % Blue for M=1 specific positive correlation or M=0 specific positive correlation
            elseif weights(w) == 1 || weights(w) == -2
                edge_colors(w, :) = [1, 0, 0];  % Red for M=1 specific negative correlation or M=0 specific negative correlation
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