%  matlab -r "reproduction_reference_data Results"
function reproduction_reference_data(output_path)
    %% Read data from two input files
    input1 = readtable("Datasets\out_CD8_exhausted.tsv", 'Delimiter', '\t', 'FileType','text');
    input2 = readtable("Datasets\out_Macrophages.tsv",'Delimiter', '\t', 'FileType','text');
    
    
    input2_Excl = input2(:, 2:end); % to remove gene column
    % store gene names
    geneNames = input1{:,1};
    % trim sample name as it is becomes too long + allows merging data
    input2_Excl.Properties.VariableNames = strrep(input2_Excl.Properties.VariableNames, "Hackathon_April2024", '');
    total_data = [input1, input2_Excl];
    %combine tables to match data matrix required by MAGIC function
    data_matrix = double(total_data{:, 2:end});
    
    %% Define GROUP
    num_vars_input1 = size(input1,2)-1; %subtract 1 for gene name column
    num_vars_input2 = size(input2_Excl, 2);
    
    % construct group vector required by MAGIC
    % assign 0 for data where modulator is not expressed, 1 for data where
    % modulator is expressed
    group = [zeros(1,num_vars_input1), ones(1,num_vars_input2)];
    
    
    
    %% Sample size average
    % since our data is separate, we can take average of samples by taking 
    % number of variables (which represent samples) of each dataset and average
    % them
    equ_sam_size = (num_vars_input1 + num_vars_input1)/2;
    
    
    %%
    data = data_matrix;
    group = group;
    bonf = 1;
    equ_sam_size = equ_sam_size;
    p_cutoff = 0.05;
    mod_score_cutoff = 0.6;
    output_filename = fullfile(output_path, 'modulation_reference');
    
    clear data_matrix mod_data_clean
    %% MAGIC algorithm
    [p1, p0, mod_score, adj_mat] = MAGIC(data, group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff, output_filename);
    
    %% Differential Regulatory Gene Network  
    % merge gene pairs
    % nodes: genes
    % edges: modulated interaction
    adj_values = unique(adj_mat);
    fprintf('Can be found in the Adjacency Matrix: %.4f\n', adj_values);
    network = graph(adj_mat); % network(i,j) = [1 OR 2 OR 0 OR -2 OR -1]
    
    % extract genes which are modulated
    gene1_nodes = network.Edges.EndNodes(:,1);
    gene2_nodes = network.Edges.EndNodes(:,2);
    
    % replace indices with gene names
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    
    % extract modulation:
    % 2: M=1 specific positive correlation; 1: M=1 specific negative correlation
    % -1: M=0 specific positive correlation; -2: M=0 specific negative correlation
    weights = network.Edges.Weight;
    
    % unique gene numbers (since we will create a new adjacency matrix)
    gene_names = unique([gene1_names; gene2_names]);
    
    % init symmetric adjacency matrix for only the modulated pairs
    num_nodes = length(gene_names);
    new_adj_mat = zeros(num_nodes);
    
    % weights
    for k = 1:length(gene1_nodes)
        l = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_adj_mat(l, j) = weights(k);
        new_adj_mat(j, l) = weights(k); 
    end
    
    
    % different edge colors for positive/negative interaction
    
    edge_colors = assignEdgeColors(weights);
    
    
    % new graph
    new_network = graph(new_adj_mat);
    
    % adjust sizes of nodes according to number of neighbors
    node_degree = degree(new_network);
    scaling_factor = 2;  % customizable to preference
    node_sizes = scaling_factor * node_degree;
    
    % plot
    figure;
    plot(new_network, 'NodeLabel', gene_names, 'MarkerSize',node_sizes, 'EdgeColor',edge_colors, 'LineWidth',1.3, 'NodeFontSize', 5);
    axis off;
    title("Differential Regulatory Gene Network of Macrophages and exhausted CD8+ T-cells")
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
    %% Finding hub genes with 3-sigma rule
    mean_degree = mean(node_degree);
    std_degree = std(node_degree);
    
    hub_gene_idx  = find(node_degree>mean_degree + 3*std_degree);
    hub_genes=gene_names(hub_gene_idx);
    
    
    %% Plot subnetworks for hub genes
    for j = 1:numel(hub_genes)
        % find hub gene index
        hub_gene_idx = find(strcmp(gene_names, hub_genes{j}));
        
        % hub gene connections
        connected_nodes = neighbors(new_network, hub_gene_idx);
        subnetwork_nodes = unique([hub_gene_idx; connected_nodes]);
        
        % select submatrix
        subnetwork_adj_mat = new_adj_mat(subnetwork_nodes, subnetwork_nodes);
        
        subnetwork = graph(subnetwork_adj_mat);
        
        % subnetwork properties
        subweights = subnetwork.Edges.Weight;
        subedge_colors = assignEdgeColors(subweights);
        subnode_degree = degree(subnetwork);
        subnode_sizes = scaling_factor * subnode_degree;
        
        % plot 
        figure;
        plot(subnetwork, 'NodeLabel', gene_names(subnetwork_nodes), 'MarkerSize', subnode_sizes, 'EdgeColor', subedge_colors, 'LineWidth', 1.3);
        axis off;
        title(['Subnetwork for Hub Gene: ', hub_genes{j}]);
    end
    
    %% Scatter Plots of Gene Pair with highest Modulation Score
    % we find the gene pair with the highest mod_score value (delta_I_adj) 
    % we costruct scatter plots for the gene expression under each condition [modulator expressed or not]
    % find raw correlation under each condition
    
    
    % to extract raw correlation under each condition separately
    out_network = readtable([output_filename, '.txt']);
    max_mod_score = max(out_network{:,6});
    max_mod_idx = find(out_network{:,6} == max_mod_score);
    
    raw_corrm1 = out_network{:,3}(max_mod_idx);
    raw_corrm0 = out_network{:,4}(max_mod_idx);
    
    mod_pair = [out_network{:,1}(max_mod_idx);out_network{:,2}(max_mod_idx)];
    max_mod_pair = geneNames(mod_pair);
    %% Scatter under M=0
    
    gene1_expression_m0_idx = strcmp(input1.Gene, max_mod_pair{1});
    gene1_expression_m0 = input1(gene1_expression_m0_idx,2:end);
    gene1_expression_m0 = gene1_expression_m0{:,:};
    
    gene2_expression_m0_idx = strcmp(input1.Gene, max_mod_pair{2});
    gene2_expression_m0 = input1(gene2_expression_m0_idx,2:end);
    gene2_expression_m0=gene2_expression_m0{:,:};
    
    figure;
    scatter(gene1_expression_m0,gene2_expression_m0, 'Marker', '.');
    
    xlabel(max_mod_pair{1});
    ylabel(max_mod_pair{2});
    title("Scatter Plot of Raw correlation under M=0");
    txt = {'Raw correlation = ', raw_corrm0, 'Modulation score = ', max_mod_score};
    annotation('textbox',[.9 .5 .1 .2], ...
        'String',txt,'EdgeColor','none')
    %% Scatter under M=1
    gene1_expression_m1_idx = strcmp(input2.Gene, max_mod_pair{1});
    gene1_expression_m1 = input2(gene1_expression_m1_idx,2:end);
    gene1_expression_m1 = gene1_expression_m1{:,:};
    
    gene2_expression_m1_idx = strcmp(input2.Gene, max_mod_pair{2});
    gene2_expression_m1 = input2(gene2_expression_m1_idx,2:end);
    gene2_expression_m1=gene2_expression_m1{:,:};
    
    figure;
    scatter(gene1_expression_m1,gene2_expression_m1, 'Marker', '.');
    
    xlabel(max_mod_pair{1});
    ylabel(max_mod_pair{2});
    title("Scatter Plot of Raw correlation under M=1");
    txt = {'Raw Correlation = ', raw_corrm1, 'Modulation score = ', max_mod_score};
    annotation('textbox',[.9 .5 .1 .2], ...
        'String',txt,'EdgeColor','none')
    
    %% Differential Regulatory Gene Network using Modulation Score as weight
    %create new adjacency matrix, where adjecency = modulation score
    
    % remove nan from modulation score matrix
    mod_no_nan = mod_score;
    mod_no_nan(isnan(mod_no_nan)) = 0;
    
    % element wise multiplication between adjacency matrix and modulation score
    % matrix
    combined_matrix = mod_no_nan .* adj_mat;
    
    % create new graph
    modulated_network = graph(combined_matrix);
    
    %% Build new matrix only with pairs
    % extract genes which are modulated
    gene1_nodes = modulated_network.Edges.EndNodes(:,1);
    gene2_nodes = modulated_network.Edges.EndNodes(:,2);
    
    % replace indices with gene names
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    
    % extract modulation:
    % 2: M=1 specific positive correlation; 1: M=1 specific negative correlation
    % -1: M=0 specific positive correlation; -2: M=0 specific negative correlation
    weights = modulated_network.Edges.Weight;
    
    % unique gene numbers (since we will create a new adjacency matrix)
    gene_names = unique([gene1_names; gene2_names]);
    
    % init symmetric adjacency matrix for only the modulated pairs
    num_nodes = length(gene_names);
    new_net = zeros(num_nodes);
    conditions = cell(size(weights));
    
    % Populate the adjacency matrix and condition information
    for k = 1:length(gene1_nodes)
        m = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_net(m, j) = weights(k);
        new_net(j, m) = weights(k); 
        
        % Determine condition based on the weight
        if weights(k) > 0
            conditions{k} = 'Modulated';
        else 
            conditions{k} = 'Not Modulated';
        end
    end
    
    
    % No graph because resulting graph is identical to original graph, with less information about type of correlation.
    
    %% Output Modulated File
    % network.tsv
        % Col 1 -> target
        % Col 2 -> regulator
        % Col 3 -> condition ???
        % Col 4 -> Weight
    
    
    output_file = fullfile(output_path,'network_modulated.tsv');
    fid = fopen(output_file, 'w');
    
    % header
    fprintf(fid, 'Target\tRegulator\tCondition\tWeight\n');
    target = gene1_names;
    regulator = gene2_names;
    w = weights; 
    
    
    % write data
    for k = 1:length(gene1_names)
        fprintf(fid, '%s\t%s\t%s\t%d\n', gene1_names{k}, gene2_names{k}, conditions{k}, weights(k));
    end
    
    
    fclose(fid);
    
    disp(['File ', output_file, ' created successfully!']);
    
    
    
    %% Functions
    
    function edge_colors = assignEdgeColors(weights)
        edge_colors = zeros(length(weights), 3);
        for wght = 1:length(weights)
            if weights(wght) == 2 || weights(wght) == -1
                edge_colors(wght, :) = [0, 0, 1];  % Blue for M=1 specific positive correlation or M=0 specific positive correlation
            elseif weights(wght) == 1 || weights(wght) == -2
                edge_colors(wght, :) = [1, 0, 0];  % Red for M=1 specific negative correlation or M=0 specific negative correlation
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