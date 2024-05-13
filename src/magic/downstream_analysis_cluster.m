function downstream_analysis_cluster(output_path)
    %% Read data from two input files%% Analysis of Algorithm
    % Comments explaining the original procedure have been left out for
    % legibility of method used for downstream analysis
    % Original construction with comments can be found:
        % reproduction_reference_data.m
    %% 1. MAGIC on CD8-exhausted and pretreated Macrophages data
    input1 = readtable("Datasets\out_CD8_exhausted.tsv", 'Delimiter', '\t', 'FileType','text');
    input2 = readtable("Datasets\out_Macrophages.tsv",'Delimiter', '\t', 'FileType','text');
    input2_Excl = input2(:, 2:end);
    geneNames = input1{:,1};
    input2_Excl.Properties.VariableNames = strrep(input2_Excl.Properties.VariableNames, "Hackathon_April2024", '');
    total_data = [input1, input2_Excl];
    data_matrix = double(total_data{:, 2:end});

    num_vars_input1 = size(input1,2)-1;
    num_vars_input2 = size(input2_Excl, 2);

    group = [zeros(1,num_vars_input1), ones(1,num_vars_input2)];

    equ_sam_size = (num_vars_input1 + num_vars_input1)/2;

    data = data_matrix;
    group = group;
    bonf = 1;
    equ_sam_size = equ_sam_size;
    p_cutoff = 0.05;
    mod_score_cutoff = 0.6;
    output_filename_original = fullfile(output_path, 'modulation');

    [p1, p0, mod_score, adj_mat] = MAGIC(data, group, bonf, equ_sam_size, p_cutoff, mod_score_cutoff, output_filename_original);

    %% Differential Regulatory Gene Network
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
        m = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_adj_mat(m, j) = weights(k);
        new_adj_mat(j, m) = weights(k); 
    end

    edge_colors = assignEdgeColors(weights);
    new_network = graph(new_adj_mat);

    node_degree = degree(new_network);
    scaling_factor = 2; 
    node_sizes = scaling_factor * node_degree;

    figure;
    plot(new_network, 'NodeLabel', gene_names, 'MarkerSize',node_sizes, 'EdgeColor',edge_colors, 'LineWidth',1.3, 'NodeFontSize', 5);
    axis off;
    title('DGRN for Macrophages and exhausted CD8+ T-Cells');

    % %% Analysis of Modules
    % 
    % % obtain GO Biological Process from all genes in dataset
    % fid = fopen("pantherGeneList.txt", 'r');
    % header = fgetl(fid);
    % variableNames = strsplit(header, '\t');
    % fclose(fid);
    % 
    % % read Gene Ontologi Biological Processes from the Panther file
    % go_bio_processes = readtable('pantherGeneList.txt','ReadVariableNames', false, 'Delimiter', '\t', 'HeaderLines',1);
    % bio_processes = cell(numel(gene_names), 1);
    % 
    % for k = 1:numel(gene_names)
    %     % find index where gene names match
    %     idx = strcmp(gene_names{k}, go_bio_processes.Var2);
    %     if any(idx)
    %         % extract biological process string
    %         bioprocess = go_bio_processes.Var4{idx};
    %         % parse into separate processes
    %         processes = strsplit(bioprocess, ';');
    %         % trim whitespace
    %         processes = strtrim(processes);
    %         bio_processes{k} = processes';
    %     else
    %         % Panther was not able to map all the genes. To perserve
    %         % dimensions:
    %         bio_processes{k} = {'no data'};
    %     end
    % end

    %% 2. Approaches to Identify Modules
    % Due to the fact that MAGIC does not generate edge weights that define
    % the strength of the interaction between two identified gene pairs,
    % three different approaches to identifying modules on the network were
    % tried:

        % a - Connected Components: This is the simplest and most
        % straightforward approach. The connected components (complete subnetworks) of the network
        % are considered modules. This is applied to the network itself.

        % b - Modularity Maximization Community Detection (Community
        % Detection Toolbox - Athanasios Kehagias): The algorithm aims to extract
        % the community structure of a network using modularity
        % optimization as a heuristic. Modularity maximization aims to find
        % the clusters whose edge density is maximally different from the expected edge density
        % in a random network with the same distribution. The function is
        % applied onai the adjacency matrix and is a direct implementation of
        % the method proposed in the paper "Fast unfolding of communities
        % in large networks"(Vincent D. Blondel, Jean-Loup Guillaume,
        % Renaud Lambiotte, Etienne Lefebvre), DOI: 10.1088/1742-5468/2008/10/P10008

        % c - Spectral Clustering: The built-in MATLAB spectralcluster was
        % applied to a new matrix that combines the adjacency matrix with
        % the matrix containing modulation scores for each pair, in an
        % attempt to come close to a weighted adjacency matrix. This is the
        % only apprach to require the number of clusters to be specified
        % and not discovered. 

    % Methods that don't require a predefined k (number of
    % clusters/classes) were preferred. Conlcusion included after the
    % implementation of all approaches. 

    %% 2a: Subgraphs with Connected components [network]
    % This was attempted as the simplest approach, motivated by the
    % topological disconnections in the network. 


    % assigns each node in the network to a bin, starting from 1
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
    
    figure;
    % build subnetworks for all the disconnected subnetworks
    for k = 1:num_clusters
        % subnetwork for current bin
        subnetwork_nodes = bin_nodes{k};
        cluster_gene_names = gene_names(subnetwork_nodes);

        %cluster_gene_bio = bio_processes{subnetwork_nodes};

        subnetwork_adj_mat = new_adj_mat(subnetwork_nodes, subnetwork_nodes);

        subnetwork = graph(subnetwork_adj_mat);
        cluster_edge_colors = assignEdgeColors(subnetwork.Edges.Weight);

        subplot(num_rows, num_cols, k);
        plot(subnetwork, 'Layout', 'force', 'NodeLabel',cluster_gene_names,'NodeFontSize',6,'EdgeColor',cluster_edge_colors);
        axis off;
        title(['Module ', num2str(k)]);

    end
    sgtitle('Modules of Differential Regulatory Gene Network of Macrophages and exhausted CD8+ T-Cells (Connected Components)');

    large_cluster_nodes = bin_nodes{1};
    large_names = gene_names(large_cluster_nodes);
    large_adj_mat = new_adj_mat(large_cluster_nodes,large_cluster_nodes);
    large_graph = graph(large_adj_mat);
    large_edge_colors = assignEdgeColors(large_graph.Edges.Weight);
    figure;
    plot(large_graph, 'Layout', 'force', 'NodeLabel',large_names,'NodeFontSize',6,'EdgeColor',large_edge_colors)
    title("Large Module");
    axis off;
    fork_nodes = bin_nodes{4};
    fork_names = gene_names(fork_nodes);
    fork_adj_mat = new_adj_mat(fork_nodes,fork_nodes);
    fork_graph = graph(fork_adj_mat);
    fork_edge_colors = assignEdgeColors(fork_graph.Edges.Weight);
    figure;
    plot(fork_graph, 'Layout', 'force', 'NodeLabel',fork_names,'NodeFontSize',6,'EdgeColor',fork_edge_colors)
    title("Fork Module");
    axis off;
    % Results: This section is brief, since the results were as expected.
    %% 2b: Modularity Maximization Community Detection [adjacency matrix]
    % Community Detection Toolbox - Athanasios Kehagias
    % This approach was chosed due to being readily available in the CM
    % Toolbox, and not requiring a predefined number of clusters. 

    % applied on the reduced adjacency matrix for a significantly improved
    % computation time
    cluster_find = GCModulMax1(new_adj_mat);
    clusters = unique(cluster_find);

    cluster_cells = cell(length(clusters), 1);
    for l = 1:length(clusters)
        cluster_cells{l} = find(cluster_find == clusters(l));
    end

    % conflicted on whether to only keep clusters with more that 3 genes
    valid_clusters = cellfun(@numel,cluster_cells)>2;
    filtered_cluster_cells = cluster_cells(valid_clusters);

    % subplot parameters
    num_clusters = length(filtered_cluster_cells);
    num_cols = ceil(sqrt(num_clusters));
    num_rows = ceil(num_clusters / num_cols);
    cluster_gene_names_all = cell(num_clusters,1);

    figure;
    % plot each cluster
    for k = 1:num_clusters
        cluster_indices = filtered_cluster_cells{k};
        cluster_gene_names = gene_names(cluster_indices);
        % extract submatrix corresponding to genes in cluster
        cluster_gene_names_all{k} = cluster_gene_names;

        sub_adj_mat = new_adj_mat(cluster_indices, cluster_indices);
        % graph
        cluster_graph = graph(sub_adj_mat);

        cluster_edge_colors = assignEdgeColors(cluster_graph.Edges.Weight);

        subplot(num_rows, num_cols, k);
        plot(cluster_graph, 'Layout', 'force', 'NodeLabel', cluster_gene_names, 'LineWidth',1.3, 'NodeFontSize',5, 'EdgeColor',cluster_edge_colors);
        title(sprintf('Module %d', k));
        axis off;

        % display number of nodes in cluster for quick check on whether the
        % clusters are of the same size across implementations (normal group vs randomized group)
        numgenes_txt = sprintf('Nodes: %d', length(cluster_indices));
        text(0.5, 0, numgenes_txt, 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');
    end
    sgtitle('Modules of Differential Regulatory Gene Network of Macrophages and exhausted CD8+ T-Cells (Modularity Maximization Community Detection)');

    % adjust to prevent empty subplot
    if num_clusters < num_cols * num_rows
        for m = num_clusters + 1 : num_cols * num_rows
            subplot(num_rows, num_cols, m);
            axis off;
        end
    end

    % Results: Notably, the dense subnetwork has been split into separate
    % clusters, and Cluster 7 (with the forked topology) is identical to Subnetwork 4 from the
    % connected components approach. 

    %% 2c: Spectral Clustering k-means [combined modulation-adjacency matrix]
    % The choice for spectral clustering with k-means was motivated by the
    % proven suitability of the algorithm for module detection in gene
    % expression data among other non-overlapping clustering methods.
    % W. Saelens, R. Cannoodt, and Y. Saeys, “A comprehensive evaluation of module detection methods for gene expression data,” 
    % Nature Communications, vol. 9, no. 1, Mar. 2018, doi: https://doi.org/10.1038/s41467-018-03424-4.‌‌


    % First, to inspect the modulation score matrix

    % larger positive elements have stronger correlation in M=1 samples compared to M=0.
    figure;
    imagesc(mod_score);
    colormap('cool');
    colorbar;
    title("Heatmap of Modulation Score between all genes in dataset");
    ylabel("Genes");
    xlabel("Genes");
    % Visible: Pink striped pattern

    % Second, create new adjacency matrix where adjecency = modulation score
    % The modulation score is being used as the similarity 

    % remove nan from modulation score matrix
    mod_no_nan = mod_score;
    mod_no_nan(isnan(mod_no_nan)) = 0;

    % element wise multiplication between adjacency matrix and modulation score
    % matrix
    % Some elements in the modulation matrix are negative, but this is
    % obviously not the case for pairs. Therefore, the element wise
    % multiplication handles this issue.
    combined_matrix = mod_no_nan .* adj_mat;

    % create new graph
    modulated_network = graph(combined_matrix);

    gene1_nodes = modulated_network.Edges.EndNodes(:,1);
    gene2_nodes = modulated_network.Edges.EndNodes(:,2);
    gene1_names = geneNames(gene1_nodes);
    gene2_names = geneNames(gene2_nodes);
    weights = modulated_network.Edges.Weight;
    gene_names = unique([gene1_names; gene2_names]);

    num_nodes = length(gene_names);
    new_mod_mat = zeros(num_nodes); 

    for k = 1:length(gene1_nodes)
        m = find(strcmp(gene_names, gene1_names{k}));
        j = find(strcmp(gene_names, gene2_names{k}));
        new_mod_mat(m, j) = weights(k);
        new_mod_mat(j, m) = weights(k); 
    end


    % Spectral Clustering

    % vector containing cluster to which that index is assigned
    cluster_find_m = spectralcluster(new_mod_mat,14,'Distance','precomputed');
    % k = 14, motivated by the number of connected components (2a) as an upper-bound approximation
    % much larger and much higher k values [2,30] were also tried. The
    % results are outlined in the Results section.
    % total clusters

    clusters_m = unique(cluster_find_m);

    % init cluster cells to hold all the genes in a cluster
    cluster_cells_m = cellfun(@(x) find(cluster_find_m == x), num2cell(clusters_m), 'UniformOutput', false);
    cluster_cells_m = cluster_cells_m(cellfun(@numel, cluster_cells_m) > 2);

    % for subplots
    num_clusters = length(cluster_cells_m);
    num_cols = ceil(sqrt(num_clusters));
    num_rows = ceil(num_clusters / num_cols);


    cluster_gene_names_all = cell(num_clusters,1);

    figure;

    for k = 1:num_clusters
        cluster_indices = cluster_cells_m{k};
        cluster_gene_names = gene_names(cluster_indices);
        % construct submatrix of combined_matrix
        sub_adj_mat_m = new_mod_mat(cluster_indices, cluster_indices);
        cluster_gene_names_all{k} = cluster_gene_names;

        cluster_graph_m = graph(sub_adj_mat_m);
        % can't assign colors to edges, since the modulation scores are
        % absolute values

        subplot(num_rows, num_cols, k);
        plot(cluster_graph_m, 'Layout', 'force', 'NodeLabel', cluster_gene_names, 'LineWidth',1.3, 'NodeFontSize',5);
        title(sprintf('Module %d', k));
        axis off;

        % display number of nodes in cluster for quick check on whether the
        % clusters are of the same size across implementations (normal group vs randomized group)
        numgenes_txt = sprintf('Nodes: %d', length(cluster_indices));
        text(0.5, 0, numgenes_txt, 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');
    end
    sgtitle('Modules of Differential Regulatory Gene Network of Macrophages and exhausted CD8+ T-Cells (Spectral Clustering using Modulation Score)');

    % prevent empty plots
    if num_clusters < num_cols * num_rows
        for l = num_clusters + 1 : num_cols * num_rows
            subplot(num_rows, num_cols, l);
            axis off;
        end
    end


    % Results: Interestingly, the clustering approach avoids separating the
    % large subnetwork even when larger k values are used. Similarly, the
    % forked cluster mentioned in the two apporaches above remains
    % consistent regardless of the k value chosen (within the specified
    % range).

    %% 2 CONCLUSION: CLUSTERING APPROACHES FOR THE DGRN BUILT USING MAGIC
    % For the validation and interpretation of the modules, Gene Set
    % Enrichment Analysis was performed. There was an unsuccesful attempt to implement
    % Gene Ontology Enrichment Analysis in MATLAB, therefore the web-based GOrilla tool was used
    % instead. GOEA was performed of each module identified by each
    % approach. Due to this, the modules in focus are the ones large enough
    % to have the analysis performed on. 

    % Notably, no enrichment was found for any of the clusters identified by the Modularity
    % Maximization approach, which succesfully split the large subnetwork.
    % The only exception was the forked module which appeared consistently
    % regardless of the clustering approach. 

    % The large module was enriched with only the GO Biological Function "cellular
    % component organization". The forked module was enriched with multiple GO cellular
    % components and GO biological processes, mainly relation to proton
    % transportation. 

    % All in all, using the connected component approach is likely to be
    % the most advantageous due to the fact that it is a straightforward
    % and computationally advantageous approach, especially if studying
    % larger gene expression datasets. 

     %% 3. Box Plot Clusters
     % For the identified clusters deemed relevant, create box plots of
     % expression data to compare. 

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
    clusters = cell(length(num_clusters));
    
    % build subnetworks for all the disconnected subnetworks
    for k = 1:num_clusters
        % subnetwork for current bin
        subnetwork_nodes = bin_nodes{k};
        cluster_gene_names = gene_names(subnetwork_nodes);
        clusters{k} = cluster_gene_names;
    
    end
    %%
    fork = clusters{4};
    matching_input1 = ismember(input1.Gene, fork);
    expression_1 = input1{matching_input1, 2:end};
    expression_1 = expression_1';
    
    matching_input2 = ismember(input2.Gene, fork);
    expression_2 = input2{matching_input2, 2:end};
    expression_2 = expression_2'; 
    
    pos1 = [1:numel(fork)];
    pos2 = [1.25:numel(fork)+1];
    figure;
    boxplot(expression_1,fork,'Positions',pos1);
    hold on;
    boxplot(expression_2, fork, 'Colors','m','Positions',pos2);
    
    
    hold off;
    title("Comparison between Expression Levels in Smaller Cluster Genes");
    xlabel('Genes In Cluster');
    ylabel('Normalized Gene Expression (TPM)');
    %%
    large = clusters{1};
    
    matching_input1 = ismember(input1.Gene, large);
    expression_1 = input1{matching_input1, 2:end};
    expression_1 = expression_1';
    
    matching_input2 = ismember(input2.Gene, large);
    expression_2 = input2{matching_input2, 2:end};
    expression_2 = expression_2'; 
    
    pos1 = [1:numel(large)];
    pos2 = [1.25:numel(large)+1];
    figure;
    boxplot(expression_1,large,'Positions',pos1);
    hold on;
    boxplot(expression_2, large, 'Colors','m','Positions',pos2);
    
    hold off;
    title("Comparison between Expression Levels in Large Cluster Genes");
    xlabel('Genes In Cluster');
    ylabel('Normalized Expression (TPM)');



    %% 4. Analyzing Condition-Specific Interactions
    % This involves inferring the Regulatory Gene Networks of the two
    % datasets from the DGRN by using the notion that positive edge weights
    % indicate correlation (or anti-correlation) in M=1 samples, and
    % negative edge weights indicated correlation in M=0 samples.

    % Initially, all the values in the adjacency matrix are checked, to
    % make sure that there are interactions specific to both conditions
    adj_values = unique(adj_mat);
    fprintf('Can be found in the Adjacency Matrix: %.4f\n', adj_values);
    % However, only positive edge weights are found.
    % This means that all the interaction pairs in the network are specific to Macrophages.

    % Trying a lower modulation score cutoff to find more viable pairs (0.5=<new_mod_score_cutoff<0.6)
    new_mod_score_cutoff = 0.5;
    output_filename_original =  fullfile(output_path, 'modulation_new_cutoff');
    
    % MAGIC algorithm
    [p1, p0, mod_score, adj_mat] = MAGIC(data, group, bonf, equ_sam_size, p_cutoff, new_mod_score_cutoff, output_filename_original);

    adj_values = unique(adj_mat);
    fprintf('Can be found in the Adjacency Matrix: %.4f\n', adj_values);

    submatrix_cd8 = zeros(size(adj_mat));

    submatrix_macrophages=zeros(size(adj_mat));
    for m = 1:size(adj_mat, 1)
        for l = 1:size(adj_mat, 2)
            if adj_mat(m, l) == 1 || adj_mat(m, l) == 2
                submatrix_macrophages(m, l) = adj_mat(m, l);
            elseif adj_mat(m, l) == -1 || adj_mat(m, l) == -2
                submatrix_cd8(m, l) = adj_mat(m, l);
            end
        end
    end

    net_macrophages = graph(submatrix_macrophages);
    gene1_nodes_mphages = net_macrophages.Edges.EndNodes(:,1);
    gene2_nodes_mphages = net_macrophages.Edges.EndNodes(:,2);   
    gene1_names_mphages = geneNames(gene1_nodes_mphages);
    gene2_names_mphages = geneNames(gene2_nodes_mphages);
    weights_mphages = net_macrophages.Edges.Weight;
    gene_names_mphages = unique([gene1_names_mphages; gene2_names_mphages]);


    num_nodes = length(gene_names_mphages);
    new_adj_mat_macrophages = zeros(num_nodes);

    for k = 1:length(gene1_nodes_mphages)
        n = find(strcmp(gene_names_mphages, gene1_names_mphages{k}));
        j = find(strcmp(gene_names_mphages, gene2_names_mphages{k}));
        new_adj_mat_macrophages(n, j) = weights_mphages(k);
        new_adj_mat_macrophages(j,n) = weights_mphages(k); 
    end

    new_net_macrophages = graph(new_adj_mat_macrophages);
    node_degree = degree(new_net_macrophages);
    scaling_factor = 1; 
    edge_colors = assignEdgeColors(weights_mphages);
    node_sizes = scaling_factor * node_degree;
    figure;
    plot(new_net_macrophages, 'NodeLabel', gene_names_mphages,'MarkerSize',node_sizes,'EdgeColor',edge_colors,'NodeFontSize',6);
    axis off;
    title('Interaction pairs that exhibit more correlation in Macrophages')
    %% Clusters from Macrophages with Connected components
    bins = conncomp(new_net_macrophages);
    bin_nodes = cell(max(bins), 1);

    % group nodes according to bins
    for l = 1:max(bins)
        bin_nodes{l} = find(bins == l);
    end

    % do not consider pairs
    valid_bins = find(cellfun(@numel, bin_nodes) > 2);
    bin_nodes = bin_nodes(valid_bins);

    num_clusters = length(bin_nodes);
    num_cols = ceil(sqrt(num_clusters));
    num_rows = ceil(num_clusters/num_cols);
    cluster_gene_names_all_mphages = cell(num_clusters,1);

    figure;
    % build subnetworks for all the disconnected subnetworks
    for j = 1:num_clusters
        % subnetwork for current bin
        subnetwork_nodes = bin_nodes{j};
        cluster_gene_names = gene_names_mphages(subnetwork_nodes);
        cluster_gene_names_all_mphages{j} = cluster_gene_names;
        subnetwork_adj_mat = new_adj_mat_macrophages(subnetwork_nodes, subnetwork_nodes);

        subnetwork = graph(subnetwork_adj_mat);

        subplot(num_rows, num_cols, j);
        plot(subnetwork, 'Layout', 'force', 'NodeLabel',cluster_gene_names,'NodeFontSize',6);
        axis off;
        title(['Module ', num2str(j)]);

    end
    sgtitle('Modules of pretreated Macrophages');

    %%
    net_cd8 = graph(submatrix_cd8);
    gene1_nodes_cd8 = net_cd8.Edges.EndNodes(:,1);
    gene2_nodes_cd8 = net_cd8.Edges.EndNodes(:,2);   
    gene1_names_cd8 = geneNames(gene1_nodes_cd8);
    gene2_names_cd8 = geneNames(gene2_nodes_cd8);
    weights_cd8 = net_cd8.Edges.Weight;

    gene_names_cd8 = unique([gene1_names_cd8; gene2_names_cd8]);

    num_nodes = length(gene_names_cd8);
    new_adj_mat_cd8 = zeros(num_nodes);

    for k = 1:length(gene1_nodes_cd8)
        m = find(strcmp(gene_names_cd8, gene1_names_cd8{k}));
        j = find(strcmp(gene_names_cd8, gene2_names_cd8{k}));
        new_adj_mat_cd8(m, j) = weights_cd8(k);
        new_adj_mat_cd8(j, m) = weights_cd8(k); 
    end

    new_net_cd8 = graph(new_adj_mat_cd8);
    node_degree = degree(new_net_cd8);
    scaling_factor = 1; 
    edge_colors = assignEdgeColors(weights_cd8);
    node_sizes = scaling_factor * node_degree;
    figure;
    plot(new_net_cd8, 'NodeLabel', gene_names_cd8, 'MarkerSize',node_sizes, 'EdgeColor',edge_colors,'NodeFontSize',6);
    axis off;
    title('Interaction pairs that exhibit more correlation in exhausted CD8+ T-cells')

    %% Clusters from exhausted CD8+ T-cells with Connected components
    bins = conncomp(new_net_cd8);
    bin_nodes = cell(max(bins), 1);

    % group nodes according to bins
    for m = 1:max(bins)
        bin_nodes{m} = find(bins == m);
    end

    % do not consider pairs
    valid_bins = find(cellfun(@numel, bin_nodes) > 2);
    bin_nodes = bin_nodes(valid_bins);

    num_clusters = length(bin_nodes);
    num_cols = ceil(sqrt(num_clusters));
    num_rows = ceil(num_clusters/num_cols);
    cluster_gene_names_all_cd8 = cell(num_clusters,1);

    figure;
    % build subnetworks for all the disconnected subnetworks
    for j = 1:num_clusters
        % subnetwork for current bin
        subnetwork_nodes = bin_nodes{j};
        cluster_gene_names = gene_names_cd8(subnetwork_nodes);
        cluster_gene_names_all_cd8{j} = cluster_gene_names;

        subnetwork_adj_mat = new_adj_mat_cd8(subnetwork_nodes, subnetwork_nodes);

        subnetwork = graph(subnetwork_adj_mat);

        subplot(num_rows, num_cols, j);
        plot(subnetwork, 'Layout', 'force', 'NodeLabel',cluster_gene_names,'NodeFontSize',6);
        axis off;
        title(['Module ', num2str(j)]);

    end
    sgtitle('Modules of exhausted CD8+ T-Cells');
    %% Box Plots of Clusters from CD8 T-cells
    
    module1 = cluster_gene_names_all_cd8{1};
    matching_input1 = ismember(input1.Gene, module1);
    expression_1 = input1{matching_input1, 2:end};
    expression_1 = expression_1';
    
    matching_input2 = ismember(input2.Gene, module1);
    expression_2 = input2{matching_input2, 2:end};
    expression_2 = expression_2'; 
    
    pos1 = [1:numel(module1)];
    pos2 = [1.25:numel(module1)+1];
    figure;
    boxplot(expression_1,module1,'Positions',pos1);
    hold on;
    boxplot(expression_2, module1, 'Colors','m','Positions',pos2);
    
    
    hold off;
    title("Comparison between Expression Levels in exhausted CD8+ T-cell Cluster Genes");
    xlabel('Genes In Cluster');
    ylabel('Normalized Gene Expression (TPM)');
    
    % RESULTS: AC090498.1 and HLA-C significantly overexpressed in CD8+
    % exhausted T-cells
    
    module1 = cluster_gene_names_all_cd8{2};
    matching_input1 = ismember(input1.Gene, module1);
    expression_1 = input1{matching_input1, 2:end};
    expression_1 = expression_1';
    
    matching_input2 = ismember(input2.Gene, module1);
    expression_2 = input2{matching_input2, 2:end};
    expression_2 = expression_2'; 
    
    pos1 = [1:numel(module1)];
    pos2 = [1.25:numel(module1)+1];
    figure;
    boxplot(expression_1,module1,'Positions',pos1);
    hold on;
    boxplot(expression_2, module1, 'Colors','m','Positions',pos2);
    
    
    hold off;
    title("Comparison between Expression Levels in exhausted CD8+ T-cell Cluster Genes");
    xlabel('Genes In Cluster');
    ylabel('Normalized Gene Expression (TPM)');

    % RESULTS: GO Enrichment Analysis highlights significance of this
    % cluster in regulation of immune response. No obvious pattern in
    % expression difference. While ACTB is more expressed in Macrophages
    % samples, it is highly expressed in both conditions
    %% Functions

    function edge_colors = assignEdgeColors(weights)
        edge_colors = zeros(length(weights), 3);
        for n = 1:length(weights)
            if weights(n) == 2 || weights(n) == -1
                edge_colors(n, :) = [0, 0, 1];  % Blue for M=1 specific positive correlation or M=0 specific positive correlation
            elseif weights(n) == 1 || weights(n) == -2
                edge_colors(n, :) = [1, 0, 0];  % Red for M=1 specific negative correlation or M=0 specific negative correlation
            end
        end
    end

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