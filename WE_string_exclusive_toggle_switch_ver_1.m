function WE_string_exclusive_toggle_switch_ver_1(number_of_string_nodes,taus,kappa,nreps,string_iter,Tmove,Tavg,timestep,endtime,save_location,parameter_set,radius,cores)

% WE_string_exclusive_toggle_switch_ver_1(number_of_string_nodes,taus,kappa,nreps,string_iter,Tmove,Tavg,timestep,endtime,save_location,gts_param,radius,cores)
% This function performs the WE-string method on an exclusive toggle switch
% BioNetGen ver 2.2.2 file and returns a matlab.mat file that returns the
% current string positions, the current replica positions/ weights, and a
% matrix containing all string positions.

% number_of_string_nodes = number of string nodes
% taus = taus parameter for string movement
% kappa = kappa parameter for string smoothening
% nreps = number replicas per bin
% string_iter = number of string movements to calculate
% Tmove = number of replica movements per string movement
% Tavg = average number of replica movements used in string movement
% timestep = timestep used by BioNetGen simulations
% endtime = total length of time to simulate each replica movement.
% save_location = filename for .mat containing data
% parameter_set = determines which parameter set (I,II,III) of the exclusive
% switch is used
% radius = radius of the hypersphere that determines the size of the
% attractor basins
% cores = number of cores matlabpool should open

% NOTE: The WE-string code written here was adapted from:
% Adelman, J. L., & Grabe, M. (2013). Simulating rare events using a weighted ensemble-based string method. The Journal of chemical physics, 138(4), 044105.
% This WE-string method script works specifically for biochemical networks written
% in BioNetGen. The WE-string method described in the above paper can be
% found at: https://simtk.org/home/westring

% Requires the Parallel Computing Toolbox and Matlab ver 2012b or higher to
% run
%
% Written by: Margaret Tse, tsemargaretld50@gmail.com

%*********************************************************************

%%INITIALIZATION OF SIMULATION
index_in_output_file_string_ind = 2; %indexer for saving string data

total_number_species = 5; %gives the number of species in the system.

%intializing the string nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_a_1 = 0;
initial_b_1 = 80;
initial_a_2 = 80;
initial_b_2 = 0;


initial_A_gene_1 = 0;
initial_A_gene_2 = 1;
initial_B_gene_1 = 1;
initial_B_gene_2 = 0;


string_initialization = linspace(0,1,number_of_string_nodes)'; %linear string for images
string_a_protein = (initial_a_2-initial_a_1)*string_initialization+initial_a_1; %the prerelaxation string
string_b_protein = (string_a_protein-initial_a_1)*(initial_b_2-initial_b_1)/(initial_a_2-initial_a_1)+initial_b_1;
string_A_gene = (string_a_protein - initial_a_1)*(initial_A_gene_2 - initial_A_gene_1)/(initial_a_2 - initial_a_1)+initial_A_gene_1;%y-coord of prerelaxation string
string_B_gene = (string_a_protein - initial_a_1)*(initial_B_gene_2 - initial_B_gene_1)/(initial_a_2 - initial_a_1)+initial_B_gene_1;

string_nodes_forwards = flipud(horzcat(string_a_protein,string_b_protein,zeros(length(string_a_protein),1),string_A_gene,string_B_gene));
string_nodes_backwards = flipud(string_nodes_forwards);

%***********************************************************************
model_name = ['exclusive_switch_' parameter_set '.net']; %name and location of the BioNetGen file. In local directory

bngsimroot = '~/BioNetGen-2.2.2-stable/bin/run_network_x86_64-linux'; %BionetGen simulation root

%%************************************************************
basinA = [80,0,0,1,0]; %Center of Basin A
basinB = [0,80,0,0,1]; %Center of Basin B
%******************************************************************

if matlabpool('size')==0 %opens 'cores' cores for parallel use if the pool isn't already open
    matlabpool(cores)
end

%%*****************************************************






%%**********************************************************************************************************

%%First replica movement



initial_reweighting_steps = 30; %Number of string movements that are used to reweight the Voronoi polyhedra
iteration = 1; %index for creating and storing BioNetGen simulations

BioNetGen_replica_folder = 'trajectories/'; %Creates folder to store BioNetGen simulations
mkdir(BioNetGen_replica_folder);

%initializing the first set of replicas with weight = delinit.
output_folder = ['trajectories/initial_conditions/exclusive_switch_' parameter_set '_init/'];
mkdir(output_folder);%location of the initial replicas


number_sim_steps = round(endtime/timestep);
initial_weight = 1/nreps;


matobj_save_info = matfile(save_location,'Writable',true); %saving the updated weights of the regions
%matobj2.weight_matrixA(1:n1,ind) = weightsA;

matobj_save_info.string_matrix_forwards = zeros(number_of_string_nodes,(string_iter+initial_reweighting_steps+1)*total_number_species);
matobj_save_info.string_matrix_backwards = zeros(number_of_string_nodes,(string_iter+initial_reweighting_steps+1)*total_number_species);


%%**************************
%Initial replica movements

%%initializing replica vectors
replicas_forwards = []; %replicas_forwards contains the replica locations at the end of simulation time.
replicas_backwards = [];

%%Initializing BioNetGen naming convention for parallelization,
%%history function and weights
rep_back_additional_info = [];
rep_for_additional_info = [];

for i = 1:number_of_string_nodes
    
    %%%******************************************************
    %%Generating initial BioNetGen .net files at the initial node positions
    
    output = [output_folder 'exclusive_switch_' parameter_set '_' num2str(i) '.net'];
    
    overwrite_file = regexp(fileread(model_name),'\n','split');
    
    a = round(string_nodes_forwards(i,1));
    b = round(string_nodes_forwards(i,2));
    o = round(string_nodes_forwards(i,3));
    oaa = round(string_nodes_forwards(i,4));
    obb = round(string_nodes_forwards(i,5));
    
    overwrite_file{25} = sprintf('    1 a()      %d',a);
    overwrite_file{26} = sprintf('    2 b()      %d',b);
    overwrite_file{27} = sprintf('    3 o()      %d',o);
    overwrite_file{28} = sprintf('    4 oaa()    %d',oaa);
    overwrite_file{29} = sprintf('    5 obb()    %d',obb);
    
    fid = fopen(output,'w');
    fprintf(fid,'%s\n',overwrite_file{:});
    fclose(fid);
    
    %%******************************************************
    tmp = (i-1)*nreps + 1;
    parfor ii = tmp:nreps-1+tmp %runs nreps parallel BioNetGen simulations and calculates the history function
        unix([bngsimroot ' -o ' BioNetGen_replica_folder 't0/newsim_' num2str(ii) ' -p ssa -h $RANDOM -a 1e-8 -r 1e-8 -e -g ' output ' ' output ' ' num2str(timestep) ' ' num2str(number_sim_steps) ' >> init.log']);
        current_replica_location = dlmread([BioNetGen_replica_folder 't0/newsim_' num2str(ii) '.gdat'],'',2,1);
        distsfrombasinB = pdist2(basinB,current_replica_location);
        distsfrombasinA = pdist2(basinA,current_replica_location);
        history_B_basin_curr = find(distsfrombasinB<radius,1,'last');
        history_A_basin_curr = find(distsfrombasinA<radius,1,'last');
        %%************* history function calculation
        if isempty(history_B_basin_curr) && isempty(history_A_basin_curr)
            history_function_A = 0; %no change
        elseif isempty(history_A_basin_curr) && not(isempty(history_B_basin_curr))
            history_function_A = 2; %basinB
        elseif isempty(history_B_basin_curr) && not(isempty(history_A_basin_curr))
            history_function_A = 1; %basinA
        else
            if history_A_basin_curr < history_B_basin_curr
                history_function_A = 2;
            elseif history_B_basin_curr < history_A_basin_curr
                history_function_A = 1;
            end
        end
        %%***************************
        if history_function_A == 2 %basin== 2 corresponds to data being in basinB
            
            
            %%creates replicas with the replica position from 1:
            replicas_backwards = vertcat(replicas_backwards,current_replica_location(end,:));
            
            %%plus a matrix with additional information, including
            %%BioNetGen file name for parallelization
            %%history function for each replica
            %%And weight of each replica
            rep_back_additional_info = vertcat(rep_back_additional_info,[ii,history_function_A,initial_weight]);
            
            
        else
            replicas_forwards = vertcat(replicas_forwards,[current_replica_location(end,:)]);
            rep_for_additional_info = vertcat(rep_for_additional_info,[ii,1,initial_weight]);
        end
    end
end

%%**************************same thing for initializing the backwards
%%string
for i = 1:number_of_string_nodes
    output = [output_folder 'exclusive_switch_' parameter_set '_' num2str(i) '.net'];
    
    
    
    tmp = (number_of_string_nodes+i-1)*nreps + 1;
    parfor jj = tmp:nreps-1+tmp
        unix([bngsimroot ' -o ' BioNetGen_replica_folder 't0/newsim_' num2str(jj) ' -p ssa -h $RANDOM -a 1e-8 -r 1e-8 -e -g ' output ' ' output ' ' num2str(timestep) ' ' num2str(number_sim_steps) ' >> init.log']);
        
        
        
        current_replica_location = dlmread([BioNetGen_replica_folder 't0/newsim_' num2str(jj) '.gdat'],'',2,1);
        distsfrombasinB = pdist2(basinB,current_replica_location);
        distsfrombasinA = pdist2(basinA,current_replica_location);
        history_B_basin_curr = find(distsfrombasinB<radius,1,'last');
        history_A_basin_curr = find(distsfrombasinA<radius,1,'last');
        if isempty(history_B_basin_curr) && isempty(history_A_basin_curr)
            history_function_B = 0; %no change
        elseif isempty(history_A_basin_curr) && not(isempty(history_B_basin_curr))
            history_function_B = 2; %basinB
        elseif isempty(history_B_basin_curr) && not(isempty(history_A_basin_curr))
            history_function_B = 1; %basinA
        else
            if history_A_basin_curr < history_B_basin_curr
                history_function_B = 2;
            elseif history_B_basin_curr < history_A_basin_curr
                history_function_B = 1;
            end
        end
        
        if history_function_B == 1 %basin== 2 corresponds to data being in basinB
            replicas_forwards = vertcat(replicas_forwards,[current_replica_location(end,:)]);
            rep_for_additional_info = vertcat(rep_for_additional_info,[jj,history_function_B,initial_weight]);
            
            
        else
            replicas_backwards = vertcat(replicas_backwards,[current_replica_location(end,:)]);
            rep_back_additional_info = vertcat(rep_back_additional_info,[jj,2,initial_weight]);
        end
        
        
    end
end

%%**************************




%% renormalizing the replica weights. Only done at the beginning.
renorm = sum(rep_for_additional_info(:,end)) + sum(rep_back_additional_info(:,end));
rep_for_additional_info(:,end) = rep_for_additional_info(:,end)./renorm;
rep_back_additional_info(:,end) = rep_back_additional_info(:,end)./renorm;


%%***************************
%%WE weight redistribution step for forwards string

bin_locations_forwards_reps = knnsearch(string_nodes_forwards,replicas_forwards); %gives the bins of the replica endpts
[temp_replica_info, weights_temp,~] = reweight([replicas_forwards,rep_for_additional_info(:,1:2)], rep_for_additional_info(:,3), bin_locations_forwards_reps, number_of_string_nodes, nreps, string_nodes_forwards);
replicas_forwards_weighted = temp_replica_info(:,1:total_number_species);
rep_for_add_info_weighted = [temp_replica_info(:,total_number_species+1:end),weights_temp];


%%WE weight redistribution step for backwards string

bin_locations_backwards_reps = knnsearch(string_nodes_backwards,replicas_backwards); %gives the bins of the replica endpts
[temp_replica_info, weights_temp,~] = reweight([replicas_backwards,rep_back_additional_info(:,1:2)], rep_back_additional_info(:,3), bin_locations_backwards_reps, number_of_string_nodes, nreps, string_nodes_backwards);
replicas_backwards_weighted = temp_replica_info(:,1:total_number_species);
rep_back_add_info_weighted = [temp_replica_info(:,total_number_species+1:end),weights_temp];


previous_BNG_dir = [BioNetGen_replica_folder '/t0'];



%%________________________________________________________________________
string_nodes_forwards_new = string_nodes_forwards;
string_nodes_backwards_new = string_nodes_backwards;



for ij = 1:initial_reweighting_steps
    
    %%****************Creating empty matrix to store replicas and bin
    %%locations for the string movement step over Tavg
    cachesize = 2*number_of_string_nodes*Tavg*nreps;
    cacherepsavg_for = zeros(cachesize,total_number_species+1);
    cacherepsavg_back = zeros(cachesize,total_number_species+1);
    cacherepbins_for = zeros(cachesize,1);
    cacherepbins_back = zeros(cachesize,1);
    
    %%********Simulates Tmove simulation and weight update steps
    for jk = 1:Tmove
        
        string_nodes_forwards = string_nodes_forwards_new; %initializing string data for the loop
        string_nodes_backwards = string_nodes_backwards_new;
        
        %%%*********Generates new folder to store new BioNetGen replicas
        newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
        mkdir(newBNGdir);
        iteration = iteration + 1;
        %%*****Merging forwards and backwards replica data for
        %%parallization purposes
        temp_replicas = [replicas_forwards_weighted; replicas_backwards_weighted];
        temp_add_info = [rep_for_add_info_weighted; rep_back_add_info_weighted];
        
        %%
        parallel_BNG_name = temp_add_info(:,1); %Ensures BNG is propogating the correct replica
        previous_basin_visited = temp_add_info(:,2); %Maintains history function
        weight = temp_add_info(:,3); %Contains replica weights
        
        new_replicas = []; %empty replica matrix for new replica positions
        new_add_rep_info = [];
        [num_replicas,~] = size(temp_replicas);
        
        parfor kk = 1:num_replicas
            unix([bngsimroot ' -o ' newBNGdir '/newsim_' num2str(kk) ' -p ssa -h $RANDOM -a 1e-8 -r 1e-8 -e -g ' previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' num2str(timestep) ' ' num2str(number_sim_steps) ' >> bng.log']);
            current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
            distsfrombasinB = pdist2(basinB,current_replica_location);
            distsfrombasinA = pdist2(basinA,current_replica_location);
            history_B_basin_curr = find(distsfrombasinB<radius,1,'last');
            history_A_basin_curr = find(distsfrombasinA<radius,1,'last');
            if isempty(history_B_basin_curr) && isempty(history_A_basin_curr)
                history_func_temp = 0; %no change
            elseif isempty(history_A_basin_curr) && not(isempty(history_B_basin_curr))
                history_func_temp = 2; %basinB
            elseif isempty(history_B_basin_curr) && not(isempty(history_A_basin_curr))
                history_func_temp = 1; %basinA
            else
                if history_A_basin_curr < history_B_basin_curr
                    history_func_temp = 2;
                elseif history_B_basin_curr < history_A_basin_curr
                    history_func_temp = 1;
                end
            end
            
            if history_func_temp == 0
                
                new_replicas = vertcat(new_replicas,[current_replica_location(end,:)]);
                new_add_rep_info = vertcat(new_add_rep_info,[kk,previous_basin_visited(kk),weight(kk)]);
            else
                new_replicas = vertcat(new_replicas,[current_replica_location(end,:)]);
                new_add_rep_info = vertcat(new_add_rep_info,[kk,history_func_temp,weight(kk)]);
            end
        end
        
        previous_BNG_dir = newBNGdir; %%resets current BNG directory
        
        %%sorts replicas back out into forwards and backwards replicas
        replicas_forwards = new_replicas(new_add_rep_info(:,2)==1,:);
        replicas_backwards = new_replicas(new_add_rep_info(:,2)==2,:);
        rep_for_additional_info = new_add_rep_info(new_add_rep_info(:,2)==1,:);
        rep_back_additional_info = new_add_rep_info(new_add_rep_info(:,2)==2,:);
        %%***************************
        %%WE weight redistribution step for forwards string
        
        bin_locations_forwards_reps = knnsearch(string_nodes_forwards,replicas_forwards); %gives the bins of the replica endpts
        [temp_replica_info, weights_temp, bins_fowards] = reweight([replicas_forwards,rep_for_additional_info(:,1:2)], rep_for_additional_info(:,3), bin_locations_forwards_reps, number_of_string_nodes, nreps, string_nodes_forwards);
        replicas_forwards_weighted = temp_replica_info(:,1:total_number_species);
        rep_for_add_info_weighted = [temp_replica_info(:,total_number_species+1:end),weights_temp];
        
        
        %%WE weight redistribution step for backwards string
        
        bin_locations_backwards_reps = knnsearch(string_nodes_backwards,replicas_backwards); %gives the bins of the replica endpts
        [temp_replica_info, weights_temp, bins_backwards] = reweight([replicas_backwards,rep_back_additional_info(:,1:2)], rep_back_additional_info(:,3), bin_locations_backwards_reps, number_of_string_nodes, nreps, string_nodes_backwards);
        replicas_backwards_weighted = temp_replica_info(:,1:total_number_species);
        rep_back_add_info_weighted = [temp_replica_info(:,total_number_species+1:end),weights_temp];
        
        
        
        %%%***********Storing replica bins and locations for string
        %%%movement
        if jk == (Tmove-Tavg)+1 %When Tavg begins, the below code stores the data
            [index_replica_avg_for,~] = size(replicas_forwards_weighted);
            [index_replica_avg_back,~] = size(replicas_backwards_weighted);
            cacherepsavg_for(1:index_replica_avg_for,:) = [replicas_forwards_weighted, rep_for_add_info_weighted(:,3)];
            cacherepsavg_back(1:index_replica_avg_back,:) = [replicas_backwards_weighted,rep_back_add_info_weighted(:,3)];
            
            cacherepbins_for(1:index_replica_avg_for) = bins_fowards;
            cacherepbins_back(1:index_replica_avg_back) = bins_backwards;
        end
        if jk>(Tmove-Tavg)+1
            [index_replica_avg_for_2,~] = size(replicas_forwards_weighted);
            [index_replica_avg_back_2,~] = size(replicas_backwards_weighted);
            cacherepsavg_for(index_replica_avg_for+1:index_replica_avg_for+index_replica_avg_for_2,:) = [replicas_forwards_weighted, rep_for_add_info_weighted(:,3)];
            cacherepsavg_back(index_replica_avg_back+1:index_replica_avg_back+index_replica_avg_back_2,:) = [replicas_backwards_weighted,rep_back_add_info_weighted(:,3)];
            
            cacherepbins_for(index_replica_avg_for+1:index_replica_avg_for+index_replica_avg_for_2) = bins_fowards;
            cacherepbins_back(index_replica_avg_back+1:index_replica_avg_back+index_replica_avg_back_2) = bins_backwards;
            index_replica_avg_for = index_replica_avg_for + index_replica_avg_for_2;
            index_replica_avg_back = index_replica_avg_back + index_replica_avg_back_2;
            
        end
        
        %%% Saving current replica locations + weights to the .mat file
        matobj_save_info.current_replicas_forwards = [replicas_forwards_weighted, rep_for_add_info_weighted(:,3)];
        matobj_save_info.current_replicas_backwards = [replicas_backwards_weighted,rep_back_add_info_weighted(:,3)];
        
    end
    
    
    cacherepsavg_for = cacherepsavg_for(any(cacherepsavg_for,2),:); %removes excess zeros from the average rep matrices
    cacherepsavg_back = cacherepsavg_back(any(cacherepsavg_back,2),:);
    cacherepbins_for = cacherepbins_for(any(cacherepbins_for,2),:);
    cacherepbins_back = cacherepbins_back(any(cacherepbins_back,2),:);
    
    %%Creates a weighted average of the cached replicas
    avg_for = average_fnc(cacherepsavg_for,number_of_string_nodes,cacherepbins_for,total_number_species,string_nodes_forwards); %caculates the average position of all stored replicas
    avg_back = average_fnc(cacherepsavg_back,number_of_string_nodes,cacherepbins_back,total_number_species,string_nodes_backwards);
    
    %%Moves the strings towards the average position in each bin and
    %%performs the string smoothing and reparameterization steps
    
    
    [string_nodes_forwards_new] = move_string(string_nodes_forwards,avg_for,number_of_string_nodes,taus,kappa);
    [string_nodes_backwards_new] = move_string(string_nodes_backwards,avg_back,number_of_string_nodes,taus,kappa);
    
    %%saves the string data to a growing matrix
    matobj_save_info.string_matrix_forwards(1:number_of_string_nodes,total_number_species*(index_in_output_file_string_ind-1)+1:total_number_species*index_in_output_file_string_ind) = string_nodes_forwards_new;
    matobj_save_info.string_matrix_backwards(1:number_of_string_nodes,total_number_species*(index_in_output_file_string_ind-1)+1:total_number_species*index_in_output_file_string_ind)= string_nodes_backwards_new;
    
    %%saves the current string information
    matobj_save_info.curr_string_forwards = string_nodes_forwards_new;
    matobj_save_info.curr_string_backwards = string_nodes_backwards_new;
    
    %%Indexes the next position in the .mat save file
    index_in_output_file_string_ind = index_in_output_file_string_ind+1;
    
end




for ij = 1:string_iter %%String movements without Tavg or Tmove
    
    
    
    string_nodes_forwards = string_nodes_forwards_new; %initializing string data for the loop
    string_nodes_backwards = string_nodes_backwards_new;
    
    %%%*********Generates new folder to store new BioNetGen replicas
    newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
    mkdir(newBNGdir);
    iteration = iteration + 1;
    %%*****Merging forwards and backwards replica data for
    %%parallization purposes
    temp_replicas = [replicas_forwards_weighted; replicas_backwards_weighted];
    temp_add_info = [rep_for_add_info_weighted; rep_back_add_info_weighted];
    
    %%
    parallel_BNG_name = temp_add_info(:,1); %Ensures BNG is propogating the correct replica
    previous_basin_visited = temp_add_info(:,2); %Maintains history function
    weight = temp_add_info(:,3); %Contains replica weights
    
    new_replicas = []; %empty replica matrix for new replica positions
    new_add_rep_info = [];
    [num_replicas,~] = size(temp_replicas);
    
    parfor kk = 1:num_replicas
        unix([bngsimroot ' -o ' newBNGdir '/newsim_' num2str(kk) ' -p ssa -h $RANDOM -a 1e-8 -r 1e-8 -e -g ' previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' num2str(timestep) ' ' num2str(number_sim_steps) ' >> bng.log']);
        current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
        distsfrombasinB = pdist2(basinB,current_replica_location);
        distsfrombasinA = pdist2(basinA,current_replica_location);
        history_B_basin_curr = find(distsfrombasinB<radius,1,'last');
        history_A_basin_curr = find(distsfrombasinA<radius,1,'last');
        if isempty(history_B_basin_curr) && isempty(history_A_basin_curr)
            history_func_temp_2 = 0; %no change
        elseif isempty(history_A_basin_curr) && not(isempty(history_B_basin_curr))
            history_func_temp_2 = 2; %basinB
        elseif isempty(history_B_basin_curr) && not(isempty(history_A_basin_curr))
            history_func_temp_2 = 1; %basinA
        else
            if history_A_basin_curr < history_B_basin_curr
                history_func_temp_2 = 2;
            elseif history_B_basin_curr < history_A_basin_curr
                history_func_temp_2 = 1;
            end
        end
        
        if history_func_temp_2 == 0
            
            new_replicas = vertcat(new_replicas,[current_replica_location(end,:)]);
            new_add_rep_info = vertcat(new_add_rep_info,[kk,previous_basin_visited(kk),weight(kk)]);
        else
            new_replicas = vertcat(new_replicas,[current_replica_location(end,:)]);
            new_add_rep_info = vertcat(new_add_rep_info,[kk,history_func_temp_2,weight(kk)]);
        end
    end
    
    previous_BNG_dir = newBNGdir; %%resets current BNG directory
    
    %%sorts replicas back out into forwards and backwards replicas
    replicas_forwards = new_replicas(new_add_rep_info(:,2)==1,:);
    replicas_backwards = new_replicas(new_add_rep_info(:,2)==2,:);
    rep_for_additional_info = new_add_rep_info(new_add_rep_info(:,2)==1,:);
    rep_back_additional_info = new_add_rep_info(new_add_rep_info(:,2)==2,:);
    %%***************************
    %%WE weight redistribution step for forwards string
    
    bin_locations_forwards_reps = knnsearch(string_nodes_forwards,replicas_forwards); %gives the bins of the replica endpts
    [temp_replica_info, weights_temp, bins_forwards] = reweight([replicas_forwards,rep_for_additional_info(:,1:2)], rep_for_additional_info(:,3), bin_locations_forwards_reps, number_of_string_nodes, nreps, string_nodes_forwards);
    replicas_forwards_weighted = temp_replica_info(:,1:total_number_species);
    rep_for_add_info_weighted = [temp_replica_info(:,total_number_species+1:end),weights_temp];
    
    
    %%WE weight redistribution step for backwards string
    
    bin_locations_backwards_reps = knnsearch(string_nodes_backwards,replicas_backwards); %gives the bins of the replica endpts
    [temp_replica_info, weights_temp, bins_backwards] = reweight([replicas_backwards,rep_back_additional_info(:,1:2)], rep_back_additional_info(:,3), bin_locations_backwards_reps, number_of_string_nodes, nreps, string_nodes_backwards);
    replicas_backwards_weighted = temp_replica_info(:,1:total_number_species);
    rep_back_add_info_weighted = [temp_replica_info(:,total_number_species+1:end),weights_temp];
    
    
    
    
    %%Saves replica information with weights
    matobj_save_info.current_replicas_forwards = [replicas_forwards_weighted, rep_for_add_info_weighted(:,3)];
    matobj_save_info.current_replicas_backwards = [replicas_backwards_weighted,rep_back_add_info_weighted(:,3)];
    
    %%calculates the average replica position per bin
    avg_for = average_fnc([replicas_forwards_weighted,rep_for_add_info_weighted(:,3)],number_of_string_nodes,bins_forwards,total_number_species,string_nodes_forwards); %caculates the average position of all stored replicas
    avg_back = average_fnc([replicas_backwards_weighted,rep_back_add_info_weighted(:,3)],number_of_string_nodes,bins_backwards,total_number_species,string_nodes_backwards);
    
    %% moves the string
    [string_nodes_forwards_new] = move_string(string_nodes_forwards,avg_for,number_of_string_nodes,taus,kappa);
    [string_nodes_backwards_new] = move_string(string_nodes_backwards,avg_back,number_of_string_nodes,taus,kappa);
    %%Saves the string matrix and current string
    matobj_save_info.string_matrix_forwards(1:number_of_string_nodes,total_number_species*(index_in_output_file_string_ind-1)+1:total_number_species*index_in_output_file_string_ind) = string_nodes_forwards_new;
    matobj_save_info.string_matrix_backwards(1:number_of_string_nodes,total_number_species*(index_in_output_file_string_ind-1)+1:total_number_species*index_in_output_file_string_ind)= string_nodes_backwards_new;
    matobj_save_info.curr_string_forwards = string_nodes_forwards_new;
    matobj_save_info.curr_string_backwards = string_nodes_backwards_new;
    
    %%indexer for matlab .mat output file
    index_in_output_file_string_ind = index_in_output_file_string_ind+1;
    
end



matlabpool close



end

function [output_replicas, output_weights, output_bins] = reweight(replicas, weights, replica_positions, total_number_of_bins, nreps,string_position)

%%******************************************************
%reweight takes in a matrix of the replica positions, assuming
%that each column is a new species, and each row is a new replica. The
%weights of each replica are assumed to be in a separate column vector.
%bin_centers is a matrix that contains a matrix of row vectors that
%determine the center of each bin.
%total_number_of_bins gives how many bins there are in bin_centeres
%nreps is the target number of replicas in each bin.
%replica_positions gives the bin associated with a given replica

%%********************************************************
%%This code performs the weighted ensemble method on the replicas, weights,
%%and bin_centers given. It ensures that there are nreps replicas in every
%%total_number_of_bins, given that there is at least one replica in the
%%bin. This code will also split and cull replicas that have weights above
%%or below a threshold.

%%*******************************************************
%outputs the redistributed weights in output_weights and the split and
%culled replicas.
%%******************************************************


output_replicas = []; %initializes the output of the code as empty matrixes
output_weights = [];
for i = 1:total_number_of_bins %iterates though each bin
    IndexesInBin = find(replica_positions == i); %Finds the indicies of the replicas that are in bin i
    
    if isempty(IndexesInBin) %if the bin is currently empty, no weight update step is performed.
    else
        ReplicasInBin = replicas(IndexesInBin,:); %Finds the replicas in the current bin
        WeightsInBin = weights(IndexesInBin,1); %'' weights
        [Current_rep_count,~] = size(ReplicasInBin); %Counts the number of replicas in the bin.
        
        [MaxWeight,idx_max] = max(WeightsInBin);
        %************Splits the large replica outliers
        while MaxWeight > 3*sum(WeightsInBin)/Current_rep_count % identifies replicas that are thrice the idea weight and splits them until they are small
            duplicates = 2;
            while MaxWeight/duplicates > 3*sum(WeightsInBin)/Current_rep_count && duplicates < Current_rep_count %this splits the replicas into m copies, but if m is larger than the current #  of reps.
                duplicates = duplicates+1;
            end
            WeightsInBin(end+1:end+duplicates-1,1) = MaxWeight/duplicates;
            ReplicasInBin(end+1:end+duplicates-1,:) = repmat(ReplicasInBin(idx_max,:),(duplicates-1),1); %creates duplicates copies of the replica
            WeightsInBin(idx_max,1) = MaxWeight/duplicates;
            [Current_rep_count,~] = size(ReplicasInBin);
            [MaxWeight,idx_max] = max(WeightsInBin);
            
        end
        
        %*********** Culls the small replica outliers if they are 1/4th the
        %ideal weight
        [SortedWeights, SortedIndex] = sort(WeightsInBin); %sorts weights from least to greatest
        MinWeight = SortedWeights(1);
        [Current_rep_count,~] = size(ReplicasInBin); %Counts the number of replicas in the bin.
        
        
        while MinWeight < MinWeight/sum(SortedWeights)/4 && Current_rep_count > 2
            index_min = 1;
            while MinWeight < MinWeight/sum(SortedWeights)/4 && index_min < Current_rep_count - 1
                
                index_min = index_min+1;
                MinWeight = sum(SortedWeights(1:index_min));
                if MinWeight > SortedWeights(1)/Current_rep_count*1.5
                    index_min = index_min-1;
                    break
                end
            end
            
            
            if index_min > 1
                weight_choice = (SortedWeights(1:index_min));
                new_pos = randsample(ReplicasInBin(SortedIndex(1:index_min))',1,true,weight_choice);
                ReplicasInBin(end+1,:) = new_pos;
                WeightsInBin(end+1,1) = sum(SortedWeights(1:index_min));
                ReplicasInBin(SortedIndex(1:index_min)) = [];
                WeightsInBin(SortedIndex(1:index_min)) = [];
                [Current_rep_count,~] = size(ReplicasInBin);
            else
                break
            end
            [SortedWeights, SortedIndex] = sort(WeightsInBin); %sorts weights from least to greatest
            MinWeight = SortedWeights(1);
        end
        
        
        
        
        
        
        
        %*********** Splits replicas
        while Current_rep_count < nreps %while the number of current replicas is less that the desired number of replicas...
            [MaxWeight,idx_max] = max(WeightsInBin); %find the max weigght
            WeightsInBin(end+1,1) = MaxWeight/2; %duplicates weight
            WeightsInBin(idx_max,1) = MaxWeight/2; %halve the weight of max weight replica
            ReplicasInBin(end+1,:) = ReplicasInBin(idx_max,:); %duplicates replica
            Current_rep_count=Current_rep_count+1; %increases count
        end
        
        %************Culls replicas
        while Current_rep_count > nreps
            [SortedWeights, SortedIndex] = sort(WeightsInBin); %sorts weights from least to greatest
            MergeWeights = SortedWeights(1:2); %finds the two smallest weights to merge
            temp = rand(1,1); %picks a random number from a uniform distribution
            if temp <= MergeWeights(1)/sum(MergeWeights) %choses one of the two replica positions
                WeightsInBin(SortedIndex(1),1) = sum(MergeWeights); %changes the weight of the chosen replica
                WeightsInBin(SortedIndex(2),:) = []; %removes the discarded replica
                ReplicasInBin(SortedIndex(2),:) = [];
            else
                WeightsInBin(SortedIndex(2),1) = sum(MergeWeights);
                WeightsInBin(SortedIndex(1),:) = [];
                ReplicasInBin(SortedIndex(1),:) = [];
                
            end
            [Current_rep_count,~] = size(ReplicasInBin); %Counts the number of replicas in the bin.
            
        end
        
        %**********Generates final replicas
        output_replicas = [output_replicas; ReplicasInBin]; %appends the new replicas to the output matrices.
        output_weights = [output_weights; WeightsInBin];
        
    end
    
    
end

[~,species] = size(string_position);
output_bins = knnsearch(string_position, output_replicas(:,1:species));


end

function [average] = average_fnc(replicas,total_number_of_bins,replica_bins,total_number_species,string_position) %finds the average weighted position of

%%******************************************************
%average_fnc takes in a matrix of the replica positions, assuming
%that each column is a new species, and each row is a new replica. The
%weights of each replica are assumed to be on the final replica column.
%replica_bins gives the bin location of each replica
%total_number_species gives the total number of species (columns) in the
%replica matrix
%string_position gives the matrix containing string nodes
%%********************************************************
%%This code finds the weighted average position of the replicas in each bin

%%*******************************************************
%outputs the matrix of the average positions in 'average'
%%******************************************************

average = zeros(total_number_of_bins,total_number_species);

weightedrep = bsxfun(@times,replicas(:,end),replicas(:,1:total_number_species));

for kk = 1:total_number_of_bins
    temp = find(replica_bins == kk);
    if isempty(temp)
        average(kk,1:total_number_species) = string_position(kk,1:total_number_species);
    else
        
        average(kk,1:total_number_species) = sum(weightedrep(temp,:),1)...
            ./(sum(replicas(temp,end),1));
    end
end
end

function [output_string] = move_string(initial_string,average_string_position,total_number_of_bins,taus,kappa)

%%******************************************************
%move_string takes in a the matrix of the current string position in
%initial_string, the average position of the replicas in
%average_string_position, the total number of string nodes in
%total_number_of_bins, the string movement parameter in taus, and the
%string smoothening parameter in kappa
%%********************************************************
%%This code moves the string towards the average position of replicas in
%%each bin, smoothes the string according to kappa and reparameterizes the
%%string such that the nodes are equidistant from each other
%%*******************************************************
%outputs the final string position in 'output_string'
%%******************************************************


%moves thes tring according to the average position and the parameters for
%movement and smoothening

%%Moves the strings to the weighted average of the replica location. The endpoints of the string are averaged against each other.
temp_string = initial_string + taus.*(average_string_position-initial_string);


%1 fixes the ends of the string
temp_string(end,:) = initial_string(end,:);
temp_string(1,:)= initial_string(1,:);

%string smoothening step
temp_string_smooth = temp_string;

for i = 2:total_number_of_bins-1
    temp_string_smooth(i,:) = temp_string(i,:) + kappa.*(temp_string(i+1,:)+temp_string(i-1,:) - 2*temp_string(i,:));
end

%%string reparameterization step***********

%%Finding the euclidian distance between points on the string.
eucdistance = zeros(total_number_of_bins,1);
for l = 2:total_number_of_bins
    eucdistance(l) = sqrt(sum((temp_string_smooth(l-1,:)-temp_string_smooth(l,:)).^2,2));
end

distot = sum(eucdistance); %total length of string
output_string = temp_string_smooth;


%%********String interpolation step
for iii = 1:(total_number_of_bins-1)
    for jjj = 1:(total_number_of_bins-1)
        jdistbefore = sum(eucdistance(1:jjj));
        jdistafter = sum(eucdistance(1:(jjj+1)));
        requirement = (iii-1)/(total_number_of_bins - 1)*distot;
        if (jdistbefore < requirement) && (requirement <jdistafter)
            reparam = (iii-1)/(total_number_of_bins - 1)*distot - jdistbefore;
            output_string(iii,:) = temp_string_smooth(jjj,:) + (temp_string_smooth(jjj+1,:) - temp_string_smooth(jjj,:)).*reparam/eucdistance(jjj+1);
            break
        end
    end
end



end


