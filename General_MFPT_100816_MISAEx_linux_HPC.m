function General_MFPT_092416_MISAEx_linux(loops2,temp_save_location)
loops2 = str2double(loops2);
tic = cputime;
save_location = ['MISAEX_MFPT_100816_circle_HPC'];
temp_dir_transition = ['Tmat_MISAEX_MFPT_100816_circle_HPC/'];

if not(exist(temp_dir_transition))
    mkdir(temp_dir_transition);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************************************************

%***************************************************************
ReplicasRequired = 200;
timestep = 250;
number_sim_steps = 2;
num_nodes = 300;
tau = timestep*number_sim_steps;
species = 8;
species_rad = 8;

matobj_save_info = matfile(save_location,'Writable',true); %saving the updated weights of the regions


load('voronoi_MISA_Ex_100716_tau_500_bin_300.mat');
BNG_root = '/dfs2/elread/rxn-share/BioNetGen-2.2.6-stable/bin/run_network';

BioNetGen_replica_folder = [temp_save_location];

iteration = 1;



poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

if poolsize == 0
    parpool('local')
end
% gives the most recent time a replica has visited a bin, starting from t = 1. Bins with time = 0 have not yet been visited.

infiletraj = 'tVoronoi_MISA_Ex_100716_tau_500_bin_300';
untar([infiletraj '.tar.gz'],[temp_save_location '/']);
previous_BNG_dir = [temp_save_location '/' infiletraj];
load('vc_Ex_MISA_081016_tau_500_bin_300.mat');

binA = [18,4,0,1,0,0,0,1];
binB = [4,18,0,0,1,0,1,0];
radius = 3;

    distsfrombasinB = pdist2(binB,replicas_forwards(:,1:species_rad));
    distsfrombasinA = pdist2(binA,replicas_forwards(:,1:species_rad));
    b=find(distsfrombasinB<=radius);
    a=find(distsfrombasinA<=radius);
    
    
    temp_switch = zeros(length(replicas_forwards(:,1)),1);
    temp_switch(a') = 1;
    temp_switch(b') = -1;


    replicas_forwards(a',species+1) = 1;
    replicas_forwards(b',species+1) = -1;
replicas_forwards = [replicas_forwards(:,1:species),temp_switch,zeros(length(temp_switch),1),replicas_forwards(:,species+1:end)];
replicas_forwards(:,end) = replicas_forwards(:,end)./sum(replicas_forwards(:,end));
    newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
    mkdir(newBNGdir);
    iteration = iteration + 1;


    rep_time = replicas_forwards(:,species+2)+1;
    weights = replicas_forwards(:,end);
    parallel_BNG_name = replicas_forwards(:,species+3);
    binsr = replicas_forwards(:,species+1);
    names_sim = dir([previous_BNG_dir '/*end.net']);
    replicas_new = [];
    parfor kk = 1:length(names_sim)
        unix([BNG_root ' -o ' newBNGdir '/newsim_' num2str(kk) ' -e -p ssa -h $RANDOM --cdat 0 --fdat 0 -g ' previous_BNG_dir '/' names_sim(kk).name ' ' previous_BNG_dir '/' names_sim(kk).name ' ' num2str(timestep) ' ' num2str(number_sim_steps) ' >> bng.log']);
        current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
        
        data_keep = [current_replica_location(end,:)];
        
        replicas_new = vertcat(replicas_new,[data_keep,binsr(kk),rep_time(kk),kk,weights(kk)]);
    end

size(replicas_new)
replicas_forwards = replicas_new;
    
        vordisstB = pdist2(binB,VoronoiLocs(:,1:species_rad));
        vorB = find(vordisstB<=radius);
bin_div = num_nodes - numel(vorB);

    distsfrombasinB = pdist2(binB,replicas_forwards(:,1:species_rad));
    distsfrombasinA = pdist2(binA,replicas_forwards(:,1:species_rad));
    b=find(distsfrombasinB<=radius);
    a=find(distsfrombasinA<=radius);
    
    
    temp_switch = zeros(length(replicas_forwards(:,1)),1);
    temp_switch(a') = 1;
    temp_switch(b') = -1;


    replicas_forwards(a',species+1) = 1;
    replicas_forwards(b',species+1) = -1;
replicas_forwards(:,end) = 1/length(replicas_forwards(:,end));


Binsprev = knnsearch(VoronoiLocs,replicas_forwards(:,1:species));
temp_weights = replicas_forwards(:,end);
weight = zeros(1,num_nodes);
for i = 1:num_nodes
    weight(i) = sum(temp_weights(Binsprev == i));
end
snumber = 4;
for i = 1:snumber
    MFPT(i).counts = zeros(loops2,1);
    MFPT(i).prob = zeros(loops2,1);
    MFPT(i).prob_avg= zeros(loops2,1);

    MFPT(i).times = [];
    MFPT(i).weights = zeros(loops2,1);
    MFPT(i).replica_number = zeros(loops2,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    currentMFPTprobalt = zeros(loops2,snumber);
    currentMFPTprobraw = zeros(loops2,snumber);
    currentMFPTcounts = zeros(loops2,snumber);
    currentMFPTavgs = zeros(loops2,snumber);
runavg  = 1;
for ijk = 1:loops2
    
    

    
    newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
    mkdir(newBNGdir);
    iteration = iteration + 1;
    
    
    
    weights = replicas_forwards(:,end);
    parallel_BNG_name = replicas_forwards(:,species+3);
    timeind = replicas_forwards(:,species+2)+1;
    binsr = replicas_forwards(:,species+1);
    replicas_new = [];
    parfor kk = 1:length(weights)
        [~,~]=unix([BNG_root ' -o ' newBNGdir '/newsim_' num2str(kk) ' -e -p ssa -h $RANDOM --cdat 0 --fdat 0 -g ' ...
            previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
            previous_BNG_dir '/newsim_' num2str(parallel_BNG_name(kk)) '_end.net ' ...
            num2str(timestep) ' ' num2str(number_sim_steps) ' >> ' BioNetGen_replica_folder '/bng.log']);
        current_replica_location = dlmread([newBNGdir '/newsim_' num2str(kk) '.gdat'],'',2,1);
        
        data_keep = [current_replica_location(end,1:species)];
        
        replicas_new = vertcat(replicas_new,[data_keep,binsr(kk),timeind(kk),kk,weights(kk)]);
    end
    
    
    if ijk > 1
        rmdir(previous_BNG_dir,'s');
    end
    previous_BNG_dir = newBNGdir;
    pause(1);

    unix(['rm -rf ' previous_BNG_dir '/*.cdat']);

    
    distsfrombasinB = pdist2(binB,replicas_new(:,1:species_rad));
    distsfrombasinA = pdist2(binA,replicas_new(:,1:species_rad));
b=find(distsfrombasinB<=radius);
a=find(distsfrombasinA<=radius);

    temp_switch = zeros(length(replicas_new(:,1)),1);
    temp_switch(a') = 1;
    temp_switch(b') = -1;
    switchA = find(temp_switch-binsr==2);
    switchB = find(temp_switch-binsr==-2);
    
    
    replicas_new(a',species+1) = 1;
    replicas_new(b',species+1) = -1;

   
    MFPT(3).weights(ijk,1) = sum(replicas_new(replicas_new(:,species+1)==1,end));
    MFPT(4).weights(ijk,1) = sum(replicas_new(replicas_new(:,species+1)==-1,end));
    MFPT(3).replica_number(ijk,1) = numel(find(replicas_new(:,species+1)==1));
    MFPT(4).replica_number(ijk,1) = numel(find(replicas_new(:,species+1)==-1));
    
            [BinsprevA] = knnsearch(VoronoiLocs,replicas_forwards(switchA,1:species));
            fromA = numel(unique(BinsprevA));
            [BinsprevB] = knnsearch(VoronoiLocs,replicas_forwards(switchB,1:species));
            fromB = numel(unique(BinsprevB));
            
            dlmwrite([temp_dir_transition 'repsnew_' num2str(ijk) '.txt'],replicas_new);
            
            
            [sa,~] = size(switchA);
            [sb,~] = size(switchB);
            
            prob_rep = replicas_new(switchA,end);
            MFPT(3).counts(ijk,1) = MFPT(3).counts(ijk,1)+sa;
            MFPT(3).prob(ijk,1) = MFPT(3).prob(ijk,1)+sum(prob_rep);
            MFPT(3).prob_avg(ijk,1) =  MFPT(3).prob(ijk,1)+sum(prob_rep)/fromA;
            
            prob_rep = replicas_new(switchB,end);
            MFPT(4).counts(ijk,1) = MFPT(4).counts(ijk,1)+sb;
            MFPT(4).prob(ijk,1) = MFPT(4).prob(ijk,1)+sum(prob_rep);
            MFPT(4).prob_avg(ijk,1) =  MFPT(4).prob(ijk,1)+sum(prob_rep)/fromB;
            
        
    
    
    replicas_forwards = WEstep_072116(replicas_new,VoronoiLocs,ReplicasRequired,species);
    
    
    [Binsprev] = knnsearch(VoronoiLocs,replicas_forwards(:,1:species));

    
    
    
    
    temp_weights = replicas_forwards(:,end);
    weight = zeros(1,num_nodes);
    for i = 1:num_nodes
        weight(i) = sum(temp_weights(Binsprev == i));
    end
    
    matobj_save_info.replicas_forwards = replicas_forwards;
    matobj_save_info.VoronoiLocs = VoronoiLocs;
    matobj_save_info.MFPT = MFPT ;
    matobj_save_info.bin_div= bin_div;

            if rem(ijk,10)==0
            save([save_location ],'VoronoiLocs','replicas_forwards','MFPT','bin_div','-v7.3')
            
        end
            
            for i = 1:snumber
                
                
                
                currentMFPTprobraw(ijk,i) = 1/(MFPT(i).prob(ijk,1)./(MFPT(i).weights(ijk,1))./(tau));
                currentMFPTcounts(ijk,i) = 1/(MFPT(i).counts(ijk,1)./(tau));
                
                currentMFPTavgs(ijk,i) = 1/(MFPT(i).prob_avg(ijk,1)./(tau));
                
            end
        
        MFa = find(currentMFPTprobraw(:,3));
        MFb = find(currentMFPTprobraw(:,4));
        if ijk > 200
            runavg = runavg+1;
        end

            display([currentMFPTprobraw(ijk,3:4);
                currentMFPTavgs(ijk,3:4);
                [mean(currentMFPTprobraw(MFa(runavg:end),3)),mean(currentMFPTprobraw(MFb(runavg:end),4))]])
            %
            %     display(ijk);
            
            if ijk == 1
                
                dlmwrite([temp_dir_transition 'MFPT_prob.txt'],currentMFPTprobraw(ijk,:) );               
                dlmwrite([temp_dir_transition 'MFPT_counts.txt'],currentMFPTcounts(ijk,:));
                dlmwrite([temp_dir_transition 'MFPT_prob_avg.txt'],currentMFPTavgs(ijk,:));
                
                dlmwrite([temp_dir_transition 'weights.txt'],weight);
            else
                dlmwrite([temp_dir_transition 'MFPT_prob.txt'],currentMFPTprobraw(ijk,:) ,'-append');
                dlmwrite([temp_dir_transition 'MFPT_counts.txt'],currentMFPTcounts(ijk,:),'-append');
                dlmwrite([temp_dir_transition 'MFPT_prob_avg.txt'],currentMFPTavgs(ijk,:),'-append');
                dlmwrite([temp_dir_transition 'weights.txt'],weight,'-append');
            end
            
        
    
    
    
    
    
    
    
end




save([save_location '_final'],'VoronoiLocs','replicas_forwards','MFPT','bin_div' ,'-v7.3')

toc = cputime-tic
end


function [newreps] = WEstep_072116(replicas,Voronoi_List,numreps,species)
newreps = [];
bins = knnsearch(Voronoi_List,replicas(:,1:species));
%This code assumes the columns of the replicas matrix are the species, and that the weights are in the last column
parfor i = 1:length(Voronoi_List(:,1))
    replicas_in_bini = replicas(bins==i,:);
    if not(isempty(replicas_in_bini))
        newreps_temp = WEstep(replicas_in_bini,numreps);
        if  abs(sum(replicas_in_bini(:,end)) -   sum(newreps_temp(:,end)) > 10^(-12))
            sum(newreps_temp(:,end))
            sum(replicas_in_bini(:,end))
            error('WEstep failure');
        end
        newreps = [newreps; newreps_temp];
        
    end
end


end



function newreps = WEstep(curr_reps,numreps)  %%%%%error in weight calc

blah = length(curr_reps(:,1));
if isempty(blah)
    error('reps_empty!')
end
%%%%%number of repilcas incorred incremented by 1
newreps = curr_reps;

curr_weights = newreps(:,end);
total_weight = sum(curr_weights);
idx_remove = [];
[SortedWeights, SortedInd] = sort(curr_weights);
min_idx = 2;
max_idx = length(SortedInd);
smallest_weight = SortedWeights(1);

%%removes smallest weight outlier
while (smallest_weight < total_weight/(3*numreps)) && (min_idx < max_idx)
    smallest_weight = sum(SortedWeights(1:min_idx));
    idx_remove = [ SortedInd(1:min_idx)];
    
    min_idx = min_idx+1;
end
if isempty(idx_remove) == 0
    new_rep_ind = randsample(SortedInd(1:min_idx),1,true,SortedWeights(1:min_idx));
    
    new_reps_temp = newreps(new_rep_ind,:);
    new_reps_temp(:,end) = smallest_weight;
    newreps = [newreps; new_reps_temp];
    newreps(idx_remove,:) = [];
end

w2 = sum(newreps(:,end));

curr_weights =  newreps(:,end);
%%removes largest weight outlier
[largest_weight, max_idx] = max(curr_weights);
tmpweight = largest_weight;
mxx = 0;
max_factor = 1;

if largest_weight > total_weight/numreps*3
    
    while largest_weight > total_weight/numreps*3
        largest_weight = tmpweight/max_factor;
        max_factor = max_factor+1;
    end
    
    new_reps_temp = repmat(newreps(max_idx,:),max_factor+mxx,1);
    new_reps_temp(:,end) = repmat(tmpweight/(max_factor+mxx),max_factor+mxx,1);
    newreps = [newreps; new_reps_temp];
    newreps(max_idx,:) = [];
end


curr_weights = newreps(:,end);
curr_len = length(curr_weights);
w3 = sum(newreps(:,end));

%%Merges weights until the number of replicas = numreps

[SortedWeights, SortedInd] = sort(curr_weights);
min_idx = 2;
max_idx = length(SortedInd);
smallest_weight = SortedWeights(1);
idx_remove = [];

while curr_len  > numreps && (min_idx < max_idx)
    smallest_weight = sum(SortedWeights(1:min_idx));
    idx_remove = [ SortedInd(1:min_idx)];
    min_idx = min_idx+1;
    curr_len = curr_len-1;
end
if isempty(idx_remove) == 0
    new_rep_ind = randsample(SortedInd(1:min_idx),1,true,SortedWeights(1:min_idx));
    
    new_reps_temp = newreps(new_rep_ind,:);
    new_reps_temp(:,end) = smallest_weight;
    newreps = [newreps; new_reps_temp];
    newreps(sort(idx_remove),:) = [];
    
end



w4 = sum(newreps(:,end));

%%splits weights until the number of replicas = numreps

curr_weights =  newreps(:,end);
curr_len = length(curr_weights);

[largest_weight, max_idx] = max(curr_weights);
tmpweight = largest_weight;
mxx_fact2 = 0;

max_factor = 1;


if curr_len < numreps
    while curr_len < numreps
        max_factor = max_factor+1;
        curr_len = curr_len+1;
    end
    new_reps_temp = repmat(newreps(max_idx,:),max_factor+mxx,1);
    a = length(new_reps_temp(:,1));
    
    new_reps_temp(:,end) = tmpweight/(max_factor+mxx_fact2).*ones(a,1);
    newreps = [newreps; new_reps_temp];
    shouldbezero = (length(newreps(:,end))-1-numreps);
    
    
    newreps(max_idx,:) = [];
    
    
    shouldbezero2 = (length(newreps(:,end))-numreps);
    
    
    if shouldbezero2~=0
        error(' replica number from WEstep after splitting is wrong');
    end
    
    
end

w5 = sum(newreps(:,end));


%%Checking for errors in the weight calculation.

if (all(abs([total_weight-w2,w2-w3,w3-w4,w4-w5]))) > 10^(-12)
    abs([total_weight-w2,w2-w3,w3-w4,w4-w5])
    error('Error is in one of the three subsections of WE step')
elseif length(newreps(:,end)) ~= numreps
    error('Final replica length is wrong');
end


end






