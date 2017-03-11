function withinBetween_module_connectivity(adjmatpath,commPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inputs:                                                %
%        adjmatpath = path to subject adjacency matrix    %
%        commPath = path to commmunity affiliation vector %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define community affiliation vector
input_commAff=load(commPath)

% Load adjacency matrix
load(adjmatpath)
A = connectivity;

% Define Modules and Nodes in network
unique_S=unique(input_commAff);
numNodes=length(A)

% Number of communities 
numComm=length(unique_S);

% Set diagonal of adjacency matrix to nan
A=A + diag(repmat(nan,[numNodes,1]));

% Define community by community matrix
comm_comm_mat=zeros(numComm,numComm);

% Define Within/Between Module connectivity matrix
comm_wb_mat=zeros(numComm,2);
wb_vec=zeros(1,2);
com1 = 1
for i=unique_S'
	com2 = 1;
	% Define index for nodes in each community
	comidx = find(input_commAff==i);
	not_comidx=find(input_commAff~=i)
	for j = unique_S'
		comidx_2= find(input_commAff==j);
		% Get mean edge weights of edges connecting each pair of communities
		% Pair-wise Between-module connectivity
		current_edges=A(comidx,comidx_2);
		mean_edgeWeight=nanmean(nanmean(current_edges));
		% Define a community X community matrix for each pair of communities
		comm_comm_mat(com1,com2)=mean_edgeWeight;
		com2= com2 + 1;
	end

	% Within module connectivity
	comm_wb_mat(i,1) = nanmean(nanmean(A(comidx,comidx)));
	% Between module connectivity
	comm_wb_mat(i,2) = nanmean(nanmean(A(comidx,not_comidx)));

	com1 = com1 + 1;

end

% Compute the overall average within- and between-module connectivity
within = logical(bsxfun(@eq,input_commAff,input_commAff'));
wb_vec(1) = nanmean(A(within));
wb_vec(2) = nanmean(A(~within));

within_between_ratio = wb_vec(1) / wb_vec(2)

% Average Within-Module Connectivity
Avg_Within_Conn=wb_vec(1)
% Average Between-Module Connectivity
Avg_Between_Conn=wb_vec(2)
