%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjmatpath = path to subject adjacency matrix    %
% commPath = path to commmunity affiliation vector %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define community affiliation vector
Ci=load(commPath)

% Load adjacency matrix
load(adjmatpath)
A = connectivity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Avg. Within and Between-module connectivity %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

withinBetween_mat=zeros(length(Ci));
numComms=numel(unique(Ci));
withinBetween_mat=~bsxfun(@eq,Ci,Ci');
withinBetween_edge_idx=squareform(withinBetween_mat)';

% Within-Module Connectivity
within_idx=find(withinBetween_edge_idx==0);
A=squareform(A);
Avg_withinModule_strength=mean(A(within_idx));

% Between-Module Connectivity
between_idx=find(withinBetween_edge_idx==1);
Avg_betweenModule_strength=mean(A(between_idx));
