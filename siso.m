% Calculate the theoretical lower bounds of the decay rate of the SIS model
%
% When using this code, please cite the following paper:
%
% Naoki Masuda, Victor M. Preciado, Masaki Ogura.
% Analysis of the susceptible-infected-susceptible epidemic dynamics in networks via the non-backtracking matrix.
% IMA Journal of Applied Mathematics, 85, 214-230 (2020).
%
% Input: edge list file, which should be in the following format:
% The first row of the input file should be "# number-of-nodes number-of-edges"
% Starting from the second row, the first two columns have two nodes to form an edge, and the third column has the edge weight.
% However, we ignore the edge weight (i.e., third column, starting from the second row).
% The node index starts from 1, not from 0. 
% The node index should be consecutive.
% If users do not like these rules on the input file, it should be easy to rewrite the read-file part of the following code.
%
% By changing the value of "data" below, one can specify which input
% network to be used.
%
% Output: bound-[input_file_name_without_extension].txt
% This is a table whose each row corresponds to an infection rate value.
% 1st column: beta (i.e., infection rate)
% 2nd column: meanfield lower bound of the decay rate
% 3rd column: lower bound of the decay rate by Ogura & Preciado (2018)
% 4th column: lower bound of the decay rate by Masuda, Preciado, & Ogura
% (2020)

% cd '/homedirectory'

data = 3; % specify the input network data
beta_samples = 30; % number of beta values for which the decay rate is computed
 
if data==1 % RRG
    edge_list = load('rrg-n100k6-decayrate.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 0.28;
elseif data==2 % BA
    edge_list = load('ba-n100m3-decayrate.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 0.2;
elseif data==3 % cycle
    edge_list = load('cycle-n100.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 1.6;
elseif data==4 % LFR
    edge_list = load('lfr-n100k6-decayrate.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 0.22;
elseif data==5 % dolphin
    edge_list = load('dolphins.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 0.35;
elseif data==6 % LCC of the network of network scientists
    edge_list = load('netscience-lcc.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 0.24;
elseif data==7 % email network
    edge_list = load('email.mat', '-ASCII');
    beta_min = 0.0;
    beta_max = 0.066;
%%% data==8 % hamsterster is somehow missing
end

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

N = max(max(edge_list(:,1:2))); % number of nodes
M = size(edge_list,1) % number of edges

if N>390
    if_skip2018 = 1; % Skip the computation of the lower bound of
    % the decay rate by Ogura and Preciado (2018) because it is
    % computationally costly for large N.
else
    if_skip2018 = 0; % no skip
end

% Input is assumed to be an undirected network

% edge_list is assumed to be an edge list for an undirected network
adj = zeros(N, N); % adjacency matrix
for i=1:M
    adj(edge_list(i,1), edge_list(i,2)) = 1;
    adj(edge_list(i,2), edge_list(i,1)) = 1;
end

beta = 0.02; % infection rate
result = []; % store decay rate for each infection rate and for each lower bound

for ind_beta = 1:beta_samples
    if beta_samples >= 2
        beta = beta_min + (ind_beta - 1) * (beta_max - beta_min) / (beta_samples - 1);
    end
    decay_rate_meanfield = bound_mf(adj,beta);
    decay_rate_nbt = bound_nbt(edge_list,beta);
    if if_skip2018==0 
%        decay_rate_Ogura2018 = bound_Ogura2018(adj,beta);
        decay_rate_Ogura2018_sparse = bound_Ogura2018_sparse(adj,beta);
        result = [result ; beta decay_rate_meanfield decay_rate_Ogura2018_sparse decay_rate_nbt];
    else
        disp([beta decay_rate_meanfield decay_rate_nbt])
        result = [result ; beta decay_rate_meanfield 0 decay_rate_nbt];
    end
end % all infection-rate values done

disp(result)
if data==1
    save('bound-rrg-n100k6.txt','result','-ascii');
elseif data==2
    save('bound-ba-n100m3.txt','result','-ascii');
elseif data==3
    save('bound-cycle-n100.txt','result','-ascii');
elseif data==4
    save('bound-lfr-n100k6.txt','result','-ascii');
elseif data==5
    save('bound-dolphin.txt','result','-ascii');
elseif data==6
    save('bound-netscience-lcc.txt','result','-ascii');
elseif data==7
    save('bound-email.txt','result','-ascii');
end  

% A new lower bound of the decay rate using the non-backtracking matrix,
% derived in Masuda, Preciado & Ogura, IMA Journal of Applied Mathematics,
% 85, 214-230 (2020)
function decay_rate = bound_nbt(edge_list, beta)
% beta: infection rate
% recovery rate is normalized to 1

    N = max(max(edge_list(:,1:2))); % number of nodes
    M = size(edge_list,1); % number of edges
    
    % incidence matrix, C = C_+ - C_-
    Cp = zeros(N,2*M); % C_+, to
    Cm = zeros(N,2*M); % C_-, from
    for i=1:M
        Cm(edge_list(i,1),2*i-1) = 1; % edge_list(i,1) -> edge_list(i,2)
        Cp(edge_list(i,2),2*i-1) = 1;
        Cm(edge_list(i,2),2*i) = 1; % edge_list(i,2) -> edge_list(i,1)
        Cp(edge_list(i,1),2*i) = 1;
    end
%    if (if_connected(Cp+Cm) == false)
%        error('Network is not a connected network');
%    end

    % non-backtracking matrix
    H = zeros(2*M,2*M);
    for i=1:2*M
        for j=[1:i-1 i+1:2*M]
            tmp_start1 = Cm(:,i);
            tmp_end1 = Cp(:,i);
            tmp_start2 = Cm(:,j);
            tmp_end2 = Cp(:,j);
            if (find(tmp_end1) == find(tmp_start2)) & (find(tmp_end2) ~= find(tmp_start1))
                H(i,j) = 1;
            end
        end
    end

    A_nbt = [-eye(N) beta*Cp ; Cm' beta*H' - (2+beta)*eye(2*M)];
    decay_rate = - eigs(A_nbt, 1, 'largestreal');

end

% Mean-field lower bound of the decay rate
function decay_rate = bound_mf(adj, beta)
% adj: adjacency matrix
% beta: infection rate
% recovery rate is normalized to 1

    N = size(adj,1); % number of nodes
    decay_rate = - eigs(beta * adj' - eye(N), 1, 'largestreal');
end

% Lower bound of the decay rate derived in Ogura and Preciado, Systems &
% Control Letters, 113, 59-64 (2018)
function decay_rate = bound_Ogura2018(adj, beta)
% adj: adjacency matrix
% beta: infection rate
% recovery rate is normalized to 1

    N = size(adj,1); % number of nodes
    
    A11 = -eye(N); % A = [A11, A12; A21, A22] is a block matrix

    A12 = [];
    for i=1:N
        A12 = blkdiag(A12, adj(i,[1:i-1 i+1:N]));
    end
    A12 = beta * A12;

    A21 = [];
    for i=1:N
        tmp_mat = eye(N);
        A21 = [A21 ; tmp_mat([1:i-1 i+1:N],:)];
    end

    A22 = [];
    for i=1:N
        Gamma_i = [];
        tmp_mat = [];
        for j= [1:i-1 i+1:N]
            Gamma_i = blkdiag(Gamma_i, 2 + adj(i,j)*beta);
            tmp_mat = [tmp_mat ; beta*adj(j,[1:i-1 i+1:N])];
        end
        A22 = blkdiag(A22, - Gamma_i + tmp_mat);
    end

    A_ogura2018 = [A11 A12; A21 A22];
    decay_rate = - eigs(A_ogura2018, 1, 'largestreal');

end

% Lower bound of the decay rate derived in Ogura and Preciado, Systems &
% Control Letters, 113, 59-64 (2018). Sparse matrix version
function decay_rate = bound_Ogura2018_sparse(adj, beta)
% adj: adjacency matrix
% beta: infection rate
% recovery rate is normalized to 1

    N = size(adj,1); % number of nodes
    
    A11 = sparse(-eye(N));

    A12 = sparse([]);
    for i=1:N
        A12 = blkdiag(A12, sparse(adj(i,[1:i-1 i+1:N])));
    end
    A12 = beta * A12;

    A21 = sparse([]);
    for i=1:N
        tmp_mat = eye(N);
        A21 = [A21 ; sparse(tmp_mat([1:i-1 i+1:N],:))];
    end

    A22 = sparse([]);
    for i=1:N
        Gamma_i = [];
        tmp_mat = sparse([]);
        for j= [1:i-1 i+1:N]
            Gamma_i = blkdiag(Gamma_i, 2 + adj(i,j)*beta);
            tmp_mat = [tmp_mat ; sparse(beta*adj(j,[1:i-1 i+1:N]))];
        end
        A22 = blkdiag(A22, sparse(- Gamma_i) + tmp_mat);
    end

    A_ogura2018 = [A11 A12; A21 A22]; % sparse matrix
    decay_rate = - eigs(A_ogura2018, 1, 'largestreal');
end

function out = if_connected(C)
%
% Check if the network is a connected network.
%
% C: incidence matrix
%
% out = 1 if connected and = 0 otherwise
%
    L = C * C'; % Laplacian matrix
    eig_ascending = sort(eig(L));
    if eig_ascending(2) > 1e-8
        out = true;
    else
        out = false;
    end
end