function binconn = binarize_conn(conn, th)
% 
% This function binarizes a weighted connectivity matrix by first selecting the top th% strongest connections 
% (based on absolute value), and then binarizes those that are positive.
% 
% Inputs :  conn    - Weighted connectivity matrix  (n x n)
%           th      - Percentage threshold for binarization (e.g., 5 for top 5%)
% Outputs:  binconn - Binarized connectivity matrix (n x n)
% 

Nroi = size(conn, 1);
binconn = zeros(Nroi, Nroi);

W = conn;

% Check if the matrix is symmetric
if isequal(W, W.')	% if symmetric matrix
    W = triu(W);      % ensure symmetry is preserved
    ud = 2;           % halve number of removed links
else
    ud = 1;
end

% Find all links and sort by absolute magnitude
ind = find(W);                            % find all links
E = sortrows([ind abs(W(ind))], -2);      % sort by absolute magnitude

% Number of links to be preserved
en = round((Nroi^2 - Nroi) * 0.01 * th / ud);

% Apply threshold: Set the top th% to their original values, others to 0
W_temp = zeros(size(W));
W_temp(E(1:en, 1)) = W(E(1:en, 1)); 

% Keep only positive connections among the top th%
W = W_temp .* (W_temp > 0);

% If symmetric, reconstruct symmetry
if ud == 2
    W = W + W.';
end

% Binarize: Set the positive connections to 1
W = double(W > 0);
binconn = W;
end
