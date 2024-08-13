function cep =  computeCEP50(C)
% cep = computeCEP50(C)
%
% Computes the radius for a CEP_50 circle from a given error covariance
% matrix C.  The CEP_50 circle is a circle that contains half of the random
% samples defined by the error covariance matrix.
%
% Calculation is extremely complex, and requires numerical integration, so
% the equation used herein is an approximation, depending on the ratio of
% the dominant to secondary eigenvalues.  If the ratio is less than 2,
% meaning that both eigenvectors contribute roughly the same amount of
% error, then we apply the approximation:
%    cep = .59*(sqrt(lamMin)+sqrt(lamMax));
% otherwise, the dominant eigenvector is responsible for the majority of
% the error, and we apply the approximation:
%    cep = sqrt(lamMax)*(.67+.8*lamMin/lamMax);    
%
% Inputs:
%
%   C       2x2 error covariance matrix (additional dimensions are assumed
%           to correspond to independent cases, and are computed in turn)
%
% Outputs:
%
%   cep     Radius of the corresponding CEP_50 circle

% Check for multiple entries
if sum(size(C)>1)>2
    % Multiple entries given
    fullDims=size(C);

    % Determine the output dimension and number of test cases input
    outDims = fullDims(3:end);
    nIters = prod(outDims);

    % Reshape the covariance matrix input
    C = reshape(C,[fullDims(1:2),nIters]);
    cep = zeros(outDims,1);
    
    % Call the function on each individual iteration
    for ii=1:nIters
        cep(ii) = utils.computeCEP50(squeeze(C(:,:,ii)));
    end
    return;
end

if any(~isfinite(C))
    cep = Inf;
    return;
end

% Eigenvector analysis to identify independent components of error
[~,Lam] = eig(C);
lam = diag(Lam); % Pull eigenvalue vector from diagonal matrix Lam

% Sort the eigenvalues
[lamSort,~] = sort(lam,'descend'); % Sort the eigenvalues

% Dominant eigenvalue
lamMax = abs(lamSort(1));
%vMax = V(:,iSort(1)); 

% Secondary eigenvalue
lamMin = abs(lamSort(2));
%vMin = V(:,iSort(2));

% Ratio of dominant to secondary eigenvalues
ratio = sqrt(lamMin/lamMax);

% Depending on the eigenvalue ratio, use the appropriate approximation
if ratio > .5
    cep = .59*(sqrt(lamMin)+sqrt(lamMax));
else
    cep = sqrt(lamMax)*(.67+.8*lamMin/lamMax);
end 