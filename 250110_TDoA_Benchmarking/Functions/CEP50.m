function cep = CEP50(C)
    % Computes the radius for a CEP_50 circle from a given error covariance
    % matrix C.  The CEP_50 circle is a circle that contains half of the random
    % samples defined by the error covariance matrix.
    % If the ratio is less than 2, meaning that both eigenvectors contribute 
    % roughly the same amount of error, then we apply the approximation:
    %    cep = .59*(sqrt(lamMin)+sqrt(lamMax));
    % otherwise, the dominant eigenvector is responsible for the majority of
    % the error, and we apply the approximation:
    % cep = sqrt(lamMax)*(.67+.8*lamMin/lamMax)
    
    % INPUTS: C: Position Error Covariance matrix (computed as:
    % posCov = (posErrCovMtx(:,:) + ((emitterPos(1:2,:) - pos_est(1:2,mc))*(emitterPos(1:2,:) - pos_est(1:2,mc)).'))/MC;


    [~,Lam] = eig(C);
    lam = diag(Lam); % Pull eigenvalue vector from diagonal matrix Lam
    
    % Sort the eigenvalues
    [lamSort,~] = sort(lam,'descend'); % Sort the eigenvalues
    
    % Dominant eigenvalue
    lamMax = abs(lamSort(1));
        
    % Secondary eigenvalue
    lamMin = abs(lamSort(2));
        
    % Ratio of dominant to secondary eigenvalues
    ratio = sqrt(lamMin/lamMax);
    
    % Depending on the eigenvalue ratio, use the appropriate approximation
    if ratio > .5
        cep = .59*(sqrt(lamMin)+sqrt(lamMax));
    else
        cep = sqrt(lamMax)*(.67+.8*lamMin/lamMax);
    end 
end