function ErrCov = fn_restructCovMtx(sensVec,refSen,timingErr)

% Abeer Chaudhry ©
% January 2025
    c = physconst('LightSpeed');

    covIn = timingErr^2.*eye(length(sensVec) + 1).*c^2; % 4x4 Error Covariance matrix in range
    %% The objective of this function is to compute the covariance 
    n_sens = length(sensVec);
    nRef = repmat(refSen,1,n_sens);
    % Precompute Row and Column Indices
    row_idx = mod((1:n_sens) - 1, n_sens) + 1; % Row indices for sensors
    col_idx = mod((1:n_sens) - 1, n_sens) + 1; % Row indices for reference sensor

    % Parse sensor indices for rows and columns
    a_Rowidx = sensVec(row_idx); % Row indices for sensors
    b_Rowidx = nRef(row_idx);  % Row indices for reference sensor
    
    a_Col_idx = sensVec(col_idx); % Column indices for sensors 
    b_Col_idx = nRef(col_idx);  % Column indices for reference sensor 
    
    
    for i = 1:n_sens
        for j = 1:n_sens
            cov_bibj = covIn(b_Rowidx(i),b_Col_idx(j));
            cov_aiaj = covIn(a_Rowidx(i),a_Col_idx(j));
            cov_aibj = covIn(a_Rowidx(i),b_Col_idx(j));
            cov_biaj = covIn(b_Rowidx(i),a_Col_idx(j));
    
            ErrCov(i,j) = cov_bibj + cov_aiaj - cov_aibj - cov_biaj;
        end
    end
end