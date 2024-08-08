function [dop] = GDOP(timing_error, grd_res,x_min,x_max,y_min,y_max)
% This function calculates and plot the geometric dilution of precion for a
% sensor-network geometry defined in the 2D plane.
% The inputs to the function are:
% 1. timing_error: Determine the timing performance of the sensors in order
% to set a covariance matrix to be used for a calculation of the CRLB
% 2. grd_res: is the grid resolution for which the source positions will be
% initialized in order to calculate the GDOP
% 3. x_min, x_max: are the min and max limits of the 2D span along x-axis
% 4. y_min, y_max: are the min and max limits of the 2D span along y-axis


%% Check input arguments and assign default values

if nargin < 1
    timing_error = 30e-9;
    grd_res = 1e2;
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 2
     grd_res = 1e2;
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 3
    
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 4
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 5
    y_min = -200e3;
    y_max = 200e3;    

elseif nargin < 6
   y_max = 200e3;
end

% Define Sensor Positions
figure
axis([x_min/1e3 x_max/1e3 y_min/1e3 y_max/1e3])
grid on
xlabel('Cross-range (km)'); ylabel('Down-range (km)'); title('Create Sensor Network Geometry')
[x,y] = ginput();
close
nSensors = length(x);
sensor_pos = zeros(2,length(x));
sensor_pos(1,:) = x*1e3;
sensor_pos(2,:) = y*1e3;
% Define Sensor Performance
timingError = timing_error;
% Create the error covariance matrix (time of arrival)
Ctoa = timingError^2*eye(nSensors);

% Define source positions
grd_size = grd_res; % Grid resolution
x_ax = linspace(x_min,x_max,grd_size);
y_ax = linspace(y_min,y_max,grd_size);
[xx,yy] = ndgrid(x_ax,y_ax);
source_posgrid = [xx(:) yy(:)]';

% Compute CRLB and CEP50
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = computeCRLB(sensor_pos,source_posgrid,Ctoa); % Ndim x Ndim x M^2
dop = reshape(computeCEP50(crlb),[grd_size,grd_size]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,2,3,5,10,15,20];
contourLevelsLabel = [.1,1,2,3,5,10,15,20];

%% Draw Figure
fig = figure();hold on;

plot(sensor_pos(1,:)/1e3,sensor_pos(2,:)/1e3,'o','DisplayName','Sensors','LineWidth',1,'MarkerFaceColor','red');
[cp,hiso2]=contour(xx/1e3,yy/1e3,dop/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);
excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;

% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
% setPlotStyle(gca,{'equal','tight'});







%% Children Functions
function cov_out = resampleCovMtx(cov, test_idx, ref_idx, test_wts, ref_wts)
% cov_out = resampleCovMtx(cov, test_idx, ref_idx, test_wts, ref_wts)
%
% Resample a covariance matrix based on a set of reference and test 
% indices.  This assumes a linear combination of the test and reference 
% vectors.  The output is an n_pair x n_pair covariance matrix for the 
% n_pair linear combinations.
%
% The measurements can be optionally weighted, to over-emphasize some
% measurements and de-emphasize others.
%
% In the resampled covariance matrix, the i,j-th entry is given
%    [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
%       where:  ai, aj are the i-th and j-th reference indices
%               bi, bj are the i-th and j-th test indices
%               C is the input covariance matrix
%
% Any indices that are NaN will be ignored,s to represent a single-sensor
% measurement, such as AoA (which does not need a reference sensor
% measurement against which to compare), or a noise-free measurement with
% no error.
%
% If the third input, ref_idx, is missing or empty, then the second input,
% test_idx, will be passed to utils.parseReferenceSensor to generate
% matching test and reference vectors.
%
% INPUTS:
%   cov         NxN covariance matrix of individual sensor measurements
%   test_idx    n_pair x 1 vector of test sensor indices
%   ref_idx     n_pair x 1 vector of reference sensor indices [Optional]
%   test_wts    Optional n_pair x 1 vector of test measurement weights
%   ref_wts     Optional n_pair x 1 vector of reference measurement weights
%
% OUTPUTS:
%   cov_out     n_pair x n_pair output covariance matrix of sensor
%               measurement pairs

%% Input handling
% Parse array sizes and indices
n_sensor = size(cov, 1);

% Handle test/reference inputs
if nargin < 3 || isempty(ref_idx)
    [test_idx, ref_idx] = utils.parseReferenceSensor(test_idx, n_sensor);
end

% Parse output size
n_test = numel(test_idx);
n_ref = numel(ref_idx);
n_out = max(n_test, n_ref);

if n_test > 1 && n_ref > 1 && n_test ~= n_ref
    error(strcat("Error calling covariance matrix resample. Reference and", ...
                 " test vectors must have the same shape."))
end

if any(test_idx > n_sensor) || any(ref_idx > n_sensor)
	error(strcat("Error calling covariance matrix resample. Indices exceed", ...
                 " the dimensions of the covariance matrix."))
end

% Parse sensor weights
do_test_wt = ~(nargin < 4 || isempty(test_wts));
do_ref_wt = ~(nargin < 5 || isempty(ref_wts));

if do_test_wt
    n_test_wt = numel(test_wts);
end

if do_ref_wt
    n_ref_wt = numel(ref_wts);
end

% Initialize output
cov_out = zeros(n_out, n_out);

a_i_wt = 1;
a_j_wt = 1;
b_i_wt = 1;
b_j_wt = 1;

% Step through reference sensors
for idx_row = 1:n_out
    % Parse sensor indices.  The mod commands seamlessly handle scalar
    % inputs
    a_i = test_idx(1+mod(idx_row-1, n_test));
    b_i = ref_idx(1+mod(idx_row-1, n_ref));

    % Parse sensor weights
    if do_test_wt
        a_i_wt = test_wts(1+mod(idx_row-1, n_test_wt));
    end

    if do_ref_wt
        b_i_wt = ref_wts(1+mod(idx_row-1, n_ref_wt));
    end

    for idx_col =1:n_out
        % Parse sensor indices.  The mod commands seamlessly handle scalar
        % inputs.
        a_j = test_idx(1+mod(idx_col-1, n_test));
        b_j = ref_idx(1+mod(idx_col-1, n_ref));

        if do_test_wt
            a_j_wt = test_wts(1+mod(idx_col-1, n_test_wt));
        end
    
        if do_ref_wt
            b_j_wt = ref_wts(1+mod(idx_col-1, n_ref_wt));
        end
        
        % Parse Input covariances
        if isnan(b_i) || isnan(b_j)
            cov_bibj = 0;
        else
            cov_bibj = cov(b_i, b_j);
        end
        if isnan(a_i) || isnan(a_j)
            cov_aiaj = 0;
        else
            cov_aiaj = cov(a_i, a_j);
        end
        if isnan(a_i) || isnan(b_j)
            cov_aibj = 0;
        else
            cov_aibj = cov(a_i, b_j);
        end
        if isnan(b_i) || isnan(a_j)
            cov_biaj = 0;
        else
            cov_biaj = cov(b_i, a_j);
        end
        
        %  [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
        cov_out(idx_row, idx_col) = b_i_wt * b_j_wt * cov_bibj + ...
                                    a_i_wt * a_j_wt * cov_aiaj - ...
                                    a_i_wt * b_j_wt * cov_aibj - ...
                                    b_i_wt * a_j_wt * cov_biaj;

    end
end
end

function [test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, num_sensors)
%[test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, num_sensors)
%
% Accepts a reference index setting (either None, a scalar integer, or a 
% 2 x N array of sensor pairs), and returns matching vectors for test 
% and reference indices.
%
% INPUTS:
%   ref_idx         Reference index setting
%   num_sensors     Number of available sensors
%
% OUTPUTS:
%   test_idx_vec    Vector of test sensor indices
%   ref_idx_vec     Vector of reference sensor indices


%% Default Behavior
if isempty(ref_idx)
    % Default behavior is to generate all possible sensor pairs
    test_idx_vec = 1:num_sensors-1;
    ref_idx_vec = num_sensors * ones(size(test_idx_vec));
elseif strcmpi(ref_idx,'full')
    % Do the full set of N(N-1)/2 pairs
    full_set = nchoosek(1:num_sensors, 2);
    test_idx_vec = full_set(:,2)';
    ref_idx_vec = full_set(:,1)';
elseif isscalar(ref_idx)
    % Scalar reference provided, use all others as test sensors
    test_idx_vec = setdiff(1:num_sensors, ref_idx);
    ref_idx_vec = ref_idx * ones(size(test_idx_vec));
else
    % Explicit sensor pairs provided, parse them
    test_idx_vec = ref_idx(1, :);
    ref_idx_vec = ref_idx(2, :);
end
        
if size(test_idx_vec) ~= size(ref_idx_vec)
    warning('utils/parseReferenceSensor.m generated unequal test and reference vectors.  Check for bugs.');
end

return
end


function J = jacobian(x_sensor, x_source, ref_idx)
% J = jacobian(x_sensor, x_source, ref_idx)
%
% Returns the Jacobian matrix for TDOA of a source at x_source 
% (nDim x nSource) from sensors at x_sensor (nDim x nSensor).
%
% INPUTS:
%   x_sensor        nDim x nSensor vector of sensor positions
%   x_source        nDim x nSource vector of source positions
%   ref_idx         Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position

% Parse inputs
[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);

if nDim1 ~= nDim2
    error('Input variables must match along first dimension.');
end
nDim = nDim1;

if nargin < 3 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor;
end

% Parse Reference Sensors
n_sensor = size(x_sensor, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);

% Compute the range from each candidate source location to each sensor
dx = reshape(x_source,nDim,1,nSource) - reshape(x_sensor,nDim,nSensor);
R = sqrt(sum(abs(dx).^2,1)); % 1 x nSensor x nSource

% Remove any zero-range samples; replace with epsilon, to avoid a divide
% by zero error
R(R<1e-10) = 1e-10;

% Compute the Jacobians
J = dx(:,test_idx_vec,:)./R(:,test_idx_vec,:) ...
    - dx(:,ref_idx_vec,:)./R(:,ref_idx_vec,:);

end

function excludeFromLegend(graphObj)
% excludeFromLegend(graphObj)
%
% Modifies the supplied graphics object such that it will not appead in a
% legend; useful for plotting annotations that do not belong in the legend,
% such as dashed lines between a label and the appropriate line or
% marker.
%
% The handles are returned from the plotting command that was used to
% generate them, such as
%   h = plot(...)
% or
%   h = text(...)
%
% Alternatively, an axis object can be searched with the findall or findobj
% commands, such as:
%   h1 = findall(gcf);
%   h2 = findall(gcf,'Type','text');
%   h3 = findobj(gcf,'Type','text');
%
% If multiple graph objects are input the command is applied to each.
%
% Inputs:
%   graphObj        Graphics object handles (scalar or array)


if numel(graphObj)>1
    arrayfun(@(x) utils.excludeFromLegend(x),graphObj);
    return;
end

% Graphics Objects have an Annotation property.
% Within the Annotation object, there is a LegendInformation Property.
% The LegendInformation object has an IconDisplayStyle property.
%
% Set that IconDisplayProperty to 'off' and a graphics object will be
% excluded from legends
graphObj.Annotation.LegendInformation.IconDisplayStyle = 'off';
end


function cov_out = ensureInvertible(cov, epsilon)
% function cov_out = ensureInvertible(cov, epsilon)
%
% Check the input matrix for invertibility by finding the eigenvalues and
% checking that they are all >= a small value (epsilon).
%
% If any of the eigenvalues are too small, then a diagonal loading term
% is applied to ensure that the matrix is positive definite (all
% eigenvalues are >= epsilon).
%
% INPUTS:
%   cov     Input square covariance matrix.  If the input has >2
%           dimensions, then the process is repeated across the extra
%           dimensions.
%
%   epsilon (Optional) Specifies the minimum eigenvalue to use when
%           checking for invertibility. Default = 1e-10
%
% OUTPUTS:
%   cov_out Output covariance matrix, guaranteed to be invertible.

% Check for epsilon input
if nargin < 2 || isempty(epsilon)
    epsilon = 1e-20;
end

% Check input dimensions
sz = size(cov);
assert(numel(sz) > 1, 'Input must have at least two dimensions.');
assert(sz(1) == sz(2), 'First two dimensions of input matrix must be equal.');
dim = sz(1);
if numel(sz) > 2
    n_matrices = prod(sz(3:end));
else
    n_matrices = 1;
end

cov_out = zeros(size(cov));
for idx_matrix = 1:n_matrices
   % Check min eigenvalue
   if min(eig(cov(:,:,idx_matrix))) < epsilon
      % We need to add a diagonal loading term, determine appropriate size
      d = epsilon;
      while min(eig(cov(:,:,idx_matrix) + d * eye(dim))) < epsilon
          d = d * 10;
      end
      
      % Apply diagonal loading
      cov_out(:,:,idx_matrix) = cov(:,:,idx_matrix) + d * eye(dim);
   else
       % No diagonal loading
       cov_out(:,:,idx_matrix) = cov(:,:,idx_matrix);
   end
end
end

function crlb = computeCRLB(x_tdoa,xs,C,ref_idx,variance_is_toa,resample_covariance)
% crlb = computeCRLB(x_tdoa,xs,C,ref_idx,variance_is_toa,resample_covariance)
%
% Computes the CRLB on position accuracy for source at location xs and
% sensors at locations in x1 (Ndim x N).  Ctdoa is an Nx1 vector of TOA
% variances at each of the N sensors.
%
% Inputs:
%   x_tdoa      (Ndim x N) array of TDOA sensor positions
%   xs          (Ndim x M) array of source positions over which to 
%               calculate CRLB
%   C           TOA covariance matrix [s^2]
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   variance_is_toa (Optional) flag indicating whether supplied variance is
%               in units of time (TRUE) or distance (FALSE). Default=TRUE.
%   resample_covariance (Optional) flag indicating whether the covariance 
%               matrix should be resamples (TRUE) to convert from sensor
%               errors to measurement errors, or whether it was directly
%               supplied as measurement errors (FALSE). Default=TRUE.
%
% Outputs:
%   crlb    Lower bound on the error covariance matrix for an unbiased
%           TDOA estimator (Ndim x Ndim)

% Parse inputs
if nargin < 6 || ~exist('resample_covariance','var')
    resample_covariance = true;
end

if nargin < 5 || ~exist('variance_is_toa','var')
    variance_is_toa = true;
end

if nargin < 4 || ~exist('ref_idx','var')
    ref_idx = [];
end

[n_dim, n_sensor] = size(x_tdoa);
n_source = size(xs,2);

% Construct Jacobian function handle
J = @(x) tdoa.jacobian(x_tdoa,x,ref_idx);

% Preprocess covariance matrix
if variance_is_toa
    C_out = C*utils.constants.c^2;
else
    C_out = C;
end

% Parse sensor pairs
if resample_covariance
    [test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
    C_tilde = utils.resampleCovMtx(C_out, test_idx_vec, ref_idx_vec);
end

% Ensure the covariance matrix is invertible
C_tilde = utils.ensureInvertible(C_tilde);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_tilde,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_tilde);
end

% Initialize output variable
crlb = zeros([n_dim,n_dim,n_source]);

warning('off','MATLAB:nearlySingularMatrix')
        
% Repeat CRLB for each of the n_source test positions
for idx =1:n_source
    this_x = xs(:,idx);
    
    % Evaluate Jacobian at x_i
    J_i = J(this_x);
    
    % Compute Fisher Information Matrix
    if do_decomp
        F = J_i/C_d*J_i'; % Ndim x Ndim
    else
        F = J_i*C_inv*J_i';
    end
    
    if any(isnan(F(:))) || any(isinf(F(:)))
        % Problem is ill defined, Fisher Information Matrix cannot be
        % inverted
        crlb(:,:,idx) = NaN;
    elseif any(diag(F)<= 1e-15)
        % Problem is ill-defined
        valid_ind = diag(F) > 1e-15;
        crlb(~valid_ind,~valid_ind,idx) = Inf;
        crlb(valid_ind,valid_ind,idx) = inv(F(valid_ind,valid_ind));
    else
        % Invert the Fisher Information Matrix to compute the CRLB
        C = inv(F);
        if any(diag(C)<0)
            % We got a negative noise term; invalid result
            crlb(:,:,idx) = NaN;
        else
            crlb(:,:,idx) = C;
        end
    end
end

warning('on','MATLAB:nearlySingularMatrix');
end


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
end


end