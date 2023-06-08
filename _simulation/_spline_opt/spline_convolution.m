function [X, basis] = spline_convolution(S, n_lags, knots, degree)
%SPLINE_CONVOLUTION convolve signal S with a set of B-splines.
%
%   Parameters:
%
%      S : input signal. This is a n_timesteps x 1 column vector.
%
%      n_lags : number of lags that the spline basis should span.
%
%      knots : if this is a number, it is taken to be the number of the
%   (multiplicity-one) knots for the spline basis. A default set of knots
%   is generated, where the knots will be equally spaced in the [0,1]
%   interval, and extra repeated knots will be added to the ends in order
%   to ensure that on the right hand side of this interval the splines are
%   constrained to go smoothly to zero (the idea is that if you're modeling
%   a kernel you may want to assume that at long enough times your kernel
%   goes to zero, and build this assumption into your spline basis). With
%   this default, the number of basis elements will be knots+degree-2, so
%   for instance a cubic spline set with 4 knots will contain 3+4-2=5
%   splines.
%
%   If knots is a vector, it is simply the full vector of knots to be used
%   (refer to the documentation of octave-bspline in that case).
%
%      degree : the degree of the splines. Default is 3 (cubic splines).
%
%
%   Returns:
%
%      X : a n_timesteps x n_basis_elements matrix. Each column of X
%   contains the convolution of the stimulus X with a different spline
%   basis element. Note that a full convolution is performed and only the
%   first n_timesteps element are kept (so the first n_lags rows of X will
%   show the effect of zero-padding).
%
%      basis : a n_lags x n_basis element matrix. Each column of basis
%   contains one of the basis elements, evaluated at n_lags equally-spaced
%   points that cover the domain of the spline basis. To get a quick sense
%   of what your basis set is, you can just visualize it with plot(basis).

if size(S,2)~=1
    error("S should be a column vector")
end
n_timesteps = size(S,1);

if nargin < 4
    degree = 3; % cubic splines
end

if isscalar(knots)
    n_knots = knots;
    knots = [zeros(1,degree), linspace(0,1,n_knots), ones(1,degree-1)];
end

xrange = linspace(knots(1),knots(end),n_lags);
basis = full(bsplinebasis(degree, knots, xrange)); % we convert to dense format as sparse bases are convenient mostly when dealing with discrete events
n_basis_elements = size(basis,2);

X = zeros(n_timesteps, n_basis_elements);
for b=1:n_basis_elements
    conv_S = conv(S,basis(:,b));
    X(:,b) = conv_S(1:n_timesteps);
end

end