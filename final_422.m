%% BASIC MODEL
% DESCRIPTIVE TEXT

% N: Number of iterations to run the simulation over
% D: Total length of fixed landscape - D depends on parameters and is
% defined below
% L: Habitat length
% n: Number of spatial sunintervals over habitat. Number must be odd
% dx: Length of spatial subinterval
% k: Carrying capacity
% init_sigma: Mean dispersal distance for establishing initial distribution
% xe: Speed at which habitat moves at each time step
% R: Net reproductive rate

L=5
n=201
N=50
dx=L/(n-1)
k=100
init_sigma=2
xe=1
R=2

%fvals: Growth function applied to the population density within the
%habitat
%Beverton-Holt Curve:

