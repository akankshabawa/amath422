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
%{
function F=fvals(R,patch_n0,k)

F=R*patch_n0/(1+(R-1)/k*patch_n0);
%}

%Define the dispersal function and relevant parameters
    %One-dimentional dispersal kernel
    %x,y: destination and source locations in one dimentional space
    %dispersal.type: functional form of dispersal kernel
    %param: controls variance of dispersal kernel

%Parameters
    %dispersal.type: reactice parameter defining functional form of
    %dispersal kernel
        %Can accomidate Gaussian, Laplace, or Cauchy
    %param: Reactive parameter defining variance of dispersal kernel
