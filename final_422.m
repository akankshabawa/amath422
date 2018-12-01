%% BASIC MODEL

clear all;

%GENERAL PARAMETERS
L = 5; %Length of habitat
n = 201; %Number of spatial subintervals over habitat
dx = L/(n-1); %Length of spatial subinterval
k = 100; %Carrying capacity
init_sigma = 2; %Mean dispersal distance for establishing initial distribution
xe = 1; %Speed at which habitat moves for each timestep
R = 2; %Net reproductive rate
N = 50; %Number of iterations to run simulation

%DISPERSAL FUNCTION PARAMETERS
dispersal_type = "gaussian"; %Type of dispersal (ie gauss, cauchy, lapalce)
var = 1; %Variance of dispersal kernel

%GET THE POPULATION IN EQUILIBRIUM
D = xe*50 + 2*L; %Total length of landscape (ensures length is long enough that pop
                 %won't crash after first iteration)
x = -D/2:dx:D/2; %vector of positions over entire landscape D
init_patch_center = -D/2 + 1.2*L; %center of initial habitat
y = -L/2+init_patch_center:dx:L/2+init_patch_center;
left_index = (D/2 - L/2 + init_patch_center)/dx + 1; %index of x where habitat starts
right_index = left_index + L/dx; %index of x where habitat ends
%size(y) = 1x201
%size(x) = 1x2401
n0 = init_dist(init_sigma, init_patch_center, x, k); %Initial distribution

%RUN DETERMINISTIC MODEL
burn = 10; %Number of iterations to get equilibrium before adding climate change
param = 1;
for i = 1:burn
    patch_n0 = n0(left_index:right_index); %previous densities in habitat
    K = outerFunc(x,y,dispersal_type, param);
    size(K);
    size(patch_n0);
    blah = bev_holt(patch_n0, R, k);
    size(blah);
    newn = dx.*K*blah.';
    newn = newn.';
    n0 = newn;
    y = x(left_index:right_index);
    
end

for t = 2:N
          patch_n0 = n0(left_index:right_index);
          K = outerFunc(x,y,dispersal_type, param);         
          newn = dx*K*bev_holt(patch_n0, R, k).';
          newn = newn.';
          %plot new population density newn in grey
          plot(x, newn, 'b'); hold on;
         
          %end simulation if population is -->0
          if(trapz(x,newn)< 5)
          break  
          end
          %If population is not approaching 0, continue simulation         
          %assign the new population to be the current population
          n0 = newn;
          %shift habitat limits according to c_left and c_right
          left_index = left_index + xe/dx;
          right_index = min(right_index+xe/dx,length(x));  
                    
          %grab the new vector of habitat locations y
          %if the habitat has shifted past the domain boundary on the right, only grab the relevant positions
          y = x(left_index:right_index);
end

xlabel('Spatial Gradient');
ylabel('Population Density');
set(gca,'FontSize',16);

function outer = outerFunc(x,y, dispersal_type, param)
    outer = zeros(length(x), length(y));
    for i = 1 : length(x)
        for j = 1 : length(y)
            outer(i,j) = kernel(x(i),y(j),dispersal_type, param);
        end
    end
end

%DISPERSAL FUNCTION
function disp = kernel(x, y, dispersal_type, param)
    if(dispersal_type == "gaussian")
        sigma = param;
        disp = (1/(sqrt(2*pi)*sigma)*exp((-((x-y).^2))/(2*sigma.^2)));
    elseif(dispersal_type == "laplace")
        alpha = param;
        disp = (alpha/2*exp(-alpha*abs(x-y)));
    elseif(dispersal_type == "cauchy")
        alpha = param;
        disp = (alpha/pi/(alpha.^2+(x-y).^2));
    end
end

%GROWTH FUNCTION (BEVERTON-HOLT HERE)
function growth = bev_holt(patch_n0, R, k)
    growth = (R.*patch_n0 ./ (1 + (R-1)/k.*patch_n0));
end

%INITIAL DISTRIBUTION FUNCTION (USES GAUSSIAN KERNEL)
function n0 = init_dist(init_sigma, init_patch_center, x, k)
    n0 = k/(sqrt(2*pi)*init_sigma)*exp(-(x-(init_patch_center+0)).^2/(2*init_sigma.^2));
end
