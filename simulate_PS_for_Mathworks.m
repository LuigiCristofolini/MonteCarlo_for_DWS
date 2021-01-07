%% This script performs a Monte Carlo simulation to determine the path length distributions P(s) for the analysis of Diffusing Wave Spectroscopy (DWS) experiments
% Luigi Cristofolini 2020
% NB: This is  a template script, and you are free to adapt it your needs!
% More details on the script, and on its usage, can be found in: Lorusso et al., Advances in Colloids and Interface Science 288 (2021) 102341 
% It is based on the following literature:
% [1] Weitz, Pine, "Diffusing wave spectroscopy", in "Dynamic Light Scattering", Oxford University Press (1993)
% [2] Durian Physical Review E 51, 3350 (1995)
% [3] Pine, Weitz et al. Phys Rev Lett 60 1134 (1988)


% MATLAB computing environment is chosen becaues  it is quite efficient by allowing implicit parallelization if the simulation is articulated in photon bunches. 
% Optimal bunch size is of the order of nParallel=10^5 photons with a typical Win10 64-bit running on 8 GB RAM. 
% This number is optimized on my PC (Matlab 202ob on a W10 machine with 8Gb RAM), but please experience other values on you machine

% Many bunches are successively and independently generated to build up a statistical ensemble. 

% The philosophy of the simulation is to follow photon propagation by a random walk process, keeping trace of the fate of each photon. 
% We choose a reference frame in which X is the horizontal direction of the impinging laser beam, Y is the other horizontal direction and Z is vertical. 
% At the initialization step, we define the cuvette size and two detection areas, one for backscattering and one for transmission experiments. 
% Initial positions for each photon are generated inside the sample, at a depth l^*, 
% with lateral distribution in the YZ plane either uniform or gaussian to account for these two possible beam profiles. 
% Then the core of the Monte Carlo algorithm consists in generating the directions of the steps uniformly distributed over the sphere, while the step lengths, 
% following the approach of Durian et al. [2], follows a lognormal distribution with average value l^*, rather than being identically fixed to l^*. 
% This choice has no effect in the results for transmission geometry but yields better reproduction of the analytically exact formulas for the correlation functions simulated for the infinite slab in backscattering geometry. 
% At each step, the fate of each photon is checked: if it is still inside the sample volume, the counter of the path length is incremented by one unit and the photon is retained for the next steps of the simulation. 
% If on the contrary the photon has come out of the sample volume, simulation for this photon stops, 
% and depending on whether this happened outside the measurement areas or inside one of them, the corresponding counter of the steps is discarded or retained in the appropriate (backscattering or transmission) pool. 
% This corresponds to the absorbing confinement condition, which we choose as the most realistic reproduction of the real experiment. 
% When enough photons have been simulated, from the pool of path lengths recorded for the backscattering and the transmission geometries, 
% the corresponding path length distributions P(s) are generated and saved to a file. 


% The experimental geometry is described, in S.I. units, by
% X is the direction of the laser (horizontal) impinging  beam
% Y is the other horizontal
% Z is vertical
% (x0, y0, z0) initial positions
% (xx, yy, zz) at each step, the instantaneous positions of the photons in their random walk

% The shape of impinging laser beam can be chosen to be gaussian or homogeneous (if sigma_beam=Inf)

%% 

clear

% simulation parameters
nPhotons=1e8; %total number of photons to be simulated
nParallel=1e5; % number of photons to be simulated in each bunch

% sample parameters, all in S.I. units 
lstar=1.0e-3;  %put here the appropriate value for which youi want the simulation. One could desire to run a number of simualtions for a set of values of l*

cuvette.size_x=1e-2;% put here the size of your cuvette. Here it is the standard 1cm sqaure cuvette containg a sample 3cm high  
cuvette.size_y=1e-2;%
cuvette.size_z=3e-2;%
cuvette.description='1cm_real'; %put anything here to recon your cuvette. It is only used in the save filename
%
% beam parameters
beam.sigma=Inf; %2e-2; %sigma of gaussian beam, set it to inf for the theoretical results obtained for uniform beam

%detector parameters: 
% here the values for our APD fiber head 5 mm diameter (here approximated as a square)
det.size_y=5e-3;   % horizontal size of collection area
det.size_z=5e-3;   % vertical size of collection area
% another example: a Linecamera with an imaging system focused on the out surface
% det.size_y=3e-3;   % horizontal size of collection area
% det.size_z=400e-6; % vertical size of collection area

% initialization of global counters
tic %for timing purposes
BS_detectable.s=[]; BS_detectable.Ps=[];
FW_detectable.s=[]; FW_detectable.Ps=[];

%main cicle
for ii=1:nPhotons/nParallel
    %re-initialize counters for each photon bunch
    nLostPhotons=0; % counter for photons lost on the y or z sides
    BS_photons=zeros(nParallel,3); % array of data for photons exiting in BS
    FW_photons=zeros(nParallel,3); % array of data for photons exiting in FW
    nBackscattering=0;
    nTransmission=0;
    Nsteps=0;
    % photon enters the cuvette from the xx=0 face.
    % in Backscattering DWS, the initial condition is that light diffuses from a position z0 ~ lstar inside the cell. (for transmission, the choice of z0 has negligible effects).

    x0=ones(nParallel,1)*lstar; %assume starting depth as in ref [1]  pg 674 675 (there it is called z0=l*)
    
    if isfinite(beam.sigma)
        y0=randn(nParallel,1)*beam.sigma;
        while sum( abs(y0)>cuvette.size_y/2 )>0
            y0(abs(y0)>cuvette.size_y/2)=randn(sum( abs(y0)>cuvette.size_y/2 ),1)*beam.sigma;
        end
        
        z0=randn(nParallel,1)*beams.igma;
        while sum( abs(z0)>cuvette.size_z/2 )>0
            z0(abs(z0)>cuvette.size_z/2)=randn(sum( abs(z0)>cuvette.size_z/2 ),1)*beam.sigma;
        end
    else
        %  uniform impinging beam:
        y0=cuvette.size_y/2*(2*rand(nParallel,1)-1);
        z0=cuvette.size_z/2*(2*rand(nParallel,1)-1);
    end
    
    xx=x0;
    yy=y0;
    zz=z0;
    
    while length(xx)>1   % repeat was long as at least one photon is inside the  cell
        nn=length(xx); %number of residual photons, possibly decreasing at each cycle
        fi=2*pi*rand(nn,1);
        th=acos(2*rand(nn,1)-1);
        xvers=sin(th).*cos(fi);
        yvers=sin(th).*sin(fi);
        zvers=cos(th);
        stepsize=-lstar*log(rand(nn,1)); % Durian PRE 51, 3350 (1995) 
        xx=xx+stepsize.*xvers;
        yy=yy+stepsize.*yvers;
        zz=zz+stepsize.*zvers;
        Nsteps=Nsteps+1;
        
        % check positions
        % 1) lost photons; drop them from simualtion
        flag_lost=find(abs(yy)>cuvette.size_y/2 | abs(zz)>cuvette.size_z/2);
        nLost=numel(flag_lost);
        xx(flag_lost)=[];     yy(flag_lost)=[];     zz(flag_lost)=[];
        nLostPhotons=nLostPhotons+nLost;
        % 2) save BS photons and drop them from simulation
        flag_BS=find(xx<=0);
        nBS=numel(flag_BS);
        if nBS>0
            BS_photons(nBackscattering+1:nBackscattering+nBS,1)=yy(flag_BS);
            BS_photons(nBackscattering+1:nBackscattering+nBS,2)=zz(flag_BS);
            BS_photons(nBackscattering+1:nBackscattering+nBS,3)=Nsteps;
            nBackscattering=nBackscattering+nBS;
            xx(flag_BS)=[];     yy(flag_BS)=[];     zz(flag_BS)=[];
        end
        % 3) save FW photons and drop them from simulation
        flag_FW=find(xx>=cuvette.size_x);
        nFW=numel(flag_FW);
        if nFW>0
            FW_photons(nTransmission+1:nTransmission+nFW,1)=yy(flag_FW);
            FW_photons(nTransmission+1:nTransmission+nFW,2)=zz(flag_FW);
            FW_photons(nTransmission+1:nTransmission+nFW,3)=Nsteps;
            nTransmission=nTransmission+nFW;
            xx(flag_FW)=[];     yy(flag_FW)=[];     zz(flag_FW)=[];
        end
        
    end
    
    %% Finding photons on detector
    BS_photons(BS_photons(:,3)==0,:)=[]; % get rid of unused part of preallocated vectors
    FW_photons(FW_photons(:,3)==0,:)=[]; % get rid of unused part of preallocated vectors
    
    BS_detectable_photons=BS_photons((abs(BS_photons(:,1))<=det.size_y/2)&(abs(BS_photons(:,2))<=det.size_z/2),:);
    FW_detectable_photons=FW_photons((abs(FW_photons(:,1))<=det.size_y/2)&(abs(FW_photons(:,2))<=det.size_z/2),:);
    
    %% calculate P(s)
    old_max_s= numel(FW_detectable.s);
    new_max_s= max(FW_photons(:,3));
    if old_max_s < new_max_s %then I need to extend the hystogram to include longer paths than previously accounted for (1.st bunch will allways do this, then depending on needs)
%         disp(['I need to extend the hystogram from ',int2str(old_max_s),' to ',int2str(new_max_s),' bins'])
        FW_detectable.s=1:new_max_s; FW_detectable.Ps=[FW_detectable.Ps,zeros(1,new_max_s-old_max_s)];
        BS_detectable.s=1:new_max_s; BS_detectable.Ps=[BS_detectable.Ps,zeros(1,new_max_s-old_max_s)];
    end
    %update histograms of counts
    BS_detectable.Ps=BS_detectable.Ps+histcounts(BS_detectable_photons(:,3),[BS_detectable.s,BS_detectable.s(end)+1]);
    FW_detectable.Ps=FW_detectable.Ps+histcounts(FW_detectable_photons(:,3),[FW_detectable.s,FW_detectable.s(end)+1]);
end

%Finally, normalize histograms transforming counts into frequencies
BS_detectable.Ps=BS_detectable.Ps / sum(BS_detectable.Ps);  
FW_detectable.Ps=FW_detectable.Ps / sum(FW_detectable.Ps); 
toc
%% plot
figure(1)
clf
subplot(211)
semilogx(BS_detectable.s,BS_detectable.Ps,'.-')
title('BS')
xlabel('log_{10}(number of steps)');
ylabel('P(s)');
subplot(212)
semilogx(FW_detectable.s,FW_detectable.Ps,'.-')
title('FW')
xlabel('log_{10}(number of steps)');
ylabel('P(s)');
%% save
fname=['TEST_MC_for_Beam_',num2str(beam.sigma*1e3),'mm_lstar_',num2str(lstar*1e3),'mm_cuvette_',cuvette.description,'_detector_',num2str(det.size_y*1e3,2),'x',num2str(det.size_z*1e3),'mm_',date,'.mat'];
save(fname,'BS_detectable','FW_detectable','lstar','cuvette','det','nPhotons');
