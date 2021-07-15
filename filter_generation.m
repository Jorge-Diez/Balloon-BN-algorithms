clear all;

for testrun = 0:10
    
    clearvars -except testrun experiment_n
    addpath ('/home/jorge/Desktop/CasualFMRI/NETSIM-Balloon');
    nma_dir = ('/home/jorge/Desktop/CasualFMRI/NETSIM-Balloon/nma_2010_Feb23cut/');
    
    base_exp_dir = "/home/jorge/Desktop/CasualFMRI/NETSIM-Balloon/BalloonSims";
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    std_BOLD_noise=1;   % added thermal (measurement) noise
    
    
    tr=1.2;               % temporal sampling of BOLD timeseries
    nsecs=600;          % session length (per "subject") in seconds
    
    Nsub=60;            % number of "subjects"; must be even: needs to be at least as large as max(N*(N-1)/2,50) to get enough samples in null
    HRF_covar=4;        % variability in HRF delay; 4 gives ~0.5s variability
    Nsubnets=1;         % how many 5-node subnets do we setup?
    
    
    sigma_prior_mean=20;    % 1/neural lag;  e.g. 20 or 10, originally was 1 in DCM/NMA - but that gives a 1s lag neurally!!
    
    
    
    
    nscans=nsecs/tr;
    res=1/200;           % 1/20    raw data resolution in seconds
    resfactor=tr/res;
    ntpts=nsecs/res;
    
    ts=[];
    poop=[];
    
    
    
    filtercut = 200;
    
    
    
    
    switch testrun
        case {0}
            exp_name = "60sub-1200ms-10min";
        case {1}
            tr=0.25;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            exp_name = "tr250ms";
        case {2}
            tr=0.75;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            exp_name = "tr750ms";
        case {3}
            tr=3.0;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            exp_name = "tr3000ms";
        case {4}
            Nsub=10;
            exp_name = "10sub";
        case {5}
            Nsub=30;
            exp_name = "30sub";
        case {6}
            Nsub=100;
            exp_name = "100sub";
        case {7}
            nsecs=300;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            exp_name = "5min";
        case {8}
            nsecs=1800;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            exp_name = "30min";
        case {9}
            nsecs=3600;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            exp_name = "60min";
        case {10}
            tr=0.72;
            nsecs = 4500;
            nscans=nsecs/tr;
            res=1/200;           % 1/20    raw data resolution in seconds
            resfactor=tr/res;
            ntpts=nsecs/res;
            Nsub = 100;
            exp_name = "excellentconditions";
            
            
    end
    
    
    filter = filtermumford(filtercut, tr, nscans);
    
    disp('Filter for ' + exp_name + ' saved');
    
    savename =  'filter_' + exp_name  +  '.mat' ;
    
    save(savename,'filter')
    
   
    
    
end
    
    
    
    
    
    
    
    
