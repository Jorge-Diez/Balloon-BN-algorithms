clear all;

for testrun = 0:10
    for experiment_n = 1:10
        clearvars -except testrun experiment_n
        addpath ('/home/jorge/Desktop/CasualFMRI/NETSIM-Balloon');
        addpath ('/home/jorge/Desktop/CasualFMRI/NETSIM_Balloon_Data');
        nma_dir = ('/home/jorge/Desktop/CasualFMRI/NETSIM-Balloon/nma_2010_Feb23cut/');
        
        base_exp_dir = "/home/jorge/Desktop/CasualFMRI/NETSIM_Balloon_Data";
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        std_BOLD_noise=1;   % added thermal (measurement) noise
        
        
        tr=1.2;               % temporal sampling of BOLD timeseries
        nsecs=600;          % session length (per "subject") in seconds
        
        Nsub=60;            % number of "subjects"; must be even: needs to be at least as large as max(N*(N-1)/2,50) to get enough samples in null
        HRF_covar=4;        % variability in HRF delay; 4 gives ~0.5s variability
        Nsubnets=1;         % how many 5-node subnets do we setup?
        
        
        sigma_prior_mean=20;    % 1/neural lag;  e.g. 20 or 10, originally was 1 in DCM/NMA - but that gives a 1s lag neurally!!
        
        
        SharedInputs=0;
        GroupTwoDiff=0;
        ConfoundGlobalMean=0;
        BackwardsConnections=0;
        Cyclic=0;
        OneExternal=0;
        StrongLinks=0;
        BadROIs=0;
        MoreLinks=0;
        NegativeConnections=0;
        Modulation=0;
        NormICOV=0;
        inputspeedup=1;
        external_A=0;
        simdesc='';
        n = 1;
        
        
        nscans=nsecs/tr;
        pre_res = 200; %used for boxcar model specifications
        res=1/pre_res;           % 1/20    raw data resolution in seconds
        resfactor=tr/res;
        ntpts=nsecs/res;
        
        ts=[];
        poop=[];
        
        
        N_nodes = 5;
        fullbold = zeros(Nsub*nscans,N_nodes);
        filtercut = 200;
        
        
        
        
        switch testrun
            case {0}
                exp_name = "60sub-1200ms-10min";
                load('../NETSIM_filters/filter_60sub-1200ms-10min.mat');
            case {1}
                tr=0.25;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                exp_name = "tr250ms";
                load('../NETSIM_filters/filter_tr250ms.mat');
            case {2}
                tr=0.75;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                exp_name = "tr750ms";
                load('../NETSIM_filters/filters/filter_tr750ms.mat');
            case {3}
                tr=3.0;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                exp_name = "tr3000ms";
                load('../NETSIM_filters/filter_tr3000ms.mat');
            case {4}
                Nsub=10;
                exp_name = "10sub";
                load('../NETSIM_filters/filter_10sub.mat');
            case {5}
                Nsub=30;
                exp_name = "30sub";
                load('../NETSIM_filters/filter_30sub.mat');
            case {6}
                Nsub=100;
                exp_name = "100sub";
                load('../NETSIM_filters/filter_100sub.mat');
            case {7}
                nsecs=300;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                exp_name = "5min";
                load('../NETSIM_filters/filter_5min.mat');
            case {8}
                nsecs=1800;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                exp_name = "30min";
                load('../NETSIM_filters/filter_30min.mat');
            case {9}
                nsecs=3600;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                exp_name = "60min";
                load('../NETSIM_filters/filter_60min.mat');
            case {10}
                tr=0.72;
                nsecs = 4500;
                nscans=nsecs/tr;
                pre_res = 200; %used for boxcar model specifications
                res=1/pre_res;           % 1/20    raw data resolution in seconds
                resfactor=tr/res;
                ntpts=nsecs/res;
                Nsub = 100;
                exp_name = "excellentconditions";
                load('../NETSIM_filters/filter_excellentconditions.mat');
                
        end
        
        
        
        
        
        experiment_dir_name = 'Cycle-Boxcar_Resting-State_' + exp_name + '/run' + string(experiment_n);
        experiment_dir = base_exp_dir + '/' + experiment_dir_name;
        
        
        
        if (isfolder(experiment_dir))
            rmdir(experiment_dir, 's');
        end
        
        
        if ~(isfolder(experiment_dir))
            mkdir(experiment_dir)
        end
        
        all_bold_subject_data_dir = experiment_dir + "/" + "BOLD_DATA/NOFILTER";
        all_bold_subject_data_dir_filtered = experiment_dir + "/" + "BOLD_DATA/FILTERED";
        
        
        if ~(isfolder(all_bold_subject_data_dir))
            mkdir(all_bold_subject_data_dir)
        end
        
        if ~(isfolder(all_bold_subject_data_dir_filtered))
            mkdir(all_bold_subject_data_dir_filtered)
        end
        
        
        
        %box car external model
        %times are specified in seconds, and later multiplied
        %This code can be moved up, need to move it up later when long
        %experiments are due
        
        task_time = ones(20 * pre_res,1);
        rest_time = zeros(10 * pre_res,1);
        
        block = vertcat(rest_time, task_time);
        
        %How many blocks of rest + task the model has
        combinations = 20;
        boxcar_stim = repmat(block, [combinations,1]);
        
        
        
        %box car external model
        %times are specified in seconds, and later multiplied
        %This code can be moved up, need to move it up later when long
        %experiments are due
        
        task_time = ones(40 * pre_res,1);
        rest_time = zeros(20 * pre_res,1);
        
        block = vertcat(rest_time, task_time);
        
        %How many blocks of rest + task the model has
        combinations = nsecs / 60;
        boxcar_stim = repmat(block, [combinations,1]);
        
        
        
        
        %addpath('/home/jorge/Desktop/CasualFMRI/NETSIM-Balloon/BalloonSims');
        
        
        addpath (experiment_dir)
        
        outname=sprintf('netsim-%d',testrun);
        
        
        
        for n = 1:Nsub/2
            
            
            
            
            for I = 1:5:5*Nsubnets
                
                pdnormal = makedist('Normal', 0.5, 0.1 );
                pd_truncated = truncate(pdnormal, 0.3, 0.7);
                
                for i=I:I+4, a(i,i)=-1; end;    % no self-connections
                
                %                 for i=I:I+3                     % 1st forward connection
                %                     a(i,i+1) = random(pd_truncated, 1);
                %                 end
                
                a(2,1) = random(pd_truncated, 1);
                a(3,2) = random(pd_truncated, 1);
                a(4,3) = random(pd_truncated, 1);
                a(4,5) = random(pd_truncated, 1);
                a(5,1) = random(pd_truncated, 1);
                
                a(2,3) = random(pd_truncated, 1);
                
                
                %Creation of the 5x6 C matrix
                c = zeros(5,6);
                for i = 1:5
                    c(i,i) = 1;
                end
                
                c(4,6) = 1;
                %adding the external box car model input
                
                if MoreLinks == 1
                    a(I+1,I+3)=random(pd_truncated, 1);
                    a(I+2,I+4)=random(pd_truncated, 1);
                end
                
                
                
                if BackwardsConnections == 1        % backward connections
                    for i=I:I+4, for j=I:I+4
                            if a(j,i)>0 && rand>=0.5
                                a(i,j)= -random(pd_truncated, 1);
                            end
                        end; end
                end
                
                %a(I,I+4)=random(pd_truncated, 1);
                if Cyclic == 1        % cyclic connections
                    %a(I,I+4)=0;
                    a(I+1,I)=random(pd_truncated, 1);
                end
                
               
                
                if StrongLinks == 1
                    a=a + 1.25 * a .* (a>0); % strong (mean 0.9)
                    if OneExternal > 0        % pretty much only one external input, on node 1
                        c(2,2)=.1; c(3,3)=.1; c(4,4)=.1; c(5,5)=.1;
                    end
                end
                
                if SharedInputs == 1        % shared inputs with chance 20%
                    for i=I:I+4, for j=I:I+4
                            if j~=i && rand>0.8
                                c(i,j)=0.3;
                            end
                        end ; end
                end
            end
            
            %%%%%% cross-subnet connections
            % Carefil. Nsubnets == 5 might be bad
            if Nsubnets==2
                a(3,8)=random(pd_truncated, 1);
            end
            
            if Nsubnets==3
                a(3,8)=random(pd_truncated, 1);
                a(3,13)=random(pd_truncated, 1);
                a(8,13)=random(pd_truncated, 1);
            end
            if Nsubnets==10
                a(3,8)=random(pd_truncated, 1);
                a(8,13)=random(pd_truncated, 1);
                a(13,18)=random(pd_truncated, 1);
                a(18,23)=random(pd_truncated, 1);
                a(3,23)=random(pd_truncated, 1);
                a(28,33)=random(pd_truncated, 1);
                a(33,38)=random(pd_truncated, 1);
                a(38,43)=random(pd_truncated, 1);
                a(43,48)=random(pd_truncated, 1);
                a(28,48)=random(pd_truncated, 1);
                a(3,28)=random(pd_truncated, 1);
            end
            
            N=size(a,1); Nreal=N;  % used to specify when there are additional nodes to be deleted later
            if Modulation > 0
                if Modulation < 3     % original modulations used in netsim paper
                    a=[a a*0]; a=[a;a*0];    a=a + 1.25 * a .* (a>0); N=N*2;  % strong (mean 0.9)
                    c=eye(N);   for j=2:N/2, c(j,j)=.3; end;
                    if Modulation == 1
                        DOMODext=zeros(N,N,N);  DOMODext(1,2,6)=-0.12; DOMODext(2,3,7)=-0.12; DOMODext(3,4,8)=-0.12; DOMODext(4,5,9)=-0.12; DOMODext(1,5,10)=-0.12;
                    end
                    
                elseif Modulation < 5    % new modulations for low-TR HCP work
                    a=a + 0.5 * a .* (a>0);  % slightly stronger links
                    c=0.6 * eye(N);          % slightly weakened external inputs
                    a(3,28)=-a(3,28);        % make link between two halves negative
                    if Modulation == 3
                        a=[a zeros(N,4); zeros(4,N+4)]; N=N+4; c=0.6*eye(N);
                        DOMODext=zeros(N,N,N); modval=-0.085;
                        DOMODext(11,12,51)=modval; DOMODext(12,13,51)=modval; DOMODext(13,18,51)=modval;
                        DOMODext(1,2,52)=modval; DOMODext(2,3,52)=modval; DOMODext(3,4,52)=modval; DOMODext(4,5,52)=modval; DOMODext(1,5,52)=modval;
                        DOMODext(46,50,53)=modval; DOMODext(47,48,53)=modval; DOMODext(48,49,53)=modval; DOMODext(49,50,53)=modval;
                        DOMODext(43,48,54)=modval; DOMODext(48,49,54)=modval;
                    end
                    
                else  % (Mod = 10)
                    DOMODint=zeros(N,N,N);  DOMODint(6,7,5)=2; DOMODint(1,2,4)=2;
                end
            end
            
            
            if NegativeConnections > 0
                a(2,3)=-a(2,3);
            end
            
            %dont fully know whiy this is here, will leave for now TODO
            aa(n,:,:)=a; cc(n,:,:)=c;
            
            %%%%%% setup second-half of the subjects
            if GroupTwoDiff==1
                a(2,3)=a(2,3)*.5;
            end
            aa(n+Nsub/2,:,:)=a; cc(n+Nsub/2,:,:)=c;
            
        end
        
        save(experiment_dir + '/' + outname + '.mat');
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CREATE DATA FOR THE BALLOON MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %HERE WE SUPPOSE THAT SUBJECT FOLERS AND AS SUCH SUBJECT DATA HAS NOT BEEN
        %GENERATED
        
        
        
        
        
        
        
        
        for subject_number = 1:Nsub
            
            
            %cd /home/jorge/Desktop/CasualFMRI/NETSIMdist/BalloonSims;
            %homedir=pwd;
            datadir=experiment_dir;
            
            
            if ~(isfolder(datadir))
                mkdir(datadir)
            end
            
            
            subjectname=sprintf('subject%d',subject_number);
            mkdir(strcat(datadir,'/',subjectname));
            
            % setup confound evs
            confounds=mwdct(nsecs/tr,1);  % 14 for 'normal real data', 1 for simulations, using little highpass filtering
            tmp=devar(confounds(:,2:end),1);
            tmp2=[confounds(:,1), tmp];
            confounds=tmp2;
            save(strcat(datadir,'/',subjectname,'/confound_evs.txt'),'confounds','-ascii');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % irrelevant stuff
            N=length(a);
            P=N+1;
            roinames=cell(1,N);
            inputevnames=cell(1,N+1);
            
            for i=1:N, roinames{i}=sprintf('%d',i); inputevnames{i}=sprintf('%d',i); end;
            model=1;
            inputevnames{N+1} = sprintf('%d',N+1);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % a(1,2) is the input of node 1 into node 2
            a=squeeze(aa(subject_number,:,:));
            
            % b(1,2,3) is stimulation 3 modulating the input of node 1 into node 2
            b=zeros(N,P,P);
            if exist('DOMODext')
                b=DOMODext;
            end
            
            % c(1,2) is the input of stimulation 2 into node 1
            c=squeeze(cc(subject_number,:,:));
            
            % d(1,2,3) is node 3 modulating the input of node 1 into node 2
            d=zeros(N,N,N);
            if exist('DOMODint')
                d=DOMODint;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %did this to avoid unnecessary stuff
            
            model_name=strcat('model_art_',num2str(1));
            auxdatadir = strcat(datadir,'/',model_name);
            
            
            
            
            
            
            % save models
            models{model}.a=a;models{model}.b=b;models{model}.c=c;models{model}.d=d;
            for m=1:length(models),
                model_name=strcat('model_art_',num2str(subject_number));
                mkdir(strcat(datadir,'/',model_name));
                
                save(strcat(datadir,'/',model_name,'/','sigmaa_prior_mean.txt'),'sigma_prior_mean','-ascii');
                tmp=100;  % originally was 100
                save(strcat(datadir,'/',model_name,'/','sigmaa_prior_var.txt'),'tmp','-ascii');
                
                for i=1:length(roinames),
                    balloon_priors{i}=balloon_defaults();
                    balloon_priors{i}.balloon_cbf_mean=mvnrnd(balloon_priors{i}.balloon_cbf_mean,balloon_priors{i}.balloon_cbf_cov)';
                    balloon_priors{i}.balloon_cbf_mean=max(balloon_priors{i}.balloon_cbf_mean,0.1);
                    balloon_priors{i}.balloon_mean=mvnrnd(balloon_priors{i}.balloon_mean,HRF_covar*balloon_priors{i}.balloon_cov)';  % factor of 4 gives an HRF delay of +- 0.5s
                    balloon_priors{i}.balloon_mean=max(balloon_priors{i}.balloon_mean,0.1);
                    balloon_priors{i}.balloon2_mean=mvnrnd(balloon_priors{i}.balloon2_mean,balloon_priors{i}.balloon2_cov)';
                    balloon_priors{i}.balloon2_mean=max(balloon_priors{i}.balloon2_mean,0.1);
                end;
                
                for i=1:length(roinames),
                    tmp=balloon_priors{i}.balloon_cbf_mean;
                    save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_cbf_mean.txt'),'tmp','-ascii');
                    tmp=balloon_priors{i}.balloon_cbf_cov;
                    save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_cbf_cov.txt'),'tmp','-ascii');
                    
                    tmp=balloon_priors{i}.balloon_mean;
                    save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_mean.txt'),'tmp','-ascii');
                    tmp=balloon_priors{i}.balloon_cov;
                    save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon_cov.txt'),'tmp','-ascii');
                    
                    tmp=balloon_priors{i}.balloon2_mean;
                    save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon2_mean.txt'),'tmp','-ascii');
                    tmp=balloon_priors{i}.balloon2_cov;
                    save(strcat(datadir,'/',model_name,'/',roinames{i},'_balloon2_cov.txt'),'tmp','-ascii');
                end;
                
                
                tmp=models{m}.a';
                save(strcat(datadir,'/model_art_',num2str(subject_number),'/matA.txt'),'tmp','-ascii','-double');
                tmp=models{m}.c;
                save(strcat(datadir,'/model_art_',num2str(subject_number),'/matC.txt'),'tmp','-ascii','-double');
                for i=1:size(models{m}.b,3),
                    tmp=models{m}.b(:,:,i);
                    save(strcat(datadir,'/model_art_',num2str(subject_number),'/matB_into_',num2str(i),'.txt'),'tmp','-ascii','-double');
                end;
                for i=1:size(models{m}.d,1),
                    tmp=squeeze(models{m}.d(:,i,:));
                    save(strcat(datadir,'/model_art_',num2str(subject_number),'/matD_into_',roinames{i},'.txt'),'tmp','-ascii','-double');
                    strcat(datadir,'/model_art_',num2str(subject_number),'/matD_into_',roinames{i},'.txt');
                end;
            end;
            
            % sim stimuli from GMM
            
            ddd=1/(res*sigma_prior_mean);    % adjust neural stim strength according to neural time-constant
            in_strength=80/sum(exp(-[0:10*ddd]/ddd));
            
            
            
            for e=1:N
                
                mix.centres=[0 in_strength];
                mix.stds=[in_strength/20 in_strength/20];
                trans_prob= res * [0.1 0.4] * inputspeedup;
                
                if e > Nreal
                    trans_prob = res * [0.05 0.03];
                end
                
                states=zeros(ntpts,1);
                stim=zeros(ntpts,1);
                stim(1)=normrnd(mix.centres(1),mix.stds(1));
                for t=2:ntpts,
                    if(states(t-1)==0),
                        if(rand<trans_prob(1)),
                            states(t)=1;
                        else,
                            states(t)=0;
                        end;
                    else,
                        if(rand<trans_prob(2)),
                            states(t)=0;
                        else,
                            states(t)=1;
                        end;
                    end;
                    
                    if(states(t)==0)
                        stim(t)=normrnd(mix.centres(1),mix.stds(1));
                    else
                        stim(t)=normrnd(mix.centres(2),mix.stds(2));
                    end;
                end;
                
                
                if inputspeedup<0
                    stim=3*in_strength*randn(ntpts,1);
                end
                
                if exist('DONEstim')  % external stimuli already setup in matrix DONEstim
                    stim=DONEstim(:,e);
                end
                
                save(strcat(datadir,'/',subjectname,'/',inputevnames{e},'.txt'),'stim','-ascii','-double');
            end;
            
            
            % save the stimulation boxcar model data
            save(strcat(datadir,'/',subjectname,'/',inputevnames{N+1},'.txt'),'boxcar_stim','-ascii','-double');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %RUN THE BALLOON MODEL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %HERE WE SUPPOSE THAT SUBJECT FOLERS AND AS SUCH SUBJECT DATA HAS NOT BEEN
            %GENERATED
            
            
            
            
            %%%%%%%
            % sim data
            truemodel=1;
            a=models{truemodel}.a;
            b=models{truemodel}.b;
            c=models{truemodel}.c;
            d=models{truemodel}.d;
            stim_list='1';   roi_list='1';
            for i=2:N, roi_list=sprintf('%s,%d',roi_list,i); end;
            for i=2:N+1, stim_list=sprintf('%s,%d',stim_list,i); end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %cd ../nma_2010_Feb23cut
            
            model_call = sprintf ('%s./nma_generate --ld=%s/tmp_model --dd=%s --md=%s/model_art_%i  --tr=%d --resfactor=%d --sub=%s --stim=%s --nn=%s --dm=0 --hm=balloon --ns=%i --sascf --to',nma_dir,datadir,datadir,datadir,subject_number,tr,resfactor,subjectname,stim_list,roi_list,nscans);
            
            [nma_status,nma_result] = system(model_call);
            if nma_status == 0
                
                for r=1:length(roinames),
                    ytrue(:,r)=load(sprintf('%s/tmp_model/signal_%s_bold.txt',datadir,roinames{r}));
                    ztrue_nma_gen(:,r)=load(sprintf('%s/tmp_model/signal_%s_z.txt',datadir,roinames{r}));
                end;
                
                % debug HRF variability
                poop=[poop ytrue];
                
                %%%%%%%% add noise
                for(i=1:size(ytrue,2))
                    boldfolder = strcat(datadir,'/',subjectname,'/','BOLD');
                    if ~(isfolder(boldfolder))
                        mkdir(boldfolder)
                    end
                    
                    %Save original bold image
                    y=ytrue(:,i);
                    save(strcat(boldfolder,'/',roinames{i},'_boldclean.txt'),'y','-ascii');
                    y2=y+randn(size(y))*std_BOLD_noise;
                    y2s(:,i)=y2;
                    
                    %Save bold with added noise
                    save(strcat(boldfolder,'/',roinames{i},'_bold.txt'),'y2','-ascii');
                    
                    
                    
                end;
                
                %save the full bold times
                writematrix(y2s, strcat(all_bold_subject_data_dir,'/','subject',num2str(subject_number),'bold.csv'))
                
                
                %filter and save the bold times
                y2sfiltered = filter * y2s;
                writematrix(y2sfiltered, strcat(all_bold_subject_data_dir_filtered,'/','subject',num2str(subject_number),'bold_filtered.csv'))
                
                
                fullbold( (((subject_number-1)*nscans)+1):(subject_number*nscans),:) = y2s;
                
                ts=[ts;y2s];
            else
                nma_result;
            end;
            
            %Delete things that we dont need to avoid hogging too much space
            %We will only delete things that take up a lot of space
            
            %Delete the tmp model folder
            rmdir(sprintf("%s/tmp_model",datadir), 's')
            
            
            %Delete the stimulus files
            for e=1:N+1
                delete((strcat(datadir,'/',subjectname,'/',inputevnames{e},'.txt')));
                
            end
            
            disp("Subject " + subject_number  + " finished of run " + experiment_n + " of boxcar cyclic smith structure " + exp_name )
            
            
            
            
        end;
        
        
        writematrix(fullbold, strcat(experiment_dir, "/bold_data.csv"));
    end
    %Average time: 560 seconds with 50 subjects
    
    
    
end


%check the A matrix
amatrix = '/home/jorge/Desktop/CasualFMRI/NETSIM_Balloon_Data/Boxcar_Resting-State_60sub-1200ms-10min/run1/model_art_3/matA.txt';

a = readmatrix(amatrix, ...
    'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join');

a = a';










