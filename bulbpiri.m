function bulbpiri(date_hour_in)

% without clear sometimes some of the global variables get carried over
% from who knows where.
if (nargin==0)
    clear
end

% if (1==2)
%     rand('twister',sum(100*clock));   % added on Juy 18, 15:30
%     randn('state',sum(100*clock));
% else
%     rand('twister',18);   % added on Juy 18, 15:30
%     randn('state',7);
% end

% need to assign job number before the parallel jobs get started and
% possibly change the job number
run_num_for_saving_figs = dlmread('fignum.txt');
%save = run_num_for_saving_figs+1;
dlmwrite('fignum.txt',run_num_for_saving_figs+1);
% determine the computer on which this is running
[~,hostname]=system('hostname');fprintf(' Running on %s \n',hostname);

% with parfor the graphics don't work, for interactive runs the loop is
% omitted.
% parallel parameters:
% 1: th,
% 2: ts,
% 3: xposition,
% 4: wt1_later,
% 5: connM
%10: random
log_parameters=10;
parameters=[1];
n_parameters=length(parameters);

if strcmpi(hostname(1:5),'qnode')
    if (n_parameters==1)
        log_parallel_ode=1; log_parallel_param=0;
        parpool(2)
    else
        log_parallel_ode=0; log_parallel_param=1;
        parpool(min(n_parameters,20));
    end
else
    log_parallel_ode=0; log_parallel_param=0;
end
fprintf(' log_parallel_ode = %g log_parallel_param = %g \n',log_parallel_ode,log_parallel_param);

if (log_parallel_param==0)
    if (nargin>0)
        for ipar=1:length(parameters)
            neurogenesis(0,parameters,ipar,log_parameters,log_parallel_ode,run_num_for_saving_figs,date_hour_in);
        end
    else
        for ipar=1:length(parameters)
            neurogenesis(0,parameters,ipar,log_parameters,log_parallel_ode,run_num_for_saving_figs);
        end
    end
    %     if (nargin>0)
    %         neurogenesis(0,parameters,1,log_parameters,log_parallel_ode,run_num_for_saving_figs,date_hour_in);
    %     else
    %         neurogenesis(0,parameters,1,log_parameters,log_parallel_ode,run_num_for_saving_figs);
    %     end
else
    % parallel jobs
    if (nargin>0)
        parfor ipar=1:length(parameters)
            neurogenesis(0,parameters,ipar,log_parameters,log_parallel_ode,run_num_for_saving_figs,date_hour_in);
        end
    else
        parfor ipar=1:length(parameters)
            neurogenesis(0,parameters,ipar,log_parameters,log_parallel_ode,run_num_for_saving_figs);
        end
    end
end

% jobs submitted on quest should exit after bulbpiri is run
if (nargin>0)
    exit
end

end% main

%Need a way to turn interaction (saving, pulling from txt files, pulling
%from odor folder, etc...)

%Need to mess with relearning

%need to mess with cell thresholds


%reduce wt1? .2 may be too large. want less g cells

%background st m cells and g cells are not negative

%since we're learning with 2*P we should always pass 2*P to the
%cal_activity1, want the m cell input

%two avenues to control g cells use th or wt1; now our m is bigger, how
%much will they be driven down?


function objfun = neurogenesis(run_num,parameters,ipar,log_parameters,log_parallel_ode,run_num_for_saving_figs,date_hour_in)
global tau_P W_PM W_PP_i W_PP c_assoc thresh steep P_amp wt1_initial wt1_later wt1_later2 CS W_PP_max
global file_type_data file_type_figure file_type_movie file_type_matlab
global time_wt1_1 time_wt1_2 W_PM_scale P_scale choose mix_factor probe wt1
global n_P Nc Wgm Wmg N_mitral learning
global file_location_figure file_tag file_location format_string
global pos fig_label
global stim_type
global prod_gp prod_gm I Wpg
global Mthresh Msteep M_amp Gthresh Gsteep G_amp
global non_lin_M non_lin_G
global num_light_cells
global read_out_weights
global fisher_plot fisher_thresh
global i_state state_probe learning_time  probe_case lxp1 lxp2 lxp3 lxp4
global randomize random_order odor_random shuffle
global c_assoc_learn c_assoc_recall
global log_probe_random
global thresh_odor steep_odor P_amp_odor
global thresh_light steep_light P_amp_light
global maxv minv ts gamma
global log_first_run P_last log_parallel
global parameter_variation

log_first_run=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_parallel=log_parallel_ode;

% to control also the random number generator randi used in randperm
if (1==2)
    rng('shuffle');
else
    if (log_parameters==6||log_parameters==10)
        seed=parameters(ipar);
    else
        seed=7;
    end
    rng(seed,'twister');
end
first_random=rand;
fprintf(' First random number = %g \n',first_random);

if (log_parameters==10)
    % to vary a number of parameters randomly the parameter parameter_variation should be set to a non-zero value
    %    parameter_variation=0.15;
    parameter_variation=0;
end

start_time=cputime;
number_arguments=7;
if (nargin>=number_arguments)
    date_hour_in=num2str(date_hour_in);
    date_hour=[date_hour_in(1:4),'_',date_hour_in(5:6),'_',date_hour_in(7:8),'_',date_hour_in(9:10),'_',date_hour_in(11:12),'_',date_hour_in(13:14)];
else
    date_hour = clock;
    %date_hour = strcat(num2str(date_hour(1)),'_',num2str(date_hour(2)),'_',num2str(date_hour(3)),'_',num2str(date_hour(4)),'_',num2str(date_hour(5)),'_',num2str(round(date_hour(6))));
    date_hour=strcat(num2str(date_hour(1)),'_',num2str(date_hour(2),'%02i'),'_',num2str(date_hour(3),'%02i'),'_',num2str(date_hour(4),'%02i'),'_',num2str(date_hour(5),'%02i'));
end
%run_num_for_saving_figs = dlmread('fignum.txt');
%save = run_num_for_saving_figs+1;
%dlmwrite('fignum.txt',run_num_for_saving_figs+1);
% label the different parallel runs differently
% file_tag_0 gives the location for information that covers all parallel runs
if (length(parameters)>1)
    file_tag = strcat('R',num2str(run_num_for_saving_figs),'_',date_hour,'_p',num2str(ipar,'%02i'));
    file_tag_0 = strcat('R',num2str(run_num_for_saving_figs),'_',date_hour,'_p',num2str(1,'%02i'));
else
    file_tag = strcat('R',num2str(run_num_for_saving_figs),'_',date_hour);
    file_tag_0=file_tag;
end
mkdir('Saved_Research',file_tag);
mkdir('Saved_Research',file_tag_0);
file_location = ['Saved_Research/',file_tag,'/'];
file_location_0 = ['Saved_Research/',file_tag_0,'/'];
file_location_data = file_location;
file_location_data_0 = file_location_0;
file_location_figure = file_location;
file_location_movie = file_location;
% file_location_data = 'Saved_Research/Data/';
% file_location_figure = 'Saved_Research/Figures/';
% file_location_movie = 'Saved_Research/Movies/';
fprintf(' This run is labeled %s \n',file_tag);
file_name = strcat(file_location,'bulbpiri_',file_tag,'.m');
if (nargin<number_arguments)
    copyfile('bulbpiri.m',file_name);
else
    copyfile(['bulbpiri','_',date_hour_in,'.m'],file_name);
end

file_type_data = '.dat';
file_type_matlab = '.mat';
file_type_figure = '.fig';
file_type_movie = '.avi';
%will need to name every plot, strcat the name with file_tag and
%file_location_figure, then save every plot to the folder

% for parallel jobs the parameter value that is run in a given run needs to
% be written to the corresponding directory
figname='parameter'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
fid_parameter=fopen(file_name,'a');
fprintf(fid_parameter,' log_parameters = %g parameters = %g \n',log_parameters,parameters(ipar));


if (1==2)
    figname='wpm_max.dat'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
    fid_wpm_max=fopen(file_name,'a');
end

last=0;
count=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the positions and sizes of the most important windows
positions

% make a set of formatting strings for outputs
format_string=make_formatting_strings(200);

% parameters
% output
TEXT_OUT = 1;
PLOT_OUT = 1; ANIME_OUT = 1;
PLOT2D_OUT = 0; ANIME2D_OUT = 0;

% network
non_lin_M = 0;   %
non_lin_G = 0;   %
if (log_parameters==5)
    connM=parameters(ipar);
else
    connM=16;%16;%16;%8;%16; % number of connections from GC to MC (originally 4)
end
connM=round(connM*random_factor);

connP=2; % number of connections from GC to PC
% now modulate CS depending on the number of GC. the value of CS here is
% the initial value used
CS=0.000001;%0.003;%0.0004;%0.001;%0.002;%0.00025;%0.0001;%0.0002;%0.001;%0.00025;%0.0005;%0.001;%0.016; .....015.001.002% coupling strength (Originally 0.002)
CS_minimum=1e-10; % if CS goes too low the code takes very long to recover.
number_GC_aim_final=500;%2000;%1000;
number_GC_switch_on=0.7*number_GC_aim_final;
tolerance_number_GC_1=1.2;
modify_CS_frac_1=1.02;
tolerance_number_GC_2=1.5;
modify_CS_frac_2=1.1;
tolerance_number_GC_3=1.8;
modify_CS_frac_3=1.2;
log_adapt_cs=1;%1; % for all-to-all coupling oscillations arise from the feedback
number_GC_max=2*number_GC_aim_final; number_GC_t=0;
time_max_adapt_cs=Inf; % turn off adaptation of CS after this time

% control adding of GC
%Na_times =[0, 100, 120, 0.23e6, 0.26e6];
%Na_times =[0, 100, 1e7, 1e7, 1e7];
Na_times =[0, 100, 5000, 10000, 15000, 20000, 25000];
if (1==2)
    % the intended number of GC can be changed over the course of the run
    number_GC_aim_values=[0.1, 0.125, 0.25, 0.5, 1]*number_GC_aim_final;
    time_double_number_GC=Na_times([3,4,5]);
else
    %    number_GC_aim_values=[1, 1, 1, 1, 1]*number_GC_aim_final;
    number_GC_aim_values=[1, 1, 1, 1, 1, 1, 1]*number_GC_aim_final;
    % times to duplicate each GC and to decrease CS by a factor of 2
    time_double_number_GC=Inf*Na_times([3,4,5]);
end
if (log_adapt_cs==1)
    Na_values=[0, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02].*number_GC_aim_values;
    %    Na_values=[0, 0.1, 0.05, 0.02, 0.01, 0.02, 0.05].*number_GC_aim_values;
else
    Na_values=[0, 20, 5, 2, 1];  % for all-to-all coupling adding more than 1 cell is pointless and can induce collapse and oscillations
end
% adaptive coupling is only turned on after an initial phase
switch_adapt_on=0;
%gamma_values=[1,2,3,5,5];
gamma_values=5*ones(1,length(Na_times));
gamma_values=gamma_values*random_factor;
%gamma_values=[1,1,1,1,1];
% to plot sigmoids
gamma=gamma_values(1);


rm = 0; rg = 0;
wt1_initial=0;
if (log_parameters==4)
    wt1_later=parameters(ipar);
else
    wt1_later=3;%.9;%2;%3;%1.5;%1.5;%;%0.9;%0.2;%0.2;
end
wt1_later=wt1_later*random_factor;
time_wt1_1=0;
time_wt1_2=2;%100;%10000 1compare these with time = i*step
%in probe_system wt1 can be modified
wt1_later2=0.0*wt1_later;%0;

% parameter scans
% the ts_list is flipped for the first th:
% to get a decreasing order use increasing values (so that the number of
% granule cells increases rather than decreases)
% (would use unnecessarily large arrays)
if (log_parameters==2)
    ts_all=[parameters(ipar)];
else
    %ts_all=[1,2,4,8,12,16,20,40,80,100,120,140];%[.1,1,2,4];%,12,16,24,32,48,64];%[6,7,8,8.5,9,9.5,10,10.5,11,11.5,12];%[.2,.5,1,2,4];%[0.2,0.5,1,2,4,8,16,32];%[1,2,6,10,16,20,24,32,40];
    %ts_all=[24,20,16,12,8,4,2,1,0.5,0.1];%[.1,1,2,4];%,12,16,24,32,48,64];%[6,7,8,8.5,9,9.5,10,10.5,11,11.5,12];%[.2,.5,1,2,4];%[0.2,0.5,1,2,4,8,16,32];%[1,2,6,10,16,20,24,32,40];
    %ts_all=[.1,.2,.5,1,2,4,8,12];
    ts_all=[3];%[3];
end
ts_all=ts_all*random_factor;
% to plot sigmoids
ts=ts_all(1);

% for parallel runs
if (log_parameters==1)
    th_all=parameters(ipar);
else
    th_all=[0];
end

% for parallel runs
if (log_parameters==3)
    delta_x_position=parameters(ipar);
else
    delta_x_position=16;
end

[thl,tsl]=meshgrid(th_all,ts_all);
th_list=thl(:)'
fliplist=1:2:size(tsl,2);
tsl(:,fliplist)=flipud(tsl(:,fliplist));
ts_list=tsl(:)'

% to include the transient period
th_list=[th_list(1),th_list];
ts_list=[ts_list(1),ts_list];

% note figures are saved automatically based on the count of the outer steps
% not the inner steps
inner_steps_init=100;   inner_steps_final=1;
outer_param_step=10000;%350;
time_extinguish_memory=Inf;%1000;

% plot diagnostics of removed GC only after time_GC_remove_plot
time_GC_remove_plot=Inf;

% some of these parameters are reset in some of the switch cases below
time_inner_steps_switch=Inf;%outer_steps; % for time_inner_steps_switch=outer_steps no switching occurs

% old cells survive
age_old=Inf;

% control plotting of effective connectivity for nonlinear GC
% 1: only the linear connectivity of the learned stimuli at probe times
% 2: linear and nonlinear connectivity of the learned stimuli at probe times
log_connectivity_nonlinear_GC=1;

% plot histogram of resilience
log_hist_Ca=0;

% plot connectivity for different resilience levels
% 1: show connectivity of each GC, sorted by resilience and unsorted
% 2: show also effective connectivities Wmm and Wmp for ranges of
% resilience
log_selected_connectivity=0;

% Spontaneous activity of MC
S_spontaneous=0.2;%.2;%1;% 0.5;% Originally 1
S_spontaneous=S_spontaneous*random_factor;

%default value of cortical inhibition
w_inhib=0.05;%0.4;%0.1;%0.1;%0.025;%0.2;%.15;% originally 0.8

% Learning sets
pirilearn_times{1}=[0];

shuffle=zeros(1,100); log_shuffle=0;
learning_case=13;
switch learning_case
    case 1 % odor-specific inhibition KatKom12
        subcase=3;
        switch subcase
            case 1
                learning_time=cell(1,3);
                iset=1;
                learning_time{iset}=0;
                stim_type_learn{iset}=1;
                %                xpositions_learn{iset}=[[11];[13];[31];[32]];
                xpositions_learn{iset}=[[11];[13];[42];[42+delta_x_position]];
                widths_learn{iset}=[[12];[12];[12];[12]];
                amplitude_learn=2-S_spontaneous;
                amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1]];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                n_stim_learn=size(xpositions_learn{1},1);
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                
                iset=2;
                learning_time{iset}=3;
                stim_type_learn{iset}=1;
                xpositions_learn{iset}=[[11];[13];[11];[13]];
                widths_learn{iset}=[[12];[12];[12];[12]];
                amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1]];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                randomize_learn{iset}=0;
                repetition_factor{iset}=2;
                
                iset=3;
                learning_time{iset}=1000;
                %learning_time{iset}=20;
                stim_type_learn{iset}=1;
                xpositions_learn{iset}=[[11];[13];[42];[42+delta_x_position]];
                widths_learn{iset}=[[12];[12];[12];[12]];
                amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1]];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                
            case 2
                learning_time=cell(1,3);
                iset=1;
                learning_time{iset}=0;
                stim_type_learn{iset}=1;
                xpositions_learn{iset}=[[11];[31];[100];[200]];
                widths_learn{iset}=[[12];[12];[12];[12]];
                amplitude_learn=2-S_spontaneous;
                amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1]];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                n_stim_learn=size(xpositions_learn{1},1);
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                
                iset=2;
                learning_time{iset}=3;
                stim_type_learn{iset}=1;
                xpositions_learn{iset}=[[11];[200];[100];[200]];
                widths_learn{iset}=[[12];[12];[12];[12]];
                amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1]];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                
                iset=3;
                learning_time{iset}=20000;
                %learning_time{iset}=20;
                stim_type_learn{iset}=1;
                xpositions_learn{iset}=[[11];[31];[100];[200]];
                widths_learn{iset}=[[12];[12];[12];[12]];
                amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1]];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                
            case 3
                
                stim_type=1;%1;
                %                xpositions_learn_all=[[11,11];[13,13];[15,15];[17,17];[30,30];[32,32];[52,52+delta_x_position];[52,52+delta_x_position]]
                xpositions_learn_all=[[10,10];[10,10];[10,10];[25,25];[54,54+delta_x_position];[54,54+delta_x_position]]
                %                widths_learn_all=[[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12]];
                widths_learn_all=[[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12]];
                amplitude_learn=2-S_spontaneous;
                %                amplitudes_learn_all=amplitude_learn*[[0,0];[0,0];[1,0];[1,0];[1,0];[1,0];[1.2,0.8];[0.8,1.2]];
                %                amplitudes_learn_all=amplitude_learn*[[0,0];[0,0];[1,0];[1,0];[1,0];[1,0];[1.0,0.0];[0.0,1.0]];
                %                amplitudes_learn_all=amplitude_learn*[[0,0];[0,0];[1,0];[1,0];[1.2,0.8];[0.8,1.2]];
                amplitudes_learn_all=amplitude_learn*[[0,0];[0,0];[1,0];[1,0];[1.0,0.0];[0.0,1.0]];
                
                MC_per_Glom=1;
                odor_normalized=0;
                log_shuffle=1;
                shuffle=[0,0,3,4,5,6];
                %                shuffle=[0,0,0,0,0,0];
                [stimuli_learn_probe] = gara3(MC_per_Glom,stim_type,xpositions_learn_all,widths_learn_all,amplitudes_learn_all,odor_normalized);
                
                learning_time=cell(1,3);
                iset=1;
                learning_time{iset}=0;
                %                stimuli_shuffle{iset}=[stimuli_learn_probe(:,3),stimuli_learn_probe(:,4),stimuli_learn_probe(:,5),stimuli_learn_probe(:,6),stimuli_learn_probe(:,7),stimuli_learn_probe(:,8)];
                stimuli_shuffle{iset}=[stimuli_learn_probe(:,3),stimuli_learn_probe(:,4),stimuli_learn_probe(:,5),stimuli_learn_probe(:,6)];
                lights1_learn{iset}=[0,0,0,0,0,0];
                lights2_learn{iset}=[0,0,0,0,0,0];
                n_stim_learn=size(stimuli_shuffle{1},2);
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                % xpositions_learn
                xpositions_learn{iset}=[1:size(stimuli_shuffle{iset},2)];
                widths_learn{iset}=xpositions_learn{iset};
                Ns_set{iset}=size(stimuli_shuffle{iset},2);
                
                iset=2;
                learning_time{iset}=3;
                %                stimuli_shuffle{iset}=[stimuli_learn_probe(:,3),stimuli_learn_probe(:,4),stimuli_learn_probe(:,5),stimuli_learn_probe(:,6)];
                stimuli_shuffle{iset}=[stimuli_learn_probe(:,3),stimuli_learn_probe(:,4),stimuli_learn_probe(:,3),stimuli_learn_probe(:,4)];
                lights1_learn{iset}=[0,0,0,0];
                lights2_learn{iset}=[0,0,0,0];
                randomize_learn{iset}=0;
                repetition_factor{iset}=2;
                xpositions_learn{iset}=[1:size(stimuli_shuffle{iset},2)];
                widths_learn{iset}=xpositions_learn{iset};
                Ns_set{iset}=size(stimuli_shuffle{iset},2);
                
                iset=3;
                learning_time{iset}=1000;
                %                stimuli_shuffle{iset}=[stimuli_learn_probe(:,3),stimuli_learn_probe(:,4),stimuli_learn_probe(:,5),stimuli_learn_probe(:,6),stimuli_learn_probe(:,7),stimuli_learn_probe(:,8)];
                stimuli_shuffle{iset}=[stimuli_learn_probe(:,3),stimuli_learn_probe(:,4),stimuli_learn_probe(:,5),stimuli_learn_probe(:,6)];
                lights1_learn{iset}=[0,0,0,0,0,0];
                lights2_learn{iset}=[0,0,0,0,0,0];
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
                xpositions_learn{iset}=[1:size(stimuli_shuffle{iset},2)];
                widths_learn{iset}=xpositions_learn{iset};
                Ns_set{iset}=size(stimuli_shuffle{iset},2);
                
        end
        % need to check in the runs that CS is not changed during the
        % evolution that is relevant
        log_adapt_cs=1;
        time_max_adapt_cs=0.75*learning_time{3};%10000;
        CS=0.0000001;%0.00025;%0.001; % need to set suitable weight to generate a reasonable # GC
        Na_times =[0, 100, 1e7, 1e7, 1e7];
        Na_values=[0, 40, 5, 2, 1];  % for all-to-all coupling adding more than 1 cell is pointless and can induce collapse and oscillations
        
    case 2 % discrimination
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        %xpositions_learn{1}=[[24];[28];[36];[40]];
        %    xpositions_learn{1}=[[16];[20];[44];[46]];
        xpositions_learn{1}=[[28];[30];[81];[82]];
        %xpositions_learn{1}=[[32];[36];[44];[46]];
        widths_learn{1}=[[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=[[amplitude_learn];[amplitude_learn];[amplitude_learn];[amplitude_learn]];
        light_switch_learn=0;
        lights1_learn{1}=light_switch_learn*[0,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 3 %background subtraction
        iset=1;
        n_stim_learn=2;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        %xpositions_learn{1}=[[24];[28];[36];[40]];
        %xpositions_learn{1}=[[16];[20];[44];[46]];
        xpositions_learn{1}=[[9];[10];[9];[10]];
        widths_learn{1}=[[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=[[amplitude_learn];[amplitude_learn];[amplitude_learn];[amplitude_learn]];
        light_switch_learn=0;
        lights1_learn{1}=light_switch_learn*[0,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,0,0];
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 4
        iset=1;
        n_stim_learn=2;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        %xpositions_learn{1}=[[24];[28];[36];[40]];
        %xpositions_learn{1}=[[16];[20];[44];[46]];
        xpositions_learn{1}=[[10];[10];[52];[52]];
        widths_learn{1}=[[12];[12];[12];[12]];
        amplitudes_learn{1}=[[1.5];[1.5];[1.5];[1.5]];
        light_switch_learn=0;
        lights1_learn{1}=light_switch_learn*[0,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,0,0];
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 5
        iset=1;
        n_stim_learn=4;
        learning_time=cell(1,3);
        learning_time{iset}=0;
        stim_type_learn{iset}=5;
        xpositions_learn{iset}=[[12];[32];[52]];
        widths_learn{iset}=[[12];[12];[12]];
        amplitudes_learn{iset}=[[1];[1];[1]];
        lights1_learn{iset}=[0,0,0];
        lights2_learn{iset}=[0,0,0];
        randomize_learn{iset}=1;
        repetition_factor{iset}=1;
        
        iset=2;
        learning_time{iset}=3;
        stim_type_learn{iset}=5;
        xpositions_learn{iset}=[[12];[12];[12];[12]];
        widths_learn{iset}=[[12];[12];[12];[12]];
        amplitudes_learn{iset}=[[1];[1];[1];[1]];
        lights1_learn{iset}=[0,0,0,0];
        lights2_learn{iset}=[0,0,0,0];
        randomize_learn{iset}=2;
        repetition_factor{iset}=1;
        
        iset=3;
        learning_time{iset}=10000;
        stim_type_learn{iset}=5;
        xpositions_learn{iset}=[[12];[12];[32];[32]];
        widths_learn{iset}=[[12];[12];[12];[12]];
        amplitudes_learn{iset}=[[1];[1];[1];[1]];
        lights1_learn{iset}=[0,0,0,0];
        lights2_learn{iset}=[0,0,0,0];
        randomize_learn{iset}=3;
        repetition_factor{iset}=1;
        
    case 6   % context and interfering stimuli
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        light_switch_learn=1;%1;
        %xpositions_learn{1}=[[24];[28];[36];[40]];
        %xpositions_learn{1}=[[16];[20];[44];[46]];
        %xpositions_learn{1}=[[34];[36];[44];[45]];
        %xpositions_learn{1}=[[25];[26];[54];[56];[64];[65]];
        if (1==1)
            if (log_parameters==3)
                subcase=1;
                xpositions_learn{1}=[[20];[21];[parameters(ipar)];[parameters(ipar)+2];[84];[85]];
                %xpositions_learn{1}=[[parameters(ipar)];[parameters(ipar)+4];[44];[46]];
            else
                subcase=2;
                switch subcase
                    case 1
                        xpositions_learn{1}=[[25];[27];[25];[27]];
                    case 2
                        xpositions_learn{1}=[[25];[27];[75];[79]];
                end
                %xpositions_learn{1}=[[20];[21];[74];[76];[84];[85]];
                %xpositions_learn{1}=[[54];[56];[44];[46]];
                %xpositions_learn{1}=[[54];[56];[44];[46]];
            end
            widths_learn{1}=[[12];[12];[12];[12]];
            amplitude_learn=2-S_spontaneous;
            amplitudes_learn{1}=amplitude_learn*[[1];[1];[1];[1]];
            switch subcase
                case 1
                    lights1_learn{1}=light_switch_learn*[1,1,0,0];
                    lights2_learn{1}=light_switch_learn*[0,0,0,0];
                case 2
                    lights1_learn{1}=light_switch_learn*[1,1,0,0];
                    lights2_learn{1}=light_switch_learn*[0,0,1,1];
            end
        else
            xpositions_learn{1}=[[54,84];[54,86];[54,85];[56,85];[20,25];[21,25]];
            widths_learn{1}=[[12,12];[12,12];[12,12];[12,12];[12,12];[12,12]];
            amplitude_learn=2-S_spontaneous;
            amplitudes_learn{1}=amplitude_learn*[[1,0];[1,0];[0,1];[0,1];[1,0];[1,0]];
            lights1_learn{1}=light_switch_learn*[1,1,0,0,0,0];
            lights2_learn{1}=light_switch_learn*[0,0,1,1,0,0];
        end
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 7   % fast learning in cortex
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        %xpositions_learn{1}=[[24];[28];[36];[40]];
        %xpositions_learn{1}=[[38];[40];[60];[61]];
        xpositions_learn{1}=[[28];[30];[81];[82]];
        widths_learn{1}=[[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=[[amplitude_learn];[amplitude_learn];[amplitude_learn];[amplitude_learn]];
        light_switch_learn=0;
        lights1_learn{1}=light_switch_learn*[0,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 8   % extinguish memory
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        %xpositions_learn{1}=[[24];[28];[36];[40]];
        %xpositions_learn{1}=[[16];[20];[44];[46]];
        xpositions_learn{1}=[[28];[30];[80];[82]];
        widths_learn{1}=[[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=[[amplitude_learn];[amplitude_learn];[amplitude_learn];[amplitude_learn]];
        light_switch_learn=0;
        lights1_learn{1}=light_switch_learn*[0,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,0,0];
        randomize_learn{iset}=0;
        time_extinguish_memory=15000;  % need to adjust outer_param_steps to get to the desired final time
        time_inner_steps_switch=time_extinguish_memory-max(50,inner_steps_init); % time at which inner step size is switched; defines the number of inner steps before switching off memory
        %time_extinguish_memory=outer_steps-950; is now set with outer_param % number of outer steps after memory is extinguished
        n_stim_learn=size(xpositions_learn{1},1);
        log_adapt_cs=0;
        CS=0.004;%0.00025;%0.001; % need to set suitable weight to generate a reasonable # GC
        Na_values=[0, 200, 5, 2, 1];  % for all-to-all coupling adding more than 1 cell is pointless and can induce collapse and oscillations
        repetition_factor{iset}=1;
        
    case 9 % learn incrementally more odors
        learning_time=cell(1,4);
        learning_duration=5000;
        iset=1;
        learning_time{iset}=0;
        stim_type_learn{iset}=1;
        xpositions_learn{iset}=[[11];[13];[11];[13];[11];[13];[11];[13]];
        widths_learn{iset}=[[12];[12];[12];[12];[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1];[1];[1];[1];[1]];
        lights1_learn{iset}=[0,0,0,0,0,0,0,0];
        lights2_learn{iset}=[0,0,0,0,0,0,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
        iset=2;
        learning_time{iset}=learning_duration;
        stim_type_learn{iset}=1;
        xpositions_learn{iset}=[[11];[13];[11];[13];[31];[33];[31];[33]];
        widths_learn{iset}=[[12];[12];[12];[12];[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1];[1];[1];[1];[1]];
        lights1_learn{iset}=[0,0,0,0,0,0,0,0];
        lights2_learn{iset}=[0,0,0,0,0,0,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
        iset=3;
        learning_time{iset}=2*learning_duration;
        stim_type_learn{iset}=1;
        xpositions_learn{iset}=[[11];[13];[31];[33];[51];[53];[11];[13]];
        widths_learn{iset}=[[12];[12];[12];[12];[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1];[1];[1];[1];[1]];
        lights1_learn{iset}=[0,0,0,0,0,0,0,0];
        lights2_learn{iset}=[0,0,0,0,0,0,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
        iset=4;
        learning_time{iset}=3*learning_duration;
        stim_type_learn{iset}=1;
        xpositions_learn{iset}=[[11];[13];[31];[33];[51];[53];[71];[73]];
        widths_learn{iset}=[[12];[12];[12];[12];[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{iset}=amplitude_learn*[[1];[1];[1];[1];[1];[1];[1];[1]];
        lights1_learn{iset}=[0,0,0,0,0,0,0,0];
        lights2_learn{iset}=[0,0,0,0,0,0,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 10   % context and interfering stimuli
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        light_switch_learn=1;
        
        xpositions_learn{1}=[[25];[27];[39];[41];[53];[55]];
        widths_learn{1}=[[12];[12];[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=amplitude_learn*[[1];[1];[1];[1];[1];[1]];
        lights1_learn{1}=light_switch_learn*[1,1,0,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,0,0,1,1];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
    case 11   % context and interfering stimuli
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        light_switch_learn=1;
        
        xpositions_learn{1}=[[20];[35];[50];[65]];
        widths_learn{1}=[[12];[12];[12];[12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=amplitude_learn*[[1];[1];[1];[1]];
        lights1_learn{1}=light_switch_learn*[1,0,0,0];
        lights2_learn{1}=light_switch_learn*[0,1,0,0];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
        w_inhib=0.03;
        
    case 12   % context and interfering stimuli
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        light_switch_learn=1;
        
        xpositions_learn{1}=[[25,45];[23,47];[65,67];[67,67]];
        widths_learn{1}=[[12,12];[12,12];[12,12];[12,12]];
        amplitude_learn=2-S_spontaneous;
        amplitudes_learn{1}=amplitude_learn*[[1,1];[1,1];[1,0];[1,0]];
        lights1_learn{1}=light_switch_learn*[1,1,0,0];
        lights2_learn{1}=light_switch_learn*[0,0,1,1];
        n_stim_learn=size(xpositions_learn{1},1);
        randomize_learn{iset}=0;
        repetition_factor{iset}=1;
        
        %w_inhib=0.07;
        
    case 13   % context and interfering stimuli
        iset=1;
        learning_time=cell(1,1);
        learning_time{1}=0;
        stim_type_learn{1}=1;
        light_switch_learn=1;
        
        subcase=1;
        switch subcase
            case 1
                xpositions_learn{1}=[[20,50];[18,52];[85,87];[87,87]];
                widths_learn{1}=[[12,12];[12,12];[12,12];[12,12]];
                amplitude_learn=2-S_spontaneous;
                amplitudes_learn{1}=amplitude_learn*[[1,1];[1,1];[1,0];[1,0]];
                lights1_learn{1}=light_switch_learn*[1,1,0,0];
                lights2_learn{1}=light_switch_learn*[0,0,1,1];
                n_stim_learn=size(xpositions_learn{1},1);
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
            case 2
                xpositions_learn{1}=[[25,45,67];[23,47,67];[25,45,67];[23,47,67]];
                widths_learn{1}=[[12,12,12];[12,12,12];[12,12,12];[12,12,12]];
                amplitude_learn=2-S_spontaneous;
                amplitudes_learn{1}=amplitude_learn*[[0,0,1];[1,1,1];[1,1,1];[1,1,1]];
                lights1_learn{1}=light_switch_learn*[1,0,0,0];
                lights2_learn{1}=light_switch_learn*[0,1,1,1];
                n_stim_learn=size(xpositions_learn{1},1);
                randomize_learn{iset}=0;
                repetition_factor{iset}=1;
        end
        
end

n_learning_set=size(learning_time,2);
% because the length is incremented it better be set to fixed value on entry (and not spill-over from previous run)
learning_time{n_learning_set+1}=Inf;

% stimuli to test for graded response
graded_sel_S=[1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% probe set
log_graded=0;
log_cortex_fastlearn=0;
log_probe_random=1;
probe_case=learning_case;
switch probe_case
    case 1  % stimulus-specific inhibition KatKom12
        subcase=3;
        switch subcase
            case 1
                iset=1;
                probe_times{1}=learning_time{3}+[-4000,-3000,-2000,-1000,-500,-200,0,2,4,6,8,10,15,20,30,40,50,100,150,200,300,500,700,1000,1500,2000,3000,4000,5000,7000,10000,12000,15000,17000,20000,25000,30000,35000,40000,50000,60000,70000,80000,100000,120000,150000,200000];
                n_probe_times{1}=length(probe_times{1});
                stim_type_probe{1}=1;
                xpositions_probe{1}=[[32];[34];[36];[40];[44];[48];[52];[56];[13];[23]];
                widths_probe{1}=    [[12];[12];[12];[12];[12];[12];[12];[12];[12];[12]];
                amplitude_probe=2-S_spontaneous;
                amplitudes_probe{1}=amplitude_probe*[[1];[1];[1];[1];[1];[1];[1];[1];[1];[1]];
                lights1_probe{1}=[0,0,0,0,0,0,0,0,0,0];
                lights2_probe{1}=[0,0,0,0,0,0,0,0,0,0];
                randomize_probe{iset}=0;
                graded_sel_probe=[1];
            case 2
                iset=1;
                probe_times{1}=learning_time{3}+[-4000,-3000,-2000,-1000,-500,-200,0,2,4,6,8,10,15,20,30,40,50,100,150,200,300,500,700,1000,1500,2000,3000,4000,5000,7000,10000,12000,15000,17000,20000,25000,30000,35000,40000,50000,60000,70000,80000,100000,120000,150000,200000];
                n_probe_times{1}=length(probe_times{1});
                stim_type_probe{1}=1;
                xpositions_probe{1}=[[11];[13];[11];[13];[42];[42+parameters(ipar)]];
                widths_probe{1}=    [[12];[12];[12];[12];[12];[12];[12];[12];[12];[12]];
                amplitude_probe=2-S_spontaneous;
                amplitudes_probe{1}=amplitude_probe*[[0];[0];[1];[1];[1];[1]];
                lights1_probe{1}=[0,0,0,0,0,0,0,0,0,0];
                lights2_probe{1}=[0,0,0,0,0,0,0,0,0,0];
                randomize_probe{iset}=0;
                graded_sel_probe=[1];
            case 3
                iset=1;
                probe_times{1}=learning_time{3}+[-4000,-3000,-2000,-1000,-500,-200,0,2,4,6,8,10,15,20,30,40,50,100,150,200,300,500,700,1000,1500,2000,3000,4000,5000,7000,10000,12000,15000,17000,20000,25000,30000,35000,40000,50000,60000,70000,80000,100000,120000,150000,200000];
                n_probe_times{1}=length(probe_times{1});
                stim_type_probe{1}=1;
                lights1_probe{1}=[0,0,0,0,0,0,0,0,0,0];
                lights2_probe{1}=[0,0,0,0,0,0,0,0,0,0];
                randomize_probe{iset}=0;
                graded_sel_probe=[1];
                probe_shuffle{iset}=[stimuli_learn_probe(:,1),stimuli_learn_probe(:,2),stimuli_learn_probe(:,5),stimuli_learn_probe(:,6)];
                xpositions_probe{1}=[1:size(probe_shuffle{iset},2)]';
                widths_probe{1}=xpositions_probe{1};
                
        end
        
    case 2
        iset=1;
        %probe_times{1}=[1000:500:4000,4200:200:5000,5010:10:5200,5250:50:5400,5600:200:6000,6500:500:10000];
        probe_times{1}=[100,1000:1000:2*length(ts_all)*length(th_all)*1000000];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        %xpositions_probe{1}=[[32];[36];[32];[36];[32];[36];[44];[46];[44];[46];[44];[46]];
        xpositions_probe{1}=[[28];[30];[81];[82]];
        widths_probe{1}=12*[[1];[1];[1];[1]];
        %widths_probe{1}=[[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12]];
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{1}=amplitude_probe*[[1];[1];[1];[1]];
        %amplitudes_probe{1}=[[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1]];
        %lights1_probe{1}=light_switch*[1,1,0,0,0,0,0,0,0,0,1,1];
        %lights2_probe{1}=light_switch*[0,0,0,0,1,1,1,1,0,0,0,0];
        light_switch_probe=light_switch_learn;
        lights1_probe{1}=light_switch_probe*[1,1,0,0];
        lights2_probe{1}=light_switch_probe*[0,0,0,0];
        randomize_probe{iset}=0;
        graded_sel_probe=[1];
        
    case 3
        iset=1;
        probe_times{1}=2000:2000:2*length(ts_all)*length(th_all)*100000;
        n_probe_times{1}=length(probe_times);
        stim_type_probe{1}=1;
        x1=10; x2=15;
        xp1=[[x1,x2];[x1,x2];[x1,x2];[x1,x2];[x1+5,x2];[x1+5,x2];[x1+10,x2];[x1+10,x2];[x1+15,x2];[x1+15,x2];[x1+20,x2];[x1+20,x2];[x1+25,x2];[x1+25,x2];[x1+35,x2];[x1+35,x2]];
        lxp1=size(xp1,1);
        x1=100; x2=95;
        xp2=[[x1,x2];[x1,x2];[x1,x2];[x1,x2];[x1-5,x2];[x1-5,x2];[x1-10,x2];[x1-10,x2];[x1-15,x2];[x1-15,x2];[x1-20,x2];[x1-20,x2];[x1-25,x2];[x1-25,x2];[x1-35,x2];[x1-35,x2]];
        lxp2=size(xp2,1);
        x1=10; x2=15;
        xp3=[[x1,x2+10];[x1,x2+10];[x1,x2+20];[x1,x2+20]];
        lxp3=size(xp3,1);
        x1=100; x2=95;
        xp4=[[x1,x2-10];[x1,x2-10];[x1,x2-20];[x1,x2-20]];
        lxp4=size(xp4,1);
        xpositions_probe{1}=[xp1;xp2;xp3;xp4];
        widths_probe{1}=repmat([12,2],size(xpositions_probe{1},1),1);
        %widths_probe{1}=[[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2];[12,2]];
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{1}=repmat([[amplitude_probe,0.5];[amplitude_probe,0]],size(xpositions_probe{1},1)/2,1);
        amplitudes_probe{1}(1,1)=0;
        amplitudes_probe{1}(2,1)=0;
        amplitudes_probe{1}(size(xp1,1)+1,1)=0;
        amplitudes_probe{1}(size(xp1,1)+2,1)=0;
        %amplitudes_probe{1}=[[0,0.5];[0,0];[amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0];...
        %    [amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0];...
        %    [0,0.5];[0,0];[amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0];[amplitude_probe,0.5];[amplitude_probe,0]];
        lights1_probe{1}=0*size(xpositions_probe{1},1);
        lights2_probe{1}=0*size(xpositions_probe{1},1);
        if (sum(size(xpositions_probe{1})==size(widths_probe{1}))+sum(size(xpositions_probe{1})==size(amplitudes_probe{1}))<4)
            size(xpositions_probe{1})
            size(widths_probe{1})
            size(amplitudes_probe{1})
            fprintf(' \n\n ******* Error in probing stimuli !! *******  \n\n')
            return
        end
        randomize_probe{iset}=0;
        graded_sel_probe=[1];
        
    case 4
        iset=1;
        probe_times{1}=[500:500:300000];
        n_probe_times{1}=length(probe_times);
        stim_type_probe{1}=1;
        xpositions_probe{1}=[[10,15];[10,15];[10,15];[10,15];[15,15];[15,15];[20,15];[20,15];[25,15];[25,15];[30,15];[30,15];[35,15];[35,15];[35,30];[35,30];[35,30];[35,30];[10,30];[10,30]];
        widths_probe{1}=[[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3];[12,3]];
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{1}=[[0,0.15];[0,0];[amplitude_probe,0.15];[amplitude_probe,0];[amplitude_probe,0.15];[amplitude_probe,0];[amplitude_probe,0.15];[amplitude_probe,0];[amplitude_probe,0.15];[amplitude_probe,0];[amplitude_probe,0.15];[amplitude_probe,0];[amplitude_probe,0.15];[amplitude_probe,0];[amplitude_probe,0.15];[amplitude_probe,0];[0,0.15];[0,0];[amplitude_probe,0.15];[amplitude_probe,0]];
        lights1_probe{1}=0*size(xpositions_probe{1},1);
        lights2_probe{1}=0*size(xpositions_probe{1},1);
        size(amplitudes_probe{1})
        if (sum(size(xpositions_probe{1})==size(widths_probe{1}))+sum(size(xpositions_probe{1})==size(amplitudes_probe{1}))<4)
            size(xpositions_probe{1})
            size(widths_probe{1})
            size(amplitudes_probe{1})
            fprintf(' \n\n ******* Error in probing stimuli !! *******  \n\n')
            return
        end
        randomize_probe{iset}=4;
        graded_sel_probe=[1];
        
    case 5 % for random patterns stored earlier
        iset=1;
        %probe_times{1}=[500:500:3000,3010:10:3200,3250:50:4000,4500:500:6000];
        probe_times{1}=[500:500:9000,9200:200:10000,10010:10:10100,10150:50:10500,10600:100:10900,11000:200:15000,15500:500:20000];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=5;
        xpositions_probe{1}=[[32];[37];[42];[47];[52];[57]];
        widths_probe{1}=[[12];[12];[12];[12];[12];[12]];
        amplitudes_probe{1}=[[1];[1];[1];[1];[1];[1]];
        lights1_probe{1}=0*size(xpositions_probe{1},1);
        lights2_probe{1}=0*size(xpositions_probe{1},1);
        randomize_probe{iset}=4;
        graded_sel_probe=[1];
        
    case 6   % context and interference
        subcase=3;
        iset=1;
        probe_times{1}=[100,2000:2000:2*length(ts_all)*length(th_all)*1e6];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        switch subcase
            case 1
                if (log_parameters==3)
                    xpos1=parameters(ipar);
                    xpos2=parameters(ipar)+2;
                    xpositions_probe{1}=[[xpos1];[xpos2];[xpos1];[xpos2];[xpos1];[xpos2];[84];[85];[84];[85];[84];[85];[20];[21]];
                    %xpositions_probe{1}=[[xpos1];[xpos2];[xpos1];[xpos2];[xpos1];[xpos2];[44];[46];[44];[46];[44];[46];[25];[26]];
                else
                    xpositions_probe{1}=[[74];[76];[74];[76];[74];[76];[84];[85];[84];[85];[84];[85];[20];[21]];
                    %xpositions_probe{1}=[[54];[56];[54];[56];[54];[56];[44];[46];[44];[46];[44];[46];[25];[26]];
                end
                widths_probe{1}=    [[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12]];
                amplitude_probe=2-S_spontaneous;
                amplitudes_probe{1}=amplitude_probe*[[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1]];
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*[1,1,0,0,0,0,0,0,1,1,0,0,0,0];
                lights2_probe{1}=light_switch_probe*[0,0,1,1,0,0,1,1,0,0,0,0,0,0];
                
            case 2
                xpositions_probe{1}=[[54,84];[54,86];[54,84];[54,86];[54,84];[54,86];[54,84];[56,86];[54,84];[56,86];[54,84];[56,86];[20,25];[22,25]];
                widths_probe{1}=[[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12];[12,12]];
                amplitude_probe=2-S_spontaneous;
                amplitudes_probe{1}=amplitude_probe*[[1,.4];[1,.4];[1,.4];[1,.4];[1,.4];[1,.4];[0,.4];[0,.4];[0,.4];[0,.4];[0,.4];[0,.4];[1,.0];[1,.0]];
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*[1,1,0,0,0,0,0,0,1,1,0,0,0,0];
                lights2_probe{1}=light_switch_probe*[0,0,1,1,0,0,1,1,0,0,0,0,0,0];
                
            case 3
                xpos_detection=[repmat([26,30],8,1);repmat([26,71],8,1);repmat([75,71],4,1);repmat([75,30],4,1)];
                xpos_discrimination=[repmat([[26,30];[26,32]],4,1);repmat([[26,75];[26,79]],4,1)];
                xpositions_probe{1}=[xpos_detection;xpos_discrimination];
                widths_probe{1}=[repmat([12,2],size(xpos_detection,1)+8,1);repmat([12,12],8,1)];
                amplitude_probe=2-S_spontaneous;
                blip=0.4;
                amplitudes_probe{1}=amplitude_probe*[repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],4,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1)];
                
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*repmat([1,1,0,0],1,size(xpositions_probe{1},1)/4);
                lights2_probe{1}=light_switch_probe*0*xpositions_probe{1};
                
            case 4
                xpos_detection=[repmat([26,30],8,1);repmat([26,71],8,1);repmat([75,71],4,1);repmat([75,30],4,1)];
                xpos_discrimination=[repmat([[26,30];[26,32]],4,1);repmat([[26,75];[26,79]],4,1);repmat([[26,48];[26,52]],4,1)];
                xpositions_probe{1}=[xpos_detection;xpos_discrimination];
                widths_probe{1}=[repmat([12,2],size(xpos_detection,1)+8,1);repmat([12,12],16,1)];
                amplitude_probe=2-S_spontaneous;
                blip=0.4;
                amplitudes_probe{1}=amplitude_probe*[repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],4,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1)];
                
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*repmat([1,1,0,0],1,size(xpositions_probe{1},1)/4);
                lights2_probe{1}=light_switch_probe*0*xpositions_probe{1};
                
        end
        
        randomize_probe{iset}=0;
        graded_sel_probe=[1,2,3];
        
    case 7   % fast learning in cortex
        log_cortex_fastlearn=1;
        iset=1;
        probe_times{1}=[100,500:500:2*length(ts_all)*length(th_all)*100000];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        %xpositions_probe{1}=[repmat([39,60.5],4,1);repmat([38,40],2,1);repmat([60,61],2,1)];
        %xpositions_probe{1}=[repmat([29,81.5],4,1);repmat([28,30],2,1);repmat([81,82],2,1)];
        xpositions_probe{1}=[repmat([29,81.5],6,1);repmat([81,82],2,1)];
        widths_probe{1}=    [repmat([12,12],size(xpositions_probe{1},1),1)];
        amplitude_probe=2-S_spontaneous;
        %amplitudes_probe{1}=[[.4*amplitude_probe,.6*amplitude_probe];[.6*amplitude_probe,.4*amplitude_probe];[.45*amplitude_probe,.55*amplitude_probe];[.55*amplitude_probe,.45*amplitude_probe];[amplitude_probe,0];[0,amplitude_probe];[amplitude_probe,0];[0,amplitude_probe]];
        %amplitudes_probe{1}=[[1.05*amplitude_probe,0.95*amplitude_probe];[0.95*amplitude_probe,1.05*amplitude_probe];[.97*amplitude_probe,1.03*amplitude_probe];[1.03*amplitude_probe,.97*amplitude_probe];[amplitude_probe,0];[0,amplitude_probe];[amplitude_probe,0];[0,amplitude_probe]];
        amplitudes_probe{1}=[[1.05*amplitude_probe,0.95*amplitude_probe];[0.95*amplitude_probe,1.05*amplitude_probe];[1.03*amplitude_probe,0.97*amplitude_probe];[0.97*amplitude_probe,1.03*amplitude_probe];[1.10*amplitude_probe,0.90*amplitude_probe];[0.9*amplitude_probe,1.1*amplitude_probe];[amplitude_probe,0];[0,amplitude_probe]];
        lights1_probe{1}=[0,0,0,0,0,0,0,0];
        lights2_probe{1}=[0,0,0,0,0,0,0,0];
        randomize_probe{iset}=0;
        graded_sel_probe=[1,3];
        
        
    case 8   % extinguish one memory
        iset=1;
        probe_times{1}=time_extinguish_memory+[-100:10:-10,0:1:9,10:5:45,50:10:90,100:50:450,500:500:2000];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        %xpositions_probe{1}=[[32];[36];[32];[36];[32];[36];[44];[46];[44];[46];[44];[46]];
        xpositions_probe{1}=xpositions_learn{1};
        widths_probe{1}=widths_learn{1};
        amplitudes_probe{1}=amplitudes_learn{1};
        light_switch_probe=light_switch_learn;
        lights1_probe{1}=light_switch_probe*[1,1,0,0];
        lights2_probe{1}=light_switch_probe*[0,0,0,0];
        randomize_probe{iset}=0;
        graded_sel_probe=[1];
        % plot diagnostics of removed GC only after time_GC_remove_plot
        time_GC_remove_plot=probe_times{1}(1);
        
    case 9    %learn multiple odors
        iset=1;
        probes=20;
        probe_times{1}=[ceil(learning_duration/probes):ceil(learning_duration/probes):4*learning_duration];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        xpositions_probe{iset}=[[11];[13];[31];[33];[51];[53];[71];[73]];
        widths_probe{iset}=[[12];[12];[12];[12];[12];[12];[12];[12]];
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{iset}=amplitude_probe*[[1];[1];[1];[1];[1];[1];[1];[1]];
        lights1_probe{iset}=[0,0,0,0,0,0,0,0];
        lights2_probe{iset}=[0,0,0,0,0,0,0,0];
        randomize_probe{iset}=0;
        graded_sel_probe=[1];
        
    case 10   % context and interference
        subcase=2;
        iset=1;
        probe_times{1}=[100,500:500:2*length(ts_all)*length(th_all)*1e6];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        switch subcase
            case 1
                if (log_parameters==3)
                    xpos1=parameters(ipar);
                    xpos2=parameters(ipar)+2;
                    xpositions_probe{1}=[[xpos1];[xpos2];[xpos1];[xpos2];[xpos1];[xpos2];[84];[85];[84];[85];[84];[85];[20];[21]];
                    %xpositions_probe{1}=[[xpos1];[xpos2];[xpos1];[xpos2];[xpos1];[xpos2];[44];[46];[44];[46];[44];[46];[25];[26]];
                else
                    xpositions_probe{1}=[[74];[76];[74];[76];[74];[76];[84];[85];[84];[85];[84];[85];[20];[21]];
                    %xpositions_probe{1}=[[54];[56];[54];[56];[54];[56];[44];[46];[44];[46];[44];[46];[25];[26]];
                end
                widths_probe{1}=    [[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12];[12]];
                amplitude_probe=2-S_spontaneous;
                amplitudes_probe{1}=amplitude_probe*[[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1];[1]];
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*[1,1,0,0,0,0,0,0,1,1,0,0,0,0];
                lights2_probe{1}=light_switch_probe*[0,0,1,1,0,0,1,1,0,0,0,0,0,0];
                
            case 2
                xpositions_probe{1}=[[27,0,47];[27,0,47];[27,0,47];[27,37,47];[27,37,47];[27,37,47]];
                widths_probe{1}=    [[12,12,12];[12,12,12];[12,12,12];[12,12,12];[12,12,12];[12,12,12]];
                amplitude_probe=2-S_spontaneous;
                amplitudes_probe{1}=amplitude_probe*[[1,0,1];[1,0,1];[1,0,1];[1,1,1];[1,1,1];[1,1,1]];
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*[0,1,1,0,1,1];
                lights2_probe{1}=light_switch_probe*[0,0,1,0,0,1];
                
            case 3
                xpos_detection=[repmat([26,30],8,1);repmat([26,71],8,1);repmat([75,71],4,1);repmat([75,30],4,1)];
                xpos_discrimination=[repmat([[26,30];[26,32]],4,1);repmat([[26,75];[26,79]],4,1)];
                xpositions_probe{1}=[xpos_detection;xpos_discrimination];
                widths_probe{1}=[repmat([12,2],size(xpos_detection,1)+8,1);repmat([12,12],8,1)];
                amplitude_probe=2-S_spontaneous;
                blip=0.4;
                amplitudes_probe{1}=amplitude_probe*[repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],4,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1)];
                
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*repmat([1,1,0,0],1,size(xpositions_probe{1},1)/4);
                lights2_probe{1}=light_switch_probe*0*xpositions_probe{1};
                
            case 4
                xpos_detection=[repmat([26,30],8,1);repmat([26,71],8,1);repmat([75,71],4,1);repmat([75,30],4,1)];
                xpos_discrimination=[repmat([[26,30];[26,32]],4,1);repmat([[26,75];[26,79]],4,1);repmat([[26,48];[26,52]],4,1)];
                xpositions_probe{1}=[xpos_detection;xpos_discrimination];
                widths_probe{1}=[repmat([12,2],size(xpos_detection,1)+8,1);repmat([12,12],16,1)];
                amplitude_probe=2-S_spontaneous;
                blip=0.4;
                amplitudes_probe{1}=amplitude_probe*[repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],2,1);
                    repmat([[0,blip];[0,.0]],2,1);
                    repmat([[1,blip];[1,0]],4,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1);
                    repmat([[1,blip];[1,blip]],2,1);
                    repmat([[0,blip];[0,blip]],2,1)];
                
                light_switch_probe=light_switch_learn;
                lights1_probe{1}=light_switch_probe*repmat([1,1,0,0],1,size(xpositions_probe{1},1)/4);
                lights2_probe{1}=light_switch_probe*0*xpositions_probe{1};
                
        end
        
        randomize_probe{iset}=0;
        graded_sel_probe=[1,2,3];
        
    case 11   % context and interference
        iset=1;
        probe_times{1}=[100,500:500:2*length(ts_all)*length(th_all)*1e6];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        
        xpositions_probe{1}=[[28,50];[28,50];[28,50];[28,50];[28,50];[28,50]];
        widths_probe{1}=    [[12,12];[12,12];[12,12];[12,12];[12,12];[12,12]];
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{1}=amplitude_probe*[[1,1];[1,1];[1,1];[1,0];[1,0];[1,0]];
        light_switch_probe=light_switch_learn;
        lights1_probe{1}=light_switch_probe*[0,1,0,0,1,0];
        lights2_probe{1}=light_switch_probe*[0,0,1,0,0,1];
        
        randomize_probe{iset}=0;
        graded_sel_probe=[1,2,3];
        
    case 12   % context and interference
        iset=1;
        probe_times{1}=[100,300,500,1000,2000:2000:2*length(ts_all)*length(th_all)*1e6];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        
        xpositions_probe{1}=[repmat([[34,65];[36,65]],6,1);repmat([35,65],12,1)];
        widths_probe{1}=    repmat([12,12],24,1);
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{1}=amplitude_probe*[repmat([1,0],6,1);repmat([1,1],6,1);repmat([[1,0];[0,0]],3,1);repmat([[1,1];[0,1]],3,1)];
        light_switch_probe=light_switch_learn;
        lights1_probe{1}=light_switch_probe*repmat([0,0,1,1,0,0],1,4);
        lights2_probe{1}=light_switch_probe*repmat([0,0,0,0,1,1],1,4);
        
        randomize_probe{iset}=0;
        graded_sel_probe=[1,2,3];
        
    case 13   % context and interference
        iset=1;
        probe_times{1}=[100,300,500,1000,2000,5000:5000:2*length(ts_all)*length(th_all)*1e6];
        n_probe_times{1}=length(probe_times{1});
        stim_type_probe{1}=1;
        
        blip=0.4;
        xpositions_probe{1}=[repmat([[33,86];[37,86]],6,1);repmat([35,86],12,1)];
        widths_probe{1}=    [repmat([12,12],12,1);repmat([2,12],12,1)];
        amplitude_probe=2-S_spontaneous;
        amplitudes_probe{1}=amplitude_probe*[repmat([0.6,0],6,1);repmat([0.6,1],6,1);repmat([[blip,0];[0,0]],3,1);repmat([[blip,1];[0,1]],3,1)];
        light_switch_probe=light_switch_learn;
        lights1_probe{1}=light_switch_probe*repmat([0,0,1,1,0,0],1,4);
        lights2_probe{1}=light_switch_probe*repmat([0,0,0,0,1,1],1,4);
        
        randomize_probe{iset}=0;
        graded_sel_probe=[1,2,3];
        
end

n_probe_set=size(stim_type_probe{1},2);
i_state=0;
state_probe=cell(n_probe_times{1},12);

if isfinite(time_extinguish_memory)
    % for extinguishing memory (case 8):
    % outer_param_step=number of outer steps
    % with inner_steps_init number of inner steps + number of outer steps with
    % inner_steps_final steps per outer step:
    outer_param_step=ceil(time_extinguish_memory/inner_steps_init)+1000;
    fprintf('Adjusting number of steps to reflect switching of inner steps: outer_param_step = %g \n',outer_param_step);
end
outer_transient=0;
outer_steps=outer_transient+outer_param_step*(length(th_list)-1);
steps_param_switch=outer_transient:outer_param_step:outer_steps;

% learn only every nskip_learn steps
%nskip_learn=ceil(.8*outer_steps);%1000000;...it gets compared to i so it needs to be less than integer
%nskip_learn=ceil(inner_steps_switch+(outer_steps-inner_steps_switch)/2)
nskip_learn=10*outer_steps;
%nskip_learn=ceil(outer_steps/2)
%1000000;...it gets compared to i so it needs to be less than integer
%that divides exp_time in the definition of step


% normalization of the odor amplitude
odor_normalized=0;

plot_cell_activity_every_time_its_calculated=1; % acutally modded via plot_cell_activity_every_time_its_calculated_mod_init
plot_cell_activity_every_time_its_calculated_mod_init=100;
plot_cell_activity_every_time_its_calculated_mod_final=1;
fisher_plot=5;
fisher_average_time=5000;   %averaging for fisher_final in terms of the time
fisher_average=ceil(fisher_average_time/plot_cell_activity_every_time_its_calculated_mod_init);

% stimulus (from old version)
sim = 40;
odor_names = 'limonene(+)_ster limonene(_)_ster propylpropionate_es3 ethylbutyrate_es3 isopropylbenzene_ModuleC1 cyclohexanone_SG18 acetone methylacetate_SG19 cycloheptanelow_cycloalk propanol_simp_2500 isoamylbutyrate_est1 butyricacid_aci1 hexanal_ald1 ethylbenzene_HC';
choose = 1:n_stim_learn;%8%In this modified version there are only four stimuli
len = 18;

% number of time steps for which the state is to be saved for each
% parameter value
number_states_saved=5;
state=cell(length(ts_list),number_states_saved,6);
%since first and last time are included in the saving
number_states_saved=number_states_saved-1;

save_video = 0;
token123 = 0;
mix_factor = .9;
counter_plot = 0;
% stepping
cont_density = 0;
% change the way the timing is defined
if (1==2)
    exp_time =100;%5000;%25000%5000
    step = round(exp_time/50);%200%50
    dt = 5;
    istep=0;
end

% controls saving of figures (this used to count outer steps, not inner steps, now it counts time)
number_save=5;
mod_save=ceil(outer_param_step/number_save)*inner_steps_init;%1000;
mod_save_pick=0; %that is the remainder for which the figures are saved
mod_save_pick=min(mod_save_pick,mod_save-1);

dt=1;
%exp_time=dt*inner_steps*outer_steps;
if (isfinite(time_inner_steps_switch))
    exp_time=dt*(time_inner_steps_switch+inner_steps_final*(outer_steps-ceil(time_inner_steps_switch/inner_steps_init)));
else
    exp_time=dt*inner_steps_init*outer_steps;
end
%step=inner_steps*dt;

istep=0;
istep_last=0;
%i goes on a for loop from 1:exp_time/step which is equal to whatever
%exp_timeis divided by above. step ==125 currently
%l_interval is compared with i*step. if mod(l_interval,i*step)==0 the memory is updated

% r_interval = 50;
% tracking vs T
tracking = 0; % 0: none, 1: in S
track_pairs = [choose]; %;5,6;7,8];


% adjust the histogram of GC activities after steps_GC_max_switch outer
% steps
steps_GC_max_switch=10;

% parameters of sigmoid in piri -> Sigmoid(p)
% for odor cells in cortex
thresh_odor=0.2;%1.7;%0.6;%1.2;%1.7;%.7
P_amp_odor=1; steep_odor=3;%8;%4;%10  %original thresh 0.7, original steep 20
% for light cells in cort5x
thresh_light=0.2;%1.7;%0.6;%1.2;%1.7;%.7
P_amp_light=1; steep_light=3;%8;%4;%10  %original thresh 0.7, original steep 20

steep_odor=steep_odor*random_factor;
steep_light=steep_odor;
thresh_odor=thresh_odor*random_factor;
thresh_light=thresh_odor;
P_amp_odor=P_amp_odor*random_factor;
P_amp_light=P_amp_odor;

Gthresh=3;%0;%10;%.7
Gthresh=Gthresh*random_factor;
G_amp=7.5; Gsteep=1;%10  %original thresh 0.7, original steep 20
%2.5
Mthresh=0;%.7
M_amp=3.5; Msteep=3;%10  %original thresh 0.7, original steep 20

% multiple random read-out weights for fisher_non_opt
n_random_sample=5000;

%W_PM_scale multiplies the connection matrix between p and m cells after
%learning has been done
% p_scale divides p cell activity calculated in cal_activity1 to scale it
% down
P_scale = 1;%1.5;

W_PM_scale=1;%10;%3;%4;%1;%1.5%4;%1.5; %4;%6;
% ...
Sstr = 1; MC_per_Glom = 1;
TYPE = 1; % 1: pearson 2: L2
% mode of making connections
% 0: random attachment, 3: a version of preferential attachment
prob_conn=3;
% preferential attachment
attach_prefer=0.0;

% threshold in Fisher discriminant (to avoid dominance by rates near 0)
fisher_thresh=0;%0.4;

perm_ratio = 0; % only for cont_density = 0
minv = 0; maxv = 1;

% PC stimulus set-up
first_time = 1;
not_first_time = 0;

fprintf(fid_parameter,'  %g %g \n',log_parameters,parameters(ipar));
fprintf(fid_parameter,'  %g %g %g %g %g %g %g %g %g %g ',th_all);
fprintf(fid_parameter,' \n');
fprintf(fid_parameter,'  %g %g %g %g %g %g %g %g %g %g ',ts_all);
fprintf(fid_parameter,' \n');
fprintf(fid_parameter,'  %g %g \n',wt1_later,Gthresh);
fprintf(fid_parameter,'  %g %g \n',learning_case,probe_case);
fprintf(fid_parameter,' --- \n');
fprintf(fid_parameter,' log_parameters = %g parameters = %g \n',log_parameters,parameters(ipar));
fprintf(fid_parameter,' th_all = %g %g %g %g %g %g %g %g %g %g \n',th_all);
fprintf(fid_parameter,' ts_all = %g %g %g %g %g %g %g %g %g %g \n',ts_all);
fprintf(fid_parameter,' wt1_later = %g Gthresh = %g \n',wt1_later,Gthresh);
fprintf(fid_parameter,' \n learning_case = %g probe_case = %g \n',learning_case,probe_case);
fprintf(fid_parameter,' seed = %g first random number = %g \n',seed,first_random);
fclose(fid_parameter);

% setup stimulus
if TEXT_OUT == 1
    fprintf('\n--- Start ---\n   - initializing odor\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%S1 and S are only different in the strength of the light signal. They use
%different weights for scaling the light signal, so the signal can be
%learned with the light signal in a higher range, and then passed to the
%rest of the code in the range that the smell signal will be suppressed to.

%S_light=back*ones(size(S_light));
num_light_cells=48;%16;%30; must be even
light_cells_1=1:num_light_cells/2; light_cells_2=num_light_cells/2+1:num_light_cells;
light_back = -1;%-.25;
light_on = 2;%1.5;

%[Sall, coord,~, name] = gara2(14, MC_per_Glom, odor_names, sim);% This, ideally, is where I would put the memory

% Stimuli
% stim_type=1: Gaussian
% stim_type=2: quadratic
% stim_type=3: Gaussian, different mixtures
% stim_type=4: Gaussian, backgrounds
S_set=cell(1,n_learning_set);
Sall=cell(1,n_learning_set);
for learning_set=1:n_learning_set
    if (log_shuffle==0)
        stim_type=stim_type_learn{learning_set};
        randomize=randomize_learn{learning_set};
        [Sall{learning_set}] = gara3(MC_per_Glom,stim_type,xpositions_learn{learning_set},widths_learn{learning_set},amplitudes_learn{learning_set},odor_normalized);% This, ideally, is where I would put the memory
    else
        Sall{learning_set}=stimuli_shuffle{learning_set};
    end
    S_smell = S_spontaneous+Sstr*(Sall{learning_set});
    
    % with the variables lights one can turn different lights on for the individual
    % stimuli
    % During learning
    %lights_1=[0,0,0,0];
    %lights_2=[0,0,0,0];
    S_light_off=light_back*ones(num_light_cells,size(S_smell,2));
    S_light_learn = S_light_off;
    S_light_learn(light_cells_1,lights1_learn{learning_set}>0)=light_on;
    S_light_learn(light_cells_2,lights2_learn{learning_set}>0)=light_on;
    S_set{learning_set}=[S_smell;S_light_learn];
end
S=S_set{1}; i_learning_set=1;

n_P=N_mitral+num_light_cells;
%if (N_mitral+num_light_cells~=n_P)
%    fprintf(' cell numbers incompatible with identity M->P N_mitral= %g num_light_cell = %g n_P = %g \n',N_mitral,num_light_cells,n_P);
%    return
%end

% Probing the Fisher discrimnant
% probe_fisher=1: wt1=0 and wt1>0
% probe_fisher=2: with light_learn and light_present1
% probe_fisher=3: with light_learn and light_present1 and light_present2
probe_fisher=1;
% I could replace odor_names with an image of a memory, like in the files.
% Better yet, the piri model could create a memory and that memory
% would replace the array or matrix that forms the odor in this code.
% Mitral cells would provide the "learning stimuli" and the memory
% would be updated in every step.

% During presentation 1
lights_1=[0,0,0,0];%[0,0,1,1];
lights_2=[0,0,0,0];%[1,1,0,0];
S_light_present_1 = S_light_off;
S_light_present_1(light_cells_1,lights_1>0)=light_on;
S_light_present_1(light_cells_2,lights_2>0)=light_on;
% During presentation 2
lights_1=[0,0,0,0];%[0,0,0,0];
lights_2=[0,0,0,0];%[0,0,0,0];
S_light_present_2 = S_light_off;
S_light_present_2(light_cells_1,lights_1>0)=light_on;
S_light_present_2(light_cells_2,lights_2>0)=light_on;

S_present_1=[S_smell;S_light_present_1];
S_present_2=[S_smell;S_light_present_2];
S_probe = S;
if (1==2)
    probe = probe+S_spontaneous;
    %probe = [probe;S_light_present(:,choose)];
    probe = [probe;S_light_learn];
end

probe_set=cell(1,n_probe_set); probe_all=cell(1,n_probe_set);
for i_probe_set=1:n_probe_set
    if (log_shuffle==0)
        stim_type=stim_type_probe{i_probe_set};
        randomize=randomize_probe{i_probe_set};
        [probe_all{i_probe_set}] = gara3(MC_per_Glom,stim_type,xpositions_probe{i_probe_set},widths_probe{i_probe_set},amplitudes_probe{i_probe_set},odor_normalized);
    else
        probe_all{i_probe_set}=probe_shuffle{i_probe_set};
    end
    % This, ideally, is where I would put the memory
    
    S_probe = S_spontaneous+Sstr*(probe_all{i_probe_set});
    
    S_light_off=light_back*ones(num_light_cells,size(S_probe,2));
    S_light_probe = S_light_off;
    S_light_probe(light_cells_1,lights1_probe{i_probe_set}>0)=light_on;
    S_light_probe(light_cells_2,lights2_probe{i_probe_set}>0)=light_on;
    probe_set{i_probe_set}=[S_probe;S_light_probe];
end

% select the odors for which the Fisher discriminant is calculated
%S_fisher=[probe(:,1),probe(:,2),probe(:,3),probe(:,4)];
S_fisher=[S_set{1}(:,1),S_set{1}(:,2),S_set{1}(:,3),S_set{1}(:,4)];
% in fisher_diagnostic the odor pairs to be compared are defined via
% compare: adjacent pairs are compared
compare=[1,2,3,4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M_learn = [ S(:,2),S(:,3)];
%Pass M_LEarn S1 because S1  is scaled to the correct size for learning. S
% is scaled to the correct size for running the rest of the code. Later
% M_learn calls might use S, if the smell signal has been suppressed.
%M_learn = [S1(:,choose)];%,S1(:,5),S1(:,6),S1(:,7),S1(:,8)];
M_learn = [S(:,1:n_stim_learn)];%,S1(:,5),S1(:,6),S1(:,7),S1(:,8)];
%M_learn = [.9*S1(:,1:2)+1.1*S1(:,3:4),1.1*S1(:,1:2)+0.9*S1(:,3:4)];%,S1(:,5),S1(:,6),S1(:,7),S1(:,8)];
% learn possibly only specific odors
learning=ones(1,size(M_learn,2));
%learning(3:4)=0;

% for testing the activity of the granule cells that have just been removed
S_GC_probe=S;

% read-out weights for fisher_non_opt
read_out_weights=2*rand(n_random_sample,N_mitral)-1;
%read_out_weights=2*ones(n_random_sample,N_mitral)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coord(2,1:size(S,1)) = 1:size(S,1);
coord(1,:) = 1;
metric =  zeros(size(S,1),size(S,1));
for i = 1:size(S,1)
    metric(i,:) = abs((1:size(S,1))-i);
end

% figures to be saved
list_save_figures=[3,15,16,17,18,19,20,21,22,30,31,32,48,49,50,51,52,80,81,82,83,84,90,91];
list_name_figures={'GC_activities','Wmm_Wmp_nonlinear','Wmg_Wmg_Wmp_Wmp','Smem_M_P_17','fisher_corr_18','fisher_corr_19','fisher_corr_20','Smem_M_P_21',...
    'Smem_M_P_22','GC_removed','response_GC_removed','response_GC_distribution',...
    'Smem_M_P_48','Smem_M_P_49','Smem_M_P_50','Smem_M_P_51','Smem_M_P_52','fisher_probe_cortex_a','fisher_probe_nocortex_a',...
    'connectivity_selective','resilience_histo','Wmg_Wpg_sorted','fisher_probe_cortex_b','fisher_probe_nocortex_b'};
fig_label(1:max(list_save_figures))=NaN;

% Nc (and with it n_P) includes mitral cells and num_light_cells `light cells', introduced via S_light_learn
Nc=N_mitral+num_light_cells; %XH
%N_mitral=size(S_smell,1); XH
%Ns = length(choose);
Ns=size(S,2);%n_stim_learn;

% with N_mitral known define now the sigmoid parameters
% P_amp=[P_amp_odor*ones(N_mitral,1);P_amp_light*ones(num_light_cells,1)];
% steep=[steep_odor*ones(N_mitral,1);steep_light*ones(num_light_cells,1)];
% thresh=[thresh_odor*ones(N_mitral,1);thresh_light*ones(num_light_cells,1)];

S_name = cell(1,size(S,2));

for i = 1:size(S,2)
    %    S_name{i} = name{choose(i)};% = name{i}
    S_name{i} = 'odor';
end

if TEXT_OUT == 1
    fprintf('       - Nc = %d\n', Nc);
    fprintf('           corr = %f\n', mean_excluNaN(uptri_1d(corr(S,TYPE))));
end

if PLOT_OUT == 1;
    if (1==2)
        setup_Pplot(N_mitral,S,corr(S_smell,TYPE),corr(S_smell',TYPE),rg,1);
        %pplot figure(1)
        drawnow;
        figname = 'P_channel_corr_pattern_corr';
        file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
        %    saveas(figure(1),file_name,file_type_figure)
        hgsave(file_name)
        %pcell activity?
    end
end


if PLOT2D_OUT == 1;
    setup_Pplot2D(S,coord,101,S_name);
    %this is figure(101)
    drawnow;
    figname = 'Smell_Pairs';
    file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
    %saveas(figure(101),file_name,file_type_figure)
    hgsave(file_name)
end


%Correlation to start
%R1 = corrcoef([S(:,1),S(:,2)]);
R1 = corr([S_smell(:,1),S_smell(:,2)],TYPE);

%Make piri learn and read stimuli
%I'll instead run piri once and save the memory to a text file

%It might be throwing errors now because W_PP doesn't expect just 4
%neurons.

%m1 = max(S1(:,2));
%m2 = max(S1(:,3));

% Parameters for piri
% The following inhibition is independent of the inhibition used when
% learning stimuli - when I used the same low stimulus for both functions,
% the decorrelation went haywire.
% log_cross_inhib regulates whether the odor and the light cortex parts
w_inhib=w_inhib*random_factor;
% inhibit each other
log_cross_inhib=1;
W_PP_i=w_inhib*ones(n_P,n_P);
W_PP_i=W_PP_i-diag(diag(W_PP_i));
if (log_cross_inhib==0)
    W_PP_i(1:n_P-num_light_cells,n_P-num_light_cells+1:n_P) = 0; %XH
    W_PP_i(n_P-num_light_cells+1:n_P,1:n_P-num_light_cells) = 0; %XH
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_PP_max=0.06;%.2;%.1;%.2;%.1;%.2;%.5;
W_PP_max=W_PP_max*random_factor;
W_PP=NaN(size(W_PP_i));
W_PM=NaN(n_P,Nc);
c_assoc_learn=0.95; c_assoc_recall=0;
tau_P = 1;

if (1==2)
    figure(111)
    m3 = max(max(S));
    imagesc(pirirecall(S/m3)); colorbar
    title(' recall activities with M=S')
    figname = 'Normalized_P_Cell_Activity_after_Learning';
    file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
    saveas(gcf,file_name,file_type_figure)
end

% dlmwrite('Winter/odortext/S.txt', S);
% m3 = max(max(S));%Normalise the recall matrix
% S_mem = pirirecall(S/m3);

%     figure(6)
%     subplot(1,2,1)
%     imagesc(S_mem)
%     subplot(1,2,2)
%     imagesc(S)

%     fprintf('   - learning odors\n')
%     S_learn = [Sall(:,11)];
%     S_recall = Sall;
%     S_mem = piri(S_learn, S_recall);
%     fprintf('   - memories generated\n')

% setup network
if TEXT_OUT == 1
    fprintf('   - initializing network\n');
end

% Xuchen 8/13: What is Isize?
Isize = 1*Nc; N = ones(Isize,1);% N is a column vector filled with ones that is Isize long
%Isize = Nc/2; N = ones(Isize,1);% N is a column vector filled with ones that is Isize long
%Nc is 26 so Isize is 52
%Is 52*26 the size of the odor?
Wmg = zeros(Nc, Isize);
Wgm = zeros(Isize, Nc);
% HR Wpg should be dimensioned with n_P not Nc
Wpg = zeros(n_P, Isize);%XH
Iage = -ones(1,Isize);% So we know that we haven't initialised any cells yet
Imark = zeros(1,Isize);
%XH 4/23
%I = NaN(Isize,length(choose));
I = NaN(Isize,Ns);

%time_axis = 0:step:exp_time;% Does this make a list?
time_axis = 0:1:exp_time;% this time_axis definition is not taking into account the different values of inner_steps. Needs to be fixed.
N_t = NaN*ones(length(time_axis),Isize);% What does all this NaN stuff do or mean?
Pcorr_t = NaN*ones(1,length(time_axis));
Tcorr_t = NaN*ones(size(track_pairs,1),length(time_axis));
Pangle_t = NaN*ones(1,length(time_axis));
Tangle_t = NaN*ones(1,length(time_axis));
Pfoc_t = NaN*ones(1,length(time_axis));
CV_t = NaN*ones(1,length(time_axis));
CVid_t = NaN*ones(Ns,length(time_axis));
F_t = NaN*ones(Ns,Ns,length(time_axis));
%G_removed=cell(step);
G_removed=cell(outer_steps);
time_istep=NaN(1,outer_steps);

%Xuchen 8/12: Not used
max_G_removed=1;
for iodor=1:Ns
    response{iodor}=NaN;
end
%%%%%%%%%%%%%%%%%%%


c_assoc=c_assoc_recall;
log_first_run=1;
[P, I, I1, I2, S_mem, ~] = cal_activity1(CS,Wpg,Wmg,Wgm,S,S,0);%Extra argument for S_mem
log_first_run=1;

if cont_density == 1
    N_t(1,:) = N;
end


if PLOT_OUT == 1;
    % This is where figure 2 comes from
    %    HP1d = setup_Pplot(P,corr(P(1:Nc/2,:),TYPE),corr(P(1:Nc/2,:)',TYPE),rg,2,Wmg,Wgm,Wpg,time_axis,N_t);
    if (1==2)
        HP1d = setup_Pplot(N_mitral,P,corr(P(1:N_mitral,:),TYPE),corr(P(1:N_mitral,:)',TYPE),rg,2,Wmg,Wgm);
        % How many G cells connect to P cells i and j? How does S_mem change?
        %PCinfo = P_connections(6,Wpg,S_mem,time_axis,N_t);
    end
    [HI, Iaxis] = setup_Iplot(cont_density,time_axis,I,I1,I2,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(1,:),3);
    VST = cell(1,size(track_pairs,1)+1);
    VST{1} = Pcorr_t;
    for i_pairs=1:size(track_pairs,1)
        VST{1+i_pairs} = Tcorr_t(i_pairs,:);
    end
    line_style = cell(1,1+size(track_pairs,1));
    line_style{1} = '-k'; line_style{2} = '-r'; line_style{3} = '-g'; line_style{4} = '-b'; line_style{5} = '-m';
    line_name = cell(1,1+size(track_pairs,1));
    line_name{1} = Pcorr_t; line_name{2} = Tcorr_t; line_name{3} = Tcorr_t;
    % This is where correlations are plotted
    if(1==2)
        Hinfo = setup_Infoplot(time_axis,VST,line_style,line_name,corr(S,TYPE),corr(P,TYPE),4);
    end
    drawnow;
end

if PLOT2D_OUT == 1;
    setup_Pplot2D(P,coord,102,S_name);
    drawnow;
end


% step
% Should include some soame way Wgm and Wmg
% are developed?
% The DiffEq for recall doesn't take too long to solve...

if TEXT_OUT == 1
    fprintf('   - running\n');
end

%Pcorr_t(1) = mean_excluNaN(uptri_1d(corr(P(1:Nc/2,:),TYPE)));
Pcorr_t(1) = mean_excluNaN(uptri_1d(corr(P(1:N_mitral,:),TYPE)));
%Pangle_t(1) = mean_excluNaN(uptri_1d(corr_angle(corr(S(1:Nc/2,:),TYPE),corr(P(1:Nc/2,:),TYPE))));
Pangle_t(1) = mean_excluNaN(uptri_1d(corr_angle(corr(S(1:N_mitral,:),TYPE),corr(P(1:N_mitral,:),TYPE))));
Pfoc_t(1) = mean_excluNaN(focality(P,metric));
CV_t(1) = std(mean(P,2))/mean(mean(P,2));
if (1==2)
    CVid_t(:,1) = std(P)./mean(P);
    F_t(:,:,1) = corr(P);
end

if tracking == 1
    %Tcorr_t(1) = mean_excluNaN(cal_track_corr(track_pairs,P));
    Tcorr_t(1:size(track_pairs,1),1) = cal_track_corr(track_pairs,P)';
    Tangle_t(1) = mean_excluNaN(corr_angle(cal_track_corr(track_pairs,S),cal_track_corr(track_pairs,P)));
end

if TEXT_OUT == 1
    fprintf('           corr = %f\n', Pcorr_t(1));
end


if save_video == 1
    
    figname = 'MC_SMEM_RATIO';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    MC_SMEM_RATIO_VIDEO = VideoWriter(file_name);
    open(MC_SMEM_RATIO_VIDEO);
    
    figname = 'Cell_Activity_with_odor_with_light';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    Cell_Activity_with_odor_with_light_video = VideoWriter(file_name);
    open(Cell_Activity_with_odor_with_light_video);
    
    figname = 'Cell_Activity_with_odor_no_light';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    Cell_Activity_with_odor_no_light_video = VideoWriter(file_name);
    open(Cell_Activity_with_odor_no_light_video);
    
    figname = 'Cell_Activity_no_odor_with_light';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    Cell_Activity_no_odor_with_light_video = VideoWriter(file_name);
    open(Cell_Activity_no_odor_with_light_video);
    
    figname = 'Cell_Activity_no_odor_no_light';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    Cell_Activity_no_odor_no_light_video = VideoWriter(file_name);
    open(Cell_Activity_no_odor_no_light_video);
    
    figname = 'Correlation_between_g_Cell_activities';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    Correlation_between_g_Cell_activities_video = VideoWriter(file_name);
    open(Correlation_between_g_Cell_activities_video);
    
    figname = 'GC_odor2_ol_nol';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    GC_odor2_ol_nol_video = VideoWriter(file_name);
    open(GC_odor2_ol_nol_video);
    
    figname = 'GC_odor6_ol_nol';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    GC_odor6_ol_nol_video = VideoWriter(file_name);
    open(GC_odor6_ol_nol_video);
    
    figname = 'new_corr_track';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    new_corr_track_video = VideoWriter(file_name);
    open(new_corr_track_video);
    
    figname = 'WmgxWgm_WmgxWgp';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    WmgxWgm_WmgxWgp_video = VideoWriter(file_name);
    open(WmgxWgm_WmgxWgp_video);
    
    figname = 'Smem_M_P';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    Smem_M_P_video = VideoWriter(file_name);
    open(Smem_M_P_video);
    
    figname = 'fisher_cortical_input';
    file_name = strcat(file_location_movie,figname,'_',file_tag,file_type_movie);
    fisher_cortical_input_video = VideoWriter(file_name);
    open(fisher_cortical_input_video);
    
end

if (1==2)
    [First1,First1_comp] = fisher_rect(P(1:N_mitral,1),P(1:N_mitral,2),fisher_thresh);
    [First2,First2_comp] = fisher_rect(P(1:N_mitral,3),P(1:N_mitral,4),fisher_thresh);
    
    fprintf('First1 =  %g  \n',First1)
    fprintf('First2 =  %g  \n',First2)
    figure(100)
    plot(1:length(First1_comp),First1_comp,1:length(First1_comp),First1_comp_thresh,1:length(First2_comp),First2_comp,1:length(First2_comp),First2_comp_thresh)
    xlabel(' MC Number')
    ylabel(' Contribution to Fisher')
end

plot_sigmoid

time=0;
i_after_switch=0;
fisher_steps=10000;
time_fisher=NaN(fisher_steps,1);
y_fisher=NaN(fisher_steps,22);
y_fisher_1=NaN(fisher_steps,12);
y_fisher_2=NaN(fisher_steps,12);
number_GC=NaN(fisher_steps,1);
number_GC_bins=40;
CS_time=NaN(fisher_steps,1);
GC_hist=zeros(number_GC_bins,size(S,2),fisher_steps);
GC_hist_cum=zeros(number_GC_bins,size(S,2),fisher_steps);
figname='GC_hist'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
fid_hist=fopen(file_name,'w');
figname='GC_hist_cum'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
fid_hist_cum=fopen(file_name,'w');

figname='parameter'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
fid_parameter=fopen(file_name,'a');
fprintf(fid_parameter,'S_spontaneous = %g Gthresh = %g ts = %g wt1_later = %g \n',S_spontaneous,Gthresh,ts,wt1_later);
fprintf(fid_parameter,' connM = %g gamma = %g W_PP_max = %g w_inhib = %g \n',connM,gamma,W_PP_max,w_inhib);
fprintf(fid_parameter,' P_amp_odor = %g steep_odor = %g thresh_odor = %g \n',P_amp_odor,steep_odor,thresh_odor);
fclose(fid_parameter);
figname='parameter_parallel'; file_name = strcat(file_location_data_0,figname,'_',file_tag_0,file_type_data);
file_name
fid_parameter=fopen(file_name,'a');
if (ipar==1)
    fprintf(fid_parameter,'S_spontaneous Gthresh ts wt1_later connM gamma W_PP_max w_inhib P_amp_odor steep_odor thresh_odor \n')
end
fprintf(fid_parameter,' %3.0f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n',ipar,S_spontaneous,Gthresh,ts,wt1_later,connM,gamma,W_PP_max,w_inhib,P_amp_odor,steep_odor,thresh_odor);
fclose(fid_parameter);

iparam=1;
ts=ts_list(iparam);
th=th_list(iparam);
log_first_run = 1;

for i = 1:outer_steps
    if (i>1)
        if (max(max(abs(P)))>1e2)
            figname='parameter'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
            fid_parameter=fopen(file_name,'a');
            fprintf(fid_parameter,' P too large %g \n\n',max(max(abs(P))));
            fprintf(fid_parameter,' %g %g %g %g %g %g %g %g %g %g \n',P)
            fclose(fid_parameter);
            figname='parameter_parallel'; file_name = strcat(file_location_data_0,figname,'_',file_tag_0,file_type_data);
            figname
            ls
            fid_parameter=fopen(file_name,'a');
            fprintf(fid_parameter,' in run %g P too large %g \n\n',ipar,max(max(abs(P))));
            fclose(fid_parameter);
            break
        end
    end
    i_after_switch=i_after_switch+1;
    
    if (time<time_inner_steps_switch)
        inner_steps=inner_steps_init;
        plot_cell_activity_every_time_its_calculated_mod=plot_cell_activity_every_time_its_calculated_mod_init;
    else
        inner_steps=inner_steps_final;
        plot_cell_activity_every_time_its_calculated_mod=plot_cell_activity_every_time_its_calculated_mod_final;
    end
    
    if (count>0)
        figname='fisher_time';
        file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
        fid_fisher=fopen(file_name,'w');
        if (probe_fisher==1 || probe_fisher==2 || probe_fisher==4)
            fprintf(fid_fisher,format_string{19},...
                [time_fisher(1:count)';y_fisher(1:count,[3,4,7:12])';y_fisher_1(1:count,[3,4,7:12])';number_GC(1:count)'/1000;CS_time(1:count)']);
        elseif (probe_fisher==3)
            fprintf(fid_fisher,format_string{27},...
                [time_fisher(1:count)';y_fisher(1:count,[3,4,7:12])';y_fisher_1(1:count,[3,4,7:12])';y_fisher_2(1:count,[3,4,7:12])';number_GC(1:count)'/1000;CS_time(1:count)']);
        end
        fclose(fid_fisher);
        
        figname='fisher_non_opt_time';
        file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
        fid_fisher=fopen(file_name,'w');
        if (probe_fisher==1 || probe_fisher==2 || probe_fisher==4)
            fprintf(fid_fisher,format_string{21},...
                [time_fisher(1:count)';y_fisher(1:count,[13:22])';y_fisher_1(1:count,[13:22])']);
        elseif (probe_fisher==3)
            fprintf(fid_fisher,format_string{31},...
                [time_fisher(1:count)';y_fisher(1:count,[13:22])';y_fisher_1(1:count,[13:22])';y_fisher_2(1:count,[13:22])']);
        end
        fclose(fid_fisher);
    end
    
    if (i==steps_param_switch(iparam))
        % output of fisher_diagnostic
        % 3,4=max and min of SigmoidM
        % 7,8=Fisher with SigmoidM;
        % 9,10=correlation SigmoidM
        % 11,12=max and min of M
        % 13,14 Fisher mean and std pair 1
        % 15,16,17 Fisher quantiles
        % 18,19 Fisher mean and std pair 2
        % 20,21,22 Fisher quantiles
        figname='fisher_final';
        file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
        fid_fisher=fopen(file_name,'a');
        averaging_range=[max(count-fisher_average,1):count];
        if (probe_fisher==1 || probe_fisher==2 || probe_fisher==4)
            fprintf(fid_fisher,format_string{16},...
                ts,th,mean(y_fisher(averaging_range,[3,4,7:10]),1),mean(y_fisher_1(averaging_range,[3,4,7:10]),1),number_GC(count)/1000,CS_time(count));
        elseif (probe_fisher==3)
            fprintf(fid_fisher,format_string{22},...
                ts,th,mean(y_fisher(averaging_range,[3,4,7:10]),1),mean(y_fisher_1(averaging_range,[3,4,7:10]),1),mean(y_fisher_2(averaging_range,[3,4,7:10]),1),number_GC(count)/1000,CS_time(count));
        end
        fclose(fid_fisher);
        figname='fisher_non_opt_time';
        file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
        fid_fisher=fopen(file_name,'w');
        if (probe_fisher==1 || probe_fisher==2 || probe_fisher==4)
            fprintf(fid_fisher,format_string{21},...
                [time_fisher(1:count)';y_fisher(1:count,[13:22])';y_fisher_1(1:count,[13:22])']);
        elseif (probe_fisher==3)
            fprintf(fid_fisher,format_string{31},...
                [time_fisher(1:count)';y_fisher(1:count,[13:22])';y_fisher_1(1:count,[13:22])';y_fisher_2(1:count,[13:22])']);
        end
        fclose(fid_fisher);
        figname='fisher_non_opt_final';
        file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
        fid_fisher=fopen(file_name,'a');
        if (probe_fisher==1 || probe_fisher==2 || probe_fisher==4)
            fprintf(fid_fisher,format_string{22},...
                ts,th,mean(y_fisher(averaging_range,[13:22]),1),mean(y_fisher_1(averaging_range,[13:22]),1));
        elseif (probe_fisher==3)
            fprintf(fid_fisher,format_string{32},...
                ts,th,mean(y_fisher(averaging_range,[13:22]),1),mean(y_fisher_1(averaging_range,[13:22]),1),mean(y_fisher_2(averaging_range,[13:22]),1));
        end
        fclose(fid_fisher);
    end
    
    state_save_end=steps_param_switch(iparam);
    state_save_begin=max(steps_param_switch(iparam)-number_states_saved,1);
    if (i>=state_save_begin && i<=state_save_end)
        state{iparam,i-state_save_begin+1,1}=[th,ts,time];
        state{iparam,i-state_save_begin+1,2}=P;
        state{iparam,i-state_save_begin+1,3}=I;
        state{iparam,i-state_save_begin+1,4}=I1;
        state{iparam,i-state_save_begin+1,5}=I2;
        state{iparam,i-state_save_begin+1,6}=P_recall;
        figname='states';
        file_name = strcat(file_location_data,figname,'_',file_tag,file_type_matlab);
        save(file_name,'state');
    end
    
    if (i>steps_param_switch(iparam))
        iparam=iparam+1;
        ts=ts_list(iparam);
        th=th_list(iparam);
        fprintf(' i = %g ts = %g th = %g \n',i,ts,th)
        if (th~=th_list(iparam-1))
            figname='fisher_final';
            file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
            fid_fisher=fopen(file_name,'a');
            fprintf(fid_fisher,' & \n');
            fclose(fid_fisher);
        end
        i_after_switch=0;
    end
    
    % bins for histograms of GC activity are adapted to the threshold th
    if (i<steps_GC_max_switch)
        GC_active_max=10;
    elseif (i==steps_GC_max_switch)
        GC_active_max=max(max(I));
        fprintf(' Adjusting GC_active_max = %g \n',GC_active_max);
    end
    GC_bins=linspace(GC_active_max/number_GC_bins/2,GC_active_max*(1-1/number_GC_bins/2),number_GC_bins);
    
    istep=istep+1;
    
    %    if i == round(exp_time/step)
    if i == outer_steps
        last = 1;
    end
    %
    %     if (time<step_Na_turn_on)
    %         Na=0;
    %     elseif (time<step_Na_switch)
    %         Na=Na_initial;
    %     else
    %         Na=Na_final;
    %     end
    
    Na_switch=time>=Na_times;
    Na_index=find(Na_switch==1);
    Na=Na_values(Na_index(end));
    gamma=gamma_values(Na_index(end));
    number_GC_aim=number_GC_aim_values(Na_index(end));
    
    %adhoc realearn
    %     if  i >1 && i<100
    %         pirilearn(M_learn,1,file_location_figure,file_tag,file_type_figure,S)
    %     end
    %     if (mod(i,nskip_learn)==0)
    %         pirilearn(M_learn,1,file_location_figure,file_tag,istep,label_connectivity)
    %     end
    
    if (ismember(time,pirilearn_times{1}))
        
        if (i==1)
            
            label_connectivity='a';
            pirilearn(M_learn,0,file_location_figure,file_tag,istep,label_connectivity)
            
            if (1==2)
                c_assoc=c_assoc_recall;
                odor_light_diagnostic(N_mitral,num_light_cells,S_smell,S_spontaneous,light_on_learn,light_back,CS,Wpg,Wmg,Wgm)
                c_assoc=c_assoc_learn;
            end
            
            W_PM(1:n_P,1:N_mitral) = W_PM_scale*W_PM(1:n_P,1:N_mitral); % to take into account that the mitral cell activities are strongly reduced with inhibition
            W_PP_all=W_PP;   % later also a modified cortical connectivity with extinguished memory may be used
            W_PP_1=W_PP; W_PM_1=W_PM;
            
            if (log_cortex_fastlearn==1)
                
                % generate second cortical connectivity for probing
                probe_learn=(probe_set{1}(:,1)+probe_set{1}(:,2));
                figure(100000)
                imagesc(probe_learn)
                colorbar
                
                label_connectivity='b';
                pirilearn(probe_learn,0,file_location_figure,file_tag,istep,label_connectivity)
                W_PM(1:n_P,1:N_mitral) = W_PM_scale*W_PM(1:n_P,1:N_mitral); % to take into account that the mitral cell activities are strongly reduced with inhibition
                W_PP_all=W_PP;   % later also a modified cortical connectivity with extinguished memory may be used
                W_PP_2=W_PP; W_PM_2=W_PM;
                
            end
            % use the first cortical connectivity for neurogenesis
            W_PP=W_PP_1; W_PM=W_PM_1;
            
        else
            
            fprintf(' piri connectivity not computed again, but read in from file \n')
            fprintf(' press any key to continue \n')
            pause
            
            % Read in Matrices and inhibition
            W_PP = dlmread('odortext/WPP.txt');
            W_PM = dlmread('odortext/WPM.txt');
            
        end
    end
    
    if (time==time_extinguish_memory)
        
        if (1==2)
            learning(3:4)=0;
            pirilearn(M_learn,1,file_location_figure,file_tag,istep,label_connectivity)
        end
        
        if (1==2)
            fprintf('Extinguishing the second memory by hand \n');
            W_PP(32:64,32:64)=0;
            figure(2000)
            imagesc(W_PP)
        end
        
        if (1==1)
            fprintf('Extinguishing the first memory by hand \n');
            figure(7000)
            subplot(2,1,1)
            imagesc(W_PP)
            xmiddle=(xpositions_learn{1}(1)+xpositions_learn{1}(2))/2;
            xdelta=0.8*widths_learn{1}(1);
            W_PP(round(xmiddle-xdelta):round(xmiddle+xdelta),round(xmiddle-xdelta):round(xmiddle+xdelta))=0;
            subplot(2,1,2)
            imagesc(W_PP)
            caxis([0,W_PP_max])
        end
    end
    
    %time = i*step;
    
    if (1==2)
        fig6=figure(6); set(fig6,'Position',pos{6});
        set(gca,'NextPlot','replaceChildren');
        axis tight
        set(gcf,'Renderer','zbuffer');
        
        h1=subplot(2,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
        imagesc(P); axis square; colorbar
        title('MC')
        h1=subplot(2,2,2); ax1 = get(h1,'position'); set(h1,'position', ax1);
        imagesc(S_mem); axis square; colorbar
        title('S_{mem}')
        h1=subplot(2,2,3); ax1 = get(h1,'position'); set(h1,'position', ax1);
        imagesc(S_mem./P); axis square;
        colorbar
        title('Ratio')
        if save_video == 1
            if last == 0
                frame =getframe(6);
                writeVideo(MC_SMEM_RATIO_VIDEO, frame);
            else
                figname = 'MC_SMEM_RATIO_VIDEO';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                saveas(gcf,file_name,file_type_figure)
                frame =getframe(6);
                writeVideo(MC_SMEM_RATIO_VIDEO, frame);
                
                close(MC_SMEM_RATIO_VIDEO);
            end
        end
    end
    if TEXT_OUT == 1
        if (count>0)
            %            fprintf('  run num = %d, time = %f fisher = %g number_GC = %g \n',run_num,i*step,y_fisher(count,7),number_GC(count));
            %            fprintf('  run num = %d, time = %f fisher = %g number_GC = %g \n',run_num,time,y_fisher(count,7),number_GC(count));
            if probe_fisher == 3
                fprintf('  time = %f i = %g fisher = %g fisher1 = %g fisher2 = %g number_GC = %g \n',time,i,y_fisher(count,7),y_fisher_1(count,7),y_fisher_2(count,7),number_GC(count));
            else
                fprintf('  time = %f i = %g fisher = %g number_GC = %g \n',time,i,y_fisher(count,7),number_GC(count));
            end
        else
            %            fprintf('  run num = %d, time = %f \n',run_num,i*step);
            fprintf('  run num = %d, time = %f \n',run_num,time);
            
        end
    end
    
    if cont_density == 1
        [~,N] = ode23(@RHS,[0 step],N_t(i,:),option);
        N_t(i+1,:) = N(end,:);
        Wmg = C*diag(N(end,:)); Wgm = C';
    else
        %        for j = 1:round(step/dt)
        for j = 1:inner_steps
            %time=istep*dt;
            time=time+dt;
            
            if (log_adapt_cs==1&&time<time_max_adapt_cs)
                number_GC_t=length(find(Iage>=0));
                if (number_GC_t>number_GC_switch_on)
                    %switch_adapt_on stays on once it has been turned on
                    %after the initial grwoing phase of the network
                    switch_adapt_on=1;
                end
                if (Na>0&&switch_adapt_on==1)
                    if (number_GC_t>number_GC_aim*tolerance_number_GC_1)
                        CS=CS*modify_CS_frac_1;
                        fprintf(' time = %g number_GC_t = %g  CS = %g  \n',time,number_GC_t,CS);
                        %                         prob_probe=survival(I);
                        %                         prob_probe_min=min(prob_probe);
                        %                         prob_probe_max=max(prob_probe);
                        %                         fprintf(' time = %g number_GC_t = %g  CS = %g prob_min=%g prob_max=%g \n',time,number_GC_t,CS,prob_probe_min,prob_probe_max);
                    end
                    if (number_GC_t<number_GC_aim/tolerance_number_GC_1)
                        CS=max(CS/modify_CS_frac_1,CS_minimum);
                        %                         prob_probe=survival(I);
                        %                         prob_probe_min=min(prob_probe);
                        %                         prob_probe_max=max(prob_probe);
                        fprintf(' time = %g number_GC_t = %g  CS = %g \n',time,number_GC_t,CS);
                    end
                    if (number_GC_t>number_GC_aim*tolerance_number_GC_2)
                        CS=CS*modify_CS_frac_2;
                        fprintf(' time = %g number_GC_t = %g  CS = %g \n',time,number_GC_t,CS);
                    end
                    if (number_GC_t<number_GC_aim/tolerance_number_GC_2)
                        CS=max(CS/modify_CS_frac_2,CS_minimum);
                        fprintf(' time = %g number_GC_t = %g  CS = %g \n',time,number_GC_t,CS);
                    end
                    if (number_GC_t>number_GC_aim*tolerance_number_GC_3)
                        CS=CS*modify_CS_frac_3;
                        fprintf(' time = %g number_GC_t = %g  CS = %g \n',time,number_GC_t,CS);
                    end
                    if (number_GC_t<number_GC_aim/tolerance_number_GC_3)
                        CS=max(CS/modify_CS_frac_3,CS_minimum);
                        fprintf(' time = %g number_GC_t = %g  CS = %g \n',time,number_GC_t,CS);
                    end
                end
            end
            
            if ismember(time,time_double_number_GC)
                double_number_GC
                CS=CS/2;
                update_Iplot(cont_density,time_axis,I,I1,I2,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(i,:),HI);
                fprintf(' Doubled GC number: CS = %g \n',CS);
            end
            
            if (number_GC_t<number_GC_max||log_adapt_cs==0)
                add_cell;
            else
                fprintf(' Too many GC: no GC added\n');
            end
            
            if (time>learning_time{i_learning_set+1})
                i_learning_set=i_learning_set+1;
                S=S_set{i_learning_set};
                Ns=Ns_set{i_learning_set};
                log_first_run=1;
                fprintf(' \n now learning_set = %g time = %g \n',i_learning_set,time);
            end
            
            c_assoc=c_assoc_recall;
            [P, I, I1, I2, S_mem, P_recall] = cal_activity1(CS,Wpg,Wmg,Wgm,S,P,time);%Extra argument for S_mem and Wpg
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            counter_plot = counter_plot+1;
            %            if plot_cell_activity_every_time_its_calculated==1  && mod(counter_plot,plot_cell_activity_every_time_its_calculated_mod) == 1
            %%%%%%%%%%
            %%%%%%%%% the Fisher diagsnotics should probably be done after remove_cell
            %%%%%%%%%%
            if plot_cell_activity_every_time_its_calculated==1  && mod(counter_plot,plot_cell_activity_every_time_its_calculated_mod) == 0
                
                if (mod(count,fisher_plot)==0)
                    plot_S_M_P(17,S_mem,P,P_recall)
                end
                
                log_first_run=1;
                [P_temp,~,~,~,S_mem_temp,P_recall_temp] = cal_activity1(CS,Wpg,Wmg,Wgm,S_fisher,P,time);%Extra argument for S_mem and Wpg
                log_first_run=1;
                
                [time_fisher,y_fisher,count]=fisher_diagnostic(time,P_temp,compare,fisher_thresh,count,time_fisher,y_fisher,18);
                number_GC(count)=length(find(Iage>=0));
                CS_time(count)=CS*1000;
                
                if (learning_case==8)
                    [GC_hist,GC_hist_cum]=plot_GC_hist(count,I,GC_bins,GC_hist,GC_hist_cum);
                end
                
                
                if (time==time_extinguish_memory+4 && count >=6 || i==outer_steps && time_extinguish_memory < exp_time)
                    for i_odor=1:size(GC_hist,2)
                        fprintf(fid_hist,' %g %g %g %g %g %g %g \n',[1:size(GC_hist,1);squeeze(GC_hist(:,i_odor,count-5))';squeeze(GC_hist(:,i_odor,count-4))';squeeze(GC_hist(:,i_odor,count-3))';squeeze(GC_hist(:,i_odor,count-2))';squeeze(GC_hist(:,i_odor,count-1))';squeeze(GC_hist(:,i_odor,count))']);
                        fprintf(fid_hist,' & odor = %g i = %g GC_active_max = %g \n',i_odor,i,GC_active_max);
                    end
                    for i_odor=1:size(GC_hist_cum,2)
                        fprintf(fid_hist_cum,' %g %g %g %g %g %g %g \n',[1:size(GC_hist_cum,1);squeeze(GC_hist_cum(:,i_odor,count-5))';squeeze(GC_hist_cum(:,i_odor,count-4))';squeeze(GC_hist_cum(:,i_odor,count-3))';squeeze(GC_hist_cum(:,i_odor,count-2))';squeeze(GC_hist_cum(:,i_odor,count-1))';squeeze(GC_hist_cum(:,i_odor,count))']);
                        fprintf(fid_hist_cum,' & odor = %g i = %g GC_active_max = %g \n',i_odor,i,GC_active_max);
                    end
                end
                
                if (probe_fisher==1)
                    wt1_initial_tmp=wt1_initial; wt1_later_tmp=wt1_later;
                    wt1_initial=0; wt1_later=0;
                    log_first_run=1;
                    [P_temp,~,~,~,~,~] = cal_activity1(CS,Wpg,Wmg,Wgm,S_fisher,P,time);%Extra argument for S_mem and Wpg
                    wt1_initial=wt1_initial_tmp; wt1_later=wt1_later_tmp;
                    log_first_run=1;
                    
                    [~,y_fisher_1,~]=fisher_diagnostic(time,P_temp,compare,fisher_thresh,-count,time_fisher,y_fisher_1,19);
                    
                elseif (probe_fisher==2)
                    
                    log_first_run=1;
                    [P_temp,~,~,~,S_mem_temp,P_recall_temp] = cal_activity1(CS,Wpg,Wmg,Wgm,S_present_1,P,time);%Extra argument for S_mem and Wpg
                    plot_S_M_P(21,S_mem_temp,P_temp,P_recall_temp)
                    log_first_run=1;
                    
                    [~,y_fisher_1,~]=fisher_diagnostic(time,P_temp,compare,fisher_thresh,-count,time_fisher,y_fisher_1,19);
                    
                elseif (probe_fisher==3)
                    
                    log_first_run=1;
                    [P_temp,~,~,~,S_mem_temp,P_recall_temp] = cal_activity1(CS,Wpg,Wmg,Wgm,S_present_1,P,time);%Extra argument for S_mem and Wpg
                    log_first_run=1;
                    plot_S_M_P(21,S_mem_temp,P_temp,P_recall_temp)
                    
                    [~,y_fisher_1,~]=fisher_diagnostic(time,P_temp,compare,fisher_thresh,-count,time_fisher,y_fisher_1,19);
                    
                    log_first_run=1;
                    [P_temp,~,~,~,S_mem_temp,P_recall_temp] = cal_activity1(CS,Wpg,Wmg,Wgm,S_present_2,P,time);%Extra argument for S_mem and Wpg
                    log_first_run=1;
                    plot_S_M_P(22,S_mem_temp,P_temp,P_recall_temp)
                    
                    [~,y_fisher_2,~]=fisher_diagnostic(time,P_temp,compare,fisher_thresh,-count,time_fisher,y_fisher_2,20);
                    
                elseif (probe_fisher==4)
                    W_PP_tmp=W_PP;
                    %fprintf('Extinguishing the first memory by hand \n');
                    W_PP(1:32,1:32)=0;
                    
                    log_first_run=1;
                    [P_temp,~,~,~,S_mem_temp,P_recall_temp] = cal_activity1(CS,Wpg,Wmg,Wgm,S_present_1,P,time);%Extra argument for S_mem and Wpg
                    log_first_run=1;
                    plot_S_M_P(21,S_mem_temp,P_temp,P_recall_temp)
                    
                    [~,y_fisher_1,~]=fisher_diagnostic(time,P_temp,compare,fisher_thresh,-count,time_fisher,y_fisher_1,19);
                    W_PP=W_PP_tmp;
                    
                end
                
                if save_video == 1
                    
                    frame =getframe(17);
                    writeVideo(Smem_M_P_video, frame);
                    
                    frame =getframe(18);
                    writeVideo(fisher_cortical_input_video, frame);
                end
                
            end
            
            remove_cell;
            
            if (1==2)
                Wmp_max=max(max(Wmg*Wpg'));
                fprintf(' time =%g Wmp_max = %g \n',time,Wmp_max);
                fprintf(fid_wpm_max,' %g %g \n',time,Wmp_max);
                if (Wmp_max<20 && time>200)
                    plot_S_M_P(17,S_mem,P,P_recall)
                    probe_system(time,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)
                    pause
                end
            end
            
            % probe response with respect to arbitrary other stimuli
            if (sum(time==ceil(probe_times{1}))>0||CS==CS_minimum&&i==steps_param_switch(iparam))
                fprintf(' Probing system at time = %g \n',time)
                label_connectivity='a';
                probe_system(time,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)
                if (log_cortex_fastlearn==1)
                    W_PP=W_PP_2; W_PM=W_PM_2;
                    fprintf(' Probing with second connectivity at time = %g \n',time)
                    label_connectivity='b';
                    probe_system(time+0.1,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)
                    W_PP=W_PP_1; W_PM=W_PM_1;
                end
                fprintf(' Saving Figures \n')
                for ifig=1:length(list_save_figures)
                    %fprintf(' figures %g %g %g \n',ifig,list_save_figures(ifig),fig_label(list_save_figures(ifig)))
                    if (ishandle(fig_label(list_save_figures(ifig))))
                        fign=figure(list_save_figures(ifig));
                        figname=list_name_figures{ifig};
                        %figname = 'fisher_corr_';
                        if (ismember(list_save_figures(ifig),[30,31,32]))
                            % do not save separate figures for different times
                            file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                        else
                            file_name = strcat(file_location_figure,figname,'_',file_tag,'_','i',num2str(istep),file_type_figure);
                        end
                        %saveas(figure(fig_number),file_name,file_type_figure)
                        hgsave(fign,file_name)
                    end
                end
                
            end
            
            % connectivity
            
            if ((sum(time==ceil(probe_times{1}))>0||CS==CS_minimum&&i==steps_param_switch(iparam)) && log_connectivity_nonlinear_GC>=1)
                
                Wmp=Wmg*Wpg';    Wmm=Wmg*Wmg';
                x1min=floor((xpositions_learn{1}(1)+xpositions_learn{1}(2))/2-widths_learn{1}(1)/2);
                x1max=ceil((xpositions_learn{1}(1)+xpositions_learn{1}(2))/2+widths_learn{1}(1)/2);
                x2min=floor((xpositions_learn{1}(3)+xpositions_learn{1}(4))/2-widths_learn{1}(2)/2);
                x2max=ceil((xpositions_learn{1}(3)+xpositions_learn{1}(4))/2+widths_learn{1}(2)/2);
                wmp_from_memory=sum(Wmp(:,x1min:x1max),2);
                wmp_memory_correct=sum(wmp_from_memory(x1min:x1max));
                %                wmp_memory_intermediate=sum(wmp_from_memory(x1max:x2min));
                wmp_memory_incorrect=sum(wmp_from_memory(x2min:x2max));
                wmp_memory_other=sum(wmp_from_memory)-wmp_memory_incorrect-wmp_memory_correct;
                wmp_ratio_incorrect=wmp_memory_incorrect/wmp_memory_correct;
                wmp_ratio_other=wmp_memory_other/wmp_memory_correct;
                wmm_from_memory=sum(Wmm(:,x1min:x1max),2);
                wmm_memory_correct=sum(wmm_from_memory(x1min:x1max));
                %                wmm_memory_intermediate=sum(wmm_from_memory(x1max:x2min));
                wmm_memory_incorrect=sum(wmm_from_memory(x2min:x2max));
                wmm_memory_other=sum(wmm_from_memory)-wmm_memory_incorrect-wmm_memory_correct;
                wmm_ratio_incorrect=wmm_memory_incorrect/wmm_memory_correct;
                wmm_ratio_other=wmm_memory_other/wmm_memory_correct;
                
                figname='fisher_learn_connectivity_0'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
                fid_fisher_learn_MC=fopen(file_name,'a');
                fprintf(fid_fisher_learn_MC,' %g %g %g %g %g %g %g %g %g %g %g \n',time,...
                    wmp_ratio_incorrect,wmp_ratio_other,wmp_memory_correct,wmp_memory_other,wmp_memory_incorrect,...
                    wmm_ratio_incorrect,wmm_ratio_other,wmm_memory_correct,wmm_memory_other,wmm_memory_incorrect);
                fclose(fid_fisher_learn_MC);
            end
            
            % determine the contributions to the P->M inhibition if the MC
            % and/or GC have nonlinear dynamics; this connectivity is
            % specific to the presented odor and does not give much
            % information about the impact for other odors
            if ((sum(time==ceil(probe_times{1}))>0||CS==CS_minimum&&i==steps_param_switch(iparam)) && log_connectivity_nonlinear_GC==2)
                [P_tmp,GC_tmp,~,~,~,P_recall_tmp] = cal_activity1(CS,Wpg,Wmg,Wgm,[S_set{1}(:,2)';S_set{1}(:,3)']',P,time);%Extra argument for S_mem and Wpg
                fig_label(15)=figure(15); set(fig_label(15),'Position',pos{15});
                % Odor 2
                [inh_incr_mc,inh_incr_pc]=inhibitory_contributions(P_tmp(:,1),GC_tmp(:,1),P_recall_tmp(:,1));
                h1=subplot_tight(2,2,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
                Wmm_eff=Wmg*inh_incr_mc;
                imagesc(Wmm_eff-diag(diag(Wmm_eff))); axis square; colorbar
                title('Wmg*Wgm Inhibition from MC to MC with odor 2')
                h1=subplot_tight(2,2,2,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
                %imagesc(Wmg*Wpg'-diag(diag(Wmg*Wpg'))); axis square; colorbar %XH
                Wmp_eff=wt1*Wmg*inh_incr_pc;
                imagesc(Wmp_eff); axis square; colorbar
                title('Wmg*Wgp Inhibition from PC to MC with odor 2')
                
                % quantify connectivity
                x1min=floor((xpositions_learn{1}(1)+xpositions_learn{1}(2))/2-widths_learn{1}(1)/2);
                x1max=ceil((xpositions_learn{1}(1)+xpositions_learn{1}(2))/2+widths_learn{1}(1)/2);
                x2min=floor((xpositions_learn{1}(3)+xpositions_learn{1}(4))/2-widths_learn{1}(2)/2);
                x2max=ceil((xpositions_learn{1}(3)+xpositions_learn{1}(4))/2+widths_learn{1}(2)/2);
                wmp_from_memory=sum(Wmp_eff(:,x1min:x1max),2);
                wmp_memory_correct=sum(wmp_from_memory(x1min:x1max));
                wmp_memory_incorrect=sum(wmp_from_memory(x2min:x2max));
                wmp_memory_other=sum(wmp_from_memory)-wmp_memory_incorrect-wmp_memory_correct;
                wmp_ratio_incorrect=wmp_memory_incorrect/wmp_memory_correct;
                wmp_ratio_other=wmp_memory_other/wmp_memory_correct;
                wmm_from_memory=sum(Wmm_eff(:,x1min:x1max),2);
                wmm_memory_correct=sum(wmm_from_memory(x1min:x1max));
                wmm_memory_incorrect=sum(wmm_from_memory(x2min:x2max));
                wmm_memory_other=sum(wmm_from_memory)-wmm_memory_incorrect-wmm_memory_correct;
                wmm_ratio_incorrect=wmm_memory_incorrect/wmm_memory_correct;
                wmm_ratio_other=wmm_memory_other/wmm_memory_correct;
                
                figname='fisher_learn_connectivity_1'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
                fid_fisher_learn_MC=fopen(file_name,'a');
                fprintf(fid_fisher_learn_MC,' %g %g %g %g %g %g %g %g %g %g %g \n',time,...
                    wmp_ratio_incorrect,wmp_ratio_other,wmp_memory_correct,wmp_memory_other,wmp_memory_incorrect,...
                    wmm_ratio_incorrect,wmm_ratio_other,wmm_memory_correct,wmm_memory_other,wmm_memory_incorrect);
                fclose(fid_fisher_learn_MC);
                
                % Odor 3
                [inh_incr_mc,inh_incr_pc]=inhibitory_contributions(P_tmp(:,2),GC_tmp(:,2),P_recall_tmp(:,2));
                h1=subplot_tight(2,2,3,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
                Wmm_eff=Wmg*inh_incr_mc;
                imagesc(Wmm_eff-diag(diag(Wmm_eff))); axis square; colorbar
                title('Wmg*Wgm Inhibition from MC to MC with odor 3')
                h1=subplot_tight(2,2,4,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
                %imagesc(Wmg*Wpg'-diag(diag(Wmg*Wpg'))); axis square; colorbar %XH
                Wmp_eff=wt1*Wmg*inh_incr_pc;
                imagesc(Wmp_eff); axis square; colorbar
                title('Wmg*Wgp Inhibition from PC to MC with odor 3')
                drawnow
                
                x1min=floor((xpositions_learn{1}(1)+xpositions_learn{1}(2))/2-widths_learn{1}(1)/2);
                x1max=ceil((xpositions_learn{1}(1)+xpositions_learn{1}(2))/2+widths_learn{1}(1)/2);
                x2min=floor((xpositions_learn{1}(3)+xpositions_learn{1}(4))/2-widths_learn{1}(2)/2);
                x2max=ceil((xpositions_learn{1}(3)+xpositions_learn{1}(4))/2+widths_learn{1}(2)/2);
                wmp_from_memory=sum(Wmp_eff(:,x2min:x2max),2);
                wmp_memory_correct=sum(wmp_from_memory(x2min:x2max));
                wmp_memory_incorrect=sum(wmp_from_memory(x1min:x1max));
                wmp_memory_other=sum(wmp_from_memory)-wmp_memory_incorrect-wmp_memory_correct;
                wmp_ratio_incorrect=wmp_memory_incorrect/wmp_memory_correct;
                wmp_ratio_other=wmp_memory_other/wmp_memory_correct;
                wmm_from_memory=sum(Wmm_eff(:,x2min:x2max),2);
                wmm_memory_correct=sum(wmm_from_memory(x2min:x2max));
                wmm_memory_incorrect=sum(wmm_from_memory(x1min:x1max));
                wmm_memory_other=sum(wmm_from_memory)-wmm_memory_incorrect-wmm_memory_correct;
                wmm_ratio_incorrect=wmm_memory_incorrect/wmm_memory_correct;
                wmm_ratio_other=wmm_memory_other/wmm_memory_correct;
                
                figname='fisher_learn_connectivity_2'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
                fid_fisher_learn_MC=fopen(file_name,'a');
                fprintf(fid_fisher_learn_MC,' %g %g %g %g %g %g %g %g %g %g %g \n',time,...
                    wmp_ratio_incorrect,wmp_ratio_other,wmp_memory_correct,wmp_memory_other,wmp_memory_incorrect,...
                    wmm_ratio_incorrect,wmm_ratio_other,wmm_memory_correct,wmm_memory_other,wmm_memory_incorrect);
                fclose(fid_fisher_learn_MC);
                
            end
            
            if (CS==CS_minimum)
                time=time+dt*(inner_steps-j);
                fprintf(' CS too small. exit loop i=%g j=%g jump to time = %g \n',i,j,time);
                CS=CS*modify_CS_frac_2;
                break
            end
            
        end   %end of inner steps loop
    end    %end of if statement for continuous vs discrete
    
    if last ==1
        %close(Smem_M_P_video);
        %           close(fisher_cortical_input_video);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ...This section is designed to check the granule cell activity with different signals
        ... after saturation. ex, light signal with no odor, odor signals with no light, etc...
        ...to remove this bit of code just delete everything within the percent signs and delete the percent sign in front of the next line down.
        %i/exp_time/step < .9 ||
    
    
    %    if i/(exp_time/step) >.5
    if (i>1000000*outer_steps)
        
        %want to probe network with wrong contextual signals
        %flip around light signal
        if token123 == 0
            S_probe(N_mitral+1:end,[1 2]) = S_probe(N_mitral+1:end,[4 3 ]);
            S=S_probe;
        end
        token123=1
        log_first_run=1;
        [P, I, I1, I2, S_mem, ~] = cal_activity1(CS,Wpg,Wmg,Wgm,S_probe,P,time);%Extra argument for S_mem and Wpg
        
        wrongsig2 = fisher_rect(P(1:N_mitral,1),P(1:N_mitral,2),fisher_thresh)
        withoutlight2 = fisher_rect(P(1:N_mitral,3),P(1:N_mitral,4),fisher_thresh)
        
    end
    
    if (i>outer_steps)
        
        stim_select=[1,2];
        [graded_stimuli]=graded(S,probe,stim_select);
        log_first_run=1;
        [P_tmp,~,~,~,S_mem,~] = cal_activity1(CS,Wpg,Wmg,Wgm,graded_stimuli,P,time);%Extra argument for S_mem and Wpg
        
        plot_S_M_P(49,graded_stimuli,P_tmp,S_mem)
        
    end
    
    % Change memory in cortex by learning different stimuli and probe its
    % effect on odor discrimination
    
    if (i>outer_steps && log_cortex_fastlearn==1)
        
        fprintf(' before learning the mixture in cortex\n');
        c_assoc=c_assoc_recall;
        probe_system(time,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)
        
        c_assoc=c_assoc_learn;
        %wt1=0;
        CS_tmp=CS; CS=0;
        W_PM_tmp=W_PM;
        %probe_learn=(probe(:,1)+probe(:,2))/2;
        probe_learn=(probe_set{1}(:,1)+probe_set{1}(:,2));
        figure(100000)
        imagesc(probe_learn)
        colorbar
        pirilearn(probe_learn,0,file_location_figure,file_tag,istep,label_connectivity)
        CS=CS_tmp;
        W_PM=W_PM_tmp;
        
        %wt1=wt1_later;
        
        fprintf(' after learning the mixture in cortex\n');
        c_assoc=c_assoc_recall;
        probe_system(time+0.1,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)
    end
    
    if (i==outer_steps-1 && log_graded==1)
        %if (i==1 && log_graded==1)
        fprintf(' i = %g istep = %g \n',i,istep)
        
        [graded_stimuli]=graded(S(:,graded_sel_S),probe_set{1}(:,graded_sel_probe));
        W_PP_tmp=W_PP; W_PP_i_tmp=W_PP_i; W_PM_tmp=W_PM;
        figure(10000)
        imagesc(graded_stimuli)
        colorbar
        
        dw_pp=W_PP_max/2; dw_inhib=w_inhib/2;
        w_pp_graded=[W_PP_max-dw_pp:dw_pp/2:W_PP_max+dw_pp];
        w_pp_i_graded=[w_inhib-dw_inhib:dw_inhib/2:w_inhib+dw_inhib];
        
        W_PM_scale_graded=[1];
        n_pp=length(w_pp_graded); n_pp_i=length(w_pp_i_graded); n_wpm_scale=length(W_PM_scale_graded);
        w_pp_graded=w_pp_graded/W_PP_max; w_pp_i_graded=w_pp_i_graded/w_inhib;
        W_PM_scale_graded=W_PM_scale_graded/W_PM_scale;
        for i_pm=1:n_wpm_scale
            for i_pp=1:n_pp
                for i_pp_i=1:n_pp_i
                    W_PP=W_PP_tmp*w_pp_graded(i_pp);
                    W_PP_i=W_PP_i_tmp*w_pp_i_graded(i_pp_i);
                    W_PM=W_PM_tmp*W_PM_scale_graded(i_pm);
                    fprintf(' W_PP = %g W_PP_i = %g W_PM = %g \n',max(max(W_PP)),max(max(W_PP_i)),max(max(W_PM)));
                    log_first_run=1;
                    [P_tmp,~,~,~,S_mem,P_recall_tmp] = cal_activity1(CS,Wpg,Wmg,Wgm,graded_stimuli,P,time);
                    
                    figure(i_pm*1000+1)
                    subplot_tight(n_pp,n_pp_i,(i_pp-1)*n_pp_i+i_pp_i,[.05,.05])
                    imagesc(P_recall_tmp)
                    colorbar
                    title(' P_{recall}')
                    
                    figure(i_pm*1000+2)
                    subplot_tight(n_pp,n_pp_i,(i_pp-1)*n_pp_i+i_pp_i,[.05,.05])
                    imagesc(P_tmp)
                    colorbar
                    title(' Mitral Cells')
                    
                    figure(i_pm*1000+3)
                    subplot_tight(n_pp,n_pp_i,(i_pp-1)*n_pp_i+i_pp_i,[.05,.05])
                    imagesc(S_mem)
                    colorbar
                    title(' S_{mem}')
                end
            end
            
            figure(i_pm*1000+1)
            file_graded=['graded_precall_',num2str(W_PM_scale_graded(i_pm)),'_',num2str(w_pp_graded(i_pp)),'_',num2str(w_pp_i_graded(i_pp_i))];
            file_name = strcat(file_location_figure,file_graded,'_',file_tag,file_type_figure);
            hgsave(file_name);
            figure(i_pm*1000+2)
            file_graded=['graded_M_',num2str(W_PM_scale_graded(i_pm)),'_',num2str(w_pp_graded(i_pp)),'_',num2str(w_pp_i_graded(i_pp_i))];
            file_name = strcat(file_location_figure,file_graded,'_',file_tag,file_type_figure);
            hgsave(file_name);
            figure(i_pm*1000+3)
            file_graded=['graded_Smem_',num2str(W_PM_scale_graded(i_pm)),'_',num2str(w_pp_graded(i_pp)),'_',num2str(w_pp_i_graded(i_pp_i))];
            file_name = strcat(file_location_figure,file_graded,'_',file_tag,file_type_figure);
            hgsave(file_name);
            
        end
        W_PP=W_PP_tmp; W_PP_i=W_PP_i_tmp; W_PM=W_PM_tmp;
        
    end
    
    if  (i<=1*outer_steps)
        c_assoc=c_assoc_recall;
        log_first_run=1;
        [P, I, I1, I2, S_mem, ~] = cal_activity1(CS,Wpg,Wmg,Wgm,S,P,time);%Extra argument for S_mem and Wpg
    else
        [r c ] = size(S);
        c_assoc=c_assoc_recall;
        odor_light = S;
        if (1==2)
            odor_nolight = [S(1:r/2,:); -1.6*ones(r/2,c)];
            noodor_light = [.4*ones(r/2,c); S(r/2+1:end,:)];
            noodor_nolight = [.4*ones(r/2,c); -1.6*ones(r/2,c)];
        else
            odor_nolight = [S(1:N_mitral,:); light_back*ones(Nc-N_mitral,c)];
            noodor_light = [S_spontaneous*ones(N_mitral,c); S(N_mitral+1:end,:)];
            noodor_nolight = [S_spontaneous*ones(N_mitral,c); light_back*ones(Nc-N_mitral,c)];
        end
        log_first_run=1;
        [P_ol,I_ol,I_ol_1,I_ol_2,S_mem_ol,P_recall_ol] = cal_activity1(CS,Wpg,Wmg,Wgm,odor_light,P,time);
        [P_onl,I_onl,I_onl_1,I_onl_2,S_mem_onl,P_recall_onl] = cal_activity1(CS,Wpg,Wmg,Wgm,odor_nolight,P,time);
        [P_nol,I_nol,I_nol_1,I_nol_2,S_mem_nol,P_recall_nol] = cal_activity1(CS,Wpg,Wmg,Wgm,noodor_light,P,time);
        [P_nonl,I_nonl,I_nonl_1,I_nonl_2,S_mem_nonl,P_recall_nonl] = cal_activity1(CS,Wpg,Wmg,Wgm,noodor_nolight,P,time);
        
        [P, I, I1, I2, S_mem,P_recall] = cal_activity1(CS,Wpg,Wmg,Wgm,S,P,time);
        
        %corr_feed = [I_ol(:,[2 3 6 7]),I_nol(:,[2 3 6 7]),I_onl(:,[2 3 6 7]),I_nonl(:,[2 3 6 7])];
        corr_feed = [I_ol(Iage>0,choose),I_nol(Iage>0,choose),I_onl(Iage>0,choose),I_nonl(Iage>0,choose)];
        
        corr_out = corrcoef(corr_feed);
        
        figure(7)
        subplot(2,3,1); imagesc(I_ol); axis square; colorbar; title('g cell odor/light');
        subplot(2,3,2); imagesc(I_ol_1); axis square; colorbar; title('2: g cell odor/light');
        subplot(2,3,3); imagesc(I_ol_2); axis square; colorbar; title('1: g cell odor/light');
        subplot(2,3,4); imagesc(P_ol); axis square; colorbar; title('m cell odor/light');
        subplot(2,3,5); imagesc(S_mem_ol); axis square; colorbar; title('p cell odor/light');
        subplot(2,3,6); imagesc(odor_light); axis square; colorbar; title('signal');
        
        figure(8)
        subplot(2,3,1); imagesc(I_onl); axis square; colorbar; title('g cell odor/no light');
        subplot(2,3,2); imagesc(I_onl_1); axis square; colorbar; title('1: g cell odor/no light');
        subplot(2,3,3); imagesc(I_onl_2); axis square; colorbar; title('2: g cell odor/no light');
        subplot(2,3,4); imagesc(P_onl); axis square; colorbar; title('m cell odor/no light');
        subplot(2,3,5); imagesc(S_mem_onl);axis square; colorbar; title('p cell odor/no light');
        subplot(2,3,6); imagesc(odor_nolight); axis square; colorbar; title('signal');
        
        figure(9)
        subplot(2,3,1); imagesc(I_nol); axis square; colorbar; title('g cell no odor/with light');
        subplot(2,3,2); imagesc(I_nol_1); axis square; colorbar; title('1: g cell no odor/with light');
        subplot(2,3,3); imagesc(I_nol_2); axis square; colorbar; title('2: g cell no odor/with light');
        subplot(2,3,4); imagesc(P_nol); axis square; colorbar; title('m cell no odor/with light');
        subplot(2,3,5); imagesc(S_mem_nol); axis square; colorbar; title('p cell no odor/with light');
        subplot(2,3,6); imagesc(noodor_light); axis square; colorbar; title('signal');
        
        figure(10)
        subplot(2,3,1); imagesc(I_nonl); axis square; colorbar; title('g cell no odor/no light');
        subplot(2,3,2); imagesc(I_nonl_1); axis square; colorbar; title('1: g cell no odor/no light');
        subplot(2,3,3); imagesc(I_nonl_2); axis square; colorbar; title('2: g cell no odor/no light');
        subplot(2,3,4); imagesc(P_nonl); axis square; colorbar; title('m cell no odor/no light');
        subplot(2,3,5); imagesc(S_mem_nonl); axis square; colorbar; title('p cell no odor/no light');
        subplot(2,3,6); imagesc(noodor_nolight); axis square; colorbar; title('signal');
        
        figure(11)
        imagesc(corr_out); title('GC-corr of two pairs of signals, 1-4 odor light, 5-8 spontaneous odor with light, 9-12 odor no light, 13-16 spontaneous odor no light signal'); axis square; colorbar
        
        figure(12)
        subplot(4,1,1)
        plot([I_ol(Iage>0,1),I_nol(Iage>0,1),I_onl(Iage>0,1)]); legend('ol1','nol1','onl1');
        subplot(4,1,2)
        plot([I_ol(Iage>0,2),I_nol(Iage>0,2),I_onl(Iage>0,2)]); legend('ol2','nol2','onl2');
        subplot(4,1,3)
        plot([I_ol(Iage>0,3),I_nol(Iage>0,3),I_onl(Iage>0,3)]);  legend('ol3','nol3','onl3');
        subplot(4,1,4)
        plot([I_ol(Iage>0,4),I_nol(Iage>0,4),I_onl(Iage>0,4)]); legend('ol4','nol4','onl4');
        
        if save_video == 1
            if last == 0
                frame =getframe(7);
                writeVideo(Cell_Activity_with_odor_with_light_video, frame);
                
                frame =getframe(8);
                writeVideo(Cell_Activity_with_odor_no_light_video, frame);
                
                frame =getframe(9);
                writeVideo(Cell_Activity_no_odor_with_light_video, frame);
                
                frame =getframe(10);
                writeVideo(Cell_Activity_no_odor_no_light_video, frame);
                
                %         frame =getframe(11);
                %         writeVideo(Correlation_between_g_Cell_activities_video, frame);
                %
                %         frame =getframe(12);
                %         writeVideo(GC_odor2_ol_nol_video, frame);
                %
                %         frame =getframe(13);
                %         writeVideo(GC_odor6_ol_nol_video, frame);
                
            else
                figname = 'Cell_Activity_with_odor_with_light';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                saveas(7,file_name,file_type_figure)
                frame =getframe(7);
                writeVideo(Cell_Activity_with_odor_with_light_video, frame);
                close(Cell_Activity_with_odor_with_light_video);
                
                figname = 'Cell_Activity_with_odor_no_light';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                saveas(8,file_name,file_type_figure)
                frame =getframe(8);
                writeVideo(Cell_Activity_with_odor_no_light_video, frame);
                close(Cell_Activity_with_odor_no_light_video);
                
                figname = 'Cell_Activity_no_odor_with_light';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                saveas(9,file_name,file_type_figure)
                frame =getframe(9);
                writeVideo(Cell_Activity_no_odor_with_light_video, frame);
                close(Cell_Activity_no_odor_with_light_video);
                
                figname = 'Cell_Activity_no_odor_no_light';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                saveas(10,file_name,file_type_figure)
                frame =getframe(10);
                writeVideo(Cell_Activity_no_odor_no_light_video, frame);
                close(Cell_Activity_no_odor_no_light_video);
                
                %             figname = 'Correlation_between_g_Cell_activities';
                %             file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                %             saveas(11,file_name,file_type_figure)
                %         frame =getframe(11);
                %         writeVideo(Correlation_between_g_Cell_activities_video, frame);
                %         close(Correlation_between_g_Cell_activities_video);
                %
                %             figname = 'GC_odor2_ol_nol';
                %             file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                %             saveas(12,file_name,file_type_figure)
                %         frame =getframe(12);
                %         writeVideo(GC_odor2_ol_nol_video, frame);
                %         close(GC_odor2_ol_nol_video);
                %
                %             figname = 'GC_odor6_ol_nol';
                %             file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                %             saveas(13,file_name,file_type_figure)
                %         frame =getframe(13);
                %         writeVideo(GC_odor6_ol_nol_video, frame);
                %         close(GC_odor6_ol_nol_video);
                
            end
        end
    end
    
    if (1==2)
        figure(14)
        title('P');
        imagesc(corrcoef(P(1:N_mitral,:)));
        caxis([-1 1])
        axis square
        colorbar;
        if last == 1
            figname = 'corr_P';
            file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
            saveas(gcf,file_name,file_type_figure)
        end
    end
    
    if (1==2)
        
        temp1 = corrcoef(P(1:N_mitral,:));
        y12(i) = temp1(1,2);
        y22(i) = temp1(3,4);
        % y32(i) = temp1(5,6);
        % y42(i) = temp1(7,8);
        
        x2(i) = i;
        figure(15)
        plot(x2,y12,'blue',x2,y22,'magenta');%,x2,y32,'green',x2,y42,'red')
        legend('pair 1', 'pair 2');%, 'pair 3', 'pair 4')
        title(' Correlation of pairs')
        if save_video == 1
            if last == 0
                frame =getframe(15);
                writeVideo(new_corr_track_video, frame);
            else
                figname = 'new_corr_track';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                saveas(15,file_name,file_type_figure)
                frame =getframe(15);
                writeVideo(new_corr_track_video, frame);
                close(new_corr_track_video);
            end
        end
    end
    
    if( 1==2)
        %    Pcorr_t(i+1) = mean_excluNaN(uptri_1d(corr(P(1:Nc/2,:),TYPE)));
        Pcorr_t(i+1) = mean_excluNaN(uptri_1d(corr(P(1:N_mitral,:),TYPE)));
        %    Pangle_t(i+1) = mean_excluNaN(uptri_1d(corr_angle(corr(S,TYPE),corr(P(1:Nc/2,:),TYPE))));
        Pangle_t(i+1) = mean_excluNaN(uptri_1d(corr_angle(corr(S,TYPE),corr(P(1:N_mitral,:),TYPE))));
        %    Pfoc_t(i+1) = mean_excluNaN(focality(P(1:Nc/2,:),metric));
        Pfoc_t(i+1) = mean_excluNaN(focality(P(1:N_mitral,:),metric));
        %    CV_t(i+1) = std(mean(P(1:Nc/2,:),2))/mean(mean(P(1:Nc/2,:),2));
        CV_t(i+1) = std(mean(P(1:N_mitral,:),2))/mean(mean(P(1:N_mitral,:),2));
        %   CVid_t(:,i+1) = std(P(1:Nc/2,:))./mean(P(1:Nc/2,:));
        if (1==2)
            CVid_t(:,i+1) = std(P(1:N_mitral,:))./mean(P(1:N_mitral,:));
        end
        %    F_t(:,:,i+1) = corr(P(1:Nc/2,:));
        F_t(:,:,i+1) = corr(P(1:N_mitral,:));
        if tracking == 1
            %Tcorr_t(i+1) = mean_excluNaN(cal_track_corr(track_pairs,P));
            Tcorr_t(1:size(track_pairs,1),i+1) = cal_track_corr(track_pairs,P)';
            Tangle_t(i+1) = mean_excluNaN(corr_angle(cal_track_corr(track_pairs,S),cal_track_corr(track_pairs,P)));
        end
    end
    
    if TEXT_OUT == 1
        %    fprintf('           corr = %f\n', Pcorr_t(i+1));
    end
    
    if ANIME_OUT == 1;
        %        update_Pplot(P,corr(P(1:Nc/2,:),TYPE),corr(P(1:Nc/2,:)',TYPE),rg,HP1d,Wmg,Wgm,N_t);
        if (1==2)
            update_Pplot(N_mitral,P,corr(P(1:N_mitral,:),TYPE),corr(P(1:N_mitral,:)',TYPE),rg,HP1d,Wmg,Wgm);
        end
        update_Iplot(cont_density,time_axis,I,I1,I2,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(i,:),HI);
        if (1==2)
            VST = cell(1,2);
            VST{1} = Pcorr_t;
            for i_pairs=1:size(track_pairs,1)
                VST{1+i_pairs} = Tcorr_t(i_pairs,:);
            end
            %        update_Infoplot(VST,corr(P(1:Nc/2,:),TYPE),Hinfo);
            if (1==2)
                update_Infoplot(VST,corr(P(1:N_mitral,:),TYPE),Hinfo);
            end
            bargraph(Wmg,'Wmg',Wpg,'Wpg')
        end
        
        drawnow;
        
        if last == 1
            
            if (1==2)
                figname = 'PPLOT';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                %            saveas(figure(2),file_name,file_type_figure)
                figure(2)
                hgsave(file_name);
                
                figname = 'G_Cell_Activity_Pattern_Corr';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                %saveas(figure(3),file_name,file_type_figure)
                figure(3)
                hgsave(file_name);
                
                figname = 'Correlation_Pairs';
                file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                %saveas(figure(4),file_name,file_type_figure)
                figure(4)
                hgsave(file_name);
            end
            
        end
        
    end
    
    if ANIME2D_OUT == 1;
        setup_Pplot2D(P,coord,102,S_name);
        drawnow;
        figname = 'Pplot 2D';
        file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
        %saveas(figure(101),file_name,file_type_figure)
        figure(101)
        hgsave(file_name);
    end
    
    %At the end because we don't want to use a different S_mem immediately
    %Let's try the beginning again
    % G cell i is connected to P cell j
    % figure(24)
    % subplot(2,1,1)
    % imagesc(Wmg)
    % title('Wmg')
    % subplot(2,1,2)
    % imagesc(Wpg)
    % title('Wpg')
    
    if (log_hist_Ca==1)
        fig_label(83)=figure(83);
        Ca_max=max(max(Ca_),1);
        nbin_Ca=200;
        hist(Ca_,nbin_Ca)
        xlim([Ca_max/nbin_Ca 1.1*Ca_max])
    end
    
    if (log_selected_connectivity==2)
        fig_label(82)=figure(82);
        res_frac=[.0,1,2,3,4,Inf];
        for i_res=1:length(res_frac)-1
            resilience_low=res_frac(i_res)*ts; resilience_high=res_frac(i_res+1)*ts;
            Wmg_select=Wmg;
            Wmg_select(:,Ca_<resilience_low | Ca_>resilience_high)=0;
            h1=subplot(length(res_frac)-1,2,(2*i_res)-1);  ax1 = get(h1,'position'); set(h1,'position', ax1);
            imagesc(Wmg_select*Wgm-diag(diag(Wmg_select*Wgm))); axis square; colorbar
            title('selected GC: Wmg*Wgm Inhibition from MC to MC')
            h1=subplot(length(res_frac)-1,2,2*i_res); ax1 = get(h1,'position'); set(h1,'position', ax1);
            %imagesc(Wmg*Wpg'-diag(diag(Wmg*Wpg'))); axis square; colorbar %XH
            imagesc(Wmg_select*Wpg'); axis square; colorbar %XH
            title('selected GC: Wmg*Wgp Inhibition from PC to MC')
        end
    end
    if (log_selected_connectivity>=1)
        [Ca_sorted,Ca_sort_indices]=sort(Ca_,'descend');
        Wmg_sorted=Wmg(1:N_mitral,Ca_sort_indices(Ca_sorted>0));
        Wpg_sorted=Wpg(1:N_mitral,Ca_sort_indices(Ca_sorted>0));
        fig_label(84)=figure(84);
        subplot(5,1,1)
        imagesc(1-Wmg_sorted);
        colormap('bone');
        title('Wmg sorted')
        subplot(5,1,2)
        imagesc(1-Wpg_sorted);
        title('Wpg_sorted');
        colormap('bone');
        subplot(5,1,3)
        plot(Ca_sorted(Ca_sorted>0));
        xlim([0 max(1,length(Ca_sorted(Ca_sorted>0)))]);
        subplot(5,1,4)
        imagesc(1-Wmg(1:N_mitral,:));
        title('Wmg');
        colormap('bone');
        subplot(5,1,5)
        imagesc(1-Wpg(1:N_mitral,:));
        title('Wpg');
        colormap('bone');
    end
    
    fig_label(16)=figure(16); set(fig_label(16),'Position',pos{16});
    % note that Fisher is calculated before GC are removed, but the
    % connectivity is shown after they are removed.
    %        h1=subplot(3,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
    h1=subplot_tight(1,2,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    Wmm_plot=Wmg*Wgm;
    imagesc(Wmm_plot(1:N_mitral,1:N_mitral)-diag(diag(Wmm_plot(1:N_mitral,1:N_mitral)))); axis square; colorbar
    title('Wmg*Wgm Inhibition from MC to MC')
    %        h1=subplot(3,2,2); ax1 = get(h1,'position'); set(h1,'position', ax1);
    h1=subplot_tight(1,2,2,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    %imagesc(Wmg*Wpg'-diag(diag(Wmg*Wpg'))); axis square; colorbar %XH
    Wmp_plot=Wmg*Wpg';
    imagesc(Wmp_plot(1:N_mitral,1:N_mitral)); axis square; colorbar %XH
    title(['Wmg*Wgc Inhibition from CC to MC time = ',num2str(time)]);
    %wpm_max2=max(max(Wmg*Wpg'));
    drawnow
    %fprintf(' time = %g wpm_max2 = %g \n',time,wpm_max2);
    %         if (time>200 && wpm_max2<20)
    %                 plot_S_M_P(17,S_mem,P,P_recall)
    %                 probe_system(time,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)
    %             pause
    %         end
    %title('Wmg*Wgp Inhibition from PC to MC')
    
    if (1==2)
        subplot_tight(2,2,1,[.05,.05])
        imagesc(Wmg*Wgm-diag(diag(Wmg*Wgm))); axis square; colorbar
        title('Wmg*Wgm Inhibition from MC to MC')
        subplot_tight(2,2,2,[.05,.05])
        imagesc(MCdiag); axis square; colorbar
        title('MCdiag Diagnostic - OK')
        subplot_tight(2,2,3,[.05,.05])
        imagesc(Wmg*Wpg'-diag(diag(Wmg*Wpg'))); axis square; colorbar
        title('Wmg*Wgp Inhibition from PC to MC')
        subplot_tight(2,2,4,[.05,.05])
        imagesc(Wpg*Wpg'-diag(diag(Wpg*Wpg'))); axis square; colorbar
        title('Wpg*Wgp Which combinations of PC connect to the same GC?')
    end
    figname = 'Wmg_Wmg_Wmp_Wmp';
    file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
    %saveas(figure(16),file_name,file_type_figure)
    hgsave(fig_label(16),file_name)
    if save_video == 1
        if last == 0
            frame =getframe(16);
            writeVideo(WmgxWgm_WmgxWgp_video, frame);
        else
            figname = 'WmgxWgm_WmgxWgp';
            file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
            %saveas(figure(16),file_name,file_type_figure)
            frame =getframe(16);
            writeVideo(WmgxWgm_WmgxWgp_video, frame);
            
            close(WmgxWgm_WmgxWgp_video);
        end
    end
    
    % Plot response spectrum of the removed GC
    nG_removed=size(G_removed{istep},1);
    if (nG_removed>0)
        max_G_removed=max(max_G_removed,nG_removed);
        for iodor=1:size(S,2)
            response_prev=response{iodor};
            response{iodor}=NaN(istep,max_G_removed);
            response{iodor}(1:size(response_prev,1),1:size(response_prev,2))=response_prev;
            if (nG_removed>0)
                response{iodor}(istep,1:nG_removed)=G_removed{istep}(1:nG_removed,iodor);
            end
        end
        if (time>=time_GC_remove_plot)
            fig_label(30)=figure(30); set(fig_label(30),'Position',pos{30});
            for iodor=1:size(S,2)
                %        subplot(size(S,2),1,2*iodor-1)
                h1=subplot_tight(size(S,2),1,iodor,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
                imagesc(response{iodor})
                colorbar
            end
            fig_label(31)=figure(31); set(fig_label(31),'Position',pos{31});
            for iodor=1:size(S,2)
                h1=subplot_tight(size(S,2),1,iodor,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
                plot(time_istep(1:istep),nanmean(response{iodor},2))
                xlabel(' Time')
                ylabel(' Mean Response')
            end
        end
    end
    
    
    % save relevant figures
    %    if (mod(istep,mod_save)==0||last==1)
    if (mod(time,mod_save)==mod_save_pick||last==1)
        fprintf(' Saving Figures \n')
        for ifig=1:length(list_save_figures)
            %fprintf(' figures %g %g %g \n',ifig,list_save_figures(ifig),fig_label(list_save_figures(ifig)))
            if (ishandle(fig_label(list_save_figures(ifig))))
                fign=figure(list_save_figures(ifig));
                figname=list_name_figures{ifig};
                %figname = 'fisher_corr_';
                if (ismember(list_save_figures(ifig),[30,31,32]))
                    % do not save separate figures for different times
                    file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
                else
                    file_name = strcat(file_location_figure,figname,'_',file_tag,'_','i',num2str(istep),file_type_figure);
                end
                %saveas(figure(fig_number),file_name,file_type_figure)
                hgsave(fign,file_name)
            end
        end
    end
    
end %step

figname='fisher_all_save';
file_name = strcat(file_location_data,figname,'_',file_tag);

if (probe_fisher==1 || probe_fisher==2 || probe_fisher==4)
    save(file_name,'ts_list','th_list','y_fisher','y_fisher_1');
elseif (probe_fisher==3)
    save(file_name,'ts_list','th_list','y_fisher','y_fisher_1','y_fisher_2');
end

%Correlation again

%R2 = corrcoef([P(:,1),P(:,2)]);
R2 = corr([P(:,1),P(:,2)],TYPE);

fprintf('\n   - Stimuli 1 and 2 corr\n')
fprintf('      - Start %f\n      - End %f\n', R1, R2)

% Plot diagnostics

if (1==2)
    % % G cell i is connected to P cell j
    fig24=figure(24);  set(fig24,'Position',pos{24});
    h1=subplot(2,1,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
    imagesc(Wmg)
    title('Wmg')
    h1=subplot(2,1,2); ax1 = get(h1,'position'); set(h1,'position', ax1);
    imagesc(Wpg)
    title('Wpg')
    figname = 'Wpg_Wmg';
    file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
    %saveas(figure(24),file_name,file_type_figure)
    hgsave(fig24,file_name)
end

% How many G cells connect to P cell i?
bargraph(Wmg,'Wmg',Wpg,'Wpg')
figname = 'Sum_of_Connections_Wpg_Wmg';
file_name = strcat(file_location_figure,figname,'_',file_tag,file_type_figure);
%saveas(gcf,file_name,file_type_figure)
hgsave(file_name);
% How many G cells connect to MC i and j, and to PC k and l?
% conntrack(Wmg, Wpg)

% PC = P_connections(6, Wpg, S_mem, time_axis,N_t)

% end

if PLOT_OUT == 1;
    %    update_Pplot(P,corr(P(1:Nc/2,:),TYPE),corr(P(1:Nc/2,:)',TYPE),rg,HP1d,Wmg,Wgm,N_t);
    if (1==2)
        update_Pplot(N_mitral,P,corr(P(1:N_mitral,:),TYPE),corr(P(1:N_mitral,:)',TYPE),rg,HP1d,Wmg,Wgm);
    end
    update_Iplot(cont_density,time_axis,I,I1,I2,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(i,:),HI);
    VST = cell(1,2);
    VST{1} = Pcorr_t;
    for i_pairs=1:size(track_pairs,1)
        VST{1+i_pairs} = Tcorr_t(i_pairs,:);
    end
    %    update_Infoplot(VST,corr(P(1:Nc/2,:),TYPE),Hinfo);
    if (1==2)
        update_Infoplot(VST,corr(P(1:N_mitral,:),TYPE),Hinfo);
    end
    drawnow;
end

if PLOT2D_OUT == 1;
    setup_Pplot2D(P,coord,102,S_name);
    drawnow;
end

objfun = return_val;
if TEXT_OUT == 1
    fprintf('---- End of run %s  -----\n',file_tag);
end

% nested function definition

    function prob_ = survival(G_)
        %        Ca_ = sum(rec(SigmoidG(G_),th,0),2); % XH 4/24 CHECK
        %        Ca_ = sum(rec(G_,th,0),2); % XH 4/24 CHECK
        % repetition_factor can reduce the resilience if an odor is presented
        % multiple times just to keep the total number of learning odors the same
        Ca_ = sum(rec(max(G_-Gthresh,0),th,0),2)/repetition_factor{i_learning_set}; % XH 4/24 CHECK
        prob_ = (tanh((Ca_-ts)*gamma)+1)*(maxv-minv)/2+minv;
    end % survival

    function val = remove_cell
        prob_ = survival(I);
        surv_ = floor(prob_ + rand(size(prob_)));
        IX_ = find((surv_==0)&(Iage'>=0)&(Iage'<age_old));
        if (probe_case==8)
            if (istep~=istep_last)
                W_PP_tmp=W_PP;
                W_PP=W_PP_all;
                log_first_run=1;
                [~,I_tmp,~,~,~,~] = cal_activity1(CS,Wpg,Wmg,Wgm,S_GC_probe,P,time);%Extra argument for S_mem and Wpg
                W_PP=W_PP_tmp;
                G_removed{istep}=I_tmp(IX_,:);
                time_istep(istep)=time;
                istep_last=istep;
            end
        end
        Iage(IX_) = -1; % the age of non-existing cells is labeled negative
        Imark(IX_) = 0;
        Wmg(:,IX_) = 0;
        Wpg(:,IX_) = 0;
        Wgm(IX_,:) = 0;
        val = length(IX_);
    end % remove_cell

    function val = add_cell
        %        N_mitral_add=dt*Na*Nc;
        N_mitral_add=dt*Na;
        IX_ = find(Iage>=0);
        Iage(IX_) = Iage(IX_)+dt;
        IX_ = find(Iage<0);
        if length(IX_) < round(N_mitral_add)
            add_space(round(N_mitral_add)-length(IX_));
            IX_ = find(Iage<0);
        end
        if prob_conn == 0 %Bring in the new matrix Wgp here, similarly to Wgm
            %Make another temp for connections from PC to GC
            %But temp_ is not necessarily for Wgm, Wmg, so I can keep it
            %like that.
            %Should we have different numbers of connections amongst MC and
            %GC and amongst GC and PC?
            %temp_ = [ones(conn,1); zeros(Nc-conn,1)];
            %temp_ = [ones(conn,1); zeros(Nc/2-conn,1)];
            % Note: Nc is the number of mitral cells + number of `light'
            % cells and the number of cells in PC is the same
            tempM_ = [ones(connM,1); zeros(N_mitral-connM,1)];
            %tempP_ = [ones(connP,1); zeros(n_P-connP,1)]; %XH
            tempP_ = [ones(connP,1); zeros(n_P-num_light_cells-connP,1)]; %HR after XH
            for i_ = 1:round(N_mitral_add)
                Iage(IX_(i_)) = 0;
                %                temp_ = temp_(randperm(Nc));
                
                %                Wpg(:, IX_(i_)) = [temp_(randperm(Nc/2));zeros(Nc/2,1)];%Random connectivity in Wpg
                Wpg(:, IX_(i_)) = [tempP_(randperm(n_P-num_light_cells));zeros(num_light_cells,1)];%Random connectivity in Wpg %XH NEED TO WORK ON THIS LATER
                
                %                temp_ = temp_(randperm(Nc));
                %                Wmg(:, IX_(i_)) = [temp_(randperm(Nc/2));zeros(Nc/2,1)];
                Wmg(:, IX_(i_)) = [tempM_(randperm(N_mitral));zeros(Nc-N_mitral,1)];
                %                IX2_ = randperm(length(temp_));
                %                temp1_ = temp_(IX2_(1:round(perm_ratio*length(temp_))));
                %                temp2_ = temp_(IX2_(1+round(perm_ratio*length(temp_)):end));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %confused- why is Wgm not equal to Wmg'?
                %Wgm(IX_(i_),IX2_) = [temp1_(randperm(length(temp1_))); temp2_]';
                Wgm(IX_(i_),:) = Wmg(:, IX_(i_))';
            end
        elseif prob_conn == 2
            tempP_ = [ones(connP,1); zeros(n_P-num_light_cells-connP,1)]; %HR after XH
            for i_ = 1:round(N_mitral_add)
                Iage(IX_(i_)) = 0;
                %                temp_ = temp_(randperm(Nc));
                
                %                Wpg(:, IX_(i_)) = [temp_(randperm(Nc/2));zeros(Nc/2,1)];%Random connectivity in Wpg
                Wpg(:, IX_(i_)) = [tempP_(randperm(n_P-num_light_cells));zeros(num_light_cells,1)];%Random connectivity in Wpg %XH NEED TO WORK ON THIS LATER
                
                %                temp_ = temp_(randperm(Nc));
                %                Wmg(:, IX_(i_)) = [temp_(randperm(Nc/2));zeros(Nc/2,1)];
                %P_mean=mean(SigmoidM(P(1:N_mitral,:)),2);
                %P_mean=mean(SigmoidM(P(1:N_mitral,:)),2);
                P_mean=mean(rect(P(1:N_mitral,:)),2);
                P_prob=rand(N_mitral,1)-attach_prob*P_mean;
                P_prob_min=min(P_prob); P_prob_max=max(P_prob);
                P_thresh=P_prob_min+(P_prob_max-P_prob_min)*connM/N_mitral;
                index_connect=P_prob<=P_thresh;
                while (sum(index_connect)>connM)
                    index_connect(ceil(rand*N_mitral))=0;
                end
                Wmg(:, IX_(i_)) = [index_connect;zeros(Nc-N_mitral,1)];
                %Wmg(:, IX_(i_)) = [tempM_(randperm(N_mitral));zeros(Nc-N_mitral,1)];
                %                IX2_ = randperm(length(temp_));
                %                temp1_ = temp_(IX2_(1:round(perm_ratio*length(temp_))));
                %                temp2_ = temp_(IX2_(1+round(perm_ratio*length(temp_)):end));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %confused- why is Wgm not equal to Wmg'?
                %Wgm(IX_(i_),IX2_) = [temp1_(randperm(length(temp1_))); temp2_]';
                Wgm(IX_(i_),:) = Wmg(:, IX_(i_))';
            end
        elseif prob_conn == 3
            tempP_ = [ones(connP,1); zeros(n_P-num_light_cells-connP,1)]; %HR after XH
            for i_ = 1:round(N_mitral_add)
                Iage(IX_(i_)) = 0;
                %                temp_ = temp_(randperm(Nc));
                
                %                Wpg(:, IX_(i_)) = [temp_(randperm(Nc/2));zeros(Nc/2,1)];%Random connectivity in Wpg
                Wpg(:, IX_(i_)) = [tempP_(randperm(n_P-num_light_cells));zeros(num_light_cells,1)];%Random connectivity in Wpg %XH NEED TO WORK ON THIS LATER
                
                %                temp_ = temp_(randperm(Nc));
                %                Wmg(:, IX_(i_)) = [temp_(randperm(Nc/2));zeros(Nc/2,1)];
                %M_max=max(SigmoidM(P(1:N_mitral,:)),[],2)';
                M_max=max(rect(P(1:N_mitral,:)),[],2)';
                M_max_max=max(M_max);
                P_prob=rand(1,N_mitral)+attach_prefer*M_max/M_max_max;
                [~,sort_indices]=sort(P_prob,'descend');
                index_connect=zeros(1,N_mitral);
                index_connect(sort_indices(1:connM))=1;
                Wmg(:, IX_(i_)) = [index_connect';zeros(Nc-N_mitral,1)];
                Wgm(IX_(i_),:) = Wmg(:, IX_(i_))';
                %                 figure(10000)
                %                 subplot(2,1,1)
                %                 imagesc(SigmoidM(P(1:N_mitral,:)))
                %                 subplot(2,1,2)
                %                 imagesc(index_connect)
            end
        else
            conn_prob_ = conn/Nc;
            for i_ = 1:round(N_mitral_add)
                Iage(IX_(i_)) = 0;
                %                if (marking == 1) && (i*step >= marking_t(1)) && (i*step < marking_t(2))
                if (marking == 1) && (time >= marking_t(1)) && (time < marking_t(2))
                    Imark(IX_(i_)) = 1;
                end
                temp_ = floor(rand(Nc,1)+conn_prob_);
                Wpg(:, IX_(i_)) = temp_;%Making Wpg just like Wmg, but randomised
                temp_ = floor(rand(Nc,1)+conn_prob_);
                Wmg(:, IX_(i_)) = temp_;
                IX2_ = randperm(length(temp_));
                temp1_ = temp_(IX2_(1:round(perm_ratio*length(temp_))));
                temp2_ = temp_(IX2_(1+round(perm_ratio*length(temp_)):end));
                Wgm(IX_(i_),IX2_) = [temp1_(randperm(length(temp1_))); temp2_]';
            end
        end
        val = round(N_mitral_add);
    end % add_cell

    function add_space(short)%Copying everything for functions of Wpg
        previous_size = Isize;
        tempI = I;
        tempIage = Iage;
        tempImark = Imark;
        tempWmg = Wmg;
        tempWgm = Wgm;
        tempWpg = Wpg;
        
        %        Isize = Isize + 5*short;
        Isize = Isize + 1*short;
        I = zeros(Isize, Ns);
        Iage = -ones(1, Isize);
        Imark = zeros(1, Isize);
        Wmg = zeros(Nc, Isize);
        Wgm = zeros(Isize, Nc);
        % HR Wpg should be dimensioned with n_P not Nc
        Wpg = zeros(n_P, Isize);%XH !! Wpg should not be the same size as Wmg!
        %We should rework this in future.
        I(1:previous_size, :) = tempI;
        Iage(1, 1:previous_size) = tempIage;
        Imark(1, 1:previous_size) = tempImark;
        Wmg(:, 1:previous_size) = tempWmg;
        Wgm(1:previous_size,:) = tempWgm;
        Wpg(:, 1:previous_size) = tempWpg;
        
        if ANIME_OUT==1
            set(Iaxis, 'YLim', [1 Isize]);
        end
    end % add_space

    function double_number_GC
        
        Isize=2*Isize;
        I=[I;I];
        Iage=[Iage,Iage];
        Imark=[Imark,Imark];
        Wmg=[Wmg,Wmg];
        Wgm=[Wgm;Wgm];
        Wpg=[Wpg,Wpg];
        
        if ANIME_OUT==1
            set(Iaxis, 'YLim', [1 Isize]);
        end
    end % double_number_GC

    function val = cal_track_corr(track_pairs_,S_)
        corr_ = zeros(1,size(track_pairs_,1));
        for i_ = 1:size(track_pairs_,1)
            %            temp_ = corr(S_(1:Nc/2,track_pairs_(i_,:)),TYPE);
            temp_ = corr(S_(1:N_mitral,track_pairs_(i_,:)),TYPE);
            corr_(i_) = temp_(1,2);
        end
        val = corr_;
    end % cal_track_corr

    function val = return_val
        
        val = 0;
        
    end % return_val

end_time=cputime;

figname='parameter'; file_name = strcat(file_location_data,figname,'_',file_tag,file_type_data);
fid_parameter=fopen(file_name,'a');
fprintf(fid_parameter,' computation time = %g hours \n',(end_time-start_time)/3600);
fclose(fid_parameter);

end % neurogenesis

%Plots the figure 4 business

function H = setup_Infoplot(time_axis,VST,line_style,line_name,Scorr,Pcorr,N)
global pos
figN=figure(N); set(figN,'Position',pos{N}); clf;
%N==4 for this plot
H = cell(1,length(VST)+1);
subplot(1,2,2);
H{1} = plot(uptri_1d(Scorr), uptri_1d(Pcorr), 'x');
hold on; line([0 0],[-1 1],'color',[0 0 0]);
hold on; line([-1 1],[0 0],'color',[0 0 0]);
hold on; line([-1 1],[-1 1],'color',[0 0 0]);
hold off;
xlabel('input'); ylabel('output'); title('in/out corr');  axis square;
subplot(1,2,1);
xlim([min(time_axis) max(time_axis)]); ylim([-1 1]);
xlabel('time'); ylabel('correlation/focality'); title('correlation');  axis square;
subplot(1,2,1);
for i = 1:length(VST)
    hold on;
    H{i+1} = plot(time_axis,VST{i},line_style{i});
end
hold off;
end

%%%%%%%%%%%%%%%

function H = setup_Pplot(N_mitral,P,Pcorr,Ccorr,rg,N,Wmg,Wgm)
global pos
%N==1 for this plot
H = cell(1,4);
figN=figure(N); set(figN,'Position',pos{N}); clf; set(gcf,'DoubleBuffer','on');
subplot(2,2,1);

if (N==1)
    H{1} = imagesc(P,[min(min(P(1:N_mitral,:))) max(max(P(1:N_mitral,:)))]);
else
    H{1} = imagesc(P,[-0.2 0.7]);
end
colormap(jet);
xlabel('Stimulus'); ylabel('MC'); title('MC Activity');  axis square;
subplot(2,2,2);
if (1==2)
    H{2} = imagesc(Ccorr); colormap(jet);
    title('Channel Corr MC'); caxis([-1 1]); axis square;
else
    H{2} = plot(P(1:N_mitral,:));
    title('MC Activity');
    if (N~=1)
        ylim([0 0.5]);
    end
    xlabel('Stimulus'); ylabel('MC');
end
subplot(2,2,3);
H{3} = imagesc(Pcorr); colormap(jet);
title('Pattern Corr MC'); caxis([-1 1]); axis square;
if nargin > 6
    conn = Wmg*Wgm;
    subplot(2,2,4);
    H{4} = imagesc(conn-diag(NaN*diag(conn)));
    title('MC-MC Connections'); axis square;
end
end % setup_Pplot

%%%%%%%%%%

function update_Infoplot(VST,Pcorr,H)
set(H{1}, 'YData', uptri_1d(Pcorr));
for i = 1:length(VST)
    set(H{i+1}, 'YData', VST{i});
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function y=Sigmoid(x)
% sigmoidal threshold function for the input/output relation of the neurons

%global thresh steep P_amp
global thresh_odor steep_odor P_amp_odor
global thresh_light steep_light P_amp_light
global n_P N_mitral

if (size(x,1)==1)
    y(1:N_mitral)=P_amp_odor./(1+exp(-steep_odor*(x(1:N_mitral)-thresh_odor)));
    y(N_mitral+1:n_P)=P_amp_light./(1+exp(-steep_light*(x(N_mitral+1:n_P)-thresh_light)));
elseif (size(x,2)==1)
    y(1:N_mitral,1)=P_amp_odor./(1+exp(-steep_odor*(x(1:N_mitral,1)-thresh_odor)));
    y(N_mitral+1:n_P,1)=P_amp_light./(1+exp(-steep_light*(x(N_mitral+1:n_P,1)-thresh_light)));
else
    y(1:N_mitral,:)=P_amp_odor./(1+exp(-steep_odor*(x(1:N_mitral,:)-thresh_odor)));
    y(N_mitral+1:n_P,:)=P_amp_light./(1+exp(-steep_light*(x(N_mitral+1:n_P,:)-thresh_light)));
end

end %sigmoid


function y=SigmoidG(x)
% sigmoidal threshold function for the input/output relation of the neurons

global Gthresh Gsteep G_amp
global non_lin_G

%if (non_lin_G==1)
%   y=G_amp./(1+exp(-Gsteep*(x-Gthresh)));
%else
y=max(x-Gthresh,0);
%end
end %sigmoid


function y=SigmoidM(x)
% sigmoidal threshold function for the input/output relation of the neurons

global Mthresh Msteep M_amp
global non_lin_M

%if (non_lin_M==1)
%    y=M_amp./(1+exp(-Msteep*(x-Mthresh)));
%else
y=rect(x);
%end
end %sigmoid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dY=P_rate(time,Y,Nc,n_P,P_amp_odor,steep_odor,thresh_odor,P_amp_light,...
    steep_light,thresh_light,Wgm,wt1,Wgp,W_PM,c_assoc,W_PP,W_PP_i,tau_P,stim_select,CS,Wmg,Gthresh)
% differential equation describing the activity P of the neurons in the
% network

%[M_i arg] = cal_activity(non_lin,CS,Wpg,Wmg,Wgm,S,P,rm,rg,S_mem)

%M_i = A\(stim_select-wt1*prod_gp*P);
M=Y(1:Nc); P=Y(Nc+1:Nc+n_P);
% coding the sigmoid and rectifier explicitly speeds up the code noticeably
% (e.g. total time in p_rate goes from 153s down to 123s)
%sigP = Sigmoid(P);

N_mitral=Nc;
sigP(1:N_mitral,1)=P_amp_odor./(1+exp(-steep_odor*(P(1:N_mitral)-thresh_odor)));
sigP(N_mitral+1:n_P,1)=P_amp_light./(1+exp(-steep_light*(P(N_mitral+1:n_P)-thresh_light)));
%sigP=P_amp./(1+exp(-steep*(P-thresh)));

%sigM = SigmoidM(M);
%sigM = rect(M);
sigM = max(M,0*M);

%fprintf(' M %g %g Wgm %g %g Wgp %g %g \n',size(M),size(Wgm),size(Wgp))

% p_scale should be also used here, or rather not at all, instead WPG
% should be modified
%XH 4/23
%GC1 = Wgm*SigmoidM(M); % might be different sigmoid?
GC1 = Wgm*sigM; % might be different sigmoid?
GC2 = wt1*Wgp*sigP(1:n_P);
GC = GC1+GC2;
sigG=max(GC-Gthresh,0);

%dM = -M+stim_select-CS*(Wmg*SigmoidG(GC));
%dM = -M+stim_select-CS*(Wmg*GC);
dM = -M+stim_select-CS*(Wmg*sigG);
%dM = -M+stim_select-CS*(prod_gm*rect(M)+wt1*prod_gp*Sigmoid(P));

%fprintf(' should be Sigmoid(P) rather than P\n');
% until 5/21/14 rect was missing in M
%dP=(-P+W_PM*M+((1-c_assoc)*W_PP-W_PP_i)*Sigmoid(P))/tau_P;

%1201 W_PP_i = 0
% allow negative weights
% relearning multiple times (more similar odors)
%dP=(-P+W_PM*SigmoidM(M)+((1-c_assoc)*W_PP-W_PP_i)*sigP)/tau_P;
dP=(-P+W_PM*sigM+((1-c_assoc)*W_PP-W_PP_i)*sigP)/tau_P;
%dP=(-P+W_PM*rect(M)+((1-c_assoc)*W_PP-W_PP_i)*Sigmoid(P))/tau_P;

%add the scale parameter in from of W_PM

dY=[dM;dP];

end %differential equation


function out=rect(in)

out=max(in,0*in);
%out=in;

end



function val = focality(P, metric)
[Nc Ns] = size(P);
tempval = [];
for sti = 1:Ns
    MCact = rec(P(:,sti),0);
    MCacttemp = MCact;
    MCacttemp(find(MCact<=max(MCact)/2)) = 0;
    nor_down = 0;
    down = 0;
    for ii = 1:length(MCact)
        for jj = ii:length(MCact)
            down = down + metric(ii,jj);
            nor_down = nor_down+1;
        end
    end
    down = down/nor_down;
    nor_up = 0;
    up = 0;
    for ii = 1:length(MCact)
        for jj = ii:length(MCact)
            up = up + MCacttemp(ii)*MCacttemp(jj)*metric(ii,jj);
            nor_up = nor_up + MCacttemp(ii)*MCacttemp(jj);
        end
    end
    up = up/nor_up;
    tempval = [tempval 1-up/down];
end
val = tempval;
end % focality

function val = corr(M,type,S)
if nargin < 2
    type = 1;
end
if type == 1
    val = corrcoef(M);
elseif type == 2
    out = eye(size(M,2));
    M2 = zeros(size(M));
    for i = 1:size(M,2)
        if norm(M(:,i)) == 0
            out = NaN*out;
            return;
        end
        M2(:,i) = M(:,i)/norm(M(:,i));
    end
    val = M2'*M2;
end
end

function val = corr_angle(Scorr,Mcorr)
tempS = 1-Scorr;
tempM = 1-Mcorr;
tempS(find(tempS<=0)) = NaN;
val = atan(tempM./tempS)/pi*4-1;
end

function val = corr_diff(M)
tempM = uptri_1d(1-corrcoef(M));
temp = atan(tempM(IX)./tempS(IX))/pi*4-1;
val = mean(temp);
end

function [n,bin] = histw(varargin)

[cax,args,nargs] = axescheck(varargin{:});

y = args{1};
x = args{2};
if nargs > 2
    w = args{3};
else
    w = ones(size(y));
end

bin = sort(x);
n = zeros(size(bin));
for i = 1:length(y)
    temp = find((bin-y(i))>=0);
    if length(temp) > 0
        n(temp(1)) = n(temp(1))+w(i);
    else
        N(end) = n(end)+w(i);
    end
end

end

function val = mean_excluNaN(V)
temp = V(find(isnan(V)==0));
if length(temp)>0
    val = mean(temp);
else val = NaN;
end
end

function out = rec(V,r,smooth)
if nargin <3
    smooth = 0.1;
end

if smooth == 0
    out = (abs(V-r)+V-r)/2;
else
    x = V-r;
    y = zeros(size(x));
    b = 1/smooth;
    c = log(2*b)/b;
    
    IX = find(x<=0);
    y(IX) = exp(b*(x(IX)-c));
    IX = find(x>0);
    y(IX) = x(IX) + exp(-b*(x(IX)+c));
    out = y;
end
end % rec

function val = std_excluNaN(V)
temp = V(find(isnan(V)==0));
if length(temp)>0
    val = std(temp);
else val = NaN;
end
end

function val = uptri_1d(M)
[r, c] = size(M);
val = zeros(1, (c-1+c-r)*r/2);
count = 1;
for i = 1:r
    for j = i+1:c
        val(count) = M(i,j);
        count = count+1;
    end
end
end % uptri_1d

function update_Iplot(cont_density,t_axis,I,I1,I2,Icorr,Iage,Wmg,N_popu,H)
if length(find(Iage>=0)) <= 0
    Iage(1) = 0;
end
max_I1_I2=max(max([I1;I2]));
set(H{1}, 'CData', SigmoidG(I));
set(H{2}, 'CData', I1);
if (~isnan(max_I1_I2))
    caxis([0 max_I1_I2]);
end
set(H{3}, 'CData', I2);
if (~isnan(max_I1_I2))
    caxis([0 max_I1_I2]);
end
set(H{4}, 'CData', Iage');
set(H{5}, 'CData', Icorr);
end % update_Iplot

function update_Pplot(N_mitral,P,Pcorr,Ccorr,rg,H,Wmg,Wgm)

set(H{1}, 'CData', P);
if (1==2)
    set(H{2}, 'CData', Ccorr);
else
    for iline=1:length(H{2})
        set(H{2}(iline), 'YData', P(1:N_mitral,iline));
    end
end
set(H{3}, 'CData', Pcorr);
if nargin > 6
    conn = Wmg*Wgm;
    set(H{4}, 'CData', conn-diag(NaN*diag(conn)));
end
end % update_Pplot

function update_Pplot2D(P,coord,H,Haxis,Pmin,Pmax)
[Nc, Ns] = size(P);
Cmin = min(coord')'; Cmax = max(coord')';
c = ceil(sqrt(Ns)); r = ceil(Ns/c);
for ii = 1:r
    for jj = 1:c
        kk = (ii-1)*c+jj;
        if kk <= Ns
            % subplot(r,c,kk);
            subplot(1,r*c,kk);
            temp = NaN*zeros(ceil(Cmax-Cmin)'+1);
            for ll = 1:Nc
                temp(ceil(coord(1,ll)-Cmin(1))+1, ceil(coord(2,ll)-Cmin(2))+1) = P(ll,kk);
            end
            set(H{kk}, 'CData', temp(:,end:-1:1)');
            if nargin < 6
                mm = mean_excluNaN(temp(:,end:-1:1));
                ss = std_excluNaN(temp(:,end:-1:1));
                Pmin = mm-3*ss;
                Pmax = mm+3*ss;
            end
            set(Haxis{kk}, 'CLim', [Pmin Pmax]);
        end
    end
end
end % setup_Pplot2D

function bargraph(m1, name1, m2, name2)
global pos
[x, y] = size(m1);
v1 = zeros(x,1);
v2 = zeros(x,1);

for i=1:y
    v1 = v1+m1(:,i);
    v2 = v2+m2(:,i);
end

fig40=figure(40);  set(fig40,'Position',pos{40});
subplot(1,2,1)
bar(1:x, v1)
title(name1)
subplot(1,2,2)
bar(1:x, v2)
title(name2)
end

function PC = P_connections(N, Wpg, S_mem, time_axis,N_t)

PC = cell(1,2);
figure(N); clf; set(gcf,'DoubleBuffer','on');

subplot(1,2,1)
conn = Wpg*(Wpg');
PC{1} = imagesc(conn-diag(NaN*diag(conn)));
title('PC via GC'); axis square;

subplot(1,2,2)
PC{2} = imagesc(S_mem);
title('S_{mem}'); axis square;

end

function conntrack(Wmg, Wpg)
global pos

MC12 = [];
MC13 = [];
MC14 = [];
MC23 = [];
MC24 = [];
MC34 = [];

[m, nGC]= size(Wmg);

connarray = zeros(6);

% Quick note on indexing: 1 indicates connections from GC to 1,2. 2 means
% GC to 1,3. 3 means GC to 1,4. 4 means GC to 2,3. 5 Means GC to 2,4. 6
% means GC to 3,4.

% Of course, whether the index is in the row r column affects whether we
% mean MC or PC, but that doesn't matter so much.

for i = 1:nGC
    if Wmg(1,i)==1;
        if Wmg(2,i)==1;
            MC12 = horzcat(MC12, [i]);
        end %if
        if Wmg(3,i)==1;
            MC13 = horzcat(MC13, [i]);
        end %if
        if Wmg(4,i)==1;
            MC14 = horzcat(MC14, [i]);
        end %if
    elseif Wmg(2,i)==1; %Because adding indices again is bad!
        if Wmg(3,i)==1;
            MC23 = horzcat(MC23, [i]);
        end %if
        if Wmg(4,i)==1;
            MC24 = horzcat(MC24, [i]);
        end %if
    elseif Wmg(3,i)==1;
        if Wmg(4,i)==1;
            MC34 = horzcat(MC34, [i]);
        end %if
    end %if
end %loop

% Now I need to see which PC the GC connect to

for j = 1:length(MC12)
    % Row 1 of connarray
    i = MC12(j);
    if Wpg(1,i)==1;
        if Wpg(2,i)==1;
            connarray(1,1) = connarray(1,1)+1;
        end %if
        if Wpg(3,i)==1;
            connarray(1,2) = connarray(1,2)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(1,3) = connarray(1,3)+1;
        end %if
    elseif Wpg(2,i)==1; %Because adding indices again is bad!
        if Wpg(3,i)==1;
            connarray(1,4) = connarray(1,4)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(1,5) = connarray(1,5)+1;
        end %if
    elseif Wpg(3,i)==1;
        if Wpg(4,i)==1;
            connarray(1,6) = connarray(1,6)+1;
        end %if
    end %if
end %loop

for j = 1:length(MC13)
    % Row 2 of connarray
    i = MC13(j);
    if Wpg(1,i)==1;
        if Wpg(2,i)==1;
            connarray(2,1) = connarray(2,1)+1;
        end %if
        if Wpg(3,i)==1;
            connarray(2,2) = connarray(2,2)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(2,3) = connarray(2,3)+1;
        end %if
    elseif Wpg(2,i)==1; %Because adding indices again is bad!
        if Wpg(3,i)==1;
            connarray(2,4) = connarray(2,4)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(2,5) = connarray(2,5)+1;
        end %if
    elseif Wpg(3,i)==1;
        if Wpg(4,i)==1;
            connarray(2,6) = connarray(2,6)+1;
        end %if
    end %if
end %loop

for j = 1:length(MC14)
    % Row 3 of connarrray
    i = MC14(j);
    if Wpg(1,i)==1;
        if Wpg(2,i)==1;
            connarray(3,1) = connarray(3,1)+1;
        end %if
        if Wpg(3,i)==1;
            connarray(3,2) = connarray(3,2)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(3,3) = connarray(3,3)+1;
        end %if
    elseif Wpg(2,i)==1; %Because adding indices again is bad!
        if Wpg(3,i)==1;
            connarray(3,4) = connarray(3,4)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(3,5) = connarray(3,5)+1;
        end %if
    elseif Wpg(3,i)==1;
        if Wpg(4,i)==1;
            connarray(3,6) = connarray(3,6)+1;
        end %if
    end %if
end %loop

for j = 1:length(MC23)
    % Row 4 of connarray
    i = MC23(j);
    if Wpg(1,i)==1;
        if Wpg(2,i)==1;
            connarray(4,1) = connarray(4,1)+1;
        end %if
        if Wpg(3,i)==1;
            connarray(4,2) = connarray(4,2)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(4,3) = connarray(4,3)+1;
        end %if
    elseif Wpg(2,i)==1; %Because adding indices again is bad!
        if Wpg(3,i)==1;
            connarray(4,4) = connarray(4,4)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(4,5) = connarray(4,5)+1;
        end %if
    elseif Wpg(3,i)==1;
        if Wpg(4,i)==1;
            connarray(4,6) = connarray(4,6)+1;
        end %if
    end %if
end %loop

for j = 1:length(MC24)
    % Row 5 of connarray
    i = MC24(j);
    if Wpg(1,i)==1;
        if Wpg(2,i)==1;
            connarray(5,1) = connarray(5,1)+1;
        end %if
        if Wpg(3,i)==1;
            connarray(5,2) = connarray(5,2)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(5,3) = connarray(5,3)+1;
        end %if
    elseif Wpg(2,i)==1; %Because adding indices again is bad!
        if Wpg(3,i)==1;
            connarray(5,4) = connarray(5,4)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(5,5) = connarray(5,5)+1;
        end %if
    elseif Wpg(3,i)==1;
        if Wpg(4,i)==1;
            connarray(5,6) = connarray(5,6)+1;
        end %if
    end %if
end %loop

for j = 1:length(MC34)
    % Row 6 of connarray
    i = MC34(j);
    if Wpg(1,i)==1;
        if Wpg(2,i)==1;
            connarray(6,1) = connarray(6,1)+1;
        end %if
        if Wpg(3,i)==1;
            connarray(6,2) = connarray(6,2)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(6,3) = connarray(6,3)+1;
        end %if
    elseif Wpg(2,i)==1; %Because adding indices again is bad!
        if Wpg(3,i)==1;
            connarray(6,4) = connarray(6,4)+1;
        end %if
        if Wpg(4,i)==1;
            connarray(6,5) = connarray(6,5)+1;
        end %if
    elseif Wpg(3,i)==1;
        if Wpg(4,i)==1;
            connarray(6,6) = connarray(6,6)+1;
        end %if
    end %if
end %loop

fig100=figure(100);  set(fig100,'Position',pos{100});
imagesc(connarray); axis square; colorbar
title('Rows: MC connections; Columns: PC connections')
end % conntrack

function [H, Iaxis] = setup_Iplot(cont_density,t_axis,I,I1,I2,Icorr,Iage,Wmg,N_popu,N)
global pos fig_label
Ns = size(I,2);
H = cell(1,5);
if length(find(Iage>=0)) <= 0
    Iage(1) = 0;
end
%N==3 for this plot
fig_label(N)=figure(N); set(fig_label(N),'Position',pos{N}); clf; set(gcf,'DoubleBuffer','on');
h1=subplot(1,5,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
H{1} = imagesc(SigmoidG(I)); colormap(jet); colorbar;  %caxis([0 1]);
Iaxis(1) = gca;
xlabel('stimulus'); ylabel('GC Index'); title('sigmoid(GC)');
h1=subplot(1,5,2); ax1 = get(h1,'position'); set(h1,'position', ax1);
H{2} = imagesc(I1); colormap(jet); colorbar;
Iaxis(2) = gca;
xlabel('stimulus'); ylabel('GC Index'); title('GC Input: Bulb');
h1=subplot(1,5,3); ax1 = get(h1,'position'); set(h1,'position', ax1);
H{3} = imagesc(I2); colormap(jet); colorbar;
Iaxis(3) = gca;
xlabel('stimulus'); ylabel('GC Index'); title('GC Input: Cortex');
h1=subplot(1,5,4); ax1 = get(h1,'position'); set(h1,'position', ax1);
H{4} = imagesc(Iage'); colormap(jet); colorbar;
Iaxis(4) = gca;
ylabel('GC Index'); title('GC Age');
h1=subplot(1,5,5); ax1 = get(h1,'position'); set(h1,'position', ax1);
H{5} = imagesc(Icorr); colormap(jet);
title('Pattern Correlation of GC'); xlim([0.5 Ns+.5]); ylim([0.5 Ns+.5]); caxis([-1 1]); axis square;
end % setup_Iplot

function [H P2daxis] = setup_Pplot2D(P,coord,N,name,Pmin,Pmax)
global pos

[Nc Ns] = size(P);
H = cell(1,Ns);
P2daxis = cell(1,Ns);
figN=figure(N); set(figN,'Position',pos{N}); clf; set(gcf,'DoubleBuffer','on');
%     M = [];

%     if nargin < 6
%         Pmin = min(min(P)); Pmax = max(max(P));
%     end
Cmin = min(coord')'; Cmax = max(coord')';
c = ceil(sqrt(Ns)); r = ceil(Ns/c);
for ii = 1:r
    for jj = 1:c
        kk = (ii-1)*c+jj;
        if kk <= Ns
            subplot(r,c,kk);
            % subplot(1,r*c,kk);
            temp = NaN*zeros(ceil(Cmax-Cmin)'+1);
            for ll = 1:Nc
                temp(ceil(coord(1,ll)-Cmin(1))+1, ceil(coord(2,ll)-Cmin(2))+1) = P(ll,kk);
            end
            H{kk} = imagesc(temp(:,end:-1:1)');%This must be each
            %                 M = [M temp(:,end:-1:1)'];
            P2daxis{kk} = gca;
            colormap(jet);%Sets the colours
            temp
            if nargin < 6
                mm = mean_excluNaN(temp(:,end:-1:1));
                ss = std_excluNaN(temp(:,end:-1:1));
                Pmin = mm-3*ss;
                Pmax = mm+3*ss;
            end
            caxis([Pmin Pmax]);
            set(gca, 'color', 'black');
            axis off; axis image;
            if nargin > 3
                title(name{kk});
            end
        end
    end
end

%     dlmwrite('S.txt',M)

end % setup_Pplot2D

function [S, coord, metric, name] = gara2(Ns,MCperGlom,name_list,sim)%,range)
global choose mix_factor probe
global stim_type
%XH
global N_mitral
% read raw pattern
name_list = strrep(name_list,' ','#');
rem = name_list;
[token, rem] = strtok(rem,'#');
temp = get_odor(token);
odor = zeros(Ns,size(temp,1),size(temp,2));
odor(1,:,:) = temp;
i = 1;

while length(rem)>1
    i = i+1;
    [token, rem] = strtok(rem,'#');
    temp = get_odor(token);
    odor(i,:,:) = temp;
end

if (1==1)
    
    if (1==1)
        N_mitral=64; Ns=max(choose);
        x=[1:N_mitral];
        
        
        if (stim_type==0)
            
            odor_redu(1,:) = double(x<=16);
            odor_redu(2,:) = double(x<=32)-odor_redu(1,:);
            odor_redu(3,:) = double(x<=48)-odor_redu(1,:)-odor_redu(2,:);
            odor_redu(4,:) = double(x<=64)-odor_redu(1,:)-odor_redu(2,:)-odor_redu(3,:);
            
            
            probe = odor_redu;
            
            
            probe = probe';
            
        end
        
        if (stim_type==1)
            
            amp=1; width=12;
            xpos=16;%24;%26;%16;%26;%16;
            odor_redu(1,:)=amp*exp(-(x-xpos).^2/width^2);
            xpos=20;%28;%30;%20;%30;%20;
            odor_redu(2,:)=amp*exp(-(x-xpos).^2/width^2);
            xpos=44;%44;
            odor_redu(3,:)=amp*exp(-(x-xpos).^2/width^2);
            xpos=48;%48;
            odor_redu(4,:)=amp*exp(-(x-xpos).^2/width^2);
            
            
            probe = zeros(4,N_mitral);
            probe(1,1:N_mitral) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (2-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            probe(2,:) = (2-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            
            probe(3,:) = odor_redu(1,:)+.1*odor_redu(3,:);
            probe(4,:) = odor_redu(1,:)+.1*odor_redu(4,:);
            
            probe = probe';
            
        end
        
        if (stim_type==2)
            
            amp=1; width=12;
            xpos=16;%26;%16;
            odor_redu(1,:)=rect(amp*(1-(x-xpos).^2/width^2));
            xpos=18;%30;%20;
            odor_redu(2,:)=rect(amp*(1-(x-xpos).^2/width^2));
            xpos=44;
            odor_redu(3,:)=rect(amp*(1-(x-xpos).^2/width^2));
            xpos=45;
            odor_redu(4,:)=rect(.8*amp*(1-(x-xpos).^2/width^2));
            
            
            probe = zeros(4,N_mitral);
            probe(1,1:N_mitral) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (2-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            probe(2,:) = (2-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            
            probe(3,:) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (1-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            probe(4,:) = (1-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            
            probe = probe';
            
        end
        
        
        if (stim_type==3)
            
            xpositions=[32]; widths=[12,12]; amplitudes=[1,.25];
            odor_redu(1,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[36]; widths=[12,12]; amplitudes=[1,.25];
            odor_redu(2,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[44]; widths=[12,12]; amplitudes=[1,1];
            odor_redu(3,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[46]; widths=[12,12]; amplitudes=[1,.5];
            odor_redu(4,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            probe = zeros(4,N_mitral);
            probe(1,1:N_mitral) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (2-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            probe(2,:) = (2-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            
            probe(3,:) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (1-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            probe(4,:) = (1-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
            
            probe = probe';
            
        end
        
        if (stim_type==4)
            
            xpositions=[19]; widths=[15]; amplitudes=[1];
            odor_redu(1,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            odor_redu(3,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[46]; widths=[15]; amplitudes=[1];
            odor_redu(2,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            odor_redu(4,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            probe = zeros(4,N_mitral);
            xpositions=[19]; widths=[15]; amplitudes=[1];
            probe(1,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[19,21]; widths=[15,4]; amplitudes=[1,.2];
            probe(2,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[33]; widths=[15]; amplitudes=[1];
            probe(3,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            xpositions=[33,35]; widths=[15,4]; amplitudes=[1,.2];
            probe(4,:)=gauss_mixture(x,amplitudes,xpositions,widths);
            
            probe = probe';
            
        end
        
    end
    
end

raw_Ns = size(odor,1);
name_list = strrep(name_list,'(_)','(-)');
name_list = strrep(name_list,'_',' ');
rem = name_list;
[token,rem] = strtok(rem,'#');
max_length = length(token);
for i = 2:raw_Ns
    [token,rem] = strtok(rem,'#');
    if length(token)>max_length
        max_length = length(token);
    end
end
rem = name_list;
name = cell(1,raw_Ns);
for i = 1:raw_Ns
    [token, rem] = strtok(rem,'#');
    temp = token;
    name{i} = temp;
end


% find revalent pix
odor_sum = reshape(sum(odor_redu,1),size(odor_redu,2),size(odor_redu,3));
odor_sum_1d = reshape(odor_sum, 1, size(odor_sum,1)*size(odor_sum,2));
IX = find(odor_sum_1d>-3);%finds entries of odor_sum_1d that are larger than -3
odor_sum_1d = odor_sum_1d(IX);

% cal coord1mc
coord1mc1 = (ones(size(odor_sum,2),1)*(size(odor_sum,1):-1:1))';
coord1mc2 = ((1:size(odor_sum,2))'*ones(1,size(odor_sum,1)))';
coord1mc1 = reshape(coord1mc1,1,size(coord1mc1,1)*size(coord1mc1,2));
coord1mc2 = reshape(coord1mc2,1,size(coord1mc2,1)*size(coord1mc2,2));
coord1mc = [coord1mc2; coord1mc1];
coord1mc = coord1mc(:,IX);

odor_redu1d = odor_redu;

coord1mc = coord1mc(:,IX(end:-1:1));

% metric
metric1mc = zeros(size(odor_redu1d,2),size(odor_redu1d,2));
for i = 1:size(odor_redu1d,2)
    for j = 1:size(odor_redu1d,2)
        metric1mc(i,j) = norm(coord1mc(:,i)-coord1mc(:,j));
    end
end

% return
if (1==1)
    fprintf(' Odors normalized \n')
    S1mc = odor_redu1d'/max(max(odor_redu1d));%Matrix dimensions error. Why?
else
    fprintf(' \n Odors not normalized !!!! \n \n')
    S1mc = odor_redu1d';
end
S1mc = S1mc(:,1:Ns);

% mult. MC per glom
coord = zeros(2,MCperGlom*size(coord1mc,2));
metric = zeros(MCperGlom*size(metric1mc,1));
S = zeros(MCperGlom*size(S1mc,1),size(S1mc,2));
for i = 1:MCperGlom
    coord(:,i:MCperGlom:end) = coord1mc;
    S(i:MCperGlom:end,:) = S1mc;
    metric(i:MCperGlom:end, 1:MCperGlom:end) = metric1mc;
end

for i = 2:MCperGlom
    metric(:, i:MCperGlom:end) = metric(:, 1:MCperGlom:end);
end

dlmwrite('Spiri1.txt',S)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,probe] = gara3(MCperGlom,stim_type,xpositions_learn,widths_learn,amplitudes_learn,odor_normalized)
%%
global mix_factor probe
%%
global N_mitral
global randomize random_order odor_random shuffle

N_mitral=110;%90;%64;
%Ns=max(choose);
x=[1:N_mitral];

if (stim_type==0)
    
    odor_redu(1,:) = double(x<=16);
    odor_redu(2,:) = double(x<=32)-odor_redu(1,:);
    odor_redu(3,:) = double(x<=48)-odor_redu(1,:)-odor_redu(2,:);
    odor_redu(4,:) = double(x<=64)-odor_redu(1,:)-odor_redu(2,:)-odor_redu(3,:);
    
    probe = odor_redu;
    probe = probe';
    
end

if (stim_type==1)
    Ns=size(xpositions_learn,1);
    permutations=cell(1,Ns);
    for is=1:Ns
        permutations{is}=randperm(N_mitral);
        xpositions=xpositions_learn(is,:); widths=widths_learn(is,:); amplitudes=amplitudes_learn(is,:);
        odor_raw=gauss_mixture(x,amplitudes,xpositions,widths);
        if (shuffle(is)>0)
            odor_redu(is,:)=odor_raw(permutations{shuffle(is)});
        else
            odor_redu(is,:)=odor_raw;
        end
    end
    
end

if (stim_type==2)
    
    amp=1; width=12;
    xpos=16;%26;%16;
    odor_redu(1,:)=rect(amp*(1-(x-xpos).^2/width^2));
    xpos=18;%30;%20;
    odor_redu(2,:)=rect(amp*(1-(x-xpos).^2/width^2));
    xpos=44;
    odor_redu(3,:)=rect(amp*(1-(x-xpos).^2/width^2));
    xpos=45;
    odor_redu(4,:)=rect(.8*amp*(1-(x-xpos).^2/width^2));
    
    
    probe = zeros(4,N_mitral);
    probe(1,1:N_mitral) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (2-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    probe(2,:) = (2-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    
    probe(3,:) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (1-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    probe(4,:) = (1-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    
    probe = probe';
    
end


if (stim_type==3)
    
    xpositions=[32]; widths=[12,12]; amplitudes=[1,.25];
    odor_redu(1,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[36]; widths=[12,12]; amplitudes=[1,.25];
    odor_redu(2,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[44]; widths=[12,12]; amplitudes=[1,1];
    odor_redu(3,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[46]; widths=[12,12]; amplitudes=[1,.5];
    odor_redu(4,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    probe = zeros(4,N_mitral);
    probe(1,1:N_mitral) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (2-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    probe(2,:) = (2-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    
    probe(3,:) = mix_factor*(odor_redu(1,:)+odor_redu(2,:))/2 + (1-mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    probe(4,:) = (1-mix_factor)*(odor_redu(1,:)+odor_redu(2,:))/2 + (mix_factor)*(odor_redu(3,:)+odor_redu(4,:))/2;
    
    probe = probe';
    
end

if (stim_type==4)
    
    xpositions=[19]; widths=[15]; amplitudes=[1];
    odor_redu(1,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    odor_redu(3,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[46]; widths=[15]; amplitudes=[1];
    odor_redu(2,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    odor_redu(4,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    probe = zeros(4,N_mitral);
    xpositions=[19]; widths=[15]; amplitudes=[1];
    probe(1,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[19,21]; widths=[15,4]; amplitudes=[1,.2];
    probe(2,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[33]; widths=[15]; amplitudes=[1];
    probe(3,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    xpositions=[33,35]; widths=[15,4]; amplitudes=[1,.2];
    probe(4,:)=gauss_mixture(x,amplitudes,xpositions,widths);
    
    probe = probe';
    
end


if (stim_type==5)
    
    if (randomize==1)
        
        Ns_ensemble=3;
        for is=1:Ns_ensemble
            xpositions=xpositions_learn(is,:); widths=widths_learn(is,:); amplitudes=amplitudes_learn(is,:);
            odor_random(is,:)=gauss_mixture(x,amplitudes,xpositions,widths);
        end
        
        
        for is=1:Ns_ensemble
            odor_random(is,:)=odor_random(is,randperm(N_mitral));
        end
        
        odor_diff=odor_random(2,:)-odor_random(3,:);
        [~,sort_indices]=sort(odor_diff);
        for is=1:Ns_ensemble
            odor_random(is,:)=odor_random(is,sort_indices);
        end
        
        noise_level=0.1;
        odor_redu(1,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(2,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(3,:)=odor_random(2,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(4,:)=odor_random(2,:).*(1+noise_level*rand(size(odor_random(1,:))));
        Ns=4;
    end
    if (randomize==2)
        noise_level=0.1;
        odor_redu(1,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(2,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(3,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(4,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        Ns=4;
    end
    if (randomize==3)
        noise_level=0.1;
        odor_redu(1,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(2,:)=odor_random(1,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(3,:)=odor_random(2,:).*(1+noise_level*rand(size(odor_random(1,:))));
        odor_redu(4,:)=odor_random(2,:).*(1+noise_level*rand(size(odor_random(1,:))));
        Ns=4;
    end
    if (randomize==4)
        ninter=6;
        for inter=1:ninter
            odor_redu(inter,:)=(1-(inter-1)/(ninter-1))*odor_random(2,:)+(inter-1)/(ninter-1)*odor_random(3,:);
        end
        Ns=ninter;
    end
end

if (stim_type==6)
    Ns=size(xpositions_learn,1);
    permutations=cell(1,Ns);
    for is=1:Ns
        permutations{is}=randperm(N_mitral);
        xpositions=xpositions_learn(is,:); widths=widths_learn(is,:); amplitudes=amplitudes_learn(is,:)
        odor_raw=linear_mixture(x,amplitudes,xpositions,widths);
        if (shuffle(is)>0)
            odor_redu(is,:)=odor_raw(permutations{shuffle(is)});
        else
            odor_redu(is,:)=odor_raw;
        end
    end
    
end


if (odor_normalized==1)
    fprintf(' Odors normalized \n')
    S1mc = odor_redu'/max(max(odor_redu));%Matrix dimensions error. Why?
else
    fprintf(' \n Odors not normalized !!!! \n \n')
    S1mc = odor_redu';
end
S1mc = S1mc(:,1:Ns);

% mult. MC per glom
S = NaN(MCperGlom*size(S1mc,1),size(S1mc,2));
for i = 1:MCperGlom
    S(i:MCperGlom:end,:) = S1mc;
end

dlmwrite('Spiri1.txt',S)

end

%%%%%%%%%%%%%%%%%%%%%%%

function stimulus=gauss_mixture(x,amplitudes,xpositions,widths)
n_components=length(xpositions);
stimulus=0*x;
for icomp=1:n_components
    stimulus=stimulus+amplitudes(icomp)*exp(-(x-xpositions(icomp)).^2/widths(icomp)^2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%

function stimulus=linear_mixture(x,amplitudes,xpositions,widths)
n_components=length(xpositions);
stimulus=0*x;
exponent=0.5;
for icomp=1:n_components
    stimulus=stimulus+amplitudes(icomp)*max((1-(abs((x-xpositions(icomp))/widths(icomp))).^exponent),0);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%

function S = get_odor(name)

odor_raw = imread(strcat('odor/',name,'.png'));
odor = zeros(size(odor_raw(:,:,1)));

for i = 1:size(odor,1)
    for j = 1:size(odor,2)
        if odor_raw(i,j,1) == 255
            if odor_raw(i,j,2) == 255
                if odor_raw(i,j,3) == 255
                    odor(i,j) = -2;
                end
            end
        end
        if odor_raw(i,j,1) == 0
            if odor_raw(i,j,2) == 0
                if odor_raw(i,j,3) == 255
                    odor(i,j) = -1;
                end
            end
        end
        if odor_raw(i,j,1) == 0
            if odor_raw(i,j,2) == 255
                if odor_raw(i,j,3) == 255
                    odor(i,j) = 0;
                end
            end
        end
        if odor_raw(i,j,1) == 0
            if odor_raw(i,j,2) == 193
                if odor_raw(i,j,3) == 0
                    odor(i,j) = 1;
                end
            end
        end
        if odor_raw(i,j,1) == 255
            if odor_raw(i,j,2) == 255
                if odor_raw(i,j,3) == 0
                    odor(i,j) = 2;
                end
            end
        end
        if odor_raw(i,j,1) == 255
            if odor_raw(i,j,2) == 133
                if odor_raw(i,j,3) == 0
                    odor(i,j) = 3;
                end
            end
        end
        if odor_raw(i,j,1) == 255
            if odor_raw(i,j,2) == 0
                if odor_raw(i,j,3) == 0
                    odor(i,j) = 4;
                end
            end
        end
        if odor_raw(i,j,1) == 0
            if odor_raw(i,j,2) == 0
                if odor_raw(i,j,3) == 0
                    odor(i,j) = 5;
                end
            end
        end
    end
end

odor(1,:) = NaN;
odor(:,1) = NaN;
for i = 2:size(odor,1)
    for j = 2:size(odor,2)
        if odor(i,j) == -2
            if isnan(odor(i-1,j))
                odor(i,j) = NaN;
            elseif isnan(odor(i,j-1))
                odor(i,j) = NaN;
            end
        end
    end
end
odor(end,:) = NaN;
odor(:,end) = NaN;
for i = size(odor,1)-1:-1:1
    for j = size(odor,2)-1:-1:1
        if odor(i,j) == -2
            if isnan(odor(i+1,j))
                odor(i,j) = NaN;
            elseif isnan(odor(i,j+1))
                odor(i,j) = NaN;
            end
        end
    end
end

S = odor;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S_mem = pirirecall(S_recall)

% currently not used.
global tau_P W_PM W_PP_i W_PP M_i c_assoc thresh steep P_amp
global n_P Nc

%Need matrix W_PP
%Need matrix W_PP_i
%Need matrix W_PM

%Need n_P
%S_recall = dlmread('Summer/odortext/S.txt');

%Inputs
%[n_P, n_recall] = size(S_recall); XH
W_PP = dlmread('odortext/WPP.txt');
% W_PP_i = dlmread('Summer/odortext/WPPi.txt');
W_PM = dlmread('odortext/WPM.txt');

%w_inhib=0.4;% originally 0.8
%W_PP_i=w_inhib*ones(n_P,n_P);

% sigmoid for P activation
%thresh=0.7; P_amp=1; steep=5;%original thresh 0.7, original steep 20
%This is to be the same as in pirilearn at the moment, but might have to
%change.

%Constants
tmax = 100;

% Recall
c_assoc=0.8;
M_recall= S_recall;%Needs to be the same M as calculated in cal_activity.
n_noise=20;

for i_noise=1:1%n_noise
    
    M_pert=M_recall+(i_noise-1)/(n_noise-1);%*M_noise;
    %M_pert needs to come from cal_activity.
    for i_recall=1:n_recall
        M_i=M_pert(:,i_recall);%M_i needs to come from cal_activity
        P_init=zeros(n_P+Nc,1);
        stim_select=S_recall(:,i_recall);
        options = odeset('RelTol',1e-4,'AbsTol',1e-5);
        [time,Y]=ode45(@P_rate,[0,tmax],P_init,options);
        
        P_recall(:,i_recall)=Y(end,Nc+1:Nc+n_P)';
        %         Ptest=Y(end,:)';
        %         exc=((1-c_assoc)*W_PP)*Sigmoid(Ptest)/tau_P;
        %         inh=(-W_PP_i)*Sigmoid(Ptest)/tau_P;
        %         fprintf(' i = %g  \n',i_recall)
        %         fprintf(' exc\n');
        %         fprintf('  %g %g %g %g %g \n',exc);
        %         fprintf(' %g %g %g %g %g\n',inh);
        % end
        
        if (i_noise==1)
            cmax_M=max(max(M_pert));
            cmax_P=max(max(P_recall));
            M0=M_pert;
            P0=P_recall;
        end
        
    end
    
end

S_mem = Sigmoid(P_recall);
% S_mem = P_recall;

fprintf('Recall Complete\n');

if 1==2
    figure(13)
    imagesc(W_PP)
    
    figure(11)
    imagesc(M_pert)
    caxis=[0 cmax_M];
    % colorbar
    title(' Stimuli')
    
    figure(12)
    imagesc(S_mem)
    caxis=[0 cmax_P];
    % % colorbar
    title(' Recall')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pirilearn(S_learn,run,file_location,file_tag,istep,label_connectivity)

% program for a simple associate memory model
% the neurons (with activity P) receive input from the sensory neurons via the connectivity matrix W_PM
% and excite each other via the connectivity matrix W_PP. They also inhibit
% each other via W_PP _i
% the connectivity matrix W_PP is modified during the learning process
% (the connections are plastic)
% the level of the reward N controls whether during the learning
% connections are strengthened or weakened (in an attempt to model
% active forgetting of a memory)

global W_PM W_PP_i W_PP M_i stim_select W_PM_scale W_PP_max
global file_type_figure
global n_P Nc learning
global pos
global I Iact  num_light_cells
global c_assoc c_assoc_learn c_assoc_recall
global Gthresh

global tau_P wt1 CS
global   Wmg Wgm Wgp
%global thresh steep P_amp
global thresh_odor steep_odor P_amp_odor
global thresh_light steep_light P_amp_light

% in pirilearn associative fibers are always weak
c_assoc=c_assoc_learn;

% if (1==1) % the random number generators should not be reinitialized in the middle of the code
% if (1==2)
%     rand('twister',sum(100*clock));   % added on Juy 18, 15:30
%     randn('state',sum(100*clock));
% else
%     rand('twister',18);   % added on Juy 18, 15:30
%     randn('state',7);
% end
% end

%Number of neurons, learn stimuli and recall stimuli
[~, n_S] = size(S_learn); %XH
%n_P = n_M;

% fixed connectivity matrices
%W_PM=eye(n_P,n_M); % Matrix with ones on the main diagonal (In this case, identity matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W_PM(17:end,:) = 0;
%w_inhib=0.4;% originally 0.8
%w_inhib = inhiblist(j);
%W_PP_i=w_inhib*ones(n_P,n_P);
%W_PP_i=W_PP_i-diag(diag(W_PP_i));
%W_PP_i(1:n_P/2,1+n_P/2:end) = 0;
%W_PP_i(1+n_P/2:end,1:n_P/2) = 0;
% plastic connectivity
W_PP_min=0;
sparseness = 0;
w_init=W_PP_max;%0.1;%0.2;
wpm_init=.5;%.27;
%w_init=0.02;

if run==0
    %thresh = .025;
    if (1==2)
        MtoP = 4;
        P = rand(n_P-num_light_cells,Nc-num_light_cells);
        A = P';
        %fixed number connections;
        sortedValues = sort(A);          %# Unique sorted values
        maxValues = sortedValues(end-MtoP+1:end,:);  %# Get the 5 largest values
        maxIndex = ismember(A,maxValues)';     %# Get a logical index of all value
        str = 1;
        W_PM_Anatomy = zeros(n_P-num_light_cells,Nc-num_light_cells);
        W_PM_Anatomy(maxIndex) =str;
    end
    
    n_learn = 100;
    
    
    % random connection
    %     thresh = .1;
    %     W_PM_Anatomy = rand(n_P,Nc);
    %     W_PM_Anatomy = double(W_PM_Anatomy<thresh);
    %
    
    % 1-1 connection
    % W_PM=eye(n_P,n_M);
    W_PM = eye(n_P,Nc);
    %  W_PM(1:n_P-num_light_cells, 1:Nc-num_light_cells) = wpm_init*ones(n_P-num_light_cells,Nc-num_light_cells).*W_PM_Anatomy;
    W_PP = w_init*ceil(rand(n_P,n_P)-sparseness);
end
if run==1
    n_learn = 300;
    % W_PM = W_PM*1.7;
    % W_PM_Anatomy = dlmread('odortext/W_PM_Anatomy.txt');
    % W_PM_Anatomy = double(W_PM>0);%THIS IS WRONG!
    
    %W_PP=dlmread('odortext/WPP.txt');
end
if run==2
    %    W_PM=eye(n_P,n_M);
    W_PP=w_init*ceil(rand(n_P,n_P)-sparseness);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W_PM MAX
W_PM_max=10000;
W_PM_min=-10000;
% learning
eta_M = .005;
eta=10; %Change eta; make it run more quickly
N_0=1;  N_reward=1.8;   N_noreward=0.95;
learn_forget=0;   % 1 = learn and forget ; 0 = learn only
%c_average=ones(n_P,1)*1;
%c_0 =.75;%6;%2;%4;
%p=2;
%tau=10;
Omega_learn=0.2;%1;%.5;
Omega_forget=2;
%Omega = Omega_learn;
Omega=Omega_learn*ones(1,n_S); %1;
Omega(learning<1)=Omega_forget;
w_pp_decay=0.1;%0.0025;%0.1;%0.01;


%n_learn=5000;%1000;
tol_learn=2e-12;%1e-12;%1e-11;
fprintf(' low values for speed: n_learn= %g  and tol = %g \n',n_learn,tol_learn);

% figure(22)
% imagesc(S_learn); colorbar

% sigmoid for P activation
% thresh=1; P_amp=1; steep=10;

%thresh=0.6; P_amp=1; steep=10;

% time scales
tmax=500;

reward=ones(1,n_S);
fprintf(' Pirilearn with Omega = %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',Omega)
fprintf('\n')

% pre-allocate memory
P_all=NaN(n_P,n_S);

% figure(11)
% imagesc(S_learn); axis square;
% colorbar

fig5=figure(5);  set(fig5,'Position',pos{5});

Y_init=zeros(Nc+n_P,n_S);

for i_learn=1:n_learn % for loop iterates through n_learn times
    
    fprintf(' Learning i_learn = %g \n', i_learn);
    
    dW_PP=zeros(n_P,n_P); %Matrix filled with zeroes?
    
    for i_S=1:n_S
        
        if (1==2)
            if (learning(i_S)==0) % Does this punish neurons for not learning?
                reward(i_S)=-reward(i_S);
            end
            if (reward(i_S)==1)
                N=N_reward;
            else
                N=N_noreward;
            end
        else
            N=N_reward;
        end
        stim_select=S_learn(:,i_S);
        %P_init=zeros(Nc+n_P,1);
        options = odeset('RelTol',1e-4,'AbsTol',1e-5);
        [~,Y]=ode15s(@(t,y) P_rate(t,y,Nc,n_P,P_amp_odor,steep_odor,thresh_odor,P_amp_light,...
            steep_light,thresh_light,Wgm,wt1,Wgp,W_PM,c_assoc,W_PP,W_PP_i,tau_P,stim_select,CS,Wmg,Gthresh),[0,tmax],Y_init(:,i_S),options); %Solves an ODE!
        %         if (1==2)
        %             figure(1)
        %             subplot(n_S,1,i_S)
        %             imagesc(Y')
        %             colorbar
        %             xlabel('time steps')
        %             ylabel(' P number')
        %         end
        P=Y(end,Nc+1:Nc+n_P)';
        %fprintf(' P(end,10)=%16.10f P(end,30)=%16.10f \n',P(10),P(30))
        
        %         if mod(n_learn,500) == 0
        %             % ,kpause;
        %         end
        %        dW_PP=dW_PP+eta*(N-N_0)*(P-Omega)*Sigmoid(P');
        %dW_PP=dW_PP+eta*((N-N_0)*(P-Omega(i_S))*Sigmoid(P')-w_pp_decay);
        dW_PP=dW_PP+eta*((N-N_0)*(Sigmoid(P)-Omega(i_S))*Sigmoid(P')-w_pp_decay);
        %         if max(eta*(N-N_0)*(P-Omega)*Sigmoid(P')) < .05
        %            max(eta*(N-N_0)*(P-Omega)*Sigmoid(P'))
        %         end
        P_all(:,i_S)=P;
        Y_init(:,i_S)=Y(end,:);
        
    end
    %Shows Associations
    h1=subplot_tight(1,4,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    imagesc(W_PP);  axis square
    title(' W_{CC}')
    xlabel('CC Index')
    ylabel('CC Index');
    colorbar
    h1=subplot_tight(1,4,2,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    imagesc(P_all); axis square;
    title(' CC Activity')
    ylabel('CC Index');
    xlabel('Stimulus Index');
    colorbar
    h1=subplot_tight(1,4,3,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    imagesc(Sigmoid(P_all)); axis square;
    title(' Sigmoid(CC)')
    xlabel('Stimulus Index');
    ylabel('CC Index');
    colorbar
    h1=subplot_tight(1,4,4,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    imagesc(Y_init); axis square;
    title(' MC and CC combined')
    xlabel('Stimulus Index');
    colorbar
    drawnow
    
    W_PP_n=W_PP+dW_PP;
    W_PP_n=W_PP_n-diag(diag(W_PP_n));
    W_PP_n=min(W_PP_n,W_PP_max);
    W_PP_n=max(W_PP_n,W_PP_min);
    
    diff_dW=max(max(abs(W_PP_n-W_PP)))/n_S/n_P^2/eta;
    % max_W=max(max(W_PP));
    %    if (diff_dW<tol_learn*n_S*n_P^2 & max_W == W_PP _max)
    W_PP=W_PP_n;
    
    if (diff_dW<tol_learn*n_S*n_P^2 && i_learn>=5)
        fprintf(' learning has converged \n');
        if (learn_forget==1)
            if (learning==1)
                learning=0;
                fprintf(' now start forgetting \n');
            else
                break
            end
        else
            break
        end
    end
end

if (1==2)
    % by hand: no connections back to the light cells?
    W_PP(65:end,:) = 0;    % do we always want this? HR 2015-06-08
    fprintf(' feedback to the light cells is set to 0 !!!\n');
end

%Shows Associations
h1=subplot_tight(1,4,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
imagesc(W_PP);  axis square
title(' W_{CC}')
xlabel('CC Index')
ylabel('CC Index');
colorbar
h1=subplot_tight(1,4,2,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
imagesc(P_all); axis square;
title(' CC Activity')
ylabel('CC Index');
xlabel('Stimulus Index');
colorbar
h1=subplot_tight(1,4,3,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
imagesc(Sigmoid(P_all)); axis square;
title(' Sigmoid(CC)')
ylabel('CC Index');
xlabel('Stimulus Index');
colorbar
h1=subplot_tight(1,4,4,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
imagesc(Y_init); axis square;
title(' MC and CC combined')
xlabel('Stimulus Index');
colorbar
drawnow

if (i_learn>=n_learn)
    fprintf(' Learning has not converged \n');
end

fprintf(' diff_dW = % g i_learn = %g \n',diff_dW,i_learn);
fprintf('Learning Complete \n');

figname = 'WPP';
file_name = strcat(file_location,figname,'_',file_tag,'_i',num2str(istep),label_connectivity,'.txt');
dlmwrite(file_name,W_PP)
figname = 'WPPi';
file_name = strcat(file_location,figname,'_',file_tag,'_i',num2str(istep),label_connectivity,'.txt');
dlmwrite(file_name,W_PP_i)
figname = 'WPM';
file_name = strcat(file_location,figname,'_',file_tag,'_i',num2str(istep),label_connectivity,'.txt');
dlmwrite(file_name,W_PM)
figname = 'P_Cell_Network_and_Activity_during_Learning';
file_name = strcat(file_location,figname,'_',file_tag,'_i',num2str(istep),label_connectivity,file_type_figure);
%saveas(gcf,file_name,file_type)
hgsave(file_name);

% outside pirilearn teh associative fibers shold be strong
c_assoc=c_assoc_recall;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MC,GC,GC1,GC2,S_mem,P_recall] = cal_activity1(CS,Wpg,Wmg,Wgm,S,P,time)

global A wt1 wt1_initial wt1_later prod_gp prod_gm stim_select
global time_wt1_1 time_wt1_2 P_scale Nc
global n_P num_light_cells Wgp
global P_last log_first_run log_parallel
global Gthresh

global tau_P W_PM W_PP_i W_PP c_assoc
global P_amp_odor steep_odor thresh_odor P_amp_light steep_light thresh_light

% why is Wpg even defined? It is never needed
Wgp = Wpg';
% The product Wmg*Wpg' comes out as a 0x0 matrix. What?
% Aha! It's because I call cal_activity1 before any cells have been
% created. Let's change that.
% No! main makes Wmg and Wpg before cal_activity1 is called, so there
% shouldn't be a problem.

% Because some things cannot be global, I need to save them wip_rateP_rate

% prod_gp = Wmg*Wpg';
% prod_gm = Wmg*Wgm;
% A = eye(size(S,1))+CS*prod_gm;

if time <= time_wt1_1
    wt1=wt1_initial;
elseif time <= time_wt1_2
    wt1=wt1_initial+(wt1_later-wt1_initial)*(time-time_wt1_1)/(time_wt1_2-time_wt1_1);
else
    wt1=wt1_later;
end

% sigmoid for P activation
%thresh=0.8; P_amp=1; steep=5;%original thresh 0.7, original steep 20
%Changing these doesn't seem to have any effect on the end correlations.

% Recall
[~,n_recall] = size(S);%XH
% It would be better to pass n_M=Nc and n_P via globals. I don't think either of these variables should
% vary from one call of cal_activity to the other
% n_recall should be read off as it is size(S)
% HR the size of P_recall should not be Nc but n_P.
MC=zeros(Nc,n_recall); P_recall=zeros(n_P,n_recall);

%Constants
tmax = 500;

%        MC = A\S;%Put the recalled stimulus in here instead of S
%        MC = A\(S-Wmg*Wpg*S_mem)%According to meeting Jun 27th. Matrix
%        sizes are mismatched
%MC = A\S;%This should make a 26x26 matrix that is right-multiplied by S_mem

if log_first_run == 1
    P_last = zeros(n_P+Nc,n_recall);
    log_first_run = 0;
end

if (1==2)
    for i_recall=1:n_recall
        P_init = P_last(:,i_recall);
        
        lastP = P_init(:,1);
        stim_select=S(:,i_recall);
        
        
        options = odeset('RelTol',1e-4,'AbsTol',1e-5','OutputFcn',@out_fun2);
        %        options = odeset('RelTol',1e-10,'AbsTol',1e-11);
        [~,Y]=ode45(@(t,y) P_rate(time,y,Nc,n_P,P_amp_odor,steep_odor,thresh_odor,P_amp_light,...
            steep_light,thresh_light,Wgm,wt1,Wgp,W_PM,...
            c_assoc,W_PP,W_PP_i,tau_P,stim_select,CS,Wmg,Gthresh),[0,tmax],P_init,options);
        %[~,Y]=ode15s(@P_rate,[0,tmax],P_init,options);
        
        %        figure(1000)
        %        plot(time,Y(:,10),time,Y(:,60),time,Y(:,80),time,Y(:,100))
        P_recall(:,i_recall)=Y(end,Nc+1:Nc+n_P)'/P_scale;
        MC(:,i_recall)=Y(end,1:Nc)';
        %         figure(1002)
        %         plot(MC(:,i_recall))
        %         hold all
        %         plot(P_recall(:,i_recall))
        %         hold off
    end %for loop
end
% make local variables for parfor loop
P_last_l=P_last; Nc_l=Nc;n_P_l=n_P;P_amp_odor_l=P_amp_odor;steep_odor_l=steep_odor;
thresh_odor_l=thresh_odor;P_amp_light_l=P_amp_light;steep_light_l=steep_light;thresh_light_l=thresh_light;
wt1_l=wt1; W_PM_l=W_PM; W_PP_l=W_PP;W_PP_i_l=W_PP_i;c_assoc_l=c_assoc;tau_P_l=tau_P;Gthresh_l=Gthresh;P_scale_l=P_scale;

if (log_parallel==1)
    parfor i_recall=1:n_recall
        %        fprintf(' i = %g time = %g %g %g %g %g %g \n',i_recall, clock)
        [P_recall_i,MC_i]=cal_activity_i(S,P_last_l,i_recall,time,Nc_l,n_P_l,P_amp_odor_l,...
            steep_odor_l,thresh_odor_l,P_amp_light_l,steep_light_l,thresh_light_l,Wgm,wt1_l,Wgp,W_PM_l,c_assoc_l,W_PP_l,W_PP_i_l,tau_P_l,CS,Wmg,Gthresh_l,P_scale_l);
        
        P_recall(:,i_recall)=P_recall_i;
        MC(:,i_recall)=MC_i;
        
    end %for loop
else
    for i_recall=1:n_recall
        
        [P_recall_i,MC_i]=cal_activity_i(S,P_last_l,i_recall,time,Nc_l,n_P_l,P_amp_odor_l,...
            steep_odor_l,thresh_odor_l,P_amp_light_l,steep_light_l,thresh_light_l,Wgm,wt1_l,Wgp,W_PM_l,c_assoc_l,W_PP_l,W_PP_i_l,tau_P_l,CS,Wmg,Gthresh_l,P_scale_l);
        
        P_recall(:,i_recall)=P_recall_i;
        MC(:,i_recall)=MC_i;
        
    end %for loop
end


P_last = [MC;P_recall];

S_mem = Sigmoid(P_recall);
%    MC = A\(S-wt1*CS*prod*S_mem);

%XH 4/23
%GC1 = Wgm*SigmoidM(MC);
GC1 = Wgm*rect(MC);
GC2 = wt1*Wpg'*S_mem(1:n_P,:);
%     GC1 = Wgm*rect(MC);
%     GC2 = wt1*Wpg'*S_mem;

GC=GC1+GC2;

if (1==2)
    figure(1001)
    subplot(4,1,1)
    imagesc(GC1)
    subplot(4,1,2)
    imagesc(GC2)
    subplot(4,1,3)
    imagesc(MC)
    colorbar
    subplot(4,1,4)
    imagesc(S_mem)
    figure(1003)
    subplot(2,1,1)
    imagesc(Wgm)
    subplot(2,1,2)
    imagesc(Wpg)
end

end %cal_activity1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_recall_i,MC_i]=cal_activity_i(S,P_last,i_recall,time,Nc,n_P,P_amp_odor,...
    steep_odor,thresh_odor,P_amp_light,steep_light,thresh_light,Wgm,wt1,Wgp,W_PM,c_assoc,W_PP,W_PP_i,tau_P,CS,Wmg,Gthresh,P_scale)
% global stim_select
%
% global tau_P W_PM W_PP_i W_PP c_assoc stim_select wt1 CS
% global n_P Nc Wmg Wgm Wgp
% %global thresh steep P_amp
% global thresh_odor steep_odor P_amp_odor
% global thresh_light steep_light P_amp_light
% global Gthresh

%Constants
tmax = 500;

stim_select=S(:,i_recall);
P_init = P_last(:,i_recall);

lastP = P_init(:,1);

options = odeset('RelTol',1e-4,'AbsTol',1e-5','OutputFcn',@out_fun2);
%        options = odeset('RelTol',1e-10,'AbsTol',1e-11);
[~,Y]=ode45(@(t,y) P_rate(time,y,Nc,n_P,P_amp_odor,steep_odor,thresh_odor,P_amp_light,steep_light,...
    thresh_light,Wgm,wt1,Wgp,W_PM,c_assoc,W_PP,W_PP_i,tau_P,stim_select,CS,Wmg,Gthresh),[0,tmax],P_init,options);
%[~,Y]=ode45(@P_rate,[0,tmax],P_init,options);
%[~,Y]=ode15s(@P_rate,[0,tmax],P_init,options);

%        figure(1000)
%        plot(time,Y(:,10),time,Y(:,60),time,Y(:,80),time,Y(:,100))
P_recall_i=Y(end,Nc+1:Nc+n_P)'/P_scale;
MC_i=Y(end,1:Nc)';

    function stop_ = out_fun2(t_,P_,flag_)
        stop_ = false;
        if strcmp(flag_,'init')
        elseif strcmp(flag_,'done')
        else
            if size(P_,2) > 1
                cond_ = norm(lastP-P_(:,end));
            else
                cond_ = norm(lastP-P_);
            end
            %if cond_<1e-2
            if cond_<1e-4
                stop_ = true;
            end
            if size(P_,2) > 1
                lastP = P_(:,end);
            else
                lastP = P_;
            end
        end
    end % out_fun

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,J_thresh,J_comp,J_comp_thresh] = fischer(m1, m2)

%Sb = (m2-m1)*(m2-m1)';
thresh=0.1;

m_max=max([m1;m2]);

Sw= diag(m1+m2);
Sw_inv=diag(1./(m1+m2));

Dm = (m2-m1);
J = Dm'*Sw_inv*Dm;
J_comp=Sw_inv*Dm.^2;
indices=find(m1+m2>m_max*thresh);
J_comp_thresh=zeros(length(m1),1);
J_comp_thresh(indices)=J_comp(indices);
J_thresh=sum(J_comp_thresh);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,J_comp] = fisher_rect(m1,m2,fisher_thresh)

%m1=SigmoidM(m1); m2=SigmoidM(m2);
m1=rect(m1); m2=rect(m2);
Dm = m2-m1;
Sm = m1+m2;

m_max=max([m1;m2]);
indices=find(Sm>fisher_thresh*m_max);
J_comp=zeros(length(m1),1);
J_comp(indices)=Dm(indices).^2./Sm(indices);
J=sum(J_comp);

end

%%%%%%%%%%%%%%%%%%%%%%%%

function [fisher_mean,fisher_std,fisher_quantiles]=fisher_non_optimal(m1,m2,fisher_thresh)

global read_out_weights

%m1=SigmoidM(m1); m2=SigmoidM(m2);
m1=rect(m1); m2=rect(m2);

% indices=find(m1+m2>fisher_thresh);
% w_opt=zeros(size(m1));
% w_opt(indices)=(m1(indices)-m2(indices))./(m1(indices)+m2(indices));
% w_opt'
% read_out_weights1=[read_out_weights;w_opt'];
% read_out_weights1=w_opt'
quantiles=[.25,.5,.75];

Dm=(read_out_weights*(m1-m2)).^2;
Sm=read_out_weights.^2*(m1+m2);
indices=find(Sm>fisher_thresh);
fisher_all=zeros(size(Sm));
fisher_all(indices)=Dm(indices)./Sm(indices);
fisher_mean=mean(fisher_all);
fisher_std=std(fisher_all);
fisher_quantiles=quantile(fisher_all,quantiles);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time_fisher,y_fisher,count]=fisher_diagnostic(time,P,compare,fisher_thresh,count,time_fisher,y_fisher,fig_number)

global N_mitral
global pos fig_label
global fisher_plot

% output:
% 3,4=max and min of SigmoidM
% 7,8=Fisher with SigmoidM;
% 9,10=correlation SigmoidM
% 11,12=max and min of M
% 13,14 Fisher mean and std pair 1
% 15,16,17 Fisher quantiles
% 18,19 Fisher mean and std pair 2
% 20,21,22 Fisher quantiles

% count<0 indicates that a second measurement is done for the same time
if (count<0)
    count=-count;
else
    count=count+1;
    time_fisher(count) = time;
end

%temp1 = corrcoef(SigmoidM(P(1:N_mitral,:)));
temp1 = corrcoef(rect(P(1:N_mitral,:)));
y_fisher(count,9) = temp1(compare(1),compare(2));
y_fisher(count,10) = temp1(compare(3),compare(4));
temp1_rect = corrcoef(rect(P(1:N_mitral,:)));
y_fisher(count,1) = temp1_rect(compare(1),compare(2));
y_fisher(count,2) = temp1_rect(compare(3),compare(4));
y_fisher(count,11)=max(max(P(1:N_mitral,1:4)));
y_fisher(count,12)=min(min(P(1:N_mitral,1:4)));
%y_fisher(count,3)=max(max(SigmoidM(P(1:N_mitral,1:4))));
%y_fisher(count,4)=min(min(SigmoidM(P(1:N_mitral,1:4))));
y_fisher(count,3)=max(max(rect(P(1:N_mitral,1:4))));
y_fisher(count,4)=min(min(rect(P(1:N_mitral,1:4))));

%  y3(count) = temp1(5,6);
%  y4(count) = temp1(7,8);
%                x(count) = count-1;
%[temp3,temp3_thresh,temp3_comp,temp3_comp_thresh] = fischer(P(1:N_mitral,1),P(1:N_mitral,2));
%[temp4,temp4_thresh,temp4_comp,temp4_comp_thresh] = fischer(P(1:N_mitral,3),P(1:N_mitral,4));
[temp3_rect,temp3_comp_rect] = fisher_rect(P(1:N_mitral,compare(1)),P(1:N_mitral,compare(2)),fisher_thresh);
[temp4_rect,temp4_comp_rect] = fisher_rect(P(1:N_mitral,compare(3)),P(1:N_mitral,compare(4)),fisher_thresh);

y_fisher(count,7) = temp3_rect;
y_fisher(count,8) = temp4_rect;

if (mod(count,fisher_plot)==0)
    
    time_axis=1:count;
    fig_label(fig_number)=figure(fig_number);
    set(fig_label(fig_number),'Position',pos{fig_number});
    h1=subplot_tight(5,1,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    plot(time_fisher(time_axis),y_fisher(time_axis,1),time_fisher(time_axis),y_fisher(time_axis,2),time_fisher(time_axis),y_fisher(time_axis,9),time_fisher(time_axis),y_fisher(time_axis,10));
    legend('P1 rect', 'P2 rect','P1 Sig', 'P2 Sig ')%, 'pair 3', 'pair 4')
    if (fig_number==18)
        title(' Correlation of Pairs with Cortical Input ')
    end
    if (fig_number==19)
        title(' Correlation of Pairs without Cortical Input ')
    end
    h1=subplot_tight(5,1,2,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    plot(time_fisher(time_axis),y_fisher(time_axis,7),time_fisher(time_axis),y_fisher(time_axis,8));
    title(' Fisher discriminant of pairs')
    
    h1=subplot_tight(5,1,3,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    %plot(1:N_mitral,temp3_comp,1:N_mitral,temp3_comp_thresh,1:N_mitral,temp4_comp,1:N_mitral,temp4_comp_thresh,1:N_mitral,temp3_comp_rect,1:N_mitral,temp4_comp_rect)
    plot(1:N_mitral,temp3_comp_rect,1:N_mitral,temp4_comp_rect)
    %ylim([-0.05,0.2]);
    xlim([1,N_mitral]);
    xlabel(' MC Number'); ylabel(' Contribution to Fisher')
    %legend('P1','P1 thresh','P2','P2 thresh','P1 rect','P2 rect')
    %legend('P1 rect','P2 rect')
    
    h1=subplot_tight(5,1,4,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    %plot(1:length(P(1:N_mitral,1)),P(1:N_mitral,1),1:length(P(1:N_mitral,2)),P(1:N_mitral,2),1:length(P(1:N_mitral,3)),P(1:N_mitral,3),1:length(P(1:N_mitral,4)),P(1:N_mitral,4))
    plot(1:N_mitral,P(1:N_mitral,1),1:N_mitral,P(1:N_mitral,2),1:N_mitral,P(1:N_mitral,3),1:N_mitral,P(1:N_mitral,4),[1,N_mitral],[0,0])
    xlabel(' MC Number'); ylabel(' MC Activity')
    xlim([1,N_mitral]);
    h1=subplot_tight(5,1,5,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    %plot(1:length(P(1:N_mitral,1)),P(1:N_mitral,1),1:length(P(1:N_mitral,2)),P(1:N_mitral,2),1:length(P(1:N_mitral,3)),P(1:N_mitral,3),1:length(P(1:N_mitral,4)),P(1:N_mitral,4))
    %plot(1:N_mitral,SigmoidM(P(1:N_mitral,1)),1:N_mitral,SigmoidM(P(1:N_mitral,2)),1:N_mitral,SigmoidM(P(1:N_mitral,3)),1:N_mitral,SigmoidM(P(1:N_mitral,4)),[1,N_mitral],[0,0])
    plot(1:N_mitral,rect(P(1:N_mitral,1)),1:N_mitral,rect(P(1:N_mitral,2)),1:N_mitral,rect(P(1:N_mitral,3)),1:N_mitral,rect(P(1:N_mitral,4)),[1,N_mitral],[0,0])
    xlabel(' MC Number'); ylabel(' rect(MC)')
    xlim([1,N_mitral]);
end


[fisher_mean,fisher_std,fisher_quantiles]=fisher_non_optimal(P(1:N_mitral,compare(1)),P(1:N_mitral,compare(2)),fisher_thresh);
y_fisher(count,13)=fisher_mean;
y_fisher(count,14)=fisher_std;
y_fisher(count,15:17)=fisher_quantiles;
[fisher_mean,fisher_std,fisher_quantiles]=fisher_non_optimal(P(1:N_mitral,compare(3)),P(1:N_mitral,compare(4)),fisher_thresh);
y_fisher(count,18)=fisher_mean;
y_fisher(count,19)=fisher_std;
y_fisher(count,20:22)=fisher_quantiles;

if (mod(count,fisher_plot)==0)
    
    %fig_label(fig_number+100)=
    figure(fig_number+1000); % set(fig_label(fig_number),'Position',pos{fig_number});
    h1=subplot_tight(2,1,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    plot(time_fisher(time_axis),y_fisher(time_axis,13),time_fisher(time_axis),y_fisher(time_axis,14),...
        time_fisher(time_axis),y_fisher(time_axis,15),time_fisher(time_axis),y_fisher(time_axis,16),time_fisher(time_axis),y_fisher(time_axis,17));
    legend('mean','std','25','50','75')
    if (fig_number==1018)
        title('Fisher for Random Read-out P1 with Cortical Input')
    end
    if (fig_number==1019)
        title('Fisher for Random Read-out P1 without Cortical Input')
    end
    h1=subplot_tight(2,1,2,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    plot(time_fisher(time_axis),y_fisher(time_axis,18),time_fisher(time_axis),y_fisher(time_axis,19),...
        time_fisher(time_axis),y_fisher(time_axis,20),time_fisher(time_axis),y_fisher(time_axis,21),time_fisher(time_axis),y_fisher(time_axis,22));
    legend('mean','std','25','50','75')
    title('Fisher for random read-out P2')
end

end    % fisher_diagnostic

%%%%%%%%%%%%%%%%%%%%%

function odor_light_diagnostic(N_mitral,num_light_cells,S_smell,S_spontaneous,light_on,light_background,CS,Wpg,Wmg,Wgm)
global pos

%Nc=64; Ns=6;
Nc=N_mitral+num_light_cells;
x=[1:N_mitral];
amp=1; width=12;
%want to change xpos
rowcount = 1;
odor_redu=NaN(N_mitral,18);
Var_signal_light_off=NaN(Nc,18);
Var_signal_light_on=NaN(Nc,18);
for xpos = 17:2:52
    odor_redu(:,rowcount)=(amp*exp(-(x-xpos).^2/width^2));
    rowcount = rowcount+1;
end
S_light_off=light_background*ones(num_light_cells,18);
S_light_on=light_on*ones(num_light_cells,18);
S_smell_2=S_smell;
S_smell = S_spontaneous+odor_redu;

SS_1 = [S_smell(:,[1:6]); S_light_off(:,[1:6])];
SS_2 = [S_smell(:,[7:12]); S_light_off(:,[7:12])];
SS_3 = [S_smell(:,[13:18]); S_light_off(:,[13:18])];
%currently Ws has reduced backgroud actvitiy
WS_1 = [(S_smell(:,[1:6])-S_spontaneous)*.5+S_spontaneous;S_light_off(:,[1:6])];
WS_2 = [(S_smell(:,[7:12])-S_spontaneous)*.5+S_spontaneous;S_light_off(:,[7:12])];
WS_3 = [(S_smell(:,[13:18])-S_spontaneous)*.5+S_spontaneous;S_light_off(:,[13:18])];

SSL_1 = [S_smell(:,[1:6]); S_light_on(:,[1:6])];
SSL_2 = [S_smell(:,[7:12]); S_light_on(:,[7:12])];
SSL_3 = [S_smell(:,[13:18]); S_light_on(:,[13:18])];

WSL_1 = [(S_smell(:,[1:6])-S_spontaneous)*.5+S_spontaneous;S_light_on(:,[1:6])];
WSL_2 = [(S_smell(:,[7:12])-S_spontaneous)*.5+S_spontaneous;S_light_on(:,[7:12])];
WSL_3 = [(S_smell(:,[13:18])-S_spontaneous)*.5+S_spontaneous;S_light_on(:,[13:18])];

for i=1:18
    Var_signal_light_off(:,i) = [(S_smell_2(:,3)-S_spontaneous)*(1-.05*(i-1))+S_spontaneous;S_light_off(:,1)];
end

for i=1:18
    Var_signal_light_on(:,i) = [(S_smell_2(:,3)-S_spontaneous)*(1-.05*(i-1))+S_spontaneous;S_light_on(:,1)];
end

log_first_run=1;
fig307=figure(307); set(fig307,'Position',pos{307});
h1=subplot(1,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
[P_1, ~,~,~, S_mem_1, p_recall_1] = cal_activity1(CS,Wpg,Wmg,Wgm,SS_1,SS_1,0);%Extra argument for S_mem
[P_2, ~,~,~, S_mem_2, p_recall_2] = cal_activity1(CS,Wpg,Wmg,Wgm,SS_2,SS_2,0);%Extra argument for S_mem
[P_3, ~,~,~, S_mem_3, p_recall_3] = cal_activity1(CS,Wpg,Wmg,Wgm,SS_3,SS_3,0);%Extra argument for S_mem
subplot(3,3,1); imagesc(P_1); title('mcell activity SS'); colorbar;
subplot(3,3,2); imagesc(p_recall_1); title('pcell activity'); colorbar;
subplot(3,3,3); imagesc(S_mem_1); title('S mem'); colorbar;

subplot(3,3,4); imagesc(P_2); title('mcell activity SS'); colorbar;
subplot(3,3,5); imagesc(p_recall_2); title('pcell activity'); colorbar;
subplot(3,3,6); imagesc(S_mem_2); title('S mem'); colorbar;

subplot(3,3,7); imagesc(P_3); title('mcell activity SS'); colorbar;
subplot(3,3,8); imagesc(p_recall_3); title('pcell activity'); colorbar;
subplot(3,3,9); imagesc(S_mem_3); title('S mem'); colorbar;

fig308=figure(308); set(fig308,'Position',pos{308});
h1=subplot(1,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
[P_1, ~,~,~, S_mem_1, p_recall_1] = cal_activity1(CS,Wpg,Wmg,Wgm,WS_1,WS_1,0);%Extra argument for S_mem
[P_2, ~,~,~, S_mem_2, p_recall_2] = cal_activity1(CS,Wpg,Wmg,Wgm,WS_2,WS_2,0);%Extra argument for S_mem
[P_3, ~,~,~, S_mem_3, p_recall_3] = cal_activity1(CS,Wpg,Wmg,Wgm,WS_3,WS_3,0);%Extra argument for S_mem
subplot(3,3,1); imagesc(P_1); title('mcell activity WS'); colorbar;
subplot(3,3,2); imagesc(p_recall_1); title('pcell activity'); colorbar;
subplot(3,3,3); imagesc(S_mem_1); title('S mem'); colorbar;

subplot(3,3,4); imagesc(P_2); title('mcell activity WS'); colorbar;
subplot(3,3,5); imagesc(p_recall_2); title('pcell activity'); colorbar;
subplot(3,3,6); imagesc(S_mem_2); title('S mem'); colorbar;

subplot(3,3,7); imagesc(P_3); title('mcell activity WS'); colorbar;
subplot(3,3,8); imagesc(p_recall_3); title('pcell activity'); colorbar;
subplot(3,3,9); imagesc(S_mem_3); title('S mem'); colorbar;

fig309=figure(309); set(fig309,'Position',pos{309});
h1=subplot(1,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
[P_1, ~,~,~, S_mem_1, p_recall_1] = cal_activity1(CS,Wpg,Wmg,Wgm,SSL_1,SSL_1,0);%Extra argument for S_mem
[P_2, ~,~,~, S_mem_2, p_recall_2] = cal_activity1(CS,Wpg,Wmg,Wgm,SSL_2,SSL_2,0);%Extra argument for S_mem
[P_3, ~,~,~, S_mem_3, p_recall_3] = cal_activity1(CS,Wpg,Wmg,Wgm,SSL_3,SSL_3,0);%Extra argument for S_mem
subplot(3,3,1); imagesc(P_1); title('mcell activity SSL'); colorbar;
subplot(3,3,2); imagesc(p_recall_1); title('pcell activity'); colorbar;
subplot(3,3,3); imagesc(S_mem_1); title('S mem'); colorbar;

subplot(3,3,4); imagesc(P_2); title('mcell activity SSL'); colorbar;
subplot(3,3,5); imagesc(p_recall_2); title('pcell activity'); colorbar;
subplot(3,3,6); imagesc(S_mem_2); title('S mem'); colorbar;

subplot(3,3,7); imagesc(P_3); title('mcell activity SSL'); colorbar;
subplot(3,3,8); imagesc(p_recall_3); title('pcell activity'); colorbar;
subplot(3,3,9); imagesc(S_mem_3); title('S mem'); colorbar;

fig310=figure(310); set(fig310,'Position',pos{310});
h1=subplot(1,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
[P_1, ~,~,~, S_mem_1, p_recall_1] = cal_activity1(CS,Wpg,Wmg,Wgm,WSL_1,WSL_1,0);%Extra argument for S_mem
[P_2, ~,~,~, S_mem_2, p_recall_2] = cal_activity1(CS,Wpg,Wmg,Wgm,WSL_2,WSL_2,0);%Extra argument for S_mem
[P_3, ~,~,~, S_mem_3, p_recall_3] = cal_activity1(CS,Wpg,Wmg,Wgm,WSL_3,WSL_3,0);%Extra argument for S_mem
subplot(3,3,1); imagesc(P_1); title('mcell activity WSL'); colorbar;
subplot(3,3,2); imagesc(p_recall_1); title('pcell activity'); colorbar;
subplot(3,3,3); imagesc(S_mem_1); title('S mem'); colorbar;

subplot(3,3,4); imagesc(P_2); title('mcell activity WSL'); colorbar;
subplot(3,3,5); imagesc(p_recall_2); title('pcell activity'); colorbar;
subplot(3,3,6); imagesc(S_mem_2); title('S mem'); colorbar;

subplot(3,3,7); imagesc(P_3); title('mcell activity WSL'); colorbar;
subplot(3,3,8); imagesc(p_recall_3); title('pcell activity'); colorbar;
subplot(3,3,9); imagesc(S_mem_3); title('S mem'); colorbar;

fig311=figure(311); set(fig311,'Position',pos{311});
h1=subplot(1,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
[P_1, ~,~,~, S_mem_1, p_recall_1] = cal_activity1(CS,Wpg,Wmg,Wgm,Var_signal_light_off(:,1:6),Var_signal_light_off(:,1:6),0);%Extra argument for S_mem
[P_2, ~,~,~, S_mem_2, p_recall_2] = cal_activity1(CS,Wpg,Wmg,Wgm,Var_signal_light_off(:,7:12),Var_signal_light_off(:,7:12),0);%Extra argument for S_mem
[P_3, ~,~,~, S_mem_3, p_recall_3] = cal_activity1(CS,Wpg,Wmg,Wgm,Var_signal_light_off(:,13:18),Var_signal_light_off(:,13:18),0);%Extra argument for S_mem
subplot(3,3,1); imagesc(P_1); title('mcell activity Var light off'); colorbar;
subplot(3,3,2); imagesc(p_recall_1); title('pcell activity'); colorbar;
subplot(3,3,3); imagesc(S_mem_1); title('S mem'); colorbar;

subplot(3,3,4); imagesc(P_2); title('mcell activity Var light off'); colorbar;
subplot(3,3,5); imagesc(p_recall_2); title('pcell activity'); colorbar;
subplot(3,3,6); imagesc(S_mem_2); title('S mem'); colorbar;

subplot(3,3,7); imagesc(P_3); title('mcell activity Var light off'); colorbar;
subplot(3,3,8); imagesc(p_recall_3); title('pcell activity'); colorbar;
subplot(3,3,9); imagesc(S_mem_3); title('S mem'); colorbar;

fig312=figure(312); set(fig312,'Position',pos{312});
h1=subplot(1,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
[P_1, ~,~,~, S_mem_1, p_recall_1] = cal_activity1(CS,Wpg,Wmg,Wgm,Var_signal_light_on(:,1:6),Var_signal_light_on(:,1:6),0);%Extra argument for S_mem
[P_2, ~,~,~, S_mem_2, p_recall_2] = cal_activity1(CS,Wpg,Wmg,Wgm,Var_signal_light_on(:,7:12),Var_signal_light_on(:,7:12),0);%Extra argument for S_mem
[P_3, ~,~,~, S_mem_3, p_recall_3] = cal_activity1(CS,Wpg,Wmg,Wgm,Var_signal_light_on(:,13:18),Var_signal_light_on(:,13:18),0);%Extra argument for S_mem
subplot(3,3,1); imagesc(P_1); title('mcell activity Var light on'); colorbar;
subplot(3,3,2); imagesc(p_recall_1); title('pcell activity'); colorbar;
subplot(3,3,3); imagesc(S_mem_1); title('S mem'); colorbar;

subplot(3,3,4); imagesc(P_2); title('mcell activity Var light on'); colorbar;
subplot(3,3,5); imagesc(p_recall_2); title('pcell activity'); colorbar;
subplot(3,3,6); imagesc(S_mem_2); title('S mem'); colorbar;

subplot(3,3,7); imagesc(P_3); title('mcell activity Var light on'); colorbar;
subplot(3,3,8); imagesc(p_recall_3); title('pcell activity'); colorbar;
subplot(3,3,9); imagesc(S_mem_3); title('S mem'); colorbar;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function positions
global screen pos

set(0,'Units','pixels');
screen=get(0,'ScreenSize');
%screen=[1 1 1800 900];
pos{1}=[.9*screen(3),.85*screen(4),.1*screen(3),.1*screen(4)];
pos{2}=[.68*screen(3),.8*screen(4),.2*screen(3),.2*screen(4)];
pos{3}=[.55*screen(3),.5*screen(4),.45*screen(3),.5*screen(4)];
pos{4}=[.9*screen(3),.75*screen(4),.1*screen(3),.1*screen(4)];
pos{5}=[.0*screen(3),.83*screen(4),.49*screen(3),.17*screen(4)];
pos{6}=[.68*screen(3),.6*screen(4),.25*screen(3),.2*screen(4)];
pos{15}=[.22*screen(3),.05*screen(4),.25*screen(3),.3*screen(4)];
pos{16}=[.22*screen(3),.4*screen(4),.45*screen(3),.4*screen(4)];
pos{17}=[.0*screen(3),.53*screen(4),.22*screen(3),.22*screen(4)];
pos{18}=[.0*screen(3),.02*screen(4),.22*screen(3),.53*screen(4)];
pos{19}=[.35*screen(3),.0*screen(4),.22*screen(3),.5*screen(4)];
pos{20}=[.64*screen(3),.0*screen(4),.22*screen(3),.5*screen(4)];
pos{21}=[.0*screen(3),.27*screen(4),.22*screen(3),.22*screen(4)];
pos{22}=[.0*screen(3),.0*screen(4),.22*screen(3),.22*screen(4)];

pos{24}=[.35*screen(3),.5*screen(4),.15*screen(3),.2*screen(4)];
pos{30}=[.0*screen(3),.0*screen(4),.3*screen(3),.5*screen(4)];
pos{31}=[.52*screen(3),.0*screen(4),.15*screen(3),.33*screen(4)];
pos{32}=[.7*screen(3),.0*screen(4),.15*screen(3),.33*screen(4)];
pos{40}=[.87*screen(3),.25*screen(4),.15*screen(3),.25*screen(4)];

pos{48}=[.53*screen(3),.0*screen(4),.25*screen(3),.17*screen(4)];
pos{49}=[.73*screen(3),.0*screen(4),.25*screen(3),.17*screen(4)];
pos{50}=[.73*screen(3),.6*screen(4),.25*screen(3),.17*screen(4)];
pos{51}=[.73*screen(3),.4*screen(4),.25*screen(3),.17*screen(4)];
pos{52}=[.73*screen(3),.2*screen(4),.25*screen(3),.17*screen(4)];

pos{80}=[.5*screen(3),.2*screen(4),.25*screen(3),.3*screen(4)];
pos{81}=[.75*screen(3),.2*screen(4),.25*screen(3),.3*screen(4)];
pos{90}=[.5*screen(3),.4*screen(4),.25*screen(3),.3*screen(4)];
pos{91}=[.75*screen(3),.4*screen(4),.25*screen(3),.3*screen(4)];

pos{100}=[.75*screen(3),.7*screen(4),.15*screen(3),.25*screen(4)];

pos{307}=[.01*screen(3),.27*screen(4),.2*screen(3),.2*screen(4)];
pos{308}=[.22*screen(3),.27*screen(4),.2*screen(3),.2*screen(4)];
pos{311}=[.43*screen(3),.27*screen(4),.2*screen(3),.2*screen(4)];
pos{309}=[.01*screen(3),.01*screen(4),.2*screen(3),.2*screen(4)];
pos{310}=[.22*screen(3),.01*screen(4),.2*screen(3),.2*screen(4)];
pos{312}=[.43*screen(3),.01*screen(4),.2*screen(3),.2*screen(4)];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [GC_hist,GC_hist_cum]=plot_GC_hist(count,I,GC_bins,GC_hist,GC_hist_cum)

global pos fig_label

GC_active=hist(I,GC_bins);
GC_active_cum=flipud(cumsum(flipud(GC_active)));
GC_hist(:,:,count)=GC_active;
GC_hist_cum(:,:,count)=GC_active_cum/GC_active_cum(1);

fig_number=32;
fig_label(fig_number)=figure(fig_number); set(fig_label(fig_number),'Position',pos{fig_number});
for i_odor=1:size(I,2)
    %    h1=subplot_tight(2,size(I,2),2*i_odor-1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
%   h1=subplot_tight(2,size(I,2),i_odor,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    h1=subplot(2,size(I,2),i_odor); ax1 = get(h1,'position'); set(h1,'position', ax1);
    %gc_tmp(:,:)=GC_hist(1:end,i_odor,1:count);
    gc_tmp(:,:)=GC_hist(2:end,i_odor,1:count);
    imagesc(gc_tmp')
    title('Hist')
    %colorbar
    %    h1=subplot_tight(2,size(I,2),2*i_odor,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
%    h1=subplot_tight(2,size(I,2),size(I,2)+i_odor,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
    h1=subplot(2,size(I,2),size(I,2)+i_odor); ax1 = get(h1,'position'); set(h1,'position', ax1);
    %gc_tmp(:,:)=GC_hist_cum(1:end,i_odor,1:count);
    gc_tmp(:,:)=GC_hist_cum(2:end,i_odor,1:count);
    imagesc(gc_tmp')
    caxis([0 1]);
    title('Cum Hist')
    %colorbar
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_S_M_P(fig_number,S_mem,P,P_recall)
global N_mitral pos
global n_P fig_label
fig_label(fig_number)=figure(fig_number);
%set(fig_label(fig_number),'Position',pos{fig_number});
%h1=subplot_tight(2,2,1,[.05,.05]); ax1 = get(h1,'position'); set(h1,'position', ax1);
h1=subplot(2,2,1); ax1 = get(h1,'position'); set(h1,'position', ax1);
imagesc(P_recall);
if (ismember(fig_number,[48,49,80,81]))
    title('Sigmoid of CC Activity');
end
if (ismember(fig_number,[17,21,22]))
    title('CC Activity')
end
xlabel('Stimulus Index');
ylabel('CC Index');
axis square;
colorbar;

%h2 = subplot_tight(2,2,2,[.05,.05]); ax2 = get(h2,'position'); set(h2,'position', ax2);
h2 = subplot(2,2,2); ax2 = get(h2,'position'); set(h2,'position', ax2);
if (fig_number<50)
    imagesc(S_mem(1:n_P,:)); %XH
    if (ismember(fig_number,[48,49,80,81,90,91]))
        title('Stimuli')
    end
    if (ismember(fig_number,[17,21,22]))
        title('Sigmoid of CC Activity')
        ylabel('CC Index');
    end
    axis square; colorbar;
else
    % this is the case probing with the mixtures
    plot(S_mem(1:N_mitral,:))
    xlim([1,N_mitral])
    title('Probe Stimuli')
    ylabel('Amplitude');
end
xlabel('Cell Index');

%h2 = subplot_tight(2,2,3,[.05,.05]); ax2 = get(h2,'position'); set(h2,'position', ax2);
h2 = subplot(2,2,3); ax2 = get(h2,'position'); set(h2,'position', ax2);
imagesc(P(1:N_mitral,:));
title('MC activity')
axis square; colorbar;
xlabel('Stimulus Index');
ylabel('MC Index');

%h2 = subplot_tight(2,2,4,[.05,.05]); ax2 = get(h2,'position'); set(h2,'position', ax2);
h2 = subplot(2,2,4); ax2 = get(h2,'position'); set(h2,'position', ax2);
plot(P(1:N_mitral,:))
hold all
plot([1,N_mitral],[0,0]);
xlim([1,N_mitral])
if (ismember(fig_number,[80,81,90,91]))
    title('MC for Probe Stimuli')
end
if (ismember(fig_number,[17]))
    title('MC for Training Stimuli')
end
xlabel('MC Index');
ylabel('MC Activity');
hold off

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [graded_stimuli]=graded(S,probe)

nsteps=5;
amplitudes=[1/nsteps:1/nsteps:1];

%graded_stimuli=[kron(S(:,stim_select(1)),amplitudes),kron(probe(:,stim_select(2)),amplitudes)];
graded_stimuli=[kron(S(:,:),amplitudes),kron(probe(:,:),amplitudes)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sigmoid
global n_P N_mitral
global maxv minv ts gamma

x=linspace(-2,20,200);
cells=ones(n_P,1);
for j=1:length(x)
    xi=x(j);
    dummy=Sigmoid(xi*cells);
    y1_odor(j)=dummy(1);
    y1_light(j)=dummy(N_mitral+1);
end
hold off
%y2=SigmoidM(x);
y2=rect(x);
y3=SigmoidG(x);
prob_=(tanh((x-ts)*gamma)+1)*(maxv-minv)/2+minv;%y3=x;

figure(1000)
subplot(4,1,1)
plot(x,y1_odor,'-o',x,y1_light,'-x')
xlim([-2,10])
legend('P_{odor}','P_{light}');
title(' Sigmoid P')
subplot(4,1,2)
plot(x,y2)
xlim([-2,10])
legend('M');
title(' Sigmoid M')
subplot(4,1,3)
plot(x,y3)
xlim([-2,10])
legend('G');
title(' Sigmoid G')
subplot(4,1,4)
plot(x,prob_)
xlim([0,20])
legend('Resilience');
title(' Survival Probability')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inh_incr_mc,inh_incr_pc]=inhibitory_contributions(MC,GC,PC)

global Wgm Wpg

% the matrix mci_to_gcj gives the contribution of each mitral cell to the activity of each granule cell
mci_to_gcj = bsxfun(@times,Wgm',SigmoidM(MC))';
% the matrix inh_remove_mc gives for each granule-mitral pair the inhibitory input (SigmoidG) that the granule cell
% would generate if the input from that mitral cell was removed
gcj_remove_mci=bsxfun(@minus,GC,mci_to_gcj);
% the matrix inh_incr_mc gives for each granule-mitral pair the increment in the inhibitory input the granule cell
% gives due to the activity of that mitral cell
inh_incr_mc=bsxfun(@minus,SigmoidG(GC),SigmoidG(gcj_remove_mci));

% the matrix pci_to_gcj gives the contribution of each mitral cell to the activity of each granule cell
pci_to_gcj = bsxfun(@times,Wpg,Sigmoid(PC))';
% the matrix inh_remove_mc gives for each granule-mitral pair the inhibitory input (SigmoidG) that the granule cell
% would generate if the input from that mitral cell was removed
gcj_remove_pci=bsxfun(@minus,GC,pci_to_gcj);
% the matrix inh_incr_mc gives for each granule-mitral pair the increment in the inhibitory input the granule cell
% gives due to the activity of that mitral cell
inh_incr_pc=bsxfun(@minus,SigmoidG(GC),SigmoidG(gcj_remove_pci));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [formatting_string]=make_formatting_strings(n_max)
formatting_string=cell(n_max,1);
for n=1:n_max
    formatting_string{n}=' ';
    for i=1:n
        formatting_string{n}=[formatting_string{n},' %g'];
    end
    formatting_string{n}=[formatting_string{n},' \n'];
end
end

%%%%%%%%%%%%%%%%%%%%%%%

function probe_system(time,probe_set,xpositions_probe,widths_probe,istep,label_connectivity)

global wt1_initial wt1_later wt1_later2 CS
global Wpg Wmg Wgm N_mitral
global file_location file_tag file_type_data file_type_matlab
global format_string
global fisher_thresh
global i_state state_probe learning_time probe_case lxp1 lxp2 lxp3 lxp4
global log_probe_random
global log_first_run

i_state=i_state+1;
n_stim=size(xpositions_probe{1},1);
n_pair=n_stim/2;
if (floor(n_pair)~=n_pair)
    n_pair=floor(n_pair);
    fprintf(' probe has an odd number of entries \n round to the lower integer n_pair = %g \n ',n_pair);
end
% P should be removed from the call of cal_activity1
P=1;
log_first_run=1;
[P_tmp,I_tmp,I1_tmp,I2_tmp,S_mem,P_recall_tmp] = cal_activity1(CS,Wpg,Wmg,Wgm,probe_set{1},P,time);%Extra argument for S_mem and Wpg
if (label_connectivity=='a')
    plot_S_M_P(80,probe_set{1},P_tmp,S_mem)
end
if (label_connectivity=='b')
    plot_S_M_P(90,probe_set{1},P_tmp,S_mem)
end

% connectivity
Wmp=Wmg*Wpg';    Wmm=Wmg*Wmg';

x1min=floor((xpositions_probe{1}(1)+xpositions_probe{1}(2))/2-widths_probe{1}(1)/2);
x1max=ceil((xpositions_probe{1}(1)+xpositions_probe{1}(2))/2+widths_probe{1}(1)/2);
x2min=floor((xpositions_probe{1}(3)+xpositions_probe{1}(4))/2-widths_probe{1}(2)/2);
x2max=ceil((xpositions_probe{1}(3)+xpositions_probe{1}(4))/2+widths_probe{1}(2)/2);
wmp_from_memory=sum(Wmp(:,x1min:x1max),2);
wmp_memory_correct=sum(wmp_from_memory(x1min:x1max));
wmp_memory_intermediate=sum(wmp_from_memory(x1max:x2min));
wmp_memory_incorrect=sum(wmp_from_memory(x2min:x2max));
wmp_ratio_incorrect=wmp_memory_incorrect/wmp_memory_correct;
wmp_ratio_other=wmp_memory_intermediate/wmp_memory_correct;
wmm_from_memory=sum(Wmm(:,x1min:x1max),2);
wmm_memory_correct=sum(wmm_from_memory(x1min:x1max));
wmm_memory_intermediate=sum(wmm_from_memory(x1max:x2min));
wmm_memory_incorrect=sum(wmm_from_memory(x2min:x2max));
wmm_ratio_incorrect=wmm_memory_incorrect/wmm_memory_correct;
wmm_ratio_other=wmm_memory_intermediate/wmm_memory_correct;

figname='fisher_probe_connectivity'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
fid_fisher_probe_MC=fopen(file_name,'a');
fprintf(fid_fisher_probe_MC,' %g %g %g %g %g %g %g %g %g %g %g \n',time,...
    wmp_ratio_incorrect,wmp_ratio_other,wmp_memory_correct,wmp_memory_intermediate,wmp_memory_incorrect,...
    wmm_ratio_incorrect,wmm_ratio_other,wmm_memory_correct,wmm_memory_intermediate,wmm_memory_incorrect);
fclose(fid_fisher_probe_MC);

state_probe{i_state,1}=time;
state_probe{i_state,2}=P_tmp;
state_probe{i_state,3}=I_tmp;
state_probe{i_state,4}=I1_tmp;
state_probe{i_state,5}=I2_tmp;
state_probe{i_state,6}=P_recall_tmp;

%for i_stim=1:n_stim
figname='fisher_probe_cortex_MC'; file_name = strcat(file_location,figname,'_',file_tag,'_i',num2str(istep),file_type_data);
fid_fisher_probe_MC=fopen(file_name,'a');
fprintf(fid_fisher_probe_MC,' & t = %g \n',time);
fprintf(fid_fisher_probe_MC,format_string{1+n_stim},[[1:N_mitral]',P_tmp(1:N_mitral,1:n_stim)]');
fclose(fid_fisher_probe_MC);
%end

figname='fisher_contrib_probe_cortex'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
fid_fisher_probe_MC=fopen(file_name,'a');
fisher_probe_input=NaN(n_pair,1); fisher_probe_output_cortex=NaN(n_pair,1); fisher_probe_output_nocortex=NaN(n_pair,1);
for i_pair=1:n_pair
    fisher_input_1=probe_set{1}(1:N_mitral,2*i_pair-1);
    fisher_input_2=probe_set{1}(1:N_mitral,2*i_pair);
    fisher_probe_input(i_pair)=fisher_rect(fisher_input_1,fisher_input_2,fisher_thresh);
    fisher_output_1=P_tmp(1:N_mitral,2*i_pair-1);
    fisher_output_2=P_tmp(1:N_mitral,2*i_pair);
    [fisher_probe_output_cortex(i_pair),fisher_contrib]=fisher_rect(fisher_output_1,fisher_output_2,fisher_thresh);
    fprintf(' Fisher: input = %g  output = %g \n',fisher_probe_input(i_pair),fisher_probe_output_cortex(i_pair));
    fprintf(fid_fisher_probe_MC,' & time = %g \n',time);
    fprintf(fid_fisher_probe_MC,' %g %g \n ',[1:N_mitral;fisher_contrib']);
end
fclose(fid_fisher_probe_MC);

if (log_probe_random==1)
    figname='fisher_non_opt_probe_cortex'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
    fid_fisher_probe=fopen(file_name,'a');
    fprintf(fid_fisher_probe,' %g ',time);
    for i_pair=1:n_pair
        fisher_output_1=P_tmp(1:N_mitral,2*i_pair-1);
        fisher_output_2=P_tmp(1:N_mitral,2*i_pair);
        [fisher_mean,fisher_std,fisher_quantiles]=fisher_non_optimal(fisher_output_1,fisher_output_2,fisher_thresh);
        fprintf(' Fisher random: mean = %g  std = %g quantiles = %g %g %g \n',fisher_mean,fisher_std,fisher_quantiles);
        fprintf(fid_fisher_probe,' %g %g %g %g %g ',fisher_mean,fisher_std,fisher_quantiles);
    end
    fprintf(fid_fisher_probe,' \n');
    fclose(fid_fisher_probe);
end

fprintf(' now with modified cortical input wt1 = %g \n',wt1_later2);
wt1_initial_tmp=wt1_initial; wt1_later_tmp=wt1_later;
wt1_initial=0; wt1_later=wt1_later2;
log_first_run=1;
[P_tmp,I_tmp,I1_tmp,I2_tmp,S_mem,P_recall_tmp] = cal_activity1(CS,Wpg,Wmg,Wgm,probe_set{1},P,time);%Extra argument for S_mem and Wpg
wt1_initial=wt1_initial_tmp; wt1_later=wt1_later_tmp;

if (label_connectivity=='a')
    plot_S_M_P(81,probe_set{1},P_tmp,S_mem)
end
if (label_connectivity=='b')
    plot_S_M_P(91,probe_set{1},P_tmp,S_mem)
end

state_probe{i_state,7}=P_tmp;
state_probe{i_state,8}=I_tmp;
state_probe{i_state,9}=I1_tmp;
state_probe{i_state,10}=I2_tmp;
state_probe{i_state,11}=P_recall_tmp;
state_probe{i_state,12}=learning_time;

figname='states_probe';
file_name = strcat(file_location,figname,'_',file_tag,file_type_matlab);
save(file_name,'state_probe');

figname='fisher_probe_cortex_wt2_MC'; file_name = strcat(file_location,figname,'_',file_tag,'_i',num2str(istep),file_type_data);
fid_fisher_probe_MC=fopen(file_name,'a');
fprintf(fid_fisher_probe_MC,' & t = %g \n',time);
fprintf(fid_fisher_probe_MC,format_string{1+n_stim},[[1:N_mitral]',P_tmp(1:N_mitral,1:n_stim)]');
fclose(fid_fisher_probe_MC);

figname='fisher_contrib_probe_cortex_wt2'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
fid_fisher_probe_MC=fopen(file_name,'a');
for i_pair=1:n_pair
    fisher_output_1=P_tmp(1:N_mitral,2*i_pair-1);
    fisher_output_2=P_tmp(1:N_mitral,2*i_pair);
    [fisher_probe_output_nocortex(i_pair),fisher_contrib]=fisher_rect(fisher_output_1,fisher_output_2,fisher_thresh);
    fprintf(' Fisher: input = %g  output = %g \n',fisher_probe_input(i_pair),fisher_probe_output_nocortex(i_pair));
    fprintf(fid_fisher_probe_MC,' & time = %g \n',time);
    fprintf(fid_fisher_probe_MC,' %g %g \n ',[1:N_mitral;fisher_contrib']);
end
fclose(fid_fisher_probe_MC);

if (log_probe_random==1)
    figname='fisher_non_opt_probe_cortex_wt2'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
    fid_fisher_probe=fopen(file_name,'a');
    fprintf(fid_fisher_probe,' %g ',time);
    for i_pair=1:n_pair
        fisher_output_1=P_tmp(1:N_mitral,2*i_pair-1);
        fisher_output_2=P_tmp(1:N_mitral,2*i_pair);
        [fisher_mean,fisher_std,fisher_quantiles]=fisher_non_optimal(fisher_output_1,fisher_output_2,fisher_thresh);
        fprintf(' Fisher random: mean = %g  std = %g quantiles = %g %g %g \n',fisher_mean,fisher_std,fisher_quantiles);
        fprintf(fid_fisher_probe,' %g %g %g %g %g ',fisher_mean,fisher_std,fisher_quantiles);
    end
    fprintf(fid_fisher_probe,' \n');
    fclose(fid_fisher_probe);
end

figname='fisher_probe_time'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
fid_fisher_probe=fopen(file_name,'a');
fprintf(fid_fisher_probe,format_string{1+4*n_pair},time,...
    fisher_probe_output_cortex,fisher_probe_output_nocortex,fisher_probe_input);
fprintf(fid_fisher_probe,'\n');
fclose(fid_fisher_probe);

if (probe_case==3)
    figname='fisher_probe_bars'; file_name = strcat(file_location,figname,'_',file_tag,file_type_data);
    fid_fisher_probe=fopen(file_name,'a');
    bar_positions_0=[1:lxp1/2,lxp1/2+2:lxp1/2+lxp2/2+1,lxp1/2+lxp2/2+3:lxp1/2+lxp2/2++lxp3/2+2,lxp1/2+lxp2/2+lxp3/2+4:lxp1/2+lxp2/2+lxp3/2+lxp4/2+3];
    bar_positions=sort([bar_positions_0,bar_positions_0+0.2,bar_positions_0+0.4,bar_positions_0+0.6]);
    fprintf(fid_fisher_probe,' & time = %g \n',time);
    fisher_probe_print=reshape([fisher_probe_input';fisher_probe_output_nocortex';fisher_probe_output_cortex';0*fisher_probe_input'],4*length(fisher_probe_input),1);
    fprintf(fid_fisher_probe,' %g %g \n',[bar_positions;fisher_probe_print']);
    fclose(fid_fisher_probe);
end

% since there is again a switch to a different stimulus set
log_first_run=1;

end

%%%%%%%%%%%%%%%%%%%
function h=subplot_tight(m,n,p,margins,varargin)
%function subplot_tight(m,n,p,margins,varargin)
%
% Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the margins between
% neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the margins between
% subplots can reach 40% of figure area, which is pretty lavish.
%
% Input arguments (defaults exist):
%   margins- two elements vector [vertical,horizontal] defining the margins between neighbouring axes. Default value
%            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%            relatively large axis.
%
% Output arguments: same as subplot- none, or axes handle according to function call.
%
% Issues & Comments: Note that if additional elements are used in order to be passed to subplot, margins parameter must
%       be defined. For default margins value use empty element- [].
%
% Author and Date:  Nikolay S. 29/03/2011.
% Last update:      Nikolay S. 21/04/2011 (accourding to Alan B comment).
%
% Usage example: h=subplot_tight((2,3,1:2,[0.5,0.2])

if (nargin<4) || isempty(margins)
    margins=[0.05,0.05]; % default margins value- 1% of figure
end

if length(margins)==1
    margins(2)=margins;
end

%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);


height=(1-(m+1)*margins(1))/m; % single subplot height
width=(1-(n+1)*margins(2))/n;  % single subplot width

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot

merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width

merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
pos_vec=[merged_left merged_bottom merged_width merged_height];

% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h_subplot=subplot('Position',pos_vec,varargin{:});

if nargout~=0
    h=h_subplot;
end

end

%%

function rand_out=random_factor
global parameter_variation
rand_out=1+2*parameter_variation*(rand-0.5);
end
