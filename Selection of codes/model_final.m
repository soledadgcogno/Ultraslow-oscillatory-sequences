%Extended data Figure 14

%% Ramp

clear all
close all
number_seq=1;
factor_v=[0.01,0.1,0.5,2,4];
fk=4;
factor_t=factor_v(fk);
g = 1.5; % Used to initialize the weight matrix
N = 500; % Number of cells

%Time vectors
dtData = 0.0641*factor_t;   % dt at which the data is sampled. Using the same dt that I use for calcium data
factor_dt_model=dtData/(0.001*factor_t);
dt = 0.001*factor_t*factor_dt_model; % integration step for the model
countT=0;
sigma=0.5*117/7.73;
T=900/number_seq;
epochs = number_seq*T/dtData; %length of the trial
tData = 0:dtData:epochs*dtData; % time steps at data sampling freq
t = 0:dt:tData(end); % time steps at simulation sampling freq


%Neural activity used as target- Sequences that are not periodic
%     nseq=1/number_seq;
%     tData_min=0:dtData:epochs*nseq*dtData; % time at data sampling freq
tData_min=0:dtData:T; % time at data sampling freq
beg_aux=-10:dtData:0-dtData;
end_aux=T+dtData:dtData:T+10;
tData_min_aux=[beg_aux,tData_min,end_aux];
com=T/N:T/N:T; %mean of the bump of each neuron
sig = sigma(1); % scaled correctly in neuron space
% sig2=exp(-(com((N/2))-tData_min(1:end)).^2/(2*sig^2));
%     sig3=repmat(sig2,1,number_seq);
% xBump = zeros(N, number_seq*length(tData_min)); % number of cells x number of data points
L=length(tData_min);
for i=1:N
    xBump_aux(i, :) = exp(-(com((i))-tData_min_aux(1:end)).^2/(2*sig^2));    
end

learning_it_max = 1000; 
T = 7177;
N = 500; 
alpha = 1;
input = zeros(N,T); 
dt_new=120/(T/2);
c=0;

for factor=[1,2,4,6,8]

    c=c+1;
    input = zeros(N,T);
    for t = 1:(T-1)/factor
        input(:,t)= xBump_aux(:,factor*t);
    end
    input2=repmat(input(:,1:(T-1)/factor),1,factor);
    input=input2(:,1:(T-1)/2);
    %input = xBump_aux;
%     output_target = zeros(1, T );
    output_target = zeros(1, size(input,2) );
    for t=1:(T-1)/2; output_target(t) = t;end
    if c==1
        norm=max(output_target((T-1)/2 - floor(20/dt_new)));
        figure
        plot(((1:(T-1)/2))*dt_new,output_target(1:(T-1)/2)./norm,'LineWidth',8);
        xlabel('Time (sec)','FontSize',16); axis square, axis([1 ceil(dt_new*(T-1)/2-20) 0 1])
        ylabel('Normalized activity');
        xticks([1 floor(T*dt_new/4 -10) ceil((T-1)*dt_new/2 -20)]);
        xticklabels({'0', num2str(floor(T*dt_new/4 -10)), num2str(ceil((T-1)*dt_new/2 -20))});
        set(gca,'fontsize',24)
    end

    output_target=output_target(1:(T-1)/2)./norm;
    %figure; plot(output_target)
    w = ones(N, 1);
    error_tot = zeros(1,learning_it_max);

%     figure
%     subplot(1,2,1)
%     plot(input(100,:));
%     subplot(1,2,2)
%     plot(output_target);
%     title(factor);

    for learning_it = 1:learning_it_max
        output = w'*input;
        error = output_target - output;
        error_tot(learning_it) = sum(abs(error(200:3200)));
        w = w +alpha*mean((ones(N,1)*error).*input,2);
    end
   
    figure; 
    imagesc(flip(input(:,1:(T-1)/2))); ylabel('Neuron #'); xlabel('Time (sec)'); 
    yticks([100 400]);
    yticklabels({'400','100'});
    xticks([1 floor(T/4) floor(T/2)]);
    xticklabels({'0', num2str(ceil(dt_new*floor(T/4))), num2str(ceil(dt_new*floor(T/2)))});
    set(gca,'fontsize',20); axis square;
    colorbar
    

    figure;
    plot(((1:(T-1)/2))*dt_new,output_target,'linewidth',10); axis([1 ceil(dt_new*(T-1)/2-20) 0 1]), hold on; 
    plot(((1:(T-1)/2))*dt_new,output,'linewidth',6); if c==1; (legend('Target','Output')); end; set(gca,'fontsize',20);...
        xlabel('Time (sec)'); axis square, axis([1 ceil(dt_new*(T-1)/2-20) 0 1]),  xticks([1 floor(T*dt_new/4)-10 ceil((T-1)*dt_new/2)-20]); ...
    xticklabels({'0', num2str(floor(T*dt_new/4)-10), num2str(ceil((T-1)*dt_new/2)-20)}), yticks([0 0.5 1]);
    ylabel('Normalized activity');  set(gca,'fontsize',20);

    error_final(c)=mean(error_tot(end-100:end));
    error_final(c)=mean(error_tot(end-100:end));

    clear input2
end

factor=[1,2,4,6,8];
figure; 
plot((T./factor)*dt_new,error_final,'k-*','linewidth',4,'MarkerSize',10);
ylabel('Mean total error'); yticks([0 400 800]);
xlabel('Cycle length (s)');
set(gca,'fontsize',20);
box off

%% Ornstein-Ulenbeck process
clear all
close all
number_seq=1;
first = 1;
learn = 1;
run = 1;
factor_v=[0.01,0.1,0.5,2,4];
fk=4;
factor_t=factor_v(fk);
g = 1.5; % Used to initialize the weight matrix
nRunTot = 40; %nRunTot has to be larger than nFree
nFree = 0;
P0 = 0.05; %Used for initializing the connectivity matrix the will undergo training
ampIn = 1;
N = 500; % Number of cells
nLearn = N; % Number of cells that learn
learnList = 1:N; %randperm(N); % random permutation of N cells
cL = learnList(1:nLearn); % Cells whose weights will undergo learning
nCL = learnList(nLearn:end); % Cells whose weights will NOT undergo learning

%Time vectors
dtData = 0.0641*factor_t;   % dt at which the data is sampled. Using the same dt that I use for calcium data
factor_dt_model=dtData/(0.001*factor_t);

dt = 0.001*factor_t*factor_dt_model; % integration step for the model
tau = 0.01*factor_t*factor_dt_model; % 10ms time constant
tauWN = 1*factor_t*1*factor_dt_model; % decay time of the stimulus


countT=0;
sigma=0.5*117/7.73;
T=900/number_seq;
epochs = number_seq*T/dtData; %length of the trial
tData = 0:dtData:epochs*dtData; % time steps at data sampling freq
t = 0:dt:tData(end); % time steps at simulation sampling freq

countT=countT+1; %counters that I am not using anymore
countS=0; %counters that I am not using anymore

%Neural activity used as target- Sequences that are not periodic
%     nseq=1/number_seq;
%     tData_min=0:dtData:epochs*nseq*dtData; % time at data sampling freq
tData_min=0:dtData:T; % time at data sampling freq
beg_aux=-10:dtData:0-dtData;
end_aux=T+dtData:dtData:T+10;
tData_min_aux=[beg_aux,tData_min,end_aux];

com=T/N:T/N:T; %mean of the bump of each neuron
sig = sigma(1); % scaled correctly in neuron space

sig2=exp(-(com((N/2))-tData_min(1:end)).^2/(2*sig^2));
%     sig3=repmat(sig2,1,number_seq);

xBump = zeros(N, number_seq*length(tData_min)); % number of cells x number of data points
L=length(tData_min);
for i=1:N
    xBump_aux(i, :) = exp(-(com((i))-tData_min_aux(1:end)).^2/(2*sig^2));
    
end

learning_it_max = 1000; 
T = 7177;
N = 500; 
alpha = 1;
input = zeros(N,T); 
dt_new=120/(T/2);

%Output target - OU process
tt=tData_min_aux;
u_ou    = zeros(length(tt),1);
b_ou  = 1/200;
mu    = 1.0; 
sig   = 0.005; 
for i = 1:length(tt)-1
    u_ou(i+1) = u_ou(i)+b_ou*(mu-u_ou(i))*dt + sig*sqrt(dt)*randn();
end
hold on
plot(u_ou)
M=min(u_ou(1:(T-1)/2));
if M<0
    u_ou = u_ou + abs(M);
end

% figure; plot(u_ou);
output_target = zeros(1, size(input,2) );
for t=1:(T-1)/2; output_target(t) =  u_ou(t);end

figure
plot(((1:(T-1)/2))*dt_new,output_target(1:(T-1)/2)./max(output_target(1:(T-1)/2)),'LineWidth',10);
xlabel('Time (sec)','FontSize',16); axis square, axis([1 ceil(dt_new*(T-1)/2-20) 0 inf])
ylabel('Normalized activity');
xticks([1 floor(T*dt_new/4 -10) ceil((T-1)*dt_new/2 -20)]);
xticklabels({'0', num2str(floor(T*dt_new/4 -10)), num2str(ceil((T-1)*dt_new/2 -20))});
set(gca,'fontsize',24)

norm=max(output_target(1:(T-1)/2));
c=0;
for factor=[1,2,4,6,8]
    c=c+1;
    input = zeros(N,T);
    for t = 1:(T-1)/factor
        input(:,t)= xBump_aux(:,factor*t);
    end

    input2=repmat(input(:,1:(T-1)/factor),1,factor);

    input=input2(:,1:(T-1)/2);

    output_target = zeros(1, size(input,2) );
    for t=1:(T-1)/2; output_target(t) =  u_ou(t)./norm;end

    %figure; plot(output_target)
    w = ones(N, 1);
    error_tot = zeros(1,learning_it_max);

    for learning_it = 1:learning_it_max
        output = w'*input;
        error = output_target - output;
        error_tot(learning_it) = sum(abs(error(200:3200)));
        w = w +alpha*mean((ones(N,1)*error).*input,2);
    end

%     figure; plot(error_tot)

    
%     figure; subplot(1,2,1); imagesc(input(:,1:(T-1)/2)); ylabel('Neuron #'); xlabel('Time'); set(gca,'fontsize',20); axis square;
%     subplot(1,2,2); plot(output_target((1:(T-1)/2)),'linewidth',2); hold on; plot(output((1:(T-1)/2)),'linewidth',2); if c==1; (legend('Target','Output')); end; set(gca,'fontsize',20);...
%         xlabel('Time'); axis square, axis([1 (T-1)/2 0 inf]);

    figure; 
    imagesc(flip(input(:,1:(T-1)/2))); ylabel('Neuron #'); xlabel('Time (sec)'); 
    yticks([100 400]);
    yticklabels({'400','100'});
    xticks([1 floor(T/4) floor(T/2)]);
    xticklabels({'0', num2str(ceil(dt_new*floor(T/4))), num2str(ceil(dt_new*floor(T/2)))});
    set(gca,'fontsize',20); axis square;

    figure;
    plot(((1:(T-1)/2))*dt_new,output_target((1:(T-1)/2)),'linewidth',10); axis([1 ceil(dt_new*(T-1)/2-20) 0 1]), hold on; 
    plot(((1:(T-1)/2))*dt_new,output((1:(T-1)/2)),'linewidth',6); if c==1; (legend('Target','Output')); end; set(gca,'fontsize',20);...
        xlabel('Time (sec)'); axis square, axis([1 ceil(dt_new*(T-1)/2-20) 0 1]),  xticks([1 floor(T*dt_new/4)-10 ceil((T-1)*dt_new/2)-20]); ...
    xticklabels({'0', num2str(floor(T*dt_new/4)-10), num2str(ceil((T-1)*dt_new/2)-20)}), yticks([0 0.5 1]);
    ylabel('Normalized activity');  set(gca,'fontsize',20);

    error_final(c)=mean(error_tot(end-100:end));

    clear input2
end

factor=[1,2,4,6,8];
figure; 
plot((T./factor)*dt_new,error_final,'k-*','linewidth',4,'MarkerSize',10);
yticks([0 300 600]);
ylabel('Mean total error');
xlabel('Sequence length (s)');
set(gca,'fontsize',20);
box off


