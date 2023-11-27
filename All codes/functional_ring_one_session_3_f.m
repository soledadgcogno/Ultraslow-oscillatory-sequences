%% Correlation VS distance on tissue and distance on ring

clear all
rec_data_path = 'C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_spikes='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
dpath_sorting='C:\Users\xscogno\MATLAB\Flavio2\Waves\Sorting\';
dpath='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';
WS_path='C:\Users\xscogno\MATLAB\Flavio2\Preprocessed_data\';

% mouse_name='L08M3';%'92229';%'L08M2';%'92229';%'L08M4';
% mouse_name='92229';%'L08M2';%'92229';%'L08M4';
% mouse_name='L09M4';%'L08M2';%'92229';%'L08M4';
mouse_name='L08M2';%'L08M2';%'92229';%'L08M4';
day=16;%19;%7;%19;%17;%V1 Day = 7
s=1; % number of session out of dates.sesnum(day) serssions

if mouse_name(2)=='0'
    mouse=[mouse_name(1),mouse_name(3:5)];
else
    mouse=mouse_name;
end

clus=10;
disc_phase=10;

%% Load files

load([rec_data_path,strcat('recording_dates_',mouse,'.mat')]);
munit=dates.ses{day}(s);
file_name_spk=[dpath ['spikes_120ms_Do_THR1p5_SNRH_FoV1','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
file_name_dff=[dpath ['DFF_120ms_Do_SNRH','_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];
file_name_anat=[dpath ['Anat_',mouse,'_day',num2str(day),'_MUnit',num2str(munit),'.mat']];


load(file_name_spk,'-mat');
spikes=full(spikes_d_s);
[N,T]=size(spikes);
N_or=N;
load([WS_path ['WS_Osc_15_sf7p73II_',mouse_name,'.mat']]);
dt=floor(WS_stat.dt(day,s));

%% Calculate correlation matrices

%Calculate sorting - Sorting PCA
[sorting_ascend,sorting_descend,sorting_0]=get_sorting(spikes);
[coeff, ~] = pca(spikes');
angle_pca=atan2(coeff(:,2),coeff(:,1));
[sorting_ascend_w,sorting_descend_w,sorting_0_w]=get_sorting_smoothed(spikes,dt);

% Lagged correlation 
maxlag=480;%250 - 520;
downsampling_factor=4;
FRp = spikes_downsample(spikes,N,downsampling_factor);
sf=7.73;
sf_new=sf/4;
bin_size_new=1/sf_new;
range=maxlag*bin_size_new;
% for i=1:N
%     FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); %Smoothing 
% end
FR=FRp;
count=0;
for i=1:N
    for j=i+1:N
        [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag)); %Check whether I need to zscore
        val=zscore(val);
        [v,in]=max(val);        
        if time(in)>=0
            Corr_mat(i,j)=v;
            Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in);
            lag_mat(j,i)=-time(in);            
        else
            Corr_mat(i,j)=v;
            Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in);
            lag_mat(j,i)=-time(in);
        end        
        clear val time
    end
end
% 
% figure
% imagesc(Corr_mat(sorting_pca,sorting_pca));
% colorbar

%% Figure Time lag VS Distance on ring

% Compute distances on the ring

for i=1:N_or
    for j=1:N_or
        delta_ring_pca(i,j)=(angdiff(angle_pca(i),angle_pca(j)));
    end
end

delta_ring_pca_vec=delta_ring_pca(:);
Corr_vec=Corr_mat(:);
lag_vec=lag_mat(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_ring_pca_vec(template)=[];
Corr_vec(template)=[];
lag_vec(template)=[];


%PCA
corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    length_c(i)=length(c);
    corrm(1:length(lag_vec(c)),i)=lag_vec(c);
    clear a b c
end

edges_time_lag=-maxlag:10:maxlag;
clear y aux
figure
for i=1:length(dist_ring_edges)-1
    aux=histogram(corrm(:,i),-maxlag:10:maxlag);
    y(i,:)=aux.Values;
end

% y normalization
y=y./sum(sum(y));

figure
imagesc(flip(y(:,2:95)'))
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-\pi','0','\pi'})
% yticks([1 13 24])
yticks([1 48 94])
% yticklabels({60,0,-60})
yticklabels({248,0,-248})
xlabel('Distance on ring (rad)');
ylabel('Time lag (s)');
co=colorbar;
% colormap puor
colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;
caxis([0 0.003]);
set(gca,'fontsize',16)
xtickangle(45)
axis square



%% Only if I want to use a subset of cells - Adapted to L9M4 Day17

locked=[1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	44	45	46	47	48	49	50	51	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	77	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	116	117	118	119	120	122	123	124	125	126	127	129	130	131	133	134	135	136	137	138	139	140	141	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	252	253	254	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	272	273	274	275	276	277	278	279	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	298	299	300	301	302	303	304	305	306	307	308	310	311	312	314	315	316	317	318	319	320	321	322	323	324	325	327	328	329	330	331	332	333	334	335	336	337	338	339	340	341	342	343	344	345	346	349	350	351	353	354	355	356	357	358	359	360	361	362	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377	378	379	380	381	382	384	385	386	387	388	389	390	391	392	393	394	395	396	397	398	399	400	402	403	404	405	406	407	408	410	411	412	413	414	415	416	417	418	419	420	421	422	423	424	425	426	427	428	429	430	431	432	433	434	435	436	437	438	439	440	441	442	443	444	445	446	447	448	449	450	451	452	453	454	455	456	457	458	459	460	461	462	463	464	465	466	467	468	469	470	471	472	473	474	475	476	477	478	479	480	481	482	483	484];
MVL=[0.830744331542054	0.852209457124681	0.874824164008808	0.837765436298618	0.884886678588304	0.860619722202613	0.739869889053488	0.662900630981510	0.274821064098610	0.747881510472578	0.576942432467035	0.256558336797857	0.847700028996744	0.701476556878339	0.619583856783054	0.706816808147832	0.489644001913405	0.936772229788196	0.606655486046027	0.834405687813549	0.690322230664419	0.0910703746391142	0.889192099114679	0.268122111540964	0.744310084004424	0.523319041186695	0.719116399178592	0.879527814623064	0.817037457010656	0.915584932556957	0.703369233224379	0.731770675325321	0.570186326724460	0.910297976278798	0.353466602048200	0.684766139225376	0.711580365981290	0.549252276414743	0.506375502445383	0.758349673467154	0.163134596539894	0.841279734043603	0.0988937766605530	0.844553166781336	0.704503721571479	0.355676938245125	0.803989449391397	0.814948594802772	0.763263283527798	0.809255853915511	0.876250181497788	0.0638872321319020	0.679908684598072	0.728493455341166	0.587450097372417	0.144120685495344	0.874936047976876	0.830886518222110	0.500821257937534	0.373953495959723	0.767065048739449	0.421659840733411	0.282384946474247	0.232507392040376	0.509752333703471	0.738521284031299	0.776024647407706	0.705611838897383	0.742898654887308	0.759548757518791	0.803469842445548	0.692741868423054	0.747909451034265	0.716848233878061	0.858769002672259	0.0771312254708299	0.499677932899914	0.0849238854560144	0.586314664236362	0.760972503264081	0.545096805255594	0.555667604208842	0.157024723249489	0.755391967871348	0.748273126032816	0.579982335051342	0.696831267273908	0.600739046602391	0.641197069440937	0.818337350576800	0.714006402119664	0.651620086520324	0.325105298934427	0.599804030111554	0.750318633290914	0.876496838660473	0.628358050672619	0.914139042236291	0.752424830825177	0.826125052515289	0.689027951174200	0.428709200004196	0.722625976207858	0.804003787256749	0.780698412836843	0.623042664093102	0.474030652869173	0.761021183369938	0.710933741610569	0.489464557734485	0.126012439309288	0.401398916504038	0.727721764262478	0.737761305711141	0.0761006272331358	0.823082365703706	0.187686319155072	0.779465350439127	0.724698049549872	0.647250142480122	0.0893138397977233	0.854014298495229	0.321384429610025	0.653353294047255	0.571411057612422	0.753522927575843	0.733048902109752	0.0892685390399692	0.815912822820013	0.759806625873523	0.622037222602977	0.0756736936877613	0.659599605858825	0.715854423282716	0.730117426007284	0.730518661607490	0.666997395894106	0.380930468919473	0.454168560010444	0.730878722225693	0.691339132952491	0.103268909353072	0.282967059634867	0.866975744952261	0.850902086151339	0.220207069111904	0.573438341379181	0.340912354267864	0.303750013527033	0.725498948199164	0.315010015690344	0.592778001540295	0.764576872959728	0.528814690040342	0.717289669778438	0.727462526538660	0.657113403181190	0.584182163803016	0.584066526991397	0.831438377346353	0.114034701278757	0.803922429075537	0.841904505496588	0.787498612998624	0.860474164755806	0.790538751307452	0.533994208336623	0.226356716894999	0.646582560754908	0.685678610541150	0.708889176269037	0.563637804298152	0.874816748100823	0.802961872208370	0.271877178706923	0.160990286266290	0.709336779470496	0.526708951641005	0.799060375239126	0.751700244651839	0.767305413864931	0.696456990900866	0.792101287549120	0.902730626402770	0.751545077517659	0.650373639710775	0.550887998829797	0.342819546586164	0.764908483957312	0.323335171663024	0.488618566293458	0.616667590154718	0.842583826412268	0.872868001648241	0.698054411915764	0.467128693384566	0.616551702283054	0.889250550078720	0.823835093598023	0.860353351268943	0.778176007941045	0.581342531149795	0.920590441982497	0.814796613055331	0.273202245383036	0.798125509525804	0.143738881676214	0.468922385100577	0.722707145373701	0.568434461459575	0.0440930602140277	0.785574837649311	0.327501236875652	0.199719916852643	0.127396576300164	0.775091122127954	0.398906646982433	0.806458441097752	0.573069328020557	0.612234008404577	0.708484959274934	0.902736071661730	0.485081565234663	0.828285818150137	0.760799264365731	0.795910875511052	0.841684033112803	0.353185233055241	0.844367724692663	0.619019903843404	0.419523368815794	0.680041061410917	0.234366271672849	0.824552451844064	0.719719395279835	0.857573036633275	0.840787697091788	0.242833071800109	0.162061058760341	0.695554564955296	0.294595000134031	0.777358762735335	0.193784839760228	0.651766064643241	0.704316558509820	0.807676929814604	0.625284072373446	0.317386657197974	0.524522363195175	0.756521916080553	0.0386604590902201	0.837387380911364	0.229848032112246	0.557303341488446	0.0456024171945760	0.848688011226692	0.491092174347856	0.885006972455600	0.579856091149977	0.593307739529807	0.524016938052724	0.446806882224421	0.836779646381043	0.827596423814701	0.189293107160198	0.412251548829239	0.236669361165743	0.498157484592805	0.822195835462465	0.261329100194867	0.103007500446498	0.503493928342180	0.542601965333207	0.641706943555410	0.433962040719324	0.722716182429415	0.554121663893569	0.242770630035187	0.468432442743842	0.106348257615858	0.742608632364942	0.721802421512540	0.661734597444034	0.681007076211206	0.221591709997693	0.814473480684526	0.813034978350016	0.796348307301553	0.617555646611808	0.715312527085092	0.511471012401214	0.854149787122322	0.571805897479650	0.730056136546613	0.565828159989268	0.480703888981454	0.0967381934189103	0.372619235830857	0.802240037304585	0.441311411653018	0.759471063117490	0.541332871443997	0.272128026243314	0.821911133423903	0.574679545427760	0.550506481127534	0.800502536853697	0.667290602968618	0.126197678866941	0.737899398287172	0.294554459671000	0.404071202178457	0.109975379683491	0.639680800878755	0.869550692113539	0.550072527673397	0.664715788977491	0.197783338946134	0.779698791172448	0.239270110945999	0.281544358499381	0.687467898359421	0.808258669216794	0.817236803880080	0.857883619497663	0.110083835905504	0.139057990940360	0.642355637281012	0.795716714189821	0.579635480817315	0.650515065253536	0.732231566274654	0.845041074372503	0.456713699407979	0.511698348385479	0.642872493158078	0.591330227728860	0.470619765327577	0.732776761616847	0.468927165611129	0.437084351977042	0.599668641891683	0.526023144926582	0.804622463110524	0.653700078671279	0.265864512015032	0.0817628097138937	0.106587476975465	0.901103901668677	0.357556242327838	0.492410666693356	0.0852622310363477	0.167101798896474	0.602239279526435	0.260069668832956	0.511847127888762	0.653290062321914	0.659830928236164	0.236252124409147	0.669789415324088	0.267823150595803	0.731335379177561	0.382937924552149	0.623763347789112	0.706031356684137	0.810925948399230	0.735039497623366	0.867687009803492	0.344854566840430	0.318835830930229	0.222522836248982	0.559727487405354	0.275022249158039	0.674432057296105	0.194674488627502	0.413883121134170	0.398947264436942	0.325905801893092	0.581965219698640	0.533567402054018	0.514377677184992	0.498974586058321	0.0994687815404226	0.793401935297849	0.573611371696457	0.150501487625123	0.336268399131045	0.338177264349133	0.836378920778389	0.823323819993200	0.685694792633052	0.526766696607085	0.628668583806268	0.164384455285415	0.524099295795025	0.411669593547452	0.448451638964862	0.391020863368041	0.570440158003127	0.670896183390078	0.130558806722901	0.785312777519020	0.195808357947406	0.728230015723792	0.496880506464046	0.498557799166351	0.318680739517698	0.801856321127603	0.0418894572493843	0.744375569461042	0.701361771757641	0.746132577355599	0.883533049621288	0.620692930190498	0.747193857535707	0.663290641980133	0.742292631358789	0.842265932135705	0.336022671428996	0.728332188758720	0.684516062682346	0.729021733486952	0.628793600984198	0.828432410893856	0.171666402042286	0.457743867379416	0.592166485681550	0.673954947092147	0.318008724170495	0.704372683488001	0.192999136739158	0.563712064885717	0.701880011591779	0.764328993290828	0.658527032092422	0.791089530649031	0.431819964399812	0.597539863867974	0.189015840136306	0.575453470678585	0.644609284947505	0.360258371191187	0.720761646716014	0.633103396334117	0.663850881595684	0.672067323179163	0.595652895173284	0.452405874983338	0.179607007341082	0.541652998650072	0.408712547906922	0.629108741667154	0.401396691131606	0.643680593395035	0.414011491089052	0.721903873768828	0.429475378817845	0.651008627880600	0.666511121540235	0.684858967977121	0.631307298691966	0.896504308207694	0.463630559720713	0.628892735379834	0.599348860507539	0.560870455907645	0.732705777948138	0.639462874459642	0.687420931645880	0.674973682166407	0.802411352933304	0.260607993363177	0.678087062627100	0.908572934270052	0.440806377695242	0.770914047334443	0.447681766690461	0.834244808129714	0.553110293562001	0.255401231862458	0.660023156344657	0.182728446356273	0.515176371919849	0.140345329092547];

clear spikes_cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Drop of 40%
p80=prctile(MVL,80);
cells_locked=find(MVL>p80);
cells_nlocked=find(MVL<p80);

spikes_cells=spikes(cells_nlocked,:);
N=size(spikes_cells,1);
% Lagged correlation 
maxlag=480;
downsampling_factor=4;

clear FRp FR
FRp = spikes_downsample(spikes_cells,N,downsampling_factor);
for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); %Smoothing 
end
FR=FRp;
count=0;
for i=1:N
    for j=i+1:N
        [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag)); %Check whether I need to zscore

        [v,in]=max(val);        
        if time(in)>=0
            Corr_mat(i,j)=v;
            Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in);
            lag_mat(j,i)=-time(in);            
        else
            Corr_mat(i,j)=v;
            Corr_mat(j,i)=v;
            lag_mat(i,j)=time(in);
            lag_mat(j,i)=-time(in);
        end        
        clear val time
    end
end

[coeff, ~] = pca(spikes_cells');
angle_pca_new=atan2(coeff(:,2),coeff(:,1));

% angle_pca_new=angle_pca(cells_nlocked);
for i=1:N
    for j=1:N
        delta_ring_pca(i,j)=(angdiff(angle_pca_new(i),angle_pca_new(j)));
    end
end

delta_ring_pca_vec=delta_ring_pca(:);
Corr_vec=Corr_mat(:);
lag_vec=lag_mat(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_ring_pca_vec(template)=[];
Corr_vec(template)=[];
lag_vec(template)=[];


%PCA
corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    length_c(i)=length(c);
    corrm(1:length(lag_vec(c)),i)=lag_vec(c);
    clear a b c
end

clear y aux
figure
for i=1:length(dist_ring_edges)-1
    aux=histogram(corrm(:,i),-maxlag:10:maxlag);
    y(i,:)=aux.Values;
end

% y normalization
% y=y./sum(sum(y));

figure
imagesc((y(:,2:95)'))
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
% yticks([1 13 24])
yticks([1 48 96])
% yticklabels({60,0,-60})
yticklabels({240,0,-240})
xlabel('Distance on ring (rad)');
ylabel('Time lag (s)');
co=colorbar;
% colormap puor
colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;
% caxis([0 0.003]);
caxis([0 500]);

set(gca,'fontsize',16)
xtickangle(45)
axis square



%% Figure Pearson VS Distance on ring

% Compute distances on the ring

for i=1:N_or
    for j=1:N_or
        delta_ring_pca(i,j)=(angdiff(angle_pca(i),angle_pca(j)));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Correlation at time lag=0

% Lagged correlation 
maxlag=480;%250 - 520;
downsampling_factor=4;
FRp = spikes_downsample(spikes,N,downsampling_factor);
for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),dt,'g')); %Smoothing
end
FR=FRp;
count=0;
for i=1:N
    for j=1:N
        [val,time]=(xcorr(FR(i,:),FR(j,:),maxlag)); %Check whether I need to zscore
        val=zscore(val);
        aux=find(time==0);
        
        Corr_mat_0(i,j)=val(aux);
%         Corr_mat_0(j,i)=val(aux);
        
        clear val time
    end
end

Corr_mat_0_vec=Corr_mat_0(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

Corr_mat_0_vec(template)=[];


corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    length_c(i)=length(c);
    corrm(1:length(Corr_mat_0(c)),i)=Corr_mat_0(c);
    clear a b c
end

clear y aux
figure
for i=1:length(dist_ring_edges)-1
    aux=histogram(corrm(:,i),-5:0.05:10);
    y(i,:)=aux.Values;
end

% y normalization
y=y./sum(sum(y));

figure
imagesc(flip(y'))
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
% yticks([1 13 24])
yticks([1 20 40])
% yticklabels({60,0,-60})
yticklabels({1,0,-1})
xlabel('Distance on ring (rad)');
ylabel('Correlation lag 0');
co=colorbar;
% colormap puor
% colormap gmt_ocean;
% colormap cividis
colormap gmt_nighttime;
caxis([0 0.004]);
set(gca,'fontsize',18)
xtickangle(45)
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %Pearson

downsampling_factor=4;
X = spikes_downsample(spikes,N,downsampling_factor);
FRp = spikes_downsample(spikes,N,downsampling_factor);
for i=1:N
    FR(i,:)=full(fire_rate(FRp(i,:),29,'g')); %Smoothing in about 10 seconds
end
       
Pearson=corr(FR',FR');

delta_ring_pca_vec=delta_ring_pca(:);
Pearson_vec=Pearson(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_ring_pca_vec(template)=[];
Pearson_vec(template)=[];

corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    length_c(i)=length(c);
    corrm(1:length(Pearson_vec(c)),i)=Pearson_vec(c);
    clear a b c
end

clear y aux
figure
for i=1:length(dist_ring_edges)-1
    aux=histogram(corrm(:,i),-1:0.05:1);
    y(i,:)=aux.Values;
end

% y normalization
y=y./sum(sum(y));

figure
imagesc(flip(y'))
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
% yticks([1 13 24])
yticks([1 20 40])
% yticklabels({60,0,-60})
yticklabels({1,0,-1})
xlabel('Distance on ring (rad)');
ylabel('Pearson Correlation');
co=colorbar;
% colormap puor
% colormap gmt_ocean;
% colormap cividis
colormap gmt_nighttime;
caxis([0 0.004]);
set(gca,'fontsize',18)
xtickangle(45)
axis square


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Correlation values at time lag different than zero
% angle_pca_new=angle_pca(cells_nlocked);

Corr_vec=Corr_mat(:);
aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);
Corr_vec(template)=[];

clear corrm;
corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    length_c(i)=length(c);
    corrm(1:length(Corr_vec(c)),i)=Corr_vec(c);
    clear a b c
end
clear y aux
figure
for i=1:length(dist_ring_edges)-1
    aux=histogram(corrm(:,i),0:0.5:10);
    y(i,:)=aux.Values;
end


% y normalization
y=y./sum(sum(y));

figure
imagesc(flip(y'))
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
% yticks([1 13 24])
yticks([1 10 20])
% yticklabels({60,0,-60})
yticklabels({10,5,0})
xlabel('Distance on ring (rad)');
ylabel('Correlation');
co=colorbar;
% colormap puor
% colormap gmt_ocean;
% colormap cividis
colormap gmt_nighttime;
caxis([0 0.003]);
set(gca,'fontsize',18)
xtickangle(45)
axis square





%% Figures correlation PCA %CHECK

% Compute distances on the ring

for i=1:N_or
    for j=1:N_or
        delta_ring_pca(i,j)=(angdiff(angle_pca(i),angle_pca(j)));
%         delta_ring_umap(i,j)=(angdiff(angleumap(i),angleumap(j))); 
    end
end

delta_ring_pca_vec=delta_ring_pca(:);
% delta_ring_umap_vec=delta_ring_umap(:);
Corr_vec=Corr_mat(:);
lag_vec=lag_mat(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_ring_pca_vec(template)=[];
% delta_ring_umap_vec(template)=[];
Corr_vec(template)=[];
lag_vec(template)=[];


%PCA
corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    corrm(1:length(Corr_vec(c)),i)=Corr_vec(c);
    clear a b c
end

clear y aux
figure
for i=1:length(dist_ring_edges)-1
%     aux=histogram(corrm(:,i),-maxlag:10:maxlag);
    aux=histogram(corrm(:,i),-0:50:600);
    y(i,:)=aux.Values;
end

% y normalization
y=y./sum(sum(y));

figure
imagesc(flip(y)')
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
yticks([1 13 24])
yticklabels({60,0,-60})
xlabel('Distance on ring (rad)');
ylabel('Time lag (s)');
co=colorbar;
% colormap puor
colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;
caxis([0 0.02]);
set(gca,'fontsize',16)
xtickangle(45)
axis square
% 
%% Figures lag PCA FOR SUBSET

to_remove=setdiff(1:N_or,cells);
angle_pca(to_remove)=[];

for i=1:N
    for j=1:N
        delta_ring_pca(i,j)=(angdiff(angle_pca(i),angle_pca(j)));
%         delta_ring_umap(i,j)=(angdiff(angleumap(i),angleumap(j))); 
    end
end

delta_ring_pca_vec=delta_ring_pca(:);
% delta_ring_umap_vec=delta_ring_umap(:);
Corr_vec=Corr_mat(:);
lag_vec=lag_mat(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_ring_pca_vec(template)=[];
% delta_ring_umap_vec(template)=[];
Corr_vec(template)=[];
lag_vec(template)=[];


%PCA
corrm=nan(30000,11);
dist_ring_edges=-3.14:2*3.14/11:3.14;
for i=1:length(dist_ring_edges)-1
    a=find(delta_ring_pca_vec>dist_ring_edges(i));
    b=find(delta_ring_pca_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    corrm(1:length(lag_vec(c)),i)=lag_vec(c);
    clear a b c
end

clear y aux
figure
for i=1:length(dist_ring_edges)-1
    aux=histogram(corrm(:,i),-120:10:120);
    y(i,:)=aux.Values;
end

figure
imagesc(flip(y)')
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 6 11])
xticklabels({'-3.1 - -2.6','-0.3 - 0.3','2.6 - 3.14'})
yticks([1 13 24])
yticklabels({30,0,-30})
xlabel('Distance on ring (rad)');
ylabel('Time lag (s)');
co=colorbar;
% colormap puor
colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;
caxis([0 1500]);
set(gca,'fontsize',16)
xtickangle(45)
axis square
% 


%% Figures lag PCA and UMAP

% Compute distances on the ring

figure
scatter(coeff(:,1),coeff(:,2),'o','filled');
alpha 0.3

for i=1:N_or
    for j=1:N_or
        delta_ring_pca(i,j)=(angdiff(angle_pca(i),angle_pca(j)));
%         delta_ring_umap(i,j)=(angdiff(angleumap(i),angleumap(j))); 
    end
end

for i=1:N_or
    for j=1:N_or
        
        ini=find(sorting_pca==i);
        fin=find(sorting_pca==j);
        dist=(ini-fin);
        
        delta_raster(i,j)=(ini-fin);
        delta_raster(j,i)=-(ini-fin);
    end
end

delta_raster_vec=delta_raster(:);
delta_ring_pca_vec=delta_ring_pca(:);
Corr_vec=Corr_mat(:);
lag_vec=lag_mat(:);

aux=eye(N,N);
aux2=aux(:);
template=find(aux2==1);

delta_raster_vec(template)=[];
delta_ring_pca_vec(template)=[];
Corr_vec(template)=[];
lag_vec(template)=[];

%PCA
corrm=nan(30000,11);
dist_ring_edges=0:N/10:N;
for i=1:length(dist_ring_edges)-1
    a=find(delta_raster_vec>dist_ring_edges(i));
    b=find(delta_raster_vec<=dist_ring_edges(i+1));
    c=intersect(a,b);
    corrm(1:length(lag_vec(c)),i)=lag_vec(c);
    clear a b c
end

clear y aux
figure
for i=1:length(dist_ring_edges)-1
    figure
    aux=histogram(corrm(:,i),0:10:480);
    y(i,:)=aux.Values./sum(aux.Values);
end

% y=y./sum(sum(y));

figure
imagesc(flip(y'))
% xticklabels({'-3.1 - -2.6','-2.6 - -2','-2 - -1.4','-1.4 - -0.8','-0.8 - -0.3','-0.3 - 0.3','0.3 - 0.8','0.8 - 1.4','1.4 - 2','2 - 2.6' , '2.6 - 3.14'})
xticks([1 5 10])
xticklabels({'0- 52','208-260','468 - 520'})
yticks([1 13 24])
yticklabels({120,60,0})
xlabel('Distance between cells (cells)');
ylabel('Time lag (s)');
co=colorbar;
% colormap puor
colormap gmt_ocean;
% colormap cividis
% colormap gmt_nighttime;
caxis([0 0.1]);
set(gca,'fontsize',16)
xtickangle(45)
axis square
