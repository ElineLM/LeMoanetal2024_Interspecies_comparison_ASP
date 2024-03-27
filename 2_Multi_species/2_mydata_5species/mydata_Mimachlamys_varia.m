function [data, auxData, metaData, txtData, weights] = mydata_Mimachlamys_varia
%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Bivalvia'; 
metaData.order      = 'Ostreoida'; 
metaData.family     = 'Pectinidae';
metaData.species    = 'Mimachlamys_varia'; 
metaData.species_en = 'Variegated scallop'; 

metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MANE','MAE','MAm'};
metaData.ecoCode.habitat = {'0jMp', 'jiMb'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biPp'};
metaData.ecoCode.gender  = {'Hsb'}; %bidirectional hermaphrodite
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(16); % K. body temp
metaData.data_0     = {'ab'; 'aj'; 'tp'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li';'Ri'; 'GSI'; 'E0'};      
metaData.data_1     = {'t-L'; 't-Wd'; 'L-Wd'}; 

metaData.COMPLETE = 3; % using criteria of LikaKear2011      

metaData.author   = {'Le Moan & RegnierBrisson'};    
metaData.date_subm = [2023 01 11];              
metaData.email    = {'Laure.Regnier.Brisson@ifremer.fr'};            
metaData.address  = {'Univ Brest'};   

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'sta@akvaplan.niva.no'}; 
metaData.date_acc    = [2021 02 21]; 

%% set data
% -------------------------------------------------------------------------
% zero-variate data
% -------------------------------------------------------------------------

data.ab = 2;      units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = 'TinduffHatchery';   
  temp.ab = C2K(19.5);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
  comment.ab = 'Based on Tinduff hatchery protocol, with controlled temperature';
data.aj = 18;    units.aj = 'd';    label.aj = 'time since fertilisation at metam'; bibkey.aj = 'TinduffHatchery';
  temp.aj = C2K(18);  units.temp.aj = 'K'; label.temp.aj = 'temperature';
  comment.aj = 'Based on Tinduff hatchery protocol, with controlled temperature';
data.tp = 270;     units.tp = 'd';    label.tp = 'time since birth at puberty'; bibkey.tp = {'RegnierBrisson'};  
  temp.tp = C2K(14);  units.temp.tp = 'K'; label.temp.tp = 'temperature'; 
data.am = 7*365;   units.am = 'd';    label.am = 'life span';                bibkey.am = 'sealifebase';   
  temp.am = C2K(16);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
  comment.am = 'Agree with the original data set of AmP, from Sealifebase website';

 
data.Lb  = 0.0105;   units.Lb  = 'cm';  label.Lb  = 'shell length at birth';   bibkey.Lb  = 'TinduffHatchery';
  comment.Lb = 'Based on Tinduff hatchery protocol, with controlled temperature';
  temp.Lb = C2K(19.5); units.temp.Lb = 'K'; label.temp.Lb = 'temperature';
data.Lj  = 0.0215;   units.Lj  = 'cm';  label.Lj  = 'metamorphosis shell length';   bibkey.Lj  = 'TinduffHatchery';
  comment.Lj = 'Based on Tinduff hatchery protocol, with controlled temperature';
  temp.Lj = C2K(18); units.temp.Lj = 'K'; label.temp.Lj = 'temperature';
data.Lp  = 1.5;    units.Lp  = 'cm';  label.Lp  = 'shell height at puberty'; bibkey.Lp  = 'RegnierBrisson';
  comment.Lp = 'Suivi N20';  %ajout 24/01/2023
  
  %added on 15th May 2023
data.Li = 7; units.Li = 'cm'; label.Li = 'infinite shell height'; bibkey.Li = 'RegnierBrisson';
    comment.Li = 'guess from in situ monitoring and Perodou & Latrouite 1981';

data.R_L  = (2*4000000/365); units.R_L  = '#/d'; label.R_L  = 'reprod rate at a large size';     bibkey.R_L  = 'TinduffHatchery'; 
   temp.R_L = C2K(14); units.temp.R_L = 'K'; label.temp.R_L = 'temperature';
   comment.R_L = 'based on Tinduff Hatchery observations, for the individuals of 5-6 cm length, nb of eggs for one spawning event, and 2 events per year';
L_R.R_L = 5.5;   units.L_R.R_L = 'cm';         label.L_R.R_L = 'length at reproduction data'; 
   
 data.E_0 = 3e-3; units.E_0  = 'J';  label.E_0  = 'egg energy content'; bibkey.E_0  = 'RegnierBrisson';
  comment.E_0 = 'calculated with a 30% GSI .E0 = 1e-3J for Crassostrea gigas';  
  
data.GSI = 0.3; units.GSI  = '-';  label.GSI  = 'gonado-somatic index'; bibkey.GSI  = 'RegnierBrisson';
temp.GSI = C2K(14); units.temp.GSI = 'K'; label.temp.GSI = 'temperature';
 comment.GSI = 'calculated from 3-year individuals at spawning. We assume that they can spawn (at least) twice a year ';  
  
data.Wwp = 0.2;     units.Wwp = 'g';    label.Wwp = 'wet weight at puberty';    bibkey.Wwp = 'LetaAud1956';
    comment.Wwp = 'the weight at 1.5 cm, fixed as puberty length in our case';
data.Wdi = 1.6;     units.Wdi = 'g';    label.Wdi = 'ultimate dry weight';      bibkey.Wdi = 'RegnierBrisson';
    comment.Wdi = 'mean obs collector at 3 years';

% -------------------------------------------------------------------------
% uni-variate data
% -------------------------------------------------------------------------

%time-length larvae
% tLTinduff = [ ... % time since spawning event (d). shell length (cm)
% % 2   0.0105
% 5	0.01294
% 5	0.01326
% 5	0.01281
% 5	0.0129
% 5	0.01339
% 5	0.01251
% 5	0.01324
% 5	0.01319
% 5	0.01293
% 5	0.01316
% 5	0.01288
% 5	0.01295
% 5	0.01324
% 5	0.01298
% 5	0.01279
% 5	0.01305
% 5	0.01298
% 5	0.01263
% 5	0.01296
% 5	0.01289
% 5	0.0133
% 5	0.01331
% 5	0.01328
% 5	0.01286
% 5	0.01303
% 5	0.01316
% 5	0.01324
% 5	0.01304
% 5	0.01331
% 5	0.01301
% 7	0.01389
% 7	0.01453
% 7	0.01419
% 7	0.01419
% 7	0.01403
% 7	0.01414
% 7	0.01409
% 7	0.01381
% 7	0.01428
% 7	0.01399
% 7	0.01393
% 7	0.01456
% 7	0.0141
% 7	0.01419
% 7	0.01411
% 7	0.01423
% 7	0.0143
% 7	0.01436
% 7	0.01423
% 7	0.01443
% 7	0.01408
% 7	0.01385
% 7	0.01454
% 7	0.0139
% 7	0.01411
% 7	0.01449
% 7	0.01411
% 7	0.01423
% 7	0.01418
% 7	0.01449
% 9	0.01593
% 9	0.01565
% 9	0.01558
% 9	0.01568
% 9	0.01543
% 9	0.01556
% 9	0.01569
% 9	0.01529
% 9	0.01626
% 9	0.0159
% 9	0.01569
% 9	0.01551
% 9	0.01515
% 9	0.01558
% 9	0.01591
% 9	0.01528
% 9	0.01581
% 9	0.01544
% 9	0.01586
% 9	0.01508
% 9	0.01553
% 9	0.01561
% 9	0.01563
% 9	0.0159
% 9	0.01553
% 9	0.01578
% 9	0.01563
% 9	0.01546
% 9	0.01571
% 9	0.0158
% 11	0.01756
% 11	0.0174
% 11	0.01739
% 11	0.01714
% 11	0.01741
% 11	0.01736
% 11	0.01726
% 11	0.01771
% 11	0.01744
% 11	0.01664
% 11	0.01707
% 11	0.01747
% 11	0.01709
% 11	0.01763
% 11	0.0171
% 11	0.01741
% 11	0.01744
% 11	0.0166
% 11	0.01703
% 11	0.01699
% 11	0.01757
% 11	0.01731
% 11	0.01696
% 11	0.01737
% 11	0.01713
% 11	0.0172
% 11	0.01733
% 11	0.01673
% 11	0.01739
% 11	0.01734
% 13	0.0201
% 13	0.01933
% 13	0.01893
% 13	0.01836
% 13	0.01911
% 13	0.01823
% 13	0.01859
% 13	0.01865
% 13	0.0193
% 13	0.01788
% 13	0.01845
% 13	0.01844
% 13	0.01878
% 13	0.01803
% 13	0.01878
% 13	0.01883
% 13	0.01895
% 13	0.01868
% 13	0.01904
% 13	0.01875
% 13	0.01845
% 13	0.01875
% 13	0.0182
% 13	0.01935
% 13	0.01946
% 13	0.01925
% 13	0.0187
% 13	0.01853
% 13	0.01846
% 13	0.01851
% 15	0.01983
% 15	0.02014
% 15	0.02003
% 15	0.01994
% 15	0.0201
% 15	0.01999
% 15	0.0204
% 15	0.01974
% 15	0.02021
% 15	0.02023
% 15	0.02023
% 15	0.02006
% 15	0.01978
% 15	0.01973
% 15	0.01935
% 15	0.02019
% 15	0.0195
% 15	0.02013
% 15	0.01978
% 15	0.02015
% 15	0.02033
% 15	0.02009
% 15	0.02021
% 15	0.02004
% 15	0.01988
% 15	0.02035
% 15	0.01949
% 15	0.01939
% 15	0.02016
% 15	0.01998
% 17	0.02024
% 17	0.02052
% 17	0.01976
% 17	0.0204
% 17	0.02036
% 17	0.02008
% 17	0.02064
% 17	0.02022
% 17	0.0201
% 17	0.02033
% 17	0.0206
% 17	0.02038
% 17	0.02025
% 17	0.0197
% 17	0.02008
% 17	0.0202
% 17	0.02005
% 17	0.01973
% 17	0.0202
% 17	0.02063
% 17	0.0201
% 17	0.01957
% 17	0.0205
% 17	0.01977
% 17	0.0205
% 17	0.02025
% 17	0.0195
% 17	0.0199
% 17	0.0195
% 17	0.0207
% % 18  0.0215  %aj     Lj 
% 30	0.044
% 30	0.0561
% 30	0.0506
% 30	0.051
% 30	0.0418
% 30	0.0522
% 30	0.0468
% 30	0.0435
% 30	0.0479
% 30	0.0607
% 30	0.047
% 30	0.0502
% 30	0.0543
% 30	0.041
% 30	0.0343
% 30	0.0578
% 30	0.0569
% 30	0.0604
% 30	0.0394
% 30	0.0659
% 30	0.0524
% 30	0.0617
% 30	0.0586
% 30	0.0396
% 30	0.0593
% 30	0.0561
% 30	0.0366
% 30	0.0372
% 30	0.0444
% 30	0.0574
% 35	0.0717
% 35	0.0797
% 35	0.0684
% 35	0.0709
% 35	0.078
% 35	0.0656
% 35	0.0768
% 35	0.0877
% 35	0.0637
% 35	0.0725
% 35	0.0656
% 35	0.0779
% 35	0.0778
% 35	0.0775
% 35	0.088
% 35	0.0799
% 35	0.0682
% 35	0.0649
% 35	0.0671
% 35	0.0618
% 35	0.0847
% 35	0.0793
% 35	0.0476
% 35	0.077
% 35	0.0574
% 35	0.0704
% 35	0.0677
% 35	0.0662
% 35	0.064
% 35	0.098
%     ];
% 
% %%% differenciate larvae and juvenile from one data set
% temp_l = tLTinduff(:,2);
% % data.tLlarvae(:,1) = tLTinduff(tLTinduff(:,1) < data.aj); %only for values under the metamorphosis age
% % data.tLlarvae(:,2) = temp_l(tLTinduff(:,1) < data.aj);
% data.tLlarvae(:,1) = [data.ab ; tLTinduff(tLTinduff(:,1) < data.aj)]; %only for values under the metamorphosis age
% data.tLlarvae(:,2) = [data.Lb ; temp_l(tLTinduff(:,1) < data.aj)];
% units.tLlarvae   = {'d', 'cm'};  label.tLlarvae = {'time', 'shell length'};  %shell length, parallel to the hinge line, biggest length
% temp.tLlarvae    = C2K(18);  units.temp.tLlarvae = 'K'; label.temp.tLlarvae = 'temperature'; %for the moment, 16°C to be the mean of temperature during nursery
% bibkey.tLlarvae = 'TinduffHatchery';
% comment.tLlarvae = 'Monitoring in the Tinduff Hatchery, controlled temperatures';
% 
% %%% Relation from Laure Regnier-Brisson data: Shell height = 1.0789 * shell length + 0.01905 (in cm), R² = 0.9841, N = 3208 individuals
% % data.tLjuve(:,1) = tLTinduff(tLTinduff(:,1) >= data.aj); %only for values above the metamorphosis age
% % data.tLjuve(:,2) = temp_l(tLTinduff(:,1) >= data.aj);
% tLTinduff(:,2) = 1.1451 * tLTinduff(:,2); % 
% 
% data.tLjuve(:,1) = [data.aj ; tLTinduff(tLTinduff(:,1) >= data.aj)]; %only for values above the metamorphosis age
% data.tLjuve(:,2) = [data.Lj ; temp_l(tLTinduff(:,1) >= data.aj)];
% units.tLjuve   = {'d', 'cm'};  label.tLjuve = {'time', 'shell height'};  %shell height, perpendicular to the hinge line
% temp.tLjuve    = C2K(18);  units.temp.tLjuve = 'K'; label.temp.tLjuve = 'temperature'; %for the moment, 16°C to be the mean of temperature during nursery
% bibkey.tLjuve = 'TinduffHatchery';
% comment.tLjuve = 'Monitoring in the Tinduff Hatchery, controlled temperatures';

%%% -----------------------------------------------------------------------
%%% Time - Shell height - Ash free total dry Weight - Ash free gonad dry weight
%%% Data from Laure Régnier-Brisson thesis, in Saint-Anne then Roscanvel
%%% The temperature is similar in the 2 places, so the real temperature of
%%% Saint-Anne will be considered.
%%% The temperature is defined with Fourier parameters determined in a
%%% previous step
%%% The translation shell length VS shell height (in mm) could be defined by:
            %%%     Shell length = 0.912 * Shell height - 1.309 for all size (N = 3208) and R² = 0.984
            %%%     Shell length = 0.843 * Shell height + 0.019 for height < 15 mm (N = 99) and R² = 0.972
%modified on 31/01/2023 - add the missing values with W as NaN and change
%the date to have the 0 at 24/05/2019 date of spawning
tLW19SA = [ %time since fertilisation (days), shell height (cm), ash free total dry weight (g), ash free gonad dry weight (g),
    % at 152 days, the organisms are still juveniles, but the shape is
    % considered the same for the whole data set
152	1.743	0.012	0
152	1.414	0.007	0
152	1.118	0.006	0
152	1.348	0.005	0
152	1.389	0.004	0
152	1.351	0.008	0
152	1.224	0.006	0
152	1.436	0.011	0
152	1.36	0.007	0
152	1.764	0.009	0
152	1.49	0.01	0
152	1.218	0.006	0
152	1.315	0.009	0
152	1.139	0.002	0
152	1.439	0.006	0
152	1.502	0.008	0
152	1.596	0.005	0
152	1.354	0.005	0
152	1.4     0.004	0
152	1.376	0.005	0
152	1.116	0.001	0
152	1.273	0.004	0
152	1.391	0.005	0
152	1.307	0.003	0
152	1.3     0.004	0
152	1.211	0.003	0
152	1.216	0.001	0
152	1.433	0.003	0
152	1.216	0.003	0
152	1.446	0.005	0
152	1.544	0.008	0
152	1.087	0.004	0
152	1.239	0.003	0
152	1.253	0.003	0
152	1.283	0.003	0
182	1.205	0.014	0
182	1.42	0.014	0
182	1.257	0.01	0
182	1.436	0.025	0
182	1.339	0.019	0
182	1.383	0.012	0
182	1.922	0.047	0
182	1.106	0.004	0
182	1.362	0.01	0
182	1.666	0.026	0
182	1.936	0.042	0
182	1.319	0.011	0
182	1.396	0.019	0
182	1.425	0.017	0
182	1.608	0.028	0
182	1.224	0.013	0
182	1.15	0.008	0
182	1.343	0.01	0
182	1.346	0.013	0
182	1.453	0.019	0
182	1.49	0.021	0
182	1.71	0.032	0
182	1.484	0.019	0
182	1.146	0.009	0
182	1.18	0.007	0
397	3.023	0.272	0
397	2.949	0.211	0
397	2.661	NaN	0
397	2.185	0.083	0
397	2.617	0.134	0
397	3.065	0.268	0
397	2.469	0.146	0
397	3.412	0.424	0
397	2.909	0.233	0
397	2.328	0.104	0
397	2.031	0.08	0
397	3.089	NaN	0
397	3.631	0.369	0
397	1.933	0.063	0
397	2.394	0.159	0
397	2.607	0.147	0
397	1.572	0.028	0
397	2.642	0.158	0
397	2.553	0.184	0
397	2.89	0.236	0
397	3.404	0.371	0
397	2.538	0.14	0
397	2.272	0.102	0
397	2.547	0.127	0
397	2.76	0.169	0
397	3.317	NaN	0
397	2.44	NaN	0
397	2.991	NaN	0
397	2.319	NaN	0
397	2.357	NaN	0
423	3.209	0.386	0.057
423	2.567	0.166	0.017
423	2.795	0.28	0.034
423	3.031	0.436	0.077
423	2.168	0.121	0.011
423	3.269	0.461	0.08
423	2.683	0.266	0.029
423	3.53	0.516	0.056
423	3.192	0.412	0.046
423	2.533	0.197	0.035
423	2.637	0.228	0.028
423	2.868	0.288	0.036
423	3.124	0.364	0.05
423	2.881	0.219	0.031
423	2.994	0.311	0.035
423	2.716	0.26	0.037
423	3.484	0.522	0.076
423	3.103	0.419	0.077
423	2.954	0.299	0.033
423	3.016	0.34	0.066
423	3.174	0.405	0.04
423	2.9	0.355	0.04
423	2.254	NaN	NaN
423	2.747	0.254	0.036
423	3.122	0.393	0.053
423	3.188	0.334	0.044
423	2.988	0.33	0.033
423	3.112	0.486	0.001
423	3.016	0.4	0.056
423	3.412	0.508	0.1
451	3.057	NaN	0.064
451	2.26	NaN	0.013
451	3.417	NaN	0.087
451	3.128	NaN	0.045
451	3.009	NaN	0.05
451	3.186	0.418	0.053
451	3.34	0.415	0.056
451	3.392	0.485	0.07
451	2.881	0.27	0.022
451	3.268	0.347	0.049
451	3.126	0.39	0.05
451	3.647	0.475	0.063
451	2.86	0.259	0.025
451	3.026	0.34	0.042
451	3.875	0.726	0.089
451	3.261	0.414	0.07
451	3.098	0.366	0.058
451	3.224	0.442	0.095
451	2.905	0.333	0.048
451	2.544	0.083	0.005
451	2.92	0.243	0.019
451	3.358	0.473	0.066
451	3.372	0.359	0.047
451	3.164	0.38	0.059
451	2.657	0.208	0.036
451	3.513	NaN	0.073
451	3.095	NaN	0.041
451	3.94	NaN	0.075
451	2.944	NaN	0.068
451	3.776	NaN	0.134
479	3.676	NaN	0.054
479	3.768	NaN	0.038
479	3.99	NaN	0.075
479	3.644	NaN	0.064
479	3.512	NaN	0.036
479	3.65	0.443	0.032
479	3.577	0.461	0.028
479	3.32	0.392	0.011
479	3.524	0.449	0.0090
479	4.198	0.765	0.039
479	3.453	0.418	0.037
479	4.3	0.822	0.017
479	3.574	0.62	0.048
479	3.38	0.421	0.022
479	3.416	0.554	0.046
479	2.826	0.333	0.013
479	4.629	0.905	0.072
479	4.093	0.861	0.034
479	4.353	0.869	0.033
479	3.795	0.528	0.05
479	4.007	0.706	0.03
479	3.381	0.422	0.024
479	4.271	0.903	0.044
479	3.885	0.357	0.02
479	3.437	0.513	0.013
479	3.752	NaN	0.024
479	3.786	NaN	0.023
479	3.775	NaN	0.084
479	3.447	NaN	0.047
479	4.24	NaN	0.066
508	3.592	NaN	0.012
508	3.407	NaN	0.0090
508	3.528	NaN	0.012
508	2.799	NaN	0.011
508	2.914	NaN	0.0060
508	3.846	0.724	0.033
508	3.62	0.506	0.028
508	3.093	0.322	0.012
508	3.95	0.714	0.014
508	4.158	0.729	0.012
508	4.332	0.749	0.024
508	4.273	0.835	0.021
508	3.759	0.51	0.013
508	4.239	0.928	0.021
508	3.811	0.574	0.015
508	3.138	0.432	0.013
508	3.866	0.653	0.017
508	3.937	0.685	0.017
508	3.75	0.597	0.014
508	4.123	0.681	0.016
508	3.551	0.542	0.011
508	3.545	0.315	0.015
508	3.134	0.331	0.0090
508	3.132	0.355	0.011
508	3.033	0.266	0.005
508	3.532	NaN	0.015
508	3.369	NaN	0.0090
508	4.205	NaN	0.045
508	3.944	NaN	0.022
508	3.386	NaN	0.013
536	3.946	NaN	0.012
536	3.441	NaN	0.015
536	4.655	NaN	0.014
536	3.887	NaN	0.014
536	3.857	NaN	0.01
536	4.109	0.528	0.018
536	3.12	0.22	0.005
536	4.088	0.357	0.0070
536	4.294	0.708	0.012
536	4.222	0.497	0.0080
536	4.036	0.411	0.015
536	3.463	0.574	0.019
536	3.86	0.852	0.02
536	3.635	0.374	0.013
536	3.217	0.798	0.018
536	3.811	0.816	0.018
536	3.246	0.841	0.01
536	4.436	0.794	0.016
536	4.582	0.688	0.025
536	3.709	0.301	0.0090
536	4.098	0.722	0.014
536	3.831	0.657	0.016
536	3.703	0.882	0.017
536	4.136	0.777	0.0080
536	3.805	0.818	0.017
536	4.481	NaN	0.02
536	3.472	NaN	0.0080
536	3.27	NaN	0.017
536	4.431	NaN	0.026
536	3.752	NaN	0.024
593	3.853	NaN	0.0090
593	4.517	NaN	0.013
593	3.879	NaN	0.01
593	3.999	NaN	0.011
593	4.471	NaN	0.011
593	4.599	0.697	0.02
593	4.596	0.821	0.021
593	4.147	0.761	0.019
593	4.6	0.735	0.014
593	4.362	0.653	0.005
593	4.371	0.854	0.012
593	4.61	0.804	0.024
593	3.841	0.476	0.011
593	4.345	0.701	0.015
593	4.365	0.841	0.015
593	3.652	0.409	0.0060
593	4.545	0.748	0.012
593	4.176	0.518	0.011
593	4.009	0.577	0.015
593	3.023	0.259	0.0080
593	4.153	0.493	0.013
593	4.929	0.979	0.02
593	4.561	0.761	0.019
593	4.498	0.579	0.0060
593	4.694	0.858	0.014
593	4.266	NaN	0.014
593	3.199	NaN	0.005
593	3.509	NaN	0.012
593	4.652	NaN	0.018
593	4.597	NaN	0.014
654	4.131	NaN	0.032
654	4.562	NaN	0.056
654	3.78	NaN	0.015
654	4.338	NaN	0.018
654	4.464	NaN	0.027
654	4.447	0.984	0.03
654	4.498	0.883	0.026
654	3.934	0.393	0.014
654	4.814	1.025	0.21
654	4.238	0.662	0.015
654	3.918	0.569	0.02
654	4.356	0.91	0.049
654	4.695	0.822	0.017
654	3.979	0.52	0.021
654	5.123	1.387	0.097
654	4.114	0.488	0.018
654	3.837	0.318	0.011
654	4.251	0.715	0.019
654	4.463	0.674	0.016
654	4.849	0.848	0.025
654	4.326	0.614	0.018
654	4.142	0.595	0.028
654	5.167	1.153	0.033
654	4.694	0.731	0.027
654	4.348	0.862	0.077
654	4.118	NaN	0.056
654	3.806	NaN	0.015
654	4.488	NaN	0.019
654	4.12	NaN	0.043
654	4.259	NaN	0.018
1189	4.885	0.748829	0.275445
1189	4.186	0.665021	0.208701
1189	4.634	0.994237	0.295856
1189	5.398	1.568337	0.468056
1189	5.246	1.604477	0.393071
1189	5.183	1.305273	0.386319
1189	4.522	1.159030	0.250403
1189	5.818	1.714369	0.506603
1189	5.512	1.894193	0.554741
1189	4.378	0.806979	0.23225
1189	5.313	1.273014	0.418932
1189	5.573	1.833609	0.504891
1189	5.475	1.936551	0.517255
1189	5.26	1.307569	0.38458
1189	5.246	1.621176	0.443907
1189	5.007	1.035512	0.313676
1189	5.113	1.201195	0.380353
1189	5.986	1.894551	0.524977
1189	4.794	0.868941	0.282705
1189	5.518	2.013026	0.47769
1189	5.838	2.312711	0.575377
1189	5.702	1.882206	0.550472
1189	5.059	1.599895	0.439025
1189	5.327	1.613184	0.447415
1189	5.721	1.823543	0.549468
1189	5.447	2.076962	0.478787
1189	5.354	1.609162	0.452883
1189	5.109	1.7056      0.355435
1189	5.956	2.027093	0.554408
1189	5.413	1.837877	0.495532
];
%%% to delete NaN values, in case the matlab version does not allow to deal
%%% with them
tLW19SA(any(isnan(tLW19SA), 2), :) = [];
% Lw0.tLh_19SA = 0.0215;  units.Lw0.tLh_19SA = 'cm'; label.Lw0.tLh_19SA = 'initial physical length'; %cm, length at metamorphosis 
Lw0.tLh_19SA = 1;  units.Lw0.tLh_19SA = 'cm'; label.Lw0.tLh_19SA = 'initial physical length'; %cm, length at metamorphosis 

data.tLh_19SA(:,1) = tLW19SA(:,1); data.tLh_19SA(:,2) = tLW19SA(:,2); %creat time-length data
% data.tL19SA(:,1) = [data.aj ; tLW19SA(:,1)]; data.tL19SA(:,2) = [data.Lj ; tLW19SA(:,2)]; %creat time-length data with th initial data as the metamorphosis data
units.tLh_19SA   = {'d', 'cm'};  label.tLh_19SA = {'time', 'shell height Régnier-Brisson StAnne', 'N19 at Ste Anne'};  
% temp.tL19SA    = [  365         13.712             % Fourier parameters to be used in the derivatives function, from the excel file in Enviro/
%                     -0.76009	-3.7422             % Fourier parameters estimated for Saint-Anne, because the Roscanvel temperatures are similar starts on the 1st June 2019
%                     -0.16068	0.031801
%                     -0.056534	-0.073573];  units.temp.tL19SA = 'C'; label.temp.tL19SA = 'Fourier parameters for temperature cycle for 19SA data';
bibkey.tLh_19SA = 'RegnierBrisson';

data.tWd_19SA(:,1) = tLW19SA(:,1);  data.tWd_19SA(:,2) = tLW19SA(:,3)- tLW19SA(:,4);
% data.tW19SA(:,1) = [tLW19SA(1,1) ; tLW19SA(:,1)];  data.tW19SA(:,2) = [tLW19SA(1,3) ; tLW19SA(:,3)]; %repeat the first row to have the same dimension than tL
units.tWd_19SA   = {'d', 'g'};  label.tWd_19SA = {'time', 'somatic afdw  Régnier-Brisson StAnne', 'N19 at Ste Anne'};  
% temp.tW19SA    = [  365         13.712             % Fourier parameters to be used in the derivatives function, from the excel file in Enviro/
%                     -0.76009	-3.7422             % Fourier parameters estimated for Saint-Anne, because the Roscanvel temperatures are similar
%                     -0.16068	0.031801
%                     -0.056534	-0.073573];  units.temp.tW19SA = 'C'; label.temp.tW19SA = 'Fourier parameters for temperature cycle for 19SA data';
bibkey.tWd_19SA = 'RegnierBrisson';


data.LWd_19SA(:,1) = tLW19SA(:,2); data.LWd_19SA(:,2) = tLW19SA(:,3) - tLW19SA(:,4); %create time-length data
units.LWd_19SA   = {'cm', 'g'};  label.LWd_19SA = {'shell height','afdw w/out gonad Régnier-Brisson StAnne', 'N19 at Ste Anne'};  
temp.LWd_19SA    = C2K(14) ;  units.temp.LWd_19SA = 'C'; label.temp.LWd_19SA = 'temperature';
bibkey.LWd_19SA = 'RegnierBrisson';
comment.LWd_19SA = 'temperature between 17 and 9 C - Mean temperature (monitoring)';


temp.tLh_19SA =[... % time since 24/05/2019 (d), temperature at Ste-Anne (Bay of Brest), Laure RB monitoring (°C)
149	14.8
150	14.7
151	14.6
152	14.5
153	14.5
154	14.5
155	14.5
156	14.4
157	14.2
158	14.1
159	14.1
160	14.1
161	14.2
162	14.0
163	13.9
164	13.8
165	13.6
166	13.4
167	13.3
168	13.0
169	12.8
170	12.9
171	12.8
172	12.6
173	12.5
174	12.3
175	12.2
176	12.2
177	12.1
178	12.0
179	11.9
180	11.8
194	11.1
195	11.2
196	11.6
197	11.5
198	11.6
199	11.4
200	11.3
201	11.3
202	11.2
203	11.2
204	11.1
205	11.0
206	10.9
207	10.9
208	10.8
209	10.9
210	10.9
211	10.9
212	10.8
213	10.6
214	10.6
215	10.6
216	10.8
217	10.8
218	10.8
219	10.8
220	10.8
221	10.8
222	10.9
223	10.9
224	11.0
225	10.7
226	10.5
227	10.6
228	10.8
229	10.9
230	11.0
231	10.9
232	10.9
233	10.9
234	10.8
235	10.9
236	11.0
237	10.9
238	10.8
239	10.6
240	10.3
241	10.2
242	10.3
243	10.2
244	10.1
245	10.0
246	10.0
247	10.0
248	10.1
249	10.0
250	10.0
251	10.0
252	10.1
253	10.2
254	10.3
255	10.4
256	10.3
257	10.2
258	9.9
259	10.0
260	10.2
261	10.3
262	10.3
263	10.2
264	10.1
265	10.1
266	10.1
267	10.2
268	10.4
269	10.3
270	10.3
271	10.3
272	10.2
273	10.2
274	10.2
275	10.3
276	10.4
277	10.4
278	10.3
279	10.3
280	10.3
281	10.3
282	10.2
283	10.2
284	10.1
285	9.9
286	10.1
287	10.1
288	10.2
289	10.3
290	10.3
291	10.5
292	10.6
293	10.6
294	10.6
295	10.6
296	10.6
297	10.6
298	10.6
299	10.7
300	10.8
301	10.9
302	10.9
303	10.8
304	10.8
305	10.9
306	11.0
307	11.0
308	11.0
309	10.9
310	10.7
311	10.5
312	10.4
313	10.3
314	10.3
315	10.4
316	10.5
317	10.8
318	11.1
319	11.1
320	11.3
321	11.5
322	11.7
323	11.9
324	12.0
325	12.1
326	12.0
327	12.2
328	12.4
329	12.6
330	12.7
331	12.8
332	12.8
333	12.8
334	12.9
335	13.0
336	13.2
337	13.2
338	13.3
339	13.4
340	13.5
341	13.4
342	13.3
343	13.3
344	13.4
345	13.5
346	13.7
347	13.8
348	14.0
349	14.1
350	14.3
351	14.4
352	14.6
353	14.3
354	14.1
355	14.0
356	13.8
357	13.8
358	13.9
359	14.0
360	14.1
361	14.3
362	14.5
363	14.7
364	14.8
365	14.8
366	15.0
367	15.1
368	15.2
369	15.4
370	15.6
371	15.9
372	16.2
373	16.5
374	16.6
375	16.5
376	16.3
377	16.1
378	15.8
379	15.7
380	15.5
381	15.5
382	15.5
383	15.5
384	15.5
385	15.5
386	15.7
387	15.8
388	15.8
389	15.9
390	15.9
391	15.9
392	15.9
393	16.0
394	16.1
395	16.3
396	16.6
397	16.8
398	17.0
399	16.9
400	16.9
401	16.9
402	16.9
403	16.9
404	17.0
405	16.9
406	16.9
407	16.9
408	16.9
409	16.9
410	17.0
411	17.0
412	17.0
413	17.0
414	17.1
415	17.2
416	17.4
417	17.3
418	17.1
419	17.1
420	17.2
421	17.2
422	17.1
423	17.0
424	17.0
425	17.0
426	17.0
427	17.0
428	17.0
429	17.0
430	17.0
431	17.1
432	17.2
433	17.4
434	17.4
435	17.3
436	17.3
437	17.3
438	17.3
439	17.4
440	17.5
441	17.5
442	17.6
443	17.6
444	17.7
445	17.9
446	18.0
447	17.9
448	17.8
449	17.9
450	17.8
451	17.8
452	17.9
453	18.0
454	18.0
455	18.0
456	17.9
457	17.9
458	17.8
459	17.7
460	17.8
461	17.8
462	17.8
463	17.7
464	17.7
465	17.6
466	17.5
467	17.5
468	17.5
469	17.5
470	17.4
471	17.3
472	17.2
473	17.2
474	17.4
475	17.5
476	17.4
477	17.4
478	17.5
479	17.5
480	17.5
481	17.4
482	17.4
483	17.3
484	17.2
485	17.1
486	17.0
487	17.0
488	16.9
489	16.7
490	16.5
491	16.4
492	16.1
493	16.0
494	16.0
495	15.9
496	15.7
497	15.5
498	15.3
499	15.2
500	15.0
501	14.9
502	14.8
503	14.8
504	14.8
505	14.7
506	14.7
507	14.6
508	14.5
509	14.4
510	14.2
511	14.1
512	14.0
513	13.9
514	13.9
515	14.0
516	14.0
517	14.1
518	14.1
519	14.1
520	13.9
521	13.9
522	13.8
523	13.7
524	13.7
525	13.7
526	13.7
527	13.8
528	13.8
529	13.7
530	13.5
531	13.3
532	13.1
533	13.1
534	13.1
535	13.2
536	13.2
537	13.2
538	13.2
539	13.2
540	13.3
541	13.3
542	13.2
543	13.2
544	13.2
545	13.1
546	13.0
547	12.9
548	12.9
549	12.9
550	12.7
551	12.8
552	12.8
553	12.8
554	12.8
555	12.8
556	12.7
557	12.7
558	12.7
559	12.5
560	12.3
561	12.1
562	11.9
563	11.5
564	11.4
565	11.4
566	11.4
567	11.4
568	11.4
569	11.3
570	11.3
571	11.3
572	11.3
573	11.3
574	11.4
575	11.5
576	11.4
577	11.4
578	11.5
579	11.5
580	11.3
581	11.0
582	11.1
583	11.3
584	11.1
585	10.7
586	10.5
587	10.5
588	10.4
589	10.5
590	10.5
591	10.5
592	10.4
593	10.4
594	10.3
595	10.2
596	10.2
597	10.2
598	10.3
599	10.4
600	10.4
601	10.5
602	10.4
603	10.5
604	10.5
605	10.5
606	10.5
607	10.5
608	10.4
609	10.3
610	10.2
611	9.9
612	9.7
613	9.5
614	9.9
615	10.1
616	10.2
617	10.2
618	10.3
619	10.4
620	10.5
621	10.5
622	10.5
623	11.1
624	11.3
625	11.2
626	10.6
627	10.6
628	10.5
629	10.1
630	10.0
631	9.8
632	9.8
633	10.0
634	10.1
635	10.3
636	10.4
637	10.4
638	10.5
639	10.6
640	10.6
641	10.6
642	10.8
643	10.8
644	10.8
645	10.7
646	10.7
647	10.7
648	10.7
649	10.8
650	10.8
651	10.8
652	10.7
653	10.5
654	10.4
655	10.4
656	10.4
657	10.5
658	10.5
659	10.5
660	10.5
661	10.6
662	10.7
663	10.8
664	10.9
665	11.0
666	11.0
667	11.0
668	11.0
669	11.0
670	11.2
671	11.2
672	11.3
673	11.3
674	11.3
675	11.5
676	11.7
677	11.8
678	12.0
679	12.1
680	12.0
681	12.1
682	12.0
683	11.9
684	11.4
685	11.0
686	11.0
687	11.0
688	10.9
689	11.0
690	11.1
691	11.1
692	11.1
693	11.1
694	11.2
695	11.3
696	11.4
697	11.5
698	11.9
699	12.0
700	11.9
701	12.0
702	12.1
703	12.0
704	12.1
705	12.1
706	12.2
707	12.2
708	12.2
709	12.2
710	12.3
711	12.3
712	12.3
713	12.3
714	12.4
715	12.5
716	12.6
717	12.7
718	12.7
719	12.7
720	12.8
721	12.8
722	12.9
723	12.9
724	13.0
725	13.2
726	13.2
727	13.3
728	13.3
729	13.4
730	13.3
731	13.3
732	13.3
733	13.3
734	13.4
735	13.6
736	13.8
737	14.0
738	14.1
739	14.4
740	14.6
741	14.7
742	14.8
743	14.9
744	15.0
745	15.0
746	15.1
747	15.1
748	15.1
749	15.2
750	15.4
751	15.5
752	15.7
753	15.9
754	16.0
755	15.8
756	15.7
757	15.7
758	15.6
759	15.5
760	15.4
761	15.4
762	15.4
763	15.4
764	15.4
765	15.4
766	15.5
767	15.5
768	15.5
769	15.7
770	15.9
771	16.0
772	16.0
773	16.1
774	16.0
775	16.2
776	16.2
777	16.3
778	16.3
779	16.3
780	16.3
781	16.3
782	16.3
783	16.5
784	16.7
785	16.9
786	17.1
787	17.5
788	17.8
789	17.9
790	18.0
791	17.9
792	17.6
793	17.4
794	17.4
795	17.4
796	17.4
797	17.4
798	17.4
799	17.4
800	17.5
801	17.5
802	17.4
803	17.6
804	17.6
805	17.5
806	17.5
807	17.5
808	17.5
809	17.5
810	17.6
811	17.6
812	17.7
813	17.8
814	17.7
815	17.7
816	17.6
817	17.4
818	17.3
819	17.3
820	17.2
821	17.3
822	17.4
823	17.4
824	17.4
825	17.3
826	17.3
827	17.3
828	17.4
829	17.3
830	17.3
831	17.4
832	17.4
833	17.4
834	17.5
835	17.5
836	17.5
837	17.5
838	17.6
839	17.6
840	17.5
841	17.5
842	17.5
843	17.4
844	17.5
845	17.5
846	17.6
847	17.5
848	17.5
849	17.5
850	17.4
851	17.4
852	17.3
853	17.2
854	17.2
855	17.2
856	17.2
857	17.2
858	17.0
859	17.0
860	16.9
861	17.0
862	16.9
863	16.8
864	16.7
865	16.6
866	16.5
867	16.5
868	16.4
869	16.3
870	16.3
871	16.2
872	16.1
873	16.1
874	16.1
875	16.0
876	16.0
877	15.9
878	15.9
879	16.0
880	16.0
881	15.8
882	15.6
883	15.5
884	15.3
885	15.3
886	15.3
887	15.3
888	15.2
889	15.2
890	15.1
891	15.0
892	14.8
893	14.6
894	14.4
895	14.3
896	14.2
897	14.1
898	14.1
899	14.0
900	14.1
901	14.0
902	14.0
903	14.0
904	14.0
905	13.9
906	13.8
907	13.8
908	13.8
909	13.7
910	13.7
911	13.5
912	13.5
913	13.3
914	13.0
915	12.9
916	12.7
917	12.6
918	12.3
919	12.0
920	11.9
921	11.8
922	11.9
923	11.7
924	11.7
925	11.7
926	11.4
927	11.3
928	11.3
929	11.1
930	11.1
931	11.0
932	10.9
933	10.9
934	10.8
935	10.9
936	11.0
937	11.0
938	11.1
939	11.0
940	10.9
941	10.7
942	10.6
943	10.4
944	10.4
945	10.5
946	10.5
947	10.5
948	10.6
949	10.6
950	10.7
951	10.9
952	10.9
953	11.0
954	11.2
955	11.2
956	11.2
957	11.0
958	10.9
959	10.8
960	10.8
961	10.7
962	10.7
963	10.7
964	10.7
965	10.6
966	10.4
967	10.5
968	10.6
969	10.6
970	10.6
971	10.5
972	10.5
973	10.3
974	10.3
975	10.2
976	10.1
977	10.0
978	9.9
979	9.8
980	9.9
981	10.0
982	10.0
983	10.0
984	10.0
985	10.1
986	10.1
987	10.1
988	10.0
989	10.1
990	10.1
991	10.2
992	10.2
993	10.1
994	10.0
995	9.9
996	9.9
997	9.9
998	9.9
999	10.0
1000	10.1
1001	10.2
1002	10.0
1003	10.1
1004	10.1
1005	10.1
1006	10.1
1007	10.1
1008	10.1
1009	9.9
1010	9.9
1011	10.1
1012	10.1
1013	10.2
1014	10.3
1015	10.3
1016	10.2
1017	10.1
1018	10.0
1019	10.0
1020	10.0
1021	10.0
1022	10.0
1023	10.0
1024	10.0
1025	10.1
1026	10.2
1027	10.2
1028	10.3
1029	10.4
1030	10.4
1031	10.5
1032	10.6
1033	10.8
1034	10.9
1035	11.0
1036	11.0
1037	11.2
1038	11.2
1039	11.3
1040	11.3
1041	11.3
1042	11.3
1043	11.1
1044	11.0
1045	11.0
1046	10.9
1047	10.9
1048	11.0
1049	11.1
1050	11.1
1051	11.2
1052	11.2
1053	11.3
1054	11.5
1055	11.6
1056	11.7
1057	11.8
1058	12.0
1059	12.1
1060	12.1
1061	12.1
1062	12.1
1063	12.2
1064	12.4
1065	12.5
1066	12.7
1067	12.8
1068	12.8
1069	12.8
1070	12.8
1071	12.9
1072	12.9
1073	13.0
1074	13.0
1075	13.0
1076	13.1
1077	13.2
1078	13.4
1079	13.6
1080	13.8
1081	14.1
1082	14.2
1083	14.1
1084	14.2
1085	14.3
1086	14.4
1087	14.4
1088	14.5
1089	14.6
1090	14.6
1091	14.7
1092	14.9
1093	15.1
1094	15.1
1095	15.1
1096	15.1
1097	15.0
1098	15.0
1099	15.1
1100	15.1
1101	15.1
1102	15.2
1103	15.2
1104	15.4
1105	15.5
1106	15.7
1107	15.7
1108	15.9
1109	15.8
1110	15.8
1111	15.8
1112	16.0
1113	16.1
1114	16.2
1115	16.2
1116	16.2
1117	16.2
1118	16.2
1119	16.3
1120	16.5
1121	16.5
1122	16.0
1123	16.0
1124	16.1
1125	16.1
1126	16.2
1127	16.1
1128	16.1
1129	16.2
1130	16.3
1131	16.3
1132	16.2
1133	16.2
1134	16.3
1135	16.4
1136	16.4
1137	16.4
1138	16.5
1139	16.6
1140	16.3
1141	16.7
1142	17.0
1143	17.0
1144	17.3
1145	17.4
1146	17.5
1147	17.3
1148	17.1
1149	17.1
1150	17.2
1151	17.4
1152	17.4
1153	17.4
1154	17.4
1155	17.3
1156	17.4
1157	17.5
1158	17.4
1159	17.5
1160	17.5
1161	17.5
1162	17.5
1163	17.6
1164	17.5
1165	17.5
1166	17.5
1167	17.6
1168	17.5
1169	17.4
1170	17.3
1171	17.3
1172	17.4
1173	17.3
1174	17.4
1175	17.3
1176	17.4
1177	17.4
1178	17.2
1179	17.0
1180	16.9
1181	16.8
1182	16.8
1183	16.9
1184	16.9
1185	17.0
1186	17.1
1187	16.8
1188	17.4
1189	17.4
1190	17.2
1191	17.1
1192	17.1
1193	17.1
1194	17.1
1195	17.2	
 ];
temp.tLh_19SA    = [temp.tLh_19SA(:,1) C2K(temp.tLh_19SA(:,2))];  
units.temp.tLh_19SA = 'K'; label.temp.tLh_19SA = 'temperature';


%% set weights for all real data
weights = setweights(data, []);
% weights.tLlarvae        = 0 * weights.tLlarvae;
% weights.tLjuve          = 0 * weights.tLjuve;
weights.aj = 0 * weights.aj;
weights.Lj = 0 * weights.Lj;
% weights.E_0 = 5 * weights.E_0;
% weights.GSI = 5 * weights.GSI;
weights.R_L = 10 * weights.R_L;
% weights.Li = 5 * weights.Li;
% weights.am = 5 * weights.am;
% weights.tp = 5 * weights.tp;
% weights.Lp = 5 * weights.Lp;
% weights.Wwp = 10 * weights.Wwp;
% weights.tLh_19SA = 5 * weights.tLh_19SA;
% weights.tWd_19SA = 5 * weights.tWd_19SA;
% weights.LWd_19SA = 5* weights.LWd_19SA;
%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
% weights.psd.kap = 100 * weights.psd.kap;
%% pack auxData and txtData for output
auxData.temp = temp;        % for temperature
auxData.Lw0 = Lw0;          % for initial body size
auxData.L_R = L_R;          % for size at which reproduction data is given
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;
    
%% add ultimate length as metadata
metaData.Li = data.Li;
%% Discussion points


%% Facts
% F1 = 'Incidental hermaphroditic';
% metaData.bibkey.F1 = 'Wiki';
% F2 = 'Sports jet propulsion swimming';
% metaData.bibkey.F2 = 'Wiki';
% metaData.facts = struct('F1',F1, 'F2',F2);

%% Links
metaData.links.id_CoL = '6RKGL'; % Cat of Life
metaData.links.id_ITIS = '79628'; % ITIS
metaData.links.id_EoL = '46468063'; % Ency of Life
metaData.links.id_Wiki = 'Mimachlamys'; % Wikipedia
metaData.links.id_ADW = ''; % ADW
metaData.links.id_Taxo = '3967485'; % Taxonomicon
metaData.links.id_WoRMS = '236719'; % WoRMS
metaData.links.id_molluscabase = '236719'; % molluscabase


%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Chlamys_varia}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman. S.A.L.M.}. ' ...
'year = {2010}. ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}. ' ...
'publisher = {Cambridge Univ. Press. Cambridge}. ' ...
'pages = {Table 4.2 (page 150). 8.1 (page 300)}. ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
% bibkey = 'marlin'; type = 'Misc'; bib = ...
% 'howpublished = {\urlhttp://www.marlin.ac.uk/biotic/browse.php?sp=6255}}';
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
bibkey = 'sealifebase'; type = 'Misc'; bib = ...
'howpublished = {\url{http://www.sealifebase.org/summary/Mimachlamys-varia.html}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'TinduffHatchery'; type = 'Misc'; bib = ...
'TinduffHatchery';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'RegnierBrisson'; type = 'Misc'; bib = ...
'RegnierBrisson';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'LetaAud1956'; type = 'Article'; bib = [ ... 
'author = {Letaconnoux, R. and Audouin, J.}, ' ... 
'year = {1956}, ' ...
'title = {Contribution a l étude du pétoncle \emph{Chlamys varia} ({L}amarck}},' ...
'journal = {Revue des Travaux de l Institut des Pêches Maritimes}, ' ...
'volume = {20}, ' ...
'pages = {133-155}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
