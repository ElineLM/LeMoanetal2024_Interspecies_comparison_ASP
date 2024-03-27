function [data, auxData, metaData, txtData, weights] = mydata_Pecten_maximus

%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Bivalvia'; 
metaData.order      = 'Ostreoida'; 
metaData.family     = 'Pectinidae';
metaData.species    = 'Pecten_maximus'; 
metaData.species_en = 'Great scallop';

metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MANE'};
metaData.ecoCode.habitat = {'0jMp', 'jiMb'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biPp'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(13); % K, body temp
metaData.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li'; 'Wdb'; 'Wdj'; 'Wdp'; 'Wdi'; 'Ri'; 'E0'}; 
metaData.data_1     = {'t-L'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011


metaData.author   = {'Romain Lavaud'};    
metaData.date_subm = [2011 05 10];              
metaData.email    = {'romain.lavaud@hotmail.fr'};            
metaData.address  = {'Institut Universitaire Europeen de la Mer, Plouzane, France'};   

% modification in 2023 by Eline Le Moan
metaData.author   = {'Eline Le Moan'};    
metaData.date_subm = [2024 04 01];              
metaData.email    = {'eline.lemoan@univ-brest.fr'};            
metaData.address  = {'Institut Universitaire Europeen de la Mer, Plouzane, France'};

metaData.author_mod_1   = {'Bas Kooijman'};    
metaData.date_mod_1 = [2012 09 24];              
metaData.email_mod_1    = {'bas.kooijman@vu.nl'};            
metaData.address_mod_1  = {'VU University Amsterdam'};   

metaData.author_mod_2   = {'Starrlight Augustine'};    
metaData.date_mod_2 = [2017 05 24];              
metaData.email_mod_2    = {'sta@akvaplan.niva.no'};            
metaData.address_mod_2  = {'Akvaplan-niva'};   


metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2017 05 24]; 

%% set data
% -------------------------------------------------------------------------
% zero-variate data
% -------------------------------------------------------------------------

data.ab = 2;       units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = {'ShumPars2006','GrufBeau1972', 'TinduffHatchery'};   
  temp.ab = C2K(18);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.tj = 25;      units.tj = 'd';    label.tj = 'time since birth at metam'; bibkey.tj = {'RicoBern2010','ShumPars2006', 'TinduffHatchery'};   
  temp.tj = C2K(18);  units.temp.tj = 'K'; label.temp.tj = 'temperature'; %10/01/2023 - modif from 13°C to 18°C because value given by Tinduff hatchery
data.tp = 365;     units.tp = 'd';    label.tp = 'time since birth at puberty'; bibkey.tp = {'Lava2011','ShumPars2006', 'TinduffHatchery'};
  temp.tp = C2K(13);  units.temp.tp = 'K'; label.temp.tp = 'temperature';
data.am = 20*365;   units.am = 'd';    label.am = 'life span';                bibkey.am = 'ShumPars2006';   
  temp.am = C2K(13);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
  
data.Lb  = 0.008;  units.Lb  = 'cm';  label.Lb  = 'shell height at birth';   bibkey.Lb  = {'ShumPars2006','SamaCoch1987', 'TinduffHatchery'};
data.Lj  = 0.024;  units.Lj  = 'cm';  label.Lj  = 'shell height at metam';   bibkey.Lj  = {'BuesCoch1982','ShumPars2006', 'TinduffHatchery'};
data.Lp  = 4.0;    units.Lp  = 'cm';  label.Lp  = 'shell height at puberty'; bibkey.Lp  = {'Lava2011','ShumPars2006', 'TinduffHatchery'};
data.Li  = 12.0;   units.Li  = 'cm';  label.Li  = 'ultimate shell height';   bibkey.Li  = {'Jean2011', 'Minchin2003'};

data.Wdb = 1.04e-7;  units.Wdb = 'g'; label.Wdb = 'dry weight at birth';     bibkey.Wdb = 'Lava2011';
  comment.Wdb = 'computed from length';
data.Wdj = 3e-6;   units.Wdj = 'g'; label.Wdj = 'dry weight at metam';     bibkey.Wdj = {'Ande2011','SamaCoch1987'};
data.Wdp = 1;     units.Wdp = 'g'; label.Wdp = 'dry weight at puberty';   bibkey.Wdp = 'Lava2011';
  comment.Wdp = 'obtained from Lp via LW-regression';
data.Wdi = 20;    units.Wdi = 'g'; label.Wdi = 'ultimate dry weight';     bibkey.Wdi = 'Lava2011';
  comment.Wdi = 'obtained from Li via LW-regression';

% data.Ri  = 57534; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = {'PaulFifa1989','PaulBekh1997'};   
%   temp.Ri = C2K(13); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
%   comment.Ri = '21 million per spawning event from Paulet and Fifas 1989, 21,000,000/365 = 57,534 eggs per day';

%   data.Ri  = 30904*2; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = {'CochDeva1993', 'Paulet1998'};   
%   temp.Ri = C2K(15); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
%   comment.Ri = '5.64 million eggs per spawning, and 2 spawning events (Paulet et al., 1988)';
  
   data.R_L  = 30904*2; units.R_L  = '#/d'; label.R_L  = 'reprod rate at a given size';     bibkey.R_L  = {'CochDeva1993', 'Paulet1988'}; 
   temp.R_L = C2K(15); units.temp.R_L = 'K'; label.temp.R_L = 'temperature';
   comment.R_L = '5.64 million eggs per spawning, and 2 spawning events (Paulet et al., 1988)';
L_R.R_L = 9;  units.L_R.R_L = 'cm';         label.L_R.R_L = 'length at reproduction data'; 

data.E_0 = 4.02e-3; units.E_0  = 'J';  label.E_0  = 'egg energy content'; bibkey.E_0  = 'calculation';
  comment.E_0 = 'calculated from egg diameter (65 µm), which gives egg dry weight from M. varia and P. magellanicus data, and times 26kJ/gdw from P. magellanicus litterature';  

% -------------------------------------------------------------------------
% uni-variate data
% -------------------------------------------------------------------------

% % celtic sea, Romain Lavauds thesis, unpublished work
% data.tL1 = [ ... time (julian day), shell height (cm)
% 1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192;
% 2.003	2.007	2.011	2.016	2.020	2.026	2.033	2.037	2.043	2.049	2.056	2.061	2.067	2.075	2.082	2.091	2.100	2.113	2.125	2.137	2.150	2.160	2.167	2.176	2.187	2.200	2.212	2.224	2.231	2.239	2.245	2.254	2.265	2.273	2.279	2.286	2.292	2.298	2.303	2.309	2.318	2.325	2.333	2.343	2.351	2.363	2.375	2.387	2.398	2.411	2.423	2.436	2.450	2.463	2.475	2.487	2.499	2.512	2.527	2.542	2.560	2.576	2.594	2.610	2.627	2.644	2.664	2.681	2.699	2.714	2.728	2.747	2.764	2.780	2.798	2.814	2.832	2.853	2.872	2.889	2.910	2.930	2.947	2.968	2.988	3.007	3.024	3.041	3.061	3.078	3.097	3.119	3.143	3.162	3.179	3.194	3.215	3.235	3.256	3.276	3.297	3.321	3.342	3.358	3.374	3.394	3.417	3.439	3.458	3.477	3.496	3.518	3.541	3.563	3.584	3.606	3.629	3.650	3.668	3.689	3.708	3.725	3.743	3.760	3.774	3.786	3.802	3.814	3.830	3.842	3.855	3.867	3.878	3.889	3.901	3.910	3.922	3.933	3.943	3.956	3.967	3.980	3.992	4.003	4.014	4.027	4.042	4.058	4.072	4.087	4.098	4.106	4.114	4.123	4.134	4.143	4.158	4.169	4.180	4.190	4.202	4.214	4.227	4.235	4.243	4.249	4.260	4.269	4.277	4.287	4.297	4.304	4.311	4.320	4.329	4.341	4.351	4.362	4.372	4.379	4.385	4.392	4.405	4.411	4.418	4.425	4.433	4.443	4.456	4.466	4.474	4.479
% ]';
% data.tL1(:,1) = data.tL1(:,1); %  julian d to day since birth (but here no change!, why?)
% units.tL1   = {'d', 'cm'};  label.tL1 = {'time in julian day', 'shell height'};  
% temp.tL1    = C2K(12);  units.temp.tL1 = 'K'; label.temp.tL1 = 'temperature';
% Lw0.tL1     = 2.0; units.Lw0.tL1 = 'cm'; label.Lw0.tL1 = 'initial shell height';
% bibkey.tL1 = 'Lava2014';
% comment.tL1 = 'Celtic Sea 162 m - 1 year old, data from chap 6, temp is constant 12 deg C';

% Bay of Brest, Romain lauvaud's thesis
data.tLl_B= [103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303	304	305	306	307	308	309	310	311	312	313	314	315	316	317	318	319	320	321	322	323	324	325	326	327	328	329	330	331	332	333	334	335	336	337	338	339	340	341	342	343	344	345	346	347	348	349	350	351	352	353	354	355	356	357	358	359	360	361	362	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377	;																
3.607291	3.614155	3.621019	3.629729	3.637359	3.644477	3.652045	3.659001	3.664984	3.671305	3.677356	3.684178	3.691358	3.698385	3.706307	3.714298	3.722285	3.730339	3.740105	3.749278	3.758389	3.769195	3.778855	3.789352	3.800261	3.811397	3.822441	3.834049	3.845831	3.857219	3.869449	3.881794	3.894279	3.907035	3.920014	3.933547	3.946987	3.96101	3.97686	3.992699	4.009179	4.026606	4.044062	4.061309	4.079422	4.096921	4.114505	4.132472	4.149757	4.167662	4.185645	4.202864	4.220328	4.237861	4.254643	4.271397	4.289299	4.305843	4.323343	4.340485	4.358224	4.375413	4.393086	4.409267	4.426178	4.443144	4.459508	4.475094	4.489946	4.502323	4.514823	4.526225	4.536816	4.546764	4.556341	4.565666	4.574469	4.584289	4.593757	4.60448	4.615949	4.628371	4.64282	4.658992	4.675947	4.693839	4.711962	4.729981	4.748585	4.768249	4.788452	4.810196	4.831546	4.852992	4.875532	4.897281	4.919942	4.942086	4.965308	4.987986	5.010585	5.032892	5.054261	5.074548	5.095961	5.117288	5.138786	5.160043	5.180944	5.200722	5.221639	5.243529	5.264942	5.286674	5.307934	5.329692	5.351984	5.37343	5.395005	5.415007	5.435544	5.455825	5.476259	5.496323	5.517239	5.53683	5.556091	5.575423	5.595651	5.615968	5.636522	5.656939	5.677094	5.697319	5.716861	5.737489	5.756501	5.775869	5.795333	5.814655	5.834099	5.853117	5.872253	5.891247	5.910817	5.92824	5.946134	5.962933	5.980002	5.99763	6.01536	6.032892	6.049862	6.066584	6.084084	6.102319	6.119525	6.135028	6.151359	6.168001	6.184146	6.20068	6.216895	6.233638	6.249801	6.266877	6.284006	6.30093	6.318658	6.335684	6.352954	6.370009	6.387273	6.403545	6.419854	6.437174	6.454471	6.471438	6.488035	6.505388	6.522173	6.53929	6.556225	6.572218	6.588483	6.605028	6.621342	6.637917	6.653677	6.669081	6.684106	6.699066	6.713975	6.729927	6.745311	6.76203	6.777183	6.793142	6.809325	6.824351	6.839629	6.854289	6.868897	6.882961	6.897142	6.911495	6.925906	6.939987	6.953404	6.967087	6.981072	6.993846	7.006292	7.018342	7.030227	7.042387	7.053805	7.064569	7.075795	7.087585	7.099035	7.109643	7.120976	7.131914	7.142266	7.152125	7.161952	7.172019	7.181601	7.190911	7.200394	7.209591	7.218632	7.227595	7.236357	7.236357	7.25342	7.261391	7.269061	7.27788	7.286117	7.294031	7.301878	7.308946	7.315947	7.323503	7.331503	7.338223	7.345657	7.353363	7.35994	7.368211	7.376068	7.383311	7.390045	7.396075	7.402188	7.407384	7.413575	7.420031	7.426176	7.43282	7.439027	7.444515	7.449325	7.454559	7.458324	7.462542	7.467508	7.471855	7.477205	7.48145	7.484113	7.491567	7.495531																	
]';
units.tLl_B  = {'d', 'cm'};  label.tLl_B = {'time', 'length Lavaud, Brest'};  
Lw0.tLl_B = 3.6; units.Lw0.tLl_B = 'cm'; label.Lw0.tLl_B = 'initial physical length'; 
bibkey.tLl_B = 'Lava2014'; comment.tLl_B = 'Bay of Brest - 1 year old';

% Temperature Time Series bay of brest: day (Julian days), Temperature in deg C
temp.tLl_B = [ ... 
    1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303	304	305	306	307	308	309	310	311	312	313	314	315	316	317	318	319	320	321	322	323	324	325	326	327	328	329	330	331	332	333	334	335	336	337	338	339	340	341	342	343	344	345	346	347	348	349	350	351	352	353	354	355	356	357	358	359	360	361	362	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377;
10.64	10.58	10.52	10.46	10.4	10.34	10.28	10.22	10.16	10.1	10.04	9.98	9.92	9.86	9.8	9.74	9.68	9.62	9.56	9.5	9.44	9.38	9.32	9.25	9.19	9.13	9.07	9.01	8.95	8.89	8.83	8.77	8.71	8.65	8.59	8.53	8.47	8.52	8.56	8.61	8.66	8.71	8.75	8.8	8.85	8.89	8.94	8.99	9.04	9.08	9.13	9.16	9.18	9.21	9.24	9.27	9.29	9.32	9.36	9.4	9.44	9.48	9.52	9.56	9.6	9.64	9.68	9.73	9.77	9.81	9.85	9.89	9.93	9.97	10.01	10.05	10.09	10.04	10.11	10.17	10.24	10.3	10.37	10.44	10.5	10.57	10.64	10.7	10.77	10.83	10.9	10.94	10.97	11.01	11.05	11.08	11.12	11.19	11.14	11.1	11	10.89	10.79	10.69	10.64	10.59	10.54	10.59	10.63	10.68	10.72	10.77	10.9	11.03	11.16	11.29	11.42	11.54	11.65	11.65	11.65	11.65	11.65	11.65	11.74	11.82	11.91	12.1	12.29	12.49	12.68	12.93	13.1	13.28	13.45	13.62	13.8	13.97	14.15	14.32	14.4	14.48	14.56	14.59	14.63	14.66	14.61	14.55	14.46	14.52	14.58	14.64	14.7	14.84	14.98	15.11	15.18	15.25	15.31	15.38	15.45	15.52	15.38	15.39	15.4	15.4	15.41	15.42	15.43	15.56	15.69	15.82	15.96	16.09	16.22	16.35	16.32	16.29	16.25	16.22	16.19	16.16	16.16	16.16	16.16	16.15	16.15	16.15	16.15	16.15	16.15	16.14	16.14	16.14	16.14	16.14	16.14	16.13	16.13	16.13	16.13	16.14	16.16	16.17	16.19	16.2	16.22	16.23	16.25	16.26	16.25	16.24	16.23	16.21	16.2	16.19	16.18	16.32	16.46	16.6	16.74	16.88	17.02	16.93	16.83	16.74	16.64	16.54	16.43	16.33	16.34	16.35	16.36	16.38	16.39	16.4	16.41	16.42	16.4	16.39	16.37	16.36	16.34	16.33	16.31	16.3	16.28	16.27	16.25	16.24	16.22	16.21	16.19	16.18	16.16	16.15	16.13	16.12	16.1	15.89	15.89	15.9	15.9	15.91	15.9	15.89	15.88	15.86	15.85	15.84	15.83	15.69	15.55	15.41	15.27	15.13	14.99	14.85	14.71	14.56	14.42	14.28	14.14	14	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.86	13.78	13.7	13.62	13.54	13.45	13.37	13.29	13.21	13.13	13.04	12.96	12.88	12.76	12.65	12.53	12.42	12.3	12.19	12.07	11.96	11.84	11.72	11.61	11.49	11.38	11.26	11.15	11.03	10.91	10.92	10.93	10.93	10.94	10.94	10.95	10.96	10.96	10.97	10.97	10.98	10.98	10.99	11	10.95	10.91	10.86	10.82	10.77	10.73	10.68	10.64	10.59	10.55	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.51	10.46	10.42	10.37	10.33	10.28
]';
temp.tLl_B    = [temp.tLl_B(:,1) C2K(temp.tLl_B(:,2))];  
units.temp.tLl_B = 'K'; label.temp.tLl_B = 'temperature';

% % Traena Date (Julian day)Length (cm)	
% data.tL3 = [ ...
%     151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303	304	305	306	307	308	309	310	311	312	313	314	315	316	317	318	319	320	321	322	323	324	325	326	327	328	;																																																																																																														
% 	1.815	1.828	1.841	1.859	1.869	1.880	1.894	1.910	1.928	1.944	1.961	1.981	2.000	2.016	2.030	2.041	2.054	2.067	2.084	2.102	2.118	2.136	2.153	2.169	2.183	2.197	2.214	2.228	2.242	2.260	2.282	2.302	2.320	2.340	2.358	2.375	2.395	2.415	2.434	2.450	2.465	2.484	2.504	2.524	2.544	2.565	2.585	2.606	2.627	2.648	2.669	2.692	2.715	2.738	2.761	2.784	2.807	2.827	2.849	2.874	2.897	2.919	2.941	2.964	2.987	3.010	3.034	3.055	3.079	3.101	3.125	3.149	3.171	3.195	3.218	3.242	3.265	3.289	3.312	3.336	3.359	3.381	3.405	3.429	3.452	3.475	3.497	3.520	3.544	3.565	3.586	3.608	3.630	3.650	3.671	3.693	3.715	3.736	3.759	3.783	3.806	3.827	3.848	3.870	3.891	3.912	3.932	3.953	3.976	3.996	4.017	4.037	4.059	4.081	4.103	4.123	4.143	4.164	4.186	4.208	4.230	4.251	4.272	4.293	4.317	4.338	4.358	4.380	4.402	4.423	4.445	4.466	4.489	4.513	4.535	4.557	4.579	4.601	4.623	4.646	4.669	4.690	4.711	4.731	4.751	4.771	4.791	4.811	4.831	4.851	4.871	4.891	4.908	4.926	4.945	4.961	4.978	4.995	5.012	5.029	5.045	5.060	5.073	5.085	5.094	5.103	5.112	5.121	5.129	5.136	5.143	5.149	5.156	5.164	5.172	5.179	5.186	5.190																																																																																																															
% ]';																																																																																																																																																																																																																																																																																																	
% units.tL3  = {'d', 'cm'};  label.tL3 = {'time - julian day', 'length (Traena)'};  
% Lw0.tL3 = 1.8; units.Lw0.tL3 = 'cm'; label.Lw0.tL3 = 'initial physical length'; 
% bibkey.tL3 = 'Lava2014'; comment.tL3 = 'Traena Norway - 3 year old';
% 																																																																																																																																																																																																																																																																																																	
% % Temperature Time Series Traena Norway
% 	temp.tL3 = [ ...  Temperature deg C, 	Julian Day																																																																																																																																																																																																																																																																																																
% 	6.6	6.6	6.4	6.4	6.4	6.2	6.2	5.8	5	5	5.8	5.8	5.6	5.8	5	4.8	4.8	4.8	4.8	6	6.2	6	5.4	5.2	4.2	4.2	4.6	4.8	5	5.2	4.8	5	5	4.8	4.6	4.8	4.6	4.4	4.2	4.2	4.4	4.2	4	3.9	3.6	4	4.4	4.6	4.6	4.4	4.4	4.2	4.2	4.4	4.5	4.6	4.6	4.7	4.8	4.2	4.6	5	5	5	5	5	4.8	4.8	4.8	4.8	5	5	5.2	5.2	5.2	5.4	5.6	5.6	5.8	5.8	5.8	5.8	5.6	5.6	5.8	5.8	6	6	6	6.2	6.2	6.4	6.4	6.4	6.4	6.6	6.6	6.6	6.6	7.2	7.2	7.4	7.6	7.8	7.8	8	8	8	8.2	8.2	8.2	8.2	9	9	9	9.6	9.6	9.4	10	10	10.4	10.6	10.4	9.8	9.8	9.4	9.4	10	10	10.6	12	11.6	11.8	11.6	11.6	11.6	11.6	11.6	11.6	11.6	12.4	12.4	12.4	12.4	12.6	12.6	12.8	13.2	13.4	13.4	13.4	13	12.8	12.8	13	13	12.4	12.4	12.8	12.6	12.4	12.4	12	12.2	12.2	12	12	12.2	12	12	12	12.2	12.4	12.2	12.2	12	11.8	11.8	11.8	12	12	12	11.8	11.8	11.8	11.6	11.8	11.4	11.6	11.8	12	11.8	11.8	11.8	11.8	11.7	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.8	11.6	11.4	11.2	11	11	11	11	11	10.8	10.6	10.8	11.8	11.6	11.2	11.4	11.2	11.2	11.2	11.2	11	11	11	11	11	11.2	11	10.8	10.8	11	11	10.6	10.6	10.2	10.2	10	10.2	10.2	9.8	9.8	9.8	9.8	9	9.3	8.8	8.8	8.6	8.6	8.5	8.8	9.2	9.2	9	8.8	8	8	8	7.8	6.8	6.8	6.6	6.8	6	6.4	6	6.4	6.4	6	6	5.2	5.5	5.5	5.5	5.5	6	5.8	6.6	5.8	5.9	6.2;
% 1	2	4	5	6	7	8	11	14	15	19	20	21	22	25	26	27	28	29	32	33	35	36	37	39	40	41	43	49	51	53	54	55	56	60	61	62	63	64	65	67	68	69	70	71	74	76	77	78	79	81	82	83	86	88	89	90	91	92	93	95	96	97	98	99	100	102	103	104	105	106	107	109	110	111	112	113	114	116	117	118	119	120	121	123	124	125	126	127	128	130	131	132	133	134	135	138	139	140	141	142	145	146	147	148	149	151	152	154	155	156	157	159	160	161	162	163	164	166	167	168	169	170	171	173	174	175	176	177	178	180	181	182	183	184	186	187	188	189	190	191	193	194	195	196	197	200	201	202	203	205	206	207	208	210	211	213	214	215	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	235	236	237	238	239	240	241	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	302	305	308	309	310	311	312	313	314	315	316	318	319	321	322	323	324	326	327	328	329	330	332	335	336	337	338	340	341	343	344	346	347	348	349	350	351	352	354	355	356	357	358	360	361	362	363	364	365
% ]';
% temp.tL3    = [temp.tL3(:,2) C2K(temp.tL3(:,1))];   % julian day first column and then temp converted to kelvin in the second column
% units.temp.tL3 = 'K'; label.temp.tL3 = 'temperature';

% t-L data from Bay of Brest by Fred Jean -- EVECOS data
data.tLl_EV = [ ... % time since birth (a), shell height (cm)
1	1.8
1	1.8
1	1.9
1	2
1	2
1	2
1	2
1	2.1
1	2.1
1	2.1
1	2.2
1	2.2
1	2.3
1	2.3
1	2.4
1	2.4
1	2.4
1	2.4
1	2.4
1	2.4
1	2.5
1	2.5
1	2.5
1	2.5
1	2.5
1	2.5
1	2.5
1	2.5
1	2.55
1	2.6
1	2.6
1	2.6
1	2.6
1	2.6
1	2.6
1	2.7
1	2.7
1	2.7
1	2.7
1	2.8
1	2.8
1	2.8
1	2.8
1	2.8
1	2.8
1	2.8
1	2.85
1	2.9
1	2.9
1	2.9
1	2.9
1	2.9
1	2.9
1	2.9
1	2.9
1	2.9
1	3
1	3
1	3
1	3
1	3
1	3
1	3.1
1	3.1
1	3.1
1	3.3
1	3.4
1	4
2	3.25
2	3.4
2	3.6
2	3.7
2	3.8
2	4
2	4
2	4
2	4
2	4
2	4.2
2	4.4
2	4.4
2	4.5
2	4.5
2	4.6
2	4.9
2	4.9
2	5
2	5.3
2	5.3
2	5.6
2	5.7
2	5.7
2	5.8
2	5.8
2	5.8
2	5.9
2	6
2	6.1
2	6.1
2	6.1
2	6.1
2	6.2
2	6.2
2	6.2
2	6.2
2	6.2
2	6.2
2	6.3
2	6.3
2	6.3
2	6.3
2	6.3
2	6.3
2	6.3
2	6.4
2	6.4
2	6.5
2	6.5
2	6.5
2	6.5
2	6.5
2	6.5
2	6.5
2	6.6
2	6.6
2	6.7
2	6.8
2	6.9
2	6.9
2	6.9
2	7
2	7.1
2	7.1
2	7.4
2	7.6
3	6
3	6
3	6
3	6.2
3	6.5
3	6.5
3	6.5
3	6.5
3	6.6
3	6.65
3	6.8
3	7
3	7
3	7
3	7.1
3	7.1
3	7.1
3	7.1
3	7.2
3	7.2
3	7.3
3	7.3
3	7.4
3	7.4
3	7.4
3	7.5
3	7.5
3	7.5
3	7.7
3	7.7
3	7.7
3	7.8
3	7.8
3	7.9
3	8
3	8
3	8.1
3	8.2
3	8.2
3	8.2
3	8.3
3	8.3
3	8.4
3	8.4
3	8.5
3	8.5
3	8.5
3	8.5
3	8.5
3	8.5
3	8.5
3	8.5
3	8.6
3	8.6
3	8.7
3	8.7
3	8.7
3	8.7
3	8.8
3	8.9
3	9.1
3	9.1
3	9.1
3	9.1
3	9.2
3	9.2
3	9.2
4	7.4
4	7.4
4	7.5
4	7.5
4	7.7
4	7.8
4	7.9
4	8.4
4	8.4
4	8.5
4	8.5
4	8.5
4	8.6
4	8.6
4	8.8
4	8.8
4	8.8
4	8.8
4	8.8
4	8.8
4	8.9
4	9
4	9
4	9
4	9
4	9
4	9
4	9
4	9
4	9
4	9.1
4	9.1
4	9.15
4	9.2
4	9.2
4	9.3
4	9.3
4	9.4
4	9.4
4	9.4
4	9.4
4	9.4
4	9.4
4	9.5
4	9.6
4	9.8
4	9.9
4	9.9
4	10
4	10.1
4	10.2
4	10.2
4	10.5
5	9.2
5	9.1
5	9.1
5	9.5
5	9.5
5	9.5
5	9.6
5	9.9
5	10
5	10
5	10
5	10.3
1	0.9
1	2.4
1	4.9
2	4.3
2	6.2
2	8.3
3	6.2
3	8.3
3	9.6
4	8
4	9.5
4	11.1
5	8.9
5	10.1
5	11.7
6	9.2
6	10.6
6	11.9
7	9.8
7	10.8
7	12] ;
data.tLl_EV(:,1) = (data.tLl_EV(:,1) * 365) - (365/2); % convert yr to d
units.tLl_EV   = {'d', 'cm'};  label.tLl_EV = {'time since birth', 'shell height EVECOS'};  
% temp.tL    = C2K(16);  units.temp.tL = 'K'; label.temp.tL = 'temperature';
temp.tLl_EV    = C2K(14);  units.temp.tLl_EV = 'K'; label.temp.tLl_EV = 'temperature'; % more 14°C in Bay of Brest in average
bibkey.tLl_EV = 'Jean2011';
  
% % t-Wd data from Bay of Brest -- EVECOS data -- retrieved from graphs on
% % Lavaud's thesis
% data.tWd = [ ... % time since january, somatic dry weight (g)
%     0	6.4
%     22	6.9
%     25	6.5
%     32	5.9
%     49	7.1
%     54	5.9
%     57	6.8
%     62	7.9
%     70	6.6
%     75	6.6
%     77	6.8
%     89	6.3
%     90	7.0
%     91	7.2
%     95	6.4
%     97	6.9
%     103	5.8
%     109	7.7
%     110	6.0
%     116	5.9
%     118	6.5
%     125	7.5
%     126	5.7
%     133	6.5
%     137	7.2
%     139	6.1
%     143	6.4
%     144	7.6
%     147	7.7
%     154	8.2
%     158	6.9
%     160	7.3
%     162	8.3
%     165	7.4
%     175	8.2
%     177	8.3
%     180	6.6
%     188	8.6
%     189	7.2
%     202	8.0
%     205	8.7
%     217	9.1
%     223	9.6
%     230	8.7
%     236	11.1
%     237	8.4
%     238	9.8
%     254	9.3
%     256	10.9
%     267	9.3
%     284	9.3
%     287	11.0
%     299	8.5
%     304	10.6
%     314	10.0
%     339	9.5
%     340	9.2
%     349	10.9    ];
% units.tWd   = {'d', 'g'};  label.tWd = {'time since 01/01', 'somatic dry weight'};  
% temp.tWd    = C2K(16);  units.temp.tWd = 'K'; label.temp.tWd = 'temperature';
% Lw0.tWd = 8.2; units.Lw0.tWd = 'cm'; label.Lw0.tWd = 'initial physical length in 1998'; 
% bibkey.tWd = 'Jean2011'; comment.tWd = "average of 20 individuals per point, grouped 1998, 1999 and 2000 by month, each time zero as the 1st of January of each year";


% %%% other data, from literature
% dataSet_Galicia = [ ... % shell length (cm), shell height (cm), shell dry weight (g), somatic DW (g), gonad DW (g)
%     %%% mean of several individuals, of class 1 for the age
%     %%% monitoring from April 1990 to July 1991. With one value
%     %%% approximately for each month
%     %%% the standard deviation values are also available (table below)
%     9.46	8.58	66.793	4.467	0.801
%     9.81	9.01	76.586	5.963	1.443
%     10.2	9.2     76.499	6.674	2.397
%     10.18	9.29	77.603	7.156	1.182
%     10.3	9.5     81.285	7.783	1.2
%     10.49	9.57	84.141	8.769	1.023
%     10.9	9.9     96.31	10.127	1.216
%     11.19	10.25	101.792	10.667	1.266
%     11.66	10.63	119.469	12.305	1.402
%     11.61	10.62	121.345	11.798	1.338
%     11.47	10.58	121.421	11.661	2.237
%     11.9	10.8	128.865	9.974	2.719
%     11.7	10.67	136.913	8.884	3.538
%     11.41	10.6	112.236	7.632	1.945
%     11.18	10.21	114.965	7.734	2.031
%     11.59	10.56	123.514	8.861	1.996
%     12.52	11.37	146.349	11.593	4.963
%     12.75	11.57	149.943	13.443	2.237
% ] ;   
% data.LW_Galicia(:,1) = dataSet_Galicia(:,2);
% data.LW_Galicia(:,2) = dataSet_Galicia(:,4);
% units.LW_Galicia   = {'cm', 'g'};  
% label.LW_Galicia = {'shell height', 'somatic dry weight'}; 
% temp.LW_Galicia    = C2K(15);  units.temp.LW_Galicia = 'K'; label.temp.LW_Galicia = 'temperature'; %mean of 17°C in summer and 13°C in winter
% bibkey.LW_Galicia = 'Pazos1997';
% comment.LW_Galicia = 'Monitoring every month between April 1990 to July 1991, class I age';

% 
% %time-length larvae
% data.tLlarvae = [ ... % time since spawning event (d). shell length (cm)
% % data from Tinduff hatchery
%     3	0.0111
% 3	0.0104
% 3	0.0103
% 3	0.0109
% 3	0.0111
% 3	0.011
% 3	0.0106
% 3	0.0108
% 3	0.0111
% 3	0.0112
% 3	0.0101
% 3	0.0104
% 3	0.01
% 3	0.0107
% 3	0.0102
% 3	0.0115
% 3	0.011
% 3	0.0111
% 3	0.0099
% 3	0.0115
% 3	0.0107
% 3	0.0112
% 3	0.0108
% 3	0.011
% 3	0.0107
% 3	0.0098
% 3	0.0112
% 3	0.0111
% 3	0.0112
% 3	0.0114
% 9	0.0139
% 9	0.0149
% 9	0.0137
% 9	0.014
% 9	0.0137
% 9	0.0136
% 9	0.014
% 9	0.0145
% 9	0.014
% 9	0.0148
% 9	0.0144
% 9	0.014
% 9	0.0135
% 9	0.0136
% 9	0.0143
% 9	0.0135
% 9	0.0144
% 9	0.0135
% 9	0.0129
% 9	0.0133
% 9	0.0136
% 9	0.0142
% 9	0.0145
% 9	0.014
% 9	0.0125
% 9	0.0142
% 9	0.0138
% 9	0.0143
% 9	0.0145
% 9	0.0139
% 14	0.0175
% 14	0.0168
% 14	0.0169
% 14	0.0167
% 14	0.0159
% 14	0.0176
% 14	0.0173
% 14	0.0183
% 14	0.0163
% 14	0.0178
% 14	0.0154
% 14	0.0195
% 14	0.0164
% 14	0.0211
% 14	0.0131
% 14	0.0153
% 14	0.0153
% 14	0.0171
% 14	0.0179
% 14	0.0174
% 14	0.0185
% 14	0.0163
% 14	0.0161
% 14	0.0161
% 14	0.0204
% 14	0.0195
% 14	0.0175
% 14	0.0149
% 14	0.0159
% 14	0.0181
% 16	0.0191
% 16	0.0212
% 16	0.0192
% 16	0.0202
% 16	0.0204
% 16	0.0173
% 16	0.02
% 16	0.0184
% 16	0.0206
% 16	0.0191
% 16	0.0191
% 16	0.0171
% 16	0.0199
% 16	0.0196
% 16	0.0182
% 16	0.0182
% 16	0.0212
% 16	0.0191
% 16	0.017
% 16	0.0188
% 16	0.0181
% 16	0.0178
% 16	0.0213
% 16	0.0184
% 16	0.0192
% 16	0.0189
% 16	0.0178
% 16	0.0214
% 16	0.0173
% 16	0.0193
% 20	0.022
% 20	0.0197
% 20	0.0222
% 20	0.0233
% 20	0.0218
% 20	0.019
% 20	0.0233
% 20	0.0213
% 20	0.0221
% 20	0.022
% 20	0.0179
% 20	0.0204
% 20	0.0249
% 20	0.0198
% 20	0.0195
% 20	0.0195
% 20	0.0222
% 20	0.019
% 20	0.0241
% 20	0.02
% 20	0.0202
% 20	0.0214
% 20	0.0197
% 20	0.0224
% 20	0.0197
% 20	0.0221
% 20	0.0193
% 20	0.0238
% 20	0.0239
% 20	0.0221
% 22	0.0234
% 22	0.0254
% 22	0.0271
% 22	0.026
% 22	0.0263
% 22	0.022
% 22	0.022
% 22	0.0241
% 22	0.0239
% 22	0.0238
% 22	0.0241
% 22	0.0231
% 22	0.0255
% 22	0.0229
% 22	0.0253
% 22	0.0241
% 22	0.0218
% 22	0.0215
% 22	0.0266
% 22	0.0235
% 22	0.0235
% 22	0.0221
% 22	0.0264
% 22	0.023
% 22	0.0246
% 22	0.0254
% 22	0.026
% 22	0.0218
% 22	0.0251
% 22	0.0234
% 24	0.023
% 24	0.0234
% 24	0.0262
% 24	0.0243
% 24	0.0255
% 24	0.0254
% 24	0.0264
% 24	0.0243
% 24	0.0236
% 24	0.0272
% 24	0.0252
% 24	0.0278
% 24	0.0273
% 24	0.0262
% 24	0.0259
% 24	0.0286
% 24	0.0251
% 24	0.0269
% 24	0.0259
% 24	0.0271
% 24	0.0258
% 24	0.0244
% 24	0.0253
% 24	0.0235
% 24	0.0237
% 24	0.0253
% 24	0.0256
% 24	0.0258
% 24	0.0241
% 24	0.0227];
% 
% % differenciate larvae and juvenile from one data set
% % temp_l = tLTinduff(:,2);
% % data.tLlarvae(:,1) = [data.ab ; tLTinduff(tLTinduff(:,1) < data.aj)]; %only for values under the metamorphosis age
% % data.tLlarvae(:,2) = [data.Lb ; temp_l(tLTinduff(:,1) < data.aj)];
% units.tLlarvae   = {'d', 'cm'};  label.tLlarvae = {'time', 'shell length (larvae)'};  %shell length, parallel to the hinge line, biggest length
% temp.tLlarvae    = C2K(18);  units.temp.tLlarvae = 'K'; label.temp.tLlarvae = 'temperature'; %18°C, temperature of the water at the Tinduff hatchery
% bibkey.tLlarvae = 'TinduffHatchery';
% comment.tLlarvae = 'Monitoring in the Tinduff Hatchery, controlled temperatures';



%% set weights for all real data
weights = setweights(data, []);
% weights.Lb = 2 * weights.Lb; %
% weights.ab = 5 * weights.ab;
% weights.Wdb = 5 * weights.Wdb;
weights.Lj      = 0 * weights.Lj; %
weights.tj      = 0 * weights.tj;
weights.Wdj     = 0 * weights.Wdj;
% weights.Wdp = 5 * weights.Wdp;
% weights.Li      = 5 * weights.Li;
% weights.E_0     = 5 * weights.E_0;
weights.R_L     = 10 * weights.R_L;
% weights.Wdi     = 5 * weights.Wdi;
% weights.tLl_EV    = 5 * weights.tLl_EV;
% weights.tLl_B     = 5 * weights.tLl_B; 
%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
% weights.psd.kap = 10;
% weights.psd.T_A = 5*100;
% weights.psd.kap = 100 * weights.psd.kap;
%% pack auxData and txtData for output
auxData.temp = temp;        % temperature values
auxData.Lw0 = Lw0;          % initial length for univariate data
auxData.L_R = L_R;          % length at which reproduction effort is given

txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;


%% add ultimate length as metadata
metaData.Li = data.Li; %to be able to calculate species zoom factor for body size scaling relationships

%% Links
metaData.links.id_CoL = '767GC'; % Cat of Life
metaData.links.id_ITIS = '79683'; % ITIS
metaData.links.id_EoL = '46467904'; % Ency of Life
metaData.links.id_Wiki = 'Pecten_maximus'; % Wikipedia
metaData.links.id_ADW = 'Pecten_maximus'; % ADW
metaData.links.id_Taxo = '39421'; % Taxonomicon
metaData.links.id_WoRMS = '140712'; % WoRMS
metaData.links.id_molluscabase = '140712'; % molluscabase


%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Pecten_maximus}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S. A. L. M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'ShumPars2006'; type = 'Incollection'; bib = [ ... 
'author = {Shumway, S. E. and Parsons, G. J.}, ' ... 
'year = {2006}, ' ...
'title = {Scallops: Biology, Ecology and Aquaculture.}, ' ...
'publisher = {Elsevier}, ' ...
'booktitle = {Developments in Aquaculture and Fisheries Science}, ' ...
'volume = {35}, '...
'address = {Amsterdam}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'SamaCoch1987'; type = 'Article'; bib = [ ... 
'author = {Samain J. F. and J. C. Cochard and L. Chevaulot and J. Y. Daniel and C. Jeanthon and J. R. Le Coz and Y. Marty and J. Moal and D. Prieur and M. Sala\"{u}n}, ' ... 
'year = {1987}, ' ...
'title = {Effect of water quality on the growth of \emph{Pecten maximus} larvae in nursery: {F}irst observations}, ' ...
'journal = {Haliotis}, ' ...
'volume = {16}, ' ...
'pages = {363--381}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'RicoBern2010'; type = 'Article'; bib = [ ... 
'author = {Rico-Villa, B. and Bernard, I. and Robert, R. and Pouvreau, S.}, ' ... 
'year = {2010}, ' ...
'title = {A {D}ynamic {E}nergy {B}udget ({D}{E}{B}) growth model for pacific oyster larvae, \emph{Crassostrea gigas}}, ' ...
'journal = {Aquaculture}, ' ...
'volume = {305}, ' ...
'pages = {84--94}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'PaulBekh1997'; type = 'Article'; bib = [ ... 
'author = {Paulet, Y.-M. and Bekhadra, F. and Devauchelle, N. and Donval, A. and Dorange, G.}, ' ... 
'year = {1997}, ' ...
'title = {Seasonal cycles, reproduction and oocyte quality in \emph{Pecten maximus} from the {B}ay of {B}rest}, ' ...
'journal = {Annales de l Institut oceanographique}, ' ...
'volume = {73}, ' ...
'number = {1}, '...
'pages = {101--112}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'PaulFifa1989'; type = 'Article'; bib = [ ... 
'author = {Paulet Y.-M. and Fifas, S.}, ' ... 
'year = {1989}, ' ...
'title = {Etude de la fecondite potentielle de la coquille {S}aint-{J}acques \emph{Pecten maximus}, en {B}aie de {S}aint-{B}rieue.}, ' ...
'journal = {Haliotis}, ' ...
'volume = {19}, ' ...
'pages = {275--285}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'GrufBeau1972'; type = 'Article'; bib = [ ... 
'author = {Gruffydd, L. D. and Beaumont, A. R.}, ' ... 
'year = {1972}, ' ...
'title = {A method for rearing Pecten maximus larvae in the laboratory}, ' ...
'journal = {Marine Biology}, ' ...
'volume = {15}, ' ...
'number = {4}, '...
'pages = {350--355}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'BuesCoch1982'; type = 'Article'; bib = [ ... 
'author = {Buestel, D. and Cochard, J.-C. and Dao, J.-C. and Gerard, A.}, ' ... 
'year = {1982}, ' ...
'title = {Artificial production of great scallop \emph{Pecten maximus} ({L}.) spat. {F}irst results in the {B}ay of {B}rest}, ' ...
'journal = {Vie Marine-Annales de la Fondation oceanographique Ricard}, ' ...
'volume = {4}, ' ...
'pages = {24--28}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Jean2011'; type = 'Misc'; bib = [...
'author = {Jean, f.}, ' ...     
'year = {2011}, ' ...
'note = {Own measurements in Bay of Brest}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Lava2011'; type = 'Misc'; bib = [...
'author = {Lavaud, R.}, ' ...     
'year = {2011}, ' ...
'note = {Own measurements in Bay of Brest}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Ande2011'; type = 'Misc'; bib = [...
'author = {Sissel Andersen}, ' ...     
'year = {2011}, ' ...
'note = {Pers. com}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Lava2014'; type = 'PhdThesis'; bib = [...
'author = {Romain Lavuad}, ' ...     
'year = {2014}, ' ...
'title = {Environmental variability and energetic adaptability of the great scallop, \emph{Pecten maximus}, facing climate change}, ' ...
'school = {Vrije Universiteit Amsterdam}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'CochDeva1993'; type = 'Article'; bib = [ ... 
  'title = {Spawning, Fecundity and Larval Survival and Growth in Relation to Controlled Conditioning in Native and Transplanted Populations of {{{\emph{Pecten}}}}{\emph{ Maximus}} ({{L}}.): Evidence for the Existence of Separate Stocks},'...
  'author = {Cochard, J. C. and Devauchelle, N.},'...
  'year = {1993},'...
  'journal = {Journal of Experimental Marine Biology and Ecology},'...
  'volume = {169},'...
  'number = {1},'...
  'pages = {41--56},'...
  'doi = {10.1016/0022-0981(93)90042-M},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Paulet1988'; type = 'Article'; bib = [ ... 
'author = {Paulet, Y.M. and Lucas, A. and Gerard, A.}, ' ... 
'year = {1988}, ' ...
'title = {Reproduction and Larval Development in Two {{{\emph{Pecten}}}}{\emph{ Maximus}} ({{L}}.) Populations from {{Brittany}}}, ' ...
'journal = {Journal of Experimental Marine Biology and Ecology}, ' ...
'volume = {119}, ' ...
'number = {2},'...
'pages = {145--156}',...
'doi = {10.1016/0022-0981(88)90229-8},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'TinduffHatchery'; type = 'Misc'; bib = ...
'TinduffHatchery';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Pazos1997'; type = 'Article'; bib = [ ... 
'author = {Pazos, A. J. and Román, G. and Acosta, C. P. and Abad, M. and Sánchez, J. L.}, ' ... 
'year = {1997}, ' ...
'title = {Seasonal Changes in Condition and Biochemical Composition of the Scallop {{{\emph{Pecten}}}}{\emph{ Maximus}} {{L}}. from Suspended Culture in the {{Ria}} de {{Arousa}} ({{Galicia}}, {{N}}.{{W}}. {{Spain}}) in Relation to Environmental Conditions}, ' ...
'journal = {Journal of Experimental Marine Biology and Ecology}, ' ...
'volume = {211}, ' ...
'number = {2},'...
'pages = {169--193}',...
'doi = {10.1016/S0022-0981(96)02724-4},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Minchin2003'; type = 'Article'; bib = [ ... 
'author = {Minchin, D}, ' ... 
'year = {2003}, ' ...
'title = {Introductions: Some Biological and Ecological Characteristics of Scallops}, ' ...
'journal = {Aquatic Living Resources}, ' ...
'volume = {16}, ' ...
'number = {6},'...
'pages = {521--532}',...
'doi = {10.1016/j.aquliv.2003.07.004},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%

