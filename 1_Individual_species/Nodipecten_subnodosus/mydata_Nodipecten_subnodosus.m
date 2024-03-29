function [data, auxData, metaData, txtData, weights] = mydata_Nodipecten_subnodosus


%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Bivalvia'; 
metaData.order      = 'Ostreoida'; 
metaData.family     = 'Pectinidae';
metaData.species    = 'Nodipecten_subnodosus'; 
metaData.species_en = "Lion's Paw scallop";

metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MANE'};
metaData.ecoCode.habitat = {'0jMp', 'jiMb'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biPp'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(20); % K, body temp --> No founded in litterautre, guess from P. maximus ? 
metaData.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li'; 'Wdb'; 'Wdj'; 'Wdp'; 'Wdi'; 'Ri'; 'E0'}; 
metaData.data_1     = {'t-L'; 't-Ww'; 'L-Wd'}; 

metaData.COMPLETE = 3; % using criteria of LikaKear2011


metaData.author   = {'Leo HEYER'};    
metaData.date_subm = [2023 07 10];              
metaData.email    = {'leo.heyer@etu.univ-amu.fr'};            
metaData.address  = {'Institut Universitaire Europeen de la Mer, Plouzane, France'};   

% modification summer 2023
metaData.author2   = {'Paulo Lagos'};    
metaData.date_subm2 = [2023 09 01];              
metaData.email2   = {'paulo.lagos@univ-brest.fr'};            
metaData.address2  = {'Institut Universitaire Europeen de la Mer, Plouzane, France'};   

% modifications and final editing on October 2023
metaData.author   = {'Eline Le Moan'};    
metaData.date_subm = [2024 04 01];              
metaData.email    = {'eline.lemoan@univ-brest.fr'};            
metaData.address  = {'Institut Universitaire Europeen de la Mer, Plouzane, France'};   

% metaData.author_mod_1   = {'Laure Pecquerie'};    
% metaData.date_mod_1 = [2012 09 24];              
% metaData.email_mod_1    = {'laure.pecquerie@ird.fr'};            
% metaData.address_mod_1  = {'IRD Plouzané'};   
% 
% metaData.author_mod_2   = {'Salvador Emilio Lluch Cota'};    
% metaData.date_mod_2 = [2017 05 24];              
% metaData.email_mod_2    = {''};            
% metaData.address_mod_2  = {'CIBNOR, La Paz, Mexico'};   


metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = []; 

%% set data
% -------------------------------------------------------------------------
% zero-variate data
% -------------------------------------------------------------------------

% age
data.ab = 1;       units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = {'NavaG2022','Serrano1997'};
temp.ab = C2K(26.5);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
comment.ab = "D larvae observed 20 to 24h after fertilization";

data.tj = 14;      units.tj = 'd';    label.tj = 'time since birth at metam'; bibkey.tj = 'NavaG2022';
temp.tj = C2K(26.5);  units.temp.tj = 'K'; label.temp.tj = 'temperature';
comment.tj = "experimental, not directly oberved, measure stop and begin once individual metamorphosed, metamorphosis occur during 3 to 5 days and measure stop in a same scale";

data.tp = 406;
% data.tp = 406 + (6*30);    
units.tp = 'd';    label.tp = 'time since birth at puberty'; bibkey.tp = 'Ram2009';
temp.tp = C2K(23.80);  units.temp.tp = 'K'; label.temp.tp = 'temperature';
comment.tp = ["experimental - two phase of cultivation, first in lab, with 24°C maintain until individuals reach size of 2mm in length,"... 
    "then Bahia Loreto field experiment, with data of surface temperature integrated in the 406 days of cultivation to get temperature"...
    "but day 0 of cultivation is when scallops are 6-months old (add 6*30 to the value)"];

data.am = 7*365;   units.am = 'd';    label.am = 'life span';                bibkey.am = 'Minchin2003';
temp.am = C2K(23);  units.temp.am = 'K'; label.temp.am = 'temperature';
comment.am = "chapter 29 in the book, not found mentioned elsewhere";

% length
data.Lb  = 0.006;  units.Lb  = 'cm';  label.Lb  = 'shell length at birth';   bibkey.Lb  = {'Abasolo2009'};
comment.Lb = "Experimental, measure start with stade D of larvae which is supposed to be the start of feeding, according to Serr-Guz97 [70-90] micrometer so taken the mean";

data.Lj  = 0.0210;  units.Lj  = 'cm';  label.Lj  = 'shell length at metam';   bibkey.Lj  = 'NavaG2022';
comment.Lj = "experimental - no exact measure when metam start, but measure before and after [200-400] micro m";

data.Lp  = 6.6;    units.Lp  = 'cm';  label.Lp  = 'shell length at puberty'; bibkey.Lp  = 'Ram2009';
comment.Lp = "Length, 50% of individuals reach sexual maturity at this size temperature extract from digitize R package from the study, integration from the beginning to day 406 for the temperature";

data.Li  = 19;   units.Li  = 'cm';  label.Li  = 'ultimate shell length  ';   bibkey.Li  = 'VelAbu2006';
comment.Li = 'L infinite for the von bertalanffy model, random sampling from Bahia Los Angeles, , from January to December 2006';

% weight
data.Wdb = 2.3398e-8;  units.Wdb = 'g'; label.Wdb = 'dry weight at birth';     bibkey.Wdb = 'NavaG2022';
temp.Wdb = C2K(26.5); units.temp.Wdb = 'K'; label.temp.Wdb = 'Kelvin';
comment.Wdb = "experimental, from linear regression of Shell length/dry weight of larva througout 35 days of experimentation";

data.Wdj = 1.68169e-7;   units.Wdj = 'g'; label.Wdj = 'dry weight at metam';     bibkey.Wdj = 'NavaG2022';
temp.Wdj = C2K(26.5); units.temp.Wdj = 'K'; label.temp.Wdj = 'Kelvin';
comment.Wdj = "experimental - end of measure during metamorphosis, last measure was near 200 micro m, from linear regression of DW and shell length of larvae";

data.Wwp = 25.34;     units.Wwp = 'g'; label.Wwp = 'Mass Wet weight at puberty';   bibkey.Wwp = 'VelAbu2006';
comment.Wwp = "from the linear regression of Length/Total weight (TW = 0.2979 * L^2.9193) and Mass weight/Total weight (MW = 2.2213 + 0.3144*TW) ,  random sampling from Bahia Los Angeles, from January to December 2006";

data.Wwi = 513.5 ;    units.Wwi = 'g'; label.Wwi = 'ultimate wet total weight';     bibkey.Wwi = 'VelAbu2006';
comment.Wwi = "from the linear regression extended of Length/Total weight (TW = 0.2979 * L^2.9193) and Mass weight/Total weight (MW = 2.2213 + 0.3144*TW) ,  random sampling from Bahia Los Angeles, from January to December 2006";

% reproduction
% data.Ri  = 54795; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = {'MalA2004'};
%  temp.Ri = C2K(20); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
%  comment.Ri = "Average of the values given in this paper, with spawning induced in the lab. Hypothesis of one spawning event per year 20 million/365";


%   data.Ri  = 84247*2; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = {'MalA2004'};
%  temp.Ri = C2K(20); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
%  comment.Ri = ["Average of the values given in this paper, with spawning induced in the lab. Hypothesis of 1.5 spawning events per year 20.5 million*1.5/365"...
%      "1.5 because usually one major and maybe another one (Racotta et al., 2003)"];
 
data.R_L  = 84247*2; units.R_L  = '#/d'; label.R_L  = 'reprod rate at a large size';     bibkey.R_L  = {'MalA2004'}; 
   temp.R_L = C2K(20); units.temp.R_L = 'K'; label.temp.R_L = 'temperature';
   comment.R_L = ["Average of the values given in this paper, with spawning induced in the lab. Hypothesis of 1.5 spawning events per year 20.5 million*1.5/365"...
     "1.5 because usually one major and maybe another one (Racotta et al., 2003)"];
L_R.R_L = 10;  units.L_R.R_L = 'cm';         label.L_R.R_L = 'length at reproduction data'; 


data.E_0 = 2.12e-3; units.E_0  = 'J';  label.E_0  = 'egg energy content'; bibkey.E_0  = 'calculation';
  comment.E_0 = 'calculated from egg diameter (52.5 µm), which gives egg dry weight from M. varia and P. magellanicus data, and times 26kJ/gdw from P. magellanicus litterature';  

 %% 
% -------------------------------------------------------------------------
% uni-variate data
% -------------------------------------------------------------------------

% Racotta 2003, Bahia Magdalena 1998 - 1999, shell height 
% -----------------------------
data.tLh_RC = [ ... time (day) ~ shell height (cm)
1     2.06
30    3.24
60    4.55
120   5.56
150   6.44
210   7.58
270   8.46
300   8.76
330   9.22
360   9.73
420   10.44
450   10.91
];
data.tLh_RC(:,1) = data.tLh_RC(:,1) + 5*30; %  julian d to day since birth     
units.tLh_RC   = {'d', 'cm'};  label.tLh_RC = {'time since birth', 'shell height Racotta 2003'};  
Lw0.tLh_RC = 2; units.Lw0.tLh_RC = 'cm'; label.Lw0.tLh_RC = 'initial physical length'; 
temp.tLh_RC    = C2K(22);  units.temp.tLh_RC = 'K'; label.temp.tLh_RC = 'temperature'; 
bibkey.tLh_RC = 'Racotta2003';
comment.tLh_RC = ['bahia magdalena, spawning = nov 98, beginning of measurements 5 months later in feb 99, temperature not known yet' ...
    ', data evaluated in days by multiply by 30 the months. ' ...
    'data extracted from graphics of the study with digitize R packages'];

tT_RA = [... culture day (d) ~ temperature (°C)
    % from the graph on the paper, figure 18, plotdigitizer
    36	23.5
    84	25.8
    137	29.8
    186	29.3
    231	25.8
    282	21.0
    351	19.9
    406	23.2
    476	27.5
    ];
tT_RA(:,1) = tT_RA(:,1) + (6*30);
tT_RA(:,2) = tT_RA(:,2) + 273;


% Ramirez-Arce 2009, Bahia Loreto, morphophological data 
% -----------------------------
tLW_RA =  [ ... {1}Culture day (d), {2}Height (cm), {3}Height SD min (cm), {4}Height SD max(cm), {5}Total weight(g), {6}Width(cm), {7}Mass weight(g)
    % diploid population
% 40.00	95.00	150.00	206.00	259.00	314.00	370.00	424.00	479.00; % culture day, from digitize R packages, but not the same space between x-points
36      84      137     186     231     282     351     406     476;    % real culture day in the paper
2.57	3.71	4.03	4.20	4.83	5.50	6.21	6.78	7.03;   % shell height (mean)
2.2653	3.4082	3.7143	3.8980	4.5306	5.2041	5.9184	6.4898	6.7347; % shell height (mean - SD)
2.8367	3.9796	4.3061	4.5102	5.1224	5.7755	6.4898	7.0816	7.3061; % shell height (mean + SD)
3.28	9.47	13.48	16.03	23.31	33.87	50.26	71.39	81.95;  % total weight
0.77	1.27	1.30	1.53	1.79	2.06	2.41	2.77	2.91;   % width 
1.15	2.82	2.93	3.66	5.64	9.20	14.42	22.15	19.02   % biomass weight (tissues)
]';

% wet weight at time
data.tWw_RA = sort(tLW_RA(:,[1 7])); % time(d) ~ biomass weight(g)
data.tWw_RA(:,1) = data.tWw_RA(:,1) + (6*30); % age since fertilisation     %% modified 04/10/2023 ELM
units.tWw_RA = {'d', 'g'}; label.tWw_RA = {'Age since fertilisation', 'Flesh wet weight Ramirez-Arce 2009'};  
temp.tWw_RA = tT_RA; units.temp.tWw_RA = {'d', 'K'}; label.temp.tWw_RA = {'time','temperature'};
% food.tWw_RA = tf; units.food.tWw_RA = {'d','mg/m^3'}; label.food.tWw_RA = {'time','food density'};
bibkey.tWw_RA = 'Ram2009';
comment.tWw_RA = ['wild population on Bahia Loreto, BCS, individual from Bahia de Loreto, larval culture in lab then transported to Bahia Loreto in a suspension system', ...
    'data extracted with digitize R packages from the study, for the diploid population ', ...
    'thesis, in spanish',...
    'day 0 of culture correspond to 6 months after fertilisation'];

% % wet weight at shell height
% data.LWw_RA = sort(tLW_RA(:,[2 7])); % shell height (cm) ~ biomass weight(g ww)
% units.LWw_RA = {'cm', 'g'}; label.LWw_RA = {'Shell height', 'Flesh wet weight Ramirez-Arce'};  
% temp.LWw_RA = tT_RA; units.temp.LWw_RA = {'d', 'K'}; label.temp.LWw_RA = {'time','temperature'};
% bibkey.LWw_RA = 'Ram2009';
% comment.LWw_RA = ['wild population on Bahia Loreto, BCS, individual from Bahia de Loreto, larval culture in lab then transported to Bahia Loreto in a suspension system', ...
%     'data extracted with digitize R packages from the study, for the diploid population ', ...
%     'thesis, in spanish',...
%     'day 0 of culture correspond to 6 months after fertilisation'];

% shell height at time
data.tLh_RA = sort(tLW_RA(:,[1 2])); % time(d) ~ shell height (cm)       %% modified 04/10/2023 ELM
data.tLh_RA(:,1) = data.tLh_RA(:,1) + (6*30); % age since fertilisation     %% modified 04/10/2023 ELM
units.tLh_RA = {'d', 'cm'}; label.tLh_RA = {'Age since fertilisation', 'Shell height Ramirez-Arce 2009'};  
Lw0.tLh_RA = 1.5; units.Lw0.tLh_RA = 'cm'; label.Lw0.tLh_RA = 'initial physical length'; 
temp.tLh_RA = tT_RA; units.temp.tLh_RA = 'K'; label.temp.tLh_RA = 'temperature';
bibkey.tLh_RA = 'Ram2009';
comment.tLh_RA = ['wild population on Bahia Loreto, BCS, individual from Bahia de Loreto, larval culture in lab then transported to Bahia Loreto in a suspension system', ...
    'data extracted with digitize R packages from the study, for the diploid population ',...
    'day 0 of culture correspond to 6 months after fertilisation'];

% % Nava-Gomez 2022, Experimental, length of larvae 
% % -----------------------------
% % adults from Bahia de los Angeles
% data.NG_tLl = [... age(d), shell length (micro m)
% 1	4	6	11	13	21	27	32	36	41;
% 91.26 119.02	135.52	188.22	206.32	501.72	755.17	904.02	1085.06	1310.34]';
% data.NG_tLl(:,2) = data.NG_tLl(:,2)/10000; % convert µm to cm
% temp.NG_tLl = C2K(26.5); units.temp.NG_tLl = 'K'; label.temp.NG_tLl = 'temperature';
% units.NG_tLl = {'d','cm'}; label.NG_tLl = {'age (N. subnodosus)','shell length Nava-Gomez'};
% bibkey.NG_tLl = 'NavaG2022';
% comment.NG_tLl = 'Experimental individual, individual initially from Bahia Los Angeles, adult spawning in lab and larvae culture at 26.5°C';
% 
% % Serrano-Guzman 1997, no mention of location from where spat provid, no mention of temperature, 1996
% % from Kino Bay
% % -----------------------------
% data.SG_tLl = [... age (d), length (micro m)
% 1	2	3	4	5	6	7	8	9	10	11	12	13	14	15;
% 90.94	98.44	103.13	116.25	125.63	136.88	144.38	154.69	161.25	166.88	169.69	180	180.94	187.5	193.13
% ]';
% data.SG_tLl(:,2) = data.SG_tLl(:,2)/10000; % convert µm to cm
% units.SG_tLl = {'d','cm'}; label.SG_tLl = {'age (N. subnodosus)','shell length Serrano-Guzman'};
% bibkey.SG_tLl = 'Serrano1997';
% comment.SG_tLl = 'Experimental, not yet mention of temperature, 1996';



% Arellano-Martinez 2011, Laguna Ojo de Liebre, Guerrero Negro, Temperatures 
% -----------------------------
tT_AM11 = [...  time (months), temperature (°K)
    % from graphs figure 4, plotDigitizer online (October 2023)
    2	25
    3	24
    4	23
    5	20
    6	17
    7	16
    8	16
    9	17
    10	17
    11	19
    12	22
    13	22
    14	24
    15	23
    16	22
    17	18
    18	15
    19	15
    20	17
    21	16
    22	18
    23	19]; 
tT_AM11(:,1) = tT_AM11(:,1) * 30; % convert months to days
tT_AM11(:,2) = tT_AM11(:,2) + 273; % convert °C to K

% Arellano-Martinez 2011, Laguna Ojo de Liebre, Guerrero Negro, wild populations
% -----------------------------
data.tLh_AM = [... age (months), shell height (mm)
    % from graphs figure 2, plotDigitizer online (October 2023)
    3	18
    4	29
    5	36
    6	43
    8	48
    9	51
    10	54
    11	55
    13	60
    14	69
    15	80
    16	86
    17	92
    20	105
    21	106
    22	111
    23	119
];
data.tLh_AM(:,2) = data.tLh_AM(:,2) / 10 ; % convert mm to cm
data.tLh_AM(:,1) = data.tLh_AM(:,1) * 30; % convert months to days
units.tLh_AM = {'d','cm'}; label.tLh_AM = {'age (N. subnodosus)','shell height, Arellano-Martinez'};
Lw0.tLh_AM = 1.8; units.Lw0.tLh_AM = 'cm'; label.Lw0.tLh_AM = 'initial physical length'; 
temp.tLh_AM = tT_AM11; units.temp.tLh_AM = {'d','K'}; label.temp.tLh_AM = {'time','temperature'};
bibkey.tLh_AM = 'Arrellano2011';
comment.tLh_AM = 'Spat from same locations, 4000 lions paw, at the beginning, 3mm shel height (2 months old), spat production method Racotta 2003 and then transported in laguna Ojo de liebre at Guerrero negro';

data.CL_LWd = [... shell length (cm), tissue dry weight (g dw)
% -----------------------------
    % from figure 5C, plotDigitize online (February 2024)
56.98466466	0.553635292
56.64395438	1.245676766
58.68824204	2.283738978
57.83645335	2.076124424
59.54003073	1.107270583
60.39181942	1.522494412
61.92504166	2.283738978
59.880754	2.352942069
61.75468652	1.45329132
63.96932932	1.176473675
63.11754063	1.937718241
62.77683035	2.560556624
62.77683035	2.560556624
64.31005259	2.422145161
65.16184128	3.114186635
64.31005259	2.422145161
66.86541866	2.283738978
65.50255156	1.522494412
66.69506352	3.044983544
69.25042959	2.906577361
66.86541866	2.283738978
65.6729067	2.283738978
68.73935118	1.937718241
69.76149501	2.145332795
70.6132837	2.906577361
71.12436211	2.422145161
71.12436211	2.422145161
70.78363884	2.906577361
69.25042959	2.837374269
68.22827277	3.460207373
69.59113987	3.875436481
70.6132837	3.598618835
69.59113987	4.498266945
61.92504166	4.013842664
63.79897418	5.605537528
63.62861904	6.297579002
72.82793949	3.044983544
73.50937304	2.560556624
73.16864977	3.529415744
74.02043846	4.221454579
72.65758435	5.121110607
75.5536607	4.77508987
74.53151687	5.25951679
77.08688294	4.844287682
75.5536607	5.813152082
76.06473911	6.574391368
78.79046032	5.882349894
79.81260415	6.505193557
82.02725994	7.197235031
79.47189387	8.442909157
77.59794836	8.09688842
85.09369143	6.782008562
86.11583526	7.197235031
87.98978077	7.404846945
91.39693553	8.304495054
87.13799208	8.304495054
85.09369143	8.719724163
81.85689181	9.480971369
83.90119246	9.757786375
81.85689181	10.31142167
86.79726881	10.24221594
84.92333629	10.51903358
86.28620339	10.86505432
89.18227974	10.31142167
89.0119246	11.28027815
89.86371329	12.52595227
92.41907936	11.14186932
86.11583526	14.11764714
98.38160019	11.48789006
118.6541736	11.55709579
119.165239	14.74048552
113.3730733	15.91695655
104.344121	16.26297729
106.8994871	16.88581567
125.1277599	19.30796084
108.2623412	22.21453292
125.6388383	25.95155793
124.7870496	26.4359875
122.4020386	26.92041574
136.8824464	34.39446445
    ];
data.CL_LWd(:,1) = (data.CL_LWd(:,1) / 10) / 0.977 ; % convert mm to cm, and convert from shell length to shell height
% ratio 1:0.977 between height and length, thus shell length divided by
% 0.977 to get shell height
units.CL_LWd = {'cm','g dw'}; label.CL_LWd = {'shell height', 'dry weight, Carreno-Leon'};
%Lw0.CL_LWd = 1.8; units.Lw0.CL_LWd = 'cm'; label.Lw0.CL_LWd = 'initial physical length'; 
temp.CL_LWd = C2K(21); units.temp.CL_LWd = {'d','K'}; label.temp.CL_LWd = {'time','temperature'};
%values from 2 Bays, Ojo de Liebre and Guerrero with average temperature at
%19°C, and Bahia de LosAngeles around 22°C
bibkey.CL_LWd = 'Carreno2023';
comment.CL_LWd = '';


% Maldonado-Amparo et al., 2004
% -----------------------------
data.tLl_MA = [... date (day) as 19/04/2001 the day 0, shell length (cm)
    % data given by Ilie Racotta (March 2024)
    0	1.75
    0	2.41
    0	2.19
    0	2.4
    0	2.64
    0	2.69
    0	2.07
    0	3.25
    0	2.31
    0	2.47
    0	1.74
    0	2.93
    22	3.84
    22	3.52
    22	3.33
    22	3.03
    22	3.42
    22	3.07
    22	2.53
    22	3.29
    22	3.98
    22	2.85
    22	3.31
    22	3.27
    22	3.33
    22	2.55
    22	2.52
    22	3.3
    22	3.73
    22	3.27
    22	3.16
    22	3.05
    22	3.06
    22	2.59
    22	3.53
    22	2.81
    22	3.65
    22	2.74
    22	2.99
    22	3.07
    22	2.97
    22	3.84
    22	2.84
    22	3.46
    54	3.31724
    54	4.85394
    54	4.0767
    54	4.318
    54	4.5466
    54	3.57632
    54	4.51358
    54	3.26644
    54	5.22732
    54	4.22148
    54	4.090162
    54	4.31292
    54	3.0988
    54	3.89128
    54	4.2164
    54	4.0767
    54	4.41452
    54	3.37312
    54	3.86334
    54	3.66268
    54	3.82778
    54	4.23926
    54	4.34594
    54	4.2164
    54	4.2291
    54	3.53568
    54	3.34264
    54	4.40182
    54	4.1656
    54	3.94716
    54	4.45262
    54	4.09448
    54	4.07162
    54	4.46024
    54	4.3307
    54	4.00558
    54	4.10972
    54	3.94208
    54	4.31038
    54	2.55016
    54	3.91922
    54	4.00558
    54	4.0894
    54	3.80238
    84	4.88696
    84	4.6736
    84	4.13512
    84	4.37642
    84	5.14858
    84	4.35356
    84	4.6482
    84	4.56692
    84	3.69824
    84	4.32816
    84	5.1689
    84	3.88112
    125	4.84886
    125	5.72262
    125	4.27736
    125	5.8801
    125	5.05206
    125	4.318
    125	5.29082
    125	5.45338
    125	5.36194
    125	4.7244
    125	5.2451
    125	5.0292
    125	5.74802
    125	5.84708
    125	5.588
    125	5.64388
    125	5.50164
    125	5.52704
    125	4.85394
    125	5.25018
    125	4.8895
    125	5.2197
    125	4.63042
    125	5.3086
    125	5.05714
    125	5.24002
    125	5.37718
    125	5.02666
    125	5.0419
    125	4.02082
    125	5.59054
    125	5.48132
    125	6.48716
    125	4.98348
    125	4.89966
    125	4.5466
    125	5.22732
    125	3.97764
    125	4.03098
    125	5.1689
    125	5.58292
    125	4.63804
    125	5.17398
    125	5.334
    154	5.9309
    154	5.97408
    154	5.08
    154	6.4389
    154	5.842
    154	6.54304
    154	6.56844
    154	5.69214
    154	5.98424
    154	5.97916
    154	4.374134
    242	6.69952
    242	6.28992
    242	6.75072
    242	7.15264
    242	7.23968
    242	7.1424
    242	6.06208
    242	7.5776
    242	5.47328
    242	6.6432
    242	5.67552
    273	7.27202
    273	7.10692
    273	6.63956
    273	7.62254
    273	6.48462
    273	6.17982
    273	7.10184
    273	7.53872
    273	4.12242
    273	6.52272
    273	4.1529
    273	5.9436
    307	7.239
    307	6.8961
    307	6.80466
    307	6.8453
    307	6.6548
    307	7.05358
    307	7.03072
    307	7.61238
    307	7.17804
    307	5.48894
    307	7.1882
    307	6.85546
    307	5.45846
    307	6.93166
    307	6.16966
    307	7.03072
    307	7.55904
    307	7.30504
    307	6.25856
    307	7.3914
    307	7.29488
    307	7.66064
    307	7.23138
    307	7.84098
    307	7.53364
    307	6.32968
    307	6.92404
    307	7.0358
    307	5.63626
    307	6.72592
    307	5.2705
    307	7.81812
    307	7.1628
    307	7.5692
    307	6.42366
    307	6.55828
    307	6.30682
    307	7.6835
    307	6.6675
    307	6.61162
    307	8.42518
    307	6.12394
    335	8.84936
    335	6.77418
    335	8.2042
    335	7.74192
    335	6.8834
    335	8.24992
    335	7.15518
    335	8.6995
    335	7.5692
    335	8.77062
    335	6.4008
    335	6.57352
    363	8.6995
    363	8.8011
    363	7.0231
    363	7.0866
    363	5.83184
    363	8.09244
    363	8.32612
    363	7.6327
    363	7.2263
    363	8.7122
    363	6.59892
    363	8.01116
    363	7.71652
    363	7.7089
    363	6.0706
    363	7.7343
    363	6.7945
    363	8.09498
    363	7.6962
    363	8.05942
    363	7.6454
    363	7.8359
    363	8.80364
    363	7.46252
    363	7.7216
    363	5.67944
    363	9.07288
    363	9.6393
    363	7.47014
    363	7.8867
    363	7.77748
    363	8.51408
    363	7.65302
    363	8.4074
    363	6.7945
    363	7.0739
    363	7.6962
    363	8.28294
    363	6.16458
    363	7.34568
    363	7.5184
    363	6.5024
    363	8.0137
    363	2.98958
    405	9.13384
    405	7.239
    405	8.5217
    405	9.67232
    405	8.0264
    405	7.4676
    405	8.79348
    405	9.11098
    405	7.5438
    405	9.43356
    405	9.2837
    405	8.4582
    427	9.1821
    427	7.9883
    427	8.59536
    427	8.0645
    427	7.112
    427	8.128
    427	9.95172
    427	8.5725
    427	7.9502
    427	9.19988
    427	7.8105
    427	8.18388
    427	8.636
    427	8.41248
    427	8.44804
    427	7.6327
    427	8.78586
    427	8.34136
    427	8.21436
    427	9.35736
    427	9.03478
    427	8.28548
    427	9.69518
    427	9.398
    427	8.92048
    427	9.03986
    427	8.06958
    427	7.95528
    427	8.04418
    427	7.24916
    427	7.73176
    427	8.68426
    427	8.06196
    427	8.18896
    427	9.8806
    427	7.84352
    427	6.4008
    427	8.24484
    427	9.75106
    427	9.4107
    427	8.04418
    427	7.5692
    427	9.25068
    446	8.85698
    446	8.89254
    446	8.2042
    446	9.2075
    446	7.68604
    446	9.56564
    446	7.112
    446	9.42086
    446	9.51738
    446	9.20496
    446	7.3406
    446	8.6741
    446	11.45794
    446	8.60298
    446	9.69264
    446	9.39292
    446	8.55726
    446	7.32536
    446	8.44296
    446	8.37946
    446	9.0297
    446	8.95096
    446	9.40562
    446	8.26516
    446	10.42162
    446	8.6614
    446	8.62076
    446	7.5692
    446	8.4455
    446	8.94334
    446	9.0551
    446	8.71728
    446	9.4107
    446	6.84022
    446	9.525
    446	8.7503
    446	7.2517
    446	9.42086
    446	8.128
    446	10.9728
    446	8.8265
    446	8.9789
    490	8.6868
    490	10.0838
    490	10.795
    490	9.71296
    490	9.7409
    490	8.5725
    490	10.19302
    490	9.65708
    490	8.75792
    490	10.5029
    490	10.5537
    490	9.25068
    490	9.32942
    490	8.36168
    490	9.3218
    490	6.63956
    490	10.68832
    490	8.9154
    490	7.37362
    490	6.73862
    490	6.53542
    490	8.61822
    490	9.02462
    490	9.4361
    490	9.8806
    490	8.89508
    490	9.779
    490	8.8392
    490	10.40638
    490	9.76376
    490	10.033
    490	7.7978
    490	9.58342
    490	10.2997
    490	9.7917
    490	7.3406
    490	10.4775
    490	10.91438
    490	9.6647
    490	9.55802
    490	8.98906
    490	10.71626
    490	9.62152
    539	9.90854
    539	9.3599
    539	9.11098
    539	10.30732
    539	9.07542
    539	10.0965
    539	8.87476
    539	10.25906
    539	6.94944
    539	10.00252
    539	10.38352
    539	9.73582
    539	9.5631
    539	9.16432
    539	9.2964
    539	10.62482
    539	10.1727
    539	9.60882
    539	9.53516
    539	10.19556
    539	9.2075
    539	8.74522
    539	10.414
    539	10.3632
    539	9.08558
    539	9.06526
    539	10.3632
    539	9.08304
    539	9.8171
    539	9.54024
    539	10.2108
    539	10.0584
    539	9.13892
    539	10.2616
    539	9.86536
    539	10.06348
    539	10.00506
    539	9.779
    539	10.80262
    539	10.24128
    539	8.04672
    539	9.6647
    539	9.398
    539	10.49274
    601	8.62838
    601	9.79424
    601	8.6995
    601	9.60628
    601	10.21334
    601	11.01344
    601	11.00836
    601	9.32942
    601	8.5979
    601	8.64616
    601	8.66902
    601	10.53592
    601	9.45642
    601	7.93496
    601	8.74522
    601	10.2616
    601	9.24814
    601	9.1821
    601	9.7282
    601	9.55294
    601	9.69772
    601	8.8265
    601	10.12444
    601	9.9441
    601	9.8679
    601	9.0424
    601	9.652
    601	8.8773
    601	9.64692
    601	10.7188
    601	10.6299
    601	8.20166
    601	9.4742
    601	10.64768
    601	8.8773
    601	10.17016
    601	10.29462
    601	6.87578
    601	9.77392
    601	10.8585
    601	10.54862
    601	10.08888
    601	8.13308
    601	10.22858
];
data.tLl_MA(:,1) = data.tLl_MA(:,1) + (30*3); % around 3-months old individuals when started measurements
units.tLl_MA = {'d','cm'}; label.tLl_MA = {'age (N. subnodosus)','shell height, Maldonado-Amparo'};
Lw0.tLl_MA = 1.7; units.Lw0.tLl_MA = 'cm'; label.Lw0.tLl_MA = 'initial physical length'; 
% temp.tLl_MA = C2K(20); units.temp.tLl_MA = 'K'; label.temp.tLl_MA = 'temperature';
temp.tLl_MA = [ 365.0000   22.6693      % Fourier parameters from data, start value in April as the growth data
                0.6978    3.1546
                -0.9459    0.5506]; units.temp.tLl_MA = 'K'; label.temp.tLl_MA = 'Fourier parameters for temperature';
bibkey.tLl_MA = 'MaldAmp2004';
comment.tLl_MA = 'only control treatment taken here, and assume individuals were 3-months old when starting the experiment';


%% set weights for all real data
weights = setweights(data, []);
weights.tj = weights.tj * 0;
weights.Lj = weights.Lj * 0;
weights.Wdj = weights.Wdj * 0;
% weights.E_0 = 5 * weights.E_0;
% weights.R_L = 5 * weights.R_L;
% weights.Li = 5 * weights.Li;
% weights.Wwi = 5 * weights.Wwi;
% weights.tp = 2 * weights.tp;
weights.Lp = 2 * weights.Lp;
% weights.Wwp = 2 * weights.Wwp;
% weights.tLh_RC = 0 * weights.tLh_RC;
% weights.tWw_RA = 2 * weights.tWw_RA;
% weights.tLh_RA = 10 * weights.tLh_RA;
% weights.LWw_RA = 2 * weights.LWw_RA;
% weights.tLh_AM = 10 * weights.tLh_AM;
% weights.CL_LWd = 5 * weights.CL_LWd;
% weights.tLl_MA = 2 * weights.tLl_MA;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
% weights.psd.kap = 5;
% weights.psd.pM = 10;
% weights.psd.kap = 100 * weights.psd.kap;
%% pack auxData and txtData for output
auxData.temp = temp;
% auxData.food = food;
auxData.Lw0 = Lw0;
auxData.L_R = L_R;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% add ultimate length as metadata
metaData.Li = data.Li;

%% Links
 metaData.links.id_CoL = '47L85'; % Cat of Life
metaData.links.id_ITIS = '567974'; % ITIS
% metaData.links.id_EoL = ''; % Ency of Life
 metaData.links.id_Wiki = 'Nodipecten_subnodosus'; % Wikipedia
metaData.links.id_ADW = 'Nodipecten_subnodosus'; % ADW
% metaData.links.id_Taxo = ''; % Taxonomicon
metaData.links.id_WoRMS = '394383'; % WoRMS
metaData.links.id_molluscabase = '394383'; % molluscabase


%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Nodipecten_subnodosus}}';
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
bibkey = 'NavaG2022'; type = 'Article'; bib = [ ... 
'author = {Shumway, S. E. and Parsons, G. J.}, ' ... 
'year = {2022}, ' ...
'title = {Survival, growth and biochemical composition of larvae of the the lions paw scallop, Nodipecten subnodosus, in batch- and flow-through culture}, ' ...
'journal = {Aquaculture}, ' ...
'volume = {555}, '...
'pages = {738181}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
bibkey = 'Serrano1997'; type = 'Article'; bib = [ ... 
'author = {Serrano-Guzman, S., Robles-Mungaray, M., Velasco-Blanco, G., Voltolina, L.D. & Hoyos, C.F.J.Z.}, ' ... 
'year = {1997}, ' ...
'title = {Larval culture of Mexican pectinids Argopecten ventricosus (circularis) (Sowerby 1842) and Lyropecten subnodosus (Sowerby 1835) in a commercial hatchery}, ' ...
'journal = {Haliotis}, ' ...
'volume = {16}, ' ...
'pages = {363--381}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
bibkey = 'Ram2009'; type = 'Thesis'; bib = [...
'author = {Ramirez-Arce., J.}, ' ...     
'year = {2009}, ' ...
'title = {EVALUACIÓN DE LA VENTAJA PRODUCTIVA Y GRADO DE ESTERILIDAD EN TRIPLOIDES DE ALMEJA MANO DE LEÓN Nodipecten subnodosus (Sowerby 1835)COMO UNA ALTERNATIVA PARA EL CULTIVO EN EL PARQUE NACIONAL BAHÍA DE LORETO, GOLFO DE CALIFORNIA}, ' ...
'school = {Instituto politecnico nacional}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% 
bibkey = 'VelAbu2006'; type = 'Article'; bib = [ ... 
'author = {Velazques-Abunader, I., Lopez-Rocha, J.A., Arellano-Martinez, M., Ceballos-Velasquez, B.P. & Cabrera, M.A.}, ' ... 
'year = {2016}, ' ...
'title = {Estimation of growth parameters in a wild population of lion-paw scallop (Nodipecten subnodosus ) in Bahia de Los Angeles, Baja California, Mexico}, ' ...
'journal = {Hidribiologica}, ' ...
'volume = {26}, ' ...
'pages = {133-142}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Arrellano2011'; type = 'Article'; bib = [ ... 
'author = {}, ' ... 
'year = {2011}, ' ...
'title = {Growth and reproduction of the lion’s paw scallop Nodipecten subnodosus in a suspended culture system at Guerrero Negro lagoon, Baja California Sur, Mexico}, ' ...
'journal = {Aquaculture Research}, ' ...
'volume = {42}, ' ...'...
'pages = {571-582}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Racotta2003'; type = 'Article'; bib = [ ... 
'author = {Racotta, I.S., Ramirez, J.L., Ibarra, A.M., Rodriguez-Jaramillo, M.C., Carreño, D. & Palacios, E.}, ' ... 
'year = {2003}, ' ...
'title = {Growth and gametogenesis in the lion-paw scallop Nodipecten (Lyropecten) subnodosus}, ' ...
'journal = {Aquaculture}, ' ...
'volume = {217}, ' ...
'pages = {335-349}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
