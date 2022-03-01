function x = p0001_calib
x = struct('General', {struct('program_name', {struct('name', {['CALrp2_PU.m ' ...
]} ...
,'date', {['07-Sep-2018 ' ...
]} ...
,'version', {['rp2 v3.0 ' ...
]} ...
)} ...
,'picture_number', {1 } ...
,'date', {['07-Sep-2018 ' ...
]} ...
,'time', {['13:36:49 ' ...
]} ...
,'comment', {['No comment. ' ...
]} ...
)} ...
,'Stimuli', {struct('frqlo', {0.2 } ...
,'frqhi', {20 } ...
,'fstlin', {0 } ...
,'fstoct', {40 } ...
,'bplo', {0.125 } ...
,'bphi', {64 } ...
,'n60lo', {2 } ...
,'n60hi', {64 } ...
,'n120lo', {2 } ...
,'n120hi', {64 } ...
,'ear', {2 } ...
,'BaseAtten', {30 } ...
,'SlopeAtten', {6 } ...
,'BeginSlope', {11 } ...
,'is_dBSPL', {1 } ...
,'crit', {100 } ...
)} ...
,'Line', {[]} ...
,'CalibData', {[[0.2;0.2034959384;0.2070529848;0.2106722072;0.2143546925;0.2181015465;0.2219138944;0.225792881;0.229739671;0.2337554497; ...
0.237841423;0.2419988178;0.2462288827;0.2505328877;0.2549121255;0.2593679109;0.2639015822;0.2685145006;0.2732080514;0.277983644; ...
0.2828427125;0.287786716;0.2928171392;0.2979354926;0.3031433133;0.3084421651;0.3138336392;0.3193193545;0.3249009585;0.3305801273; ...
0.3363585661;0.3422380103;0.3482202253;0.3543070076;0.360500185;0.3668016173;0.3732131966;0.3797368484;0.3863745316;0.3931282394; ...
0.4;0.4069918768;0.4141059695;0.4213444144;0.428709385;0.4362030931;0.4438277888;0.4515857619;0.459479342;0.4675108994; ...
0.475682846;0.4839976357;0.4924577653;0.5010657754;0.5098242509;0.5187358219;0.5278031643;0.5370290011;0.5464161027;0.5559672879; ...
0.5656854249;0.575573432;0.5856342784;0.5958709852;0.6062866266;0.6168843302;0.6276672783;0.6386387091;0.6498019171;0.6611602545; ...
0.6727171322;0.6844760205;0.6964404506;0.7086140153;0.7210003701;0.7336032346;0.7464263932;0.7594736968;0.7727490631;0.7862564788; ...
0.8;0.8139837537;0.8282119391;0.8426888288;0.85741877;0.8724061861;0.8876555777;0.9031715238;0.918958684;0.9350217988; ...
0.951365692;0.9679952714;0.9849155307;1.002131551;1.019648502;1.037471644;1.055606329;1.074058002;1.092832205;1.111934576; ...
1.13137085;1.151146864;1.171268557;1.19174197;1.212573253;1.23376866;1.255334557;1.277277418;1.299603834;1.322320509; ...
1.345434264;1.368952041;1.392880901;1.417228031;1.44200074;1.467206469;1.492852786;1.518947394;1.545498126;1.572512958; ...
1.6;1.627967507;1.656423878;1.685377658;1.71483754;1.744812372;1.775311155;1.806343048;1.837917368;1.870043598; ...
1.902731384;1.935990543;1.969831061;2.004263102;2.039297004;2.074943287;2.111212657;2.148116004;2.185664411;2.223869152; ...
2.2627417;2.302293728;2.342537114;2.383483941;2.425146506;2.467537321;2.510669113;2.554554836;2.599207668;2.644641018; ...
2.690868529;2.737904082;2.785761803;2.834456061;2.88400148;2.934412938;2.985705573;3.037894787;3.090996253;3.145025915; ...
3.2;3.255935015;3.312847756;3.370755315;3.42967508;3.489624745;3.550622311;3.612686095;3.675834736;3.740087195; ...
3.805462768;3.871981086;3.939662123;4.008526204;4.078594007;4.149886575;4.222425314;4.296232009;4.371328822;4.447738303; ...
4.5254834;4.604587456;4.685074227;4.766967882;4.850293013;4.935074641;5.021338227;5.109109673;5.198415337;5.289282036; ...
5.381737058;5.475808164;5.571523605;5.668912122;5.768002961;5.868825877;5.971411146;6.075789574;6.181992505;6.290051831; ...
6.4;6.511870029;6.625695513;6.74151063;6.85935016;6.979249489;7.101244621;7.225372191;7.351669472;7.480174391; ...
7.610925536;7.743962171;7.879324245;8.017052407;8.157188015;8.29977315;8.444850629;8.592464018;8.742657643;8.895476607; ...
9.050966799;9.209174912;9.370148454;9.533935764;9.700586026;9.870149283;10.04267645;10.21821935;10.39683067;10.57856407; ...
10.76347412;10.95161633;11.14304721;11.33782424;11.53600592;11.73765175;11.94282229;12.15157915;12.36398501;12.58010366; ...
12.8;13.02374006;13.25139103;13.48302126;13.71870032;13.95849898;14.20248924;14.45074438;14.70333894;14.96034878; ...
15.22185107;15.48792434;15.75864849;16.03410481;16.31437603;16.5995463;16.88970126;17.18492804;17.48531529;17.79095321; ...
18.1019336;18.41834982;18.74029691;19.06787153;19.40117205;19.74029857; ...
] ...
[92.30588658;92.41588749;92.23239591;92.17374217;92.42914976;92.19500795;92.52381852;92.36894492;92.3440852;92.36617966; ...
92.18182921;92.15321656;92.16221422;91.98495564;91.83503666;91.7215108;91.42285291;91.36754942;90.90234153;90.86407693; ...
90.76825107;90.53100015;90.57363553;90.55632964;90.42455535;89.97373601;90.38737996;90.33287286;90.56483109;90.76857267; ...
90.64672458;90.61793937;90.62018265;90.77135221;90.52906088;90.70890738;90.89351322;90.86771609;90.80918498;90.68433942; ...
90.5525038;90.48910534;90.54366868;90.42623121;90.14142103;90.14016189;90.06278865;90.03142345;89.85001681;89.78178402; ...
89.82392234;89.79784176;89.68267461;89.78851953;89.73014568;89.6755486;89.53000945;89.49336923;89.35299523;89.35471756; ...
89.3440229;89.24394682;89.11628522;89.04071177;89.06871199;88.99907273;88.91738229;88.93397655;88.72052583;88.73558729; ...
88.77678999;88.64187988;88.59728717;88.54541028;88.44831535;88.34946738;88.3358639;88.22986901;88.24085876;88.2795177; ...
88.24784969;88.28359359;88.42965304;88.45401478;88.61162493;88.70142413;88.7032868;88.81734969;88.87205339;88.88973834; ...
89.04748443;89.09187482;89.21267655;89.29834293;89.34774163;89.48226515;89.54924995;89.75446508;89.91219451;90.10501859; ...
90.15792468;90.3953082;90.5996584;90.74397426;90.92658869;91.07710762;91.27130261;91.47635234;91.6571317;91.82217895; ...
91.95918399;92.04572752;92.17696858;92.1875661;92.24963081;92.20405013;92.20914487;92.09836203;91.97917302;92.06628312; ...
92.09931438;92.22346311;92.37310841;92.51633616;92.65840951;92.83746617;93.01223217;93.21468029;93.43541028;93.55047282; ...
93.66698961;93.72596954;93.67782061;93.49219185;93.24152948;93.02058303;92.67662684;92.20644726;91.92201659;91.63446026; ...
91.59170365;91.23374761;91.01929342;90.84371988;90.88497387;91.21758033;91.53946024;91.9130306;92.35350201;92.80683968; ...
92.98386745;92.6199235;91.36624332;89.75051265;89.53626551;91.43657025;93.59657872;95.45710187;96.897682;97.83776991; ...
98.14860629;98.08017036;97.88170637;97.68217809;97.64726712;98.10183916;98.68350286;99.32653318;99.64288671;99.45037421; ...
98.81257024;98.16716921;98.0367984;98.57753478;99.39475714;100.0331095;100.0087734;99.32396668;98.46154818;98.08367708; ...
97.95055064;97.89842768;97.90101904;97.61345012;97.22410096;97.1023103;96.90521643;96.89818627;96.68333154;96.48165167; ...
96.05830389;96.10286433;96.63537951;97.32441449;97.72443002;97.91005598;97.72725973;97.89225891;98.15849286;98.20157366; ...
98.35119065;98.29565561;98.58650609;99.25262445;99.91788854;100.2474924;100.1308158;100.2371297;100.9720273;101.8659905; ...
102.6942039;102.7019052;102.4030832;102.2558244;102.587344;102.6020315;102.0544386;101.5218045;100.8077964;99.34851072; ...
97.90911088;97.04100312;95.80439976;94.25377718;93.4567233;93.31083352;92.46865298;92.09529464;92.78252611;92.4700361; ...
92.31506368;93.48659533;94.60709387;94.76191746;95.27969041;95.60111206;95.31029746;95.39192408;94.98916727;94.57926279; ...
94.21405365;93.68447589;93.44929913;92.50385735;92.2364049;90.33076364;89.12088171;87.39179922;86.96419762;87.40892838; ...
87.92456923;89.69598219;93.18729992;96.67952552;99.19532311;97.64178323;95.0439472;90.85867033;86.40847849;81.40324601; ...
76.25881025;71.79035886;66.84959431;61.41458478;59.31698031;56.81050934; ...
] ...
[-0.5275857626;-0.7702840285;0.02763258084;-0.2066369843;0.5661726176;-0.1706203197;0.5678594583;0.3440855809;-0.9433971247;0.2952008694; ...
-0.4699211035;-0.6707291258;-0.01090986418;-0.1985977856;0.4345148023;-0.3475716563;0.2479775299;0.07929587342;0.6473471885;-0.1572497377; ...
-0.9645073378;0.225727209;-0.5863871211;0.599314972;-0.2169947151;0.9648622872;0.1450713027;-0.6761753063;0.5013220512;0.453587444; ...
-0.3572375618;0.8315625471;0.02022240412;0.3774848279;-0.4479327972;0.7270099543;-0.09742212122;-0.9209550805;0.2566944058;0.5132121347; ...
-0.3220116181;0.8445739572;0.01328730763;-0.81554254;0.3584242937;0.5118420771;-0.3254271821;0.8408750336;0.01112939979;0.9028289981; ...
0.8836507868;0.8641373463;0.8442828168;0.8240812363;0.8035265384;-0.4347748987;-0.4773340142;0.2190444249;0.8706059808;-0.5238178091; ...
-0.6378521358;-0.1046556968;0.3792903199;0.2111232769;0.6171619722;0.9704676518;0.7427154381;-0.9877975713;-0.7747576239;0.9320247777; ...
-0.9450743716;0.7237715615;0.7523981088;0.3809965698;0.3110340498;-0.1030417812;-0.276095211;-0.735393265;0.9837656453;0.4765708028; ...
-0.03948964935;-0.4748563468;0.9574931632;0.3799203174;-0.2199688136;-0.8530827164;0.4195554527;-0.2722985313;-0.9762459099;0.1131599243; ...
-0.6539685276;0.293768544;-0.5401216933;0.2595743975;-0.6448323113;0.4349522291;-0.9786823276;0.02534784888;-0.9880312205;-0.6286509676; ...
0.2766194967;0.4612025268;-0.719354858;0.07945197758;0.01155281955;0.7182575437;0.4562926286;-0.9341068849;-0.3488101539;-0.8915102226; ...
-0.4103470677;0.8315625471;-0.7969930532;0.2199841682;0.4758602135;0.701249509;-0.6078602056;-0.506292162;-0.0649374187;0.3095713504; ...
0.2450637938;0.3588573302;0.3379387635;0.1787860317;-0.1222048964;-0.4739348425;-0.9981755036;0.3256747664;-0.5063123452;0.5018439319; ...
-0.1885716748;0.5654938514;-0.8515560596;-0.4441006572;-0.2166200042;0.6958526344;0.6346476012;0.3808884985;-0.9409102876;0.4951320432; ...
-0.9645073378;0.1472354365;0.5416773931;0.9080546895;-0.4455422499;-0.2362587857;0.09919406573;0.09273845823;-0.09991186325;-0.1307288895; ...
-0.9450743716;0.298362082;0.1700795054;0.7399947227;0.1586200712;-0.4329168303;-0.1530919048;-0.776629224;-0.8433475394;0.859424817; ...
-0.3948964935;-0.4114266984;0.4782175988;-0.9272838098;0.6801871184;-0.4644611808;-0.65815235;-0.4783786382;0.7468438274;-0.3321040454; ...
0.2768251779;0.5654938514;-0.4232743769;-0.4441006572;0.9974861493;0.08720324237;0.8159754872;0.2664877959;0.6133263129;-0.8191934406; ...
-0.6205634172;-0.5427747014;-0.9362581047;-0.4949044061;0.468647889;0.09096734129;-0.7055176635;0.1561910875;0.8045831109;0.6954168884; ...
-0.6817849195;-0.9721118934;0.9245526823;-0.6080042218;-0.6096559146;0.7219992144;-0.4036059308;0.70561796;0.951296936;-0.8514379576; ...
0.7038276298;-0.2582827174;0.05526516167;0.2908647608;-0.2932917515;0.4834085041;0.6727483214;-0.8756891846;-0.2910302436;-0.8831522119; ...
-0.8926992883;0.435777865;-0.8124484769;-0.2909102969;0.5592184394;0.2616097271;-0.7307047977;-0.3301928008;0.4472717699;0.4951320432; ...
-0.4112683444;0.4027138745;0.009806445436;-0.6220275702;-0.4092112203;-0.2172964123;-0.1165529904;0.4685732626;0.005288204897;0.9123964431; ...
-0.8206941354;0.2656037213;-0.6035785415;0.7198522366;0.7337777647;0.619664179;0.8866043289;-0.07101126396;0.9837656453;-0.179446718; ...
0.5388829739;-0.2102774728;-0.7906397126;-0.2788469856;-0.5418618461;-0.891027715;0.6763973143;0.1621361834;-0.07574814197;-0.6568323633; ...
0.7682517794;-0.9520983777;-0.9147493077;0.447194742;-0.2871507718;-0.7819918941;0.1759346326;0.07537278137;0.8242392331;0.7426980648; ...
0.3549266223;0.2316150691;0.2549675813;0.528875032;-0.1254084442;-0.7618189516; ...
] ...
[Inf;0;Inf;0;0;0;0;0;Inf;0; ...
Inf;0;Inf;0;Inf;0;0;Inf;Inf;Inf; ...
0;0;0;0;0;Inf;0;0;0;0; ...
0;0;0;0;0;Inf;0;0;Inf;0; ...
0;0;0;Inf;0;Inf;Inf;0;Inf;0; ...
Inf;0;Inf;0;0;Inf;Inf;0;0;0; ...
0;Inf;0;Inf;0;0;0;0;Inf;0; ...
0;0;0;Inf;0;0;0;0;Inf;0; ...
0;0;0;Inf;0;0;0;Inf;0;0; ...
0;0;0;0;Inf;0;0;0;0;0; ...
0;Inf;0;Inf;0;Inf;0;Inf;0;0; ...
Inf;Inf;0;Inf;Inf;0;0;0;0;Inf; ...
0;Inf;0;Inf;Inf;0;0;Inf;0;0; ...
0;0;0;0;0;0;Inf;0;0;Inf; ...
Inf;0;Inf;Inf;0;Inf;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
Inf;0;0;Inf;0;0;0;0;0;0; ...
0;0;0;0;0;Inf;0;0;0;0; ...
0;0;0;Inf;Inf;0;Inf;0;Inf;0; ...
0;0;0;0;0;0;0;Inf;Inf;0; ...
0;0;Inf;0;0;0;0;0;0;0; ...
0;0;0;0;0;Inf;0;0;0;0; ...
Inf;0;0;0;0;0;0;Inf;0;0; ...
Inf;0;0;0;0;Inf;Inf;0;0;0; ...
0;0;0;0;0;0;0;Inf;Inf;0; ...
Inf;0;Inf;0;0;0;0;0;0;0; ...
0;0;Inf;0;0;0; ...
] ...
[0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0; ...
] ...
] ...
} ...
,'User', {[]} ...
,'Hardware', {struct('mic', {['ER7c ' ...
]} ...
,'MicDate', {['12-Jul-2006 ' ...
]} ...
)} ...
);
