function x = p0001_calib
x = struct('General', {struct('program_name', {struct('name', {['CALrx8_PU.m ' ...
]} ...
,'date', {['26-Jul-2019 ' ...
]} ...
,'version', {['rx8 v1.0 ' ...
]} ...
)} ...
,'picture_number', {1 } ...
,'date', {['26-Jul-2019 ' ...
]} ...
,'time', {['16:55:43 ' ...
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
[98.77009556;98.95165244;98.96154899;98.96075401;99.09449938;99.2387115;99.10198524;99.22129099;99.20373018;99.24042927; ...
99.30061024;99.18409073;99.35871126;99.24255751;99.13149274;99.33313384;99.23904016;99.26179668;99.40799854;99.19800582; ...
99.29050717;99.545908;99.56234341;99.97095715;99.76287831;99.92907225;99.94077365;100.0545709;99.78497929;99.81457041; ...
99.49596223;99.2946076;98.89641081;99.03273754;98.9263895;98.722878;98.72141538;98.75311432;98.5571294;98.65480248; ...
98.49270455;98.37968438;98.23938117;98.20280511;98.0564769;97.96780361;97.98451527;97.80624654;97.67704507;97.60633319; ...
97.52359386;97.4852737;97.47730916;97.50165612;97.34509959;97.44190793;97.37742053;97.27630556;97.21008462;97.24791229; ...
97.12201191;97.03946197;97.04111216;96.92485419;96.90159525;96.82223151;96.70389719;96.62534752;96.61251697;96.69144397; ...
96.59607417;96.60919534;96.63038347;96.5761397;96.59379293;96.58731588;96.65349348;96.61260186;96.52938665;96.4736211; ...
96.46017005;96.30438917;96.42872432;96.33785893;96.31829701;96.34818256;96.42090875;96.5043305;96.67074825;96.85561104; ...
97.05246998;97.32222746;97.51883166;97.7845053;97.95129795;98.24526694;98.39794984;98.53220107;98.68405452;98.82590117; ...
98.92992038;99.09548124;99.28005767;99.40210589;99.58136196;99.68998849;99.8919051;100.0416606;100.2321416;100.3871412; ...
100.5302998;100.666262;100.7622233;100.8472779;100.9329365;100.942821;100.9223062;100.8409865;100.8304187;100.7999369; ...
100.982121;101.1273391;101.2070498;101.3473371;101.5032726;101.6945514;101.8595826;102.1899028;102.440707;102.6597142; ...
102.7679741;102.6362129;102.2831314;101.8189053;101.5765482;101.0991543;100.4854168;100.0790895;100.0207899;99.88729764; ...
100.0628211;101.1944937;102.6074179;103.320066;103.2380338;102.5826841;101.6351846;100.4373837;98.59421215;96.86198663; ...
96.49192149;98.21475584;100.2185219;101.7973106;102.8363818;103.5471307;104.1795925;104.9566321;105.8966104;106.7994952; ...
107.4309158;107.622291;107.3903337;107.1853826;107.1381554;107.2416992;107.8773969;108.7237545;109.3591656;109.3869725; ...
109.2758874;109.0227546;108.9131214;109.0862499;109.3246128;109.566049;109.6186999;109.2644249;108.9035673;108.5492776; ...
108.5473695;108.6933485;108.6437689;108.2117679;107.4734764;106.8446702;106.2454254;106.111751;106.0160096;105.8084448; ...
105.2408479;104.9728404;105.0616174;105.3604048;105.8430809;106.0379647;105.6714749;105.7928961;106.1191703;106.9222488; ...
107.3998473;107.5387383;107.6182707;108.0409315;108.6355903;108.8930332;109.524085;109.9102026;110.4509322;111.1114345; ...
111.7659824;112.0893052;112.2975556;112.5962079;113.0064593;113.355984;113.2903664;113.316973;113.294068;112.7995304; ...
112.0205446;111.4728931;110.5587166;108.9855491;107.8309519;107.2568098;106.1509341;104.8082864;104.3978992;103.7286406; ...
102.4415396;102.0850884;101.7871111;100.8051798;100.945505;101.3819145;101.1547523;101.2957163;101.3235014;100.9452429; ...
101.0856349;100.9078501;100.9360788;100.4048269;100.3176388;98.78786945;98.03070728;97.09859518;95.68305061;93.4719959; ...
90.95549628;89.64317068;89.42303899;87.88853729;89.83699886;92.79193968;87.49435272;82.54029888;84.6370615;85.09423152; ...
76.65093054;79.81283057;78.65857524;70.72694607;73.79297833;75.39843333; ...
] ...
[0.8805153147;0.1848588052;0.9994710172;0.7821888887;0.07533798485;0.8530764224;0.1302454272;-0.5961170893;0.6740862665;-0.05904262774; ...
0.6464281097;0.4651335279;-0.855192681;0.3893615963;-0.3690117361;0.869815793;0.1059780406;-0.6603856016;0.5708699456;-0.8524861759; ...
-0.973151697;0.2518911332;-0.5248064754;0.6969314126;-0.08271258494;-0.8635489881;0.3546187297;-0.4280056706;0.02630001704;0.005221555114; ...
0.432140537;0.437914917;0.8374378615;0.04048509173;-0.7558638008;0.4486557958;-0.3456824405;0.8614044705;-0.8365482706;0.3584264595; ...
-0.444543515;0.7548596247;-0.04303581962;-0.8378907939;-0.6354654615;0.5592355613;-0.2422458902;0.9604696932;0.1677739214;-0.6199291873; ...
0.5977771074;0.8641373463;0.8442828168;0.8240812363;0.8035265384;0.7826125507;0.7613329929;-0.52063705;-0.5646970096;0.08570931452; ...
0.6897182914;0.5968962021;-0.8719354858;-0.3904657627;-0.5371285208;-0.1340807153;0.2157470561;0.009761942971;0.2751983077;0.03529299993; ...
0.2124380236;-0.06352369877;0.02125993818;-0.2930029025;-0.3048279573;-0.3813959397;-0.7726158496;-0.953033714;0.6107949397;0.3217851277; ...
-0.1620215463;-0.5645706793;0.9011700359;0.3575720634;-0.207748324;-0.8056892322;0.5859180366;-0.1527038717;-0.8196013883;0.5018439319; ...
-0.4212701012;0.8372191588;0.08274703943;0.9632494523;0.1410615352;0.8697271278;-0.02401429912;0.5459848989;-0.4233342109;-0.01912384396; ...
0.931760351;-0.8372455743;0.02941933631;0.07945197758;0.8572623265;-0.3862908234;-0.5976441354;0.06345262931;-0.3488101539;0.2117579996; ...
-0.4103470677;0.04426728689;0.4718687762;-0.4540153041;-0.1400017936;0.701249509;0.8956191558;-0.94157306;-0.4379081244;-0.403288139; ...
0.2450637938;0.3588573302;0.225292509;0.1340895238;-0.09776391716;-0.4739348425;-0.8318129196;0.5648640855;-0.193023302;0.8905279394; ...
-0.1885716748;0.5654938514;-0.8515560596;0.2595743975;0.5692738423;0.6958526344;-0.4106843703;-0.5778374014;0.1884837316;-0.2858137095; ...
0.3457743708;0.9534430322;0.03922578175;0.9080546895;-0.4455422499;-0.02716205154;0.04959703286;0.004880971486;-0.1998237265;-0.4130728889; ...
-0.6300495811;-0.850818959;0.9245526823;0.04399683364;0.9268960569;-0.4329168303;-0.1530919048;-0.776629224;-0.7784101207;0.9531416057; ...
0.411148289;-0.4114266984;0.4782175988;0.3242184121;-0.5865835029;-0.4644611808;-0.9835795321;-0.6973009343;0.4558135838;0.4452639697; ...
-0.99301746;-0.4773209053;-0.2286873269;-0.4441006572;0.9974861493;0.0879125642;-0.1867211395;-0.6415424062;-0.5878803834;-0.8191934406; ...
-0.4785927683;0.3255088515;-0.9362581047;-0.4949044061;0.468647889;0.09096734129;-0.7055176635;0.1561910875;0.6012339145;0.6954168884; ...
-0.3522311494;0.6852587377;0.9245526823;-0.7760401072;-0.9020716869;0.5366653574;-0.5149454979;0.70561796;-0.9302103989;-0.6640043802; ...
0.8617862273;-0.2582827174;-0.7906397126;0.3575720634;-0.9776391716;0.4834085041;0.7055892572;0.518912684;0.4810629645;0.9052793944; ...
-0.06222375444;0.3501483517;0.5255752077;0.7411991237;0.2821230704;-0.6061667657;0.3518692651;-0.146253205;-0.9394019172;0.1140470972; ...
-0.4112683444;0.2013569372;0.009806445436;-0.7491507344;0.7326901678;0.510442318;0.1983881315;0.1952388594;0.007050939863;0.7816675536; ...
0.8788100537;-0.9278446065;0.6803180218;-0.8481308761;-0.973806463;0.8049980361;0.997943896;-0.07101126396;0.9837656453;-0.179446718; ...
0.5388829739;-0.2102774728;-0.7906397126;-0.2788469856;0.8268329941;0.6919297644;-0.6399073856;0.4324272367;-0.6830580098;-0.4231367872; ...
-0.8926992883;-0.7808393511;-0.6419470923;-0.8882013145;0.4391253968;0.08720324237;0.3410827582;0.03768639069;-0.3515215338;-0.03894365446; ...
-0.6095807154;-0.5738126799;0.4510964901;-0.983014687;0.310253557;0.8743116832; ...
] ...
[0;0;0;Inf;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
Inf;0;0;0;0;0;Inf;0;0;0; ...
0;Inf;Inf;Inf;0;0;0;0;Inf;0; ...
0;Inf;0;0;Inf;0;Inf;0;0;0; ...
0;0;Inf;0;0;0;0;0;Inf;Inf; ...
0;0;0;Inf;0;Inf;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
Inf;0;Inf;0;0;0;0;0;Inf;Inf; ...
0;0;0;0;0;Inf;Inf;0;Inf;0; ...
0;0;0;0;0;0;0;0;Inf;0; ...
0;0;Inf;0;Inf;Inf;0;Inf;Inf;0; ...
Inf;0;0;0;0;0;0;0;0;0; ...
0;0;0;Inf;0;0;0;0;Inf;0; ...
0;0;0;0;0;Inf;0;0;0;0; ...
0;Inf;0;Inf;0;0;Inf;0;0;0; ...
0;Inf;0;0;0;0;0;Inf;0;Inf; ...
0;0;0;0;0;0;0;Inf;0;0; ...
0;0;0;0;Inf;0;Inf;0;0;0; ...
0;0;0;0;0;0;0;Inf;0;Inf; ...
0;0;Inf;0;0;0;0;Inf;0;0; ...
0;0;Inf;0;0;0;Inf;Inf;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;Inf;0;0;0;0; ...
0;0;Inf;0;Inf;0;0;0;0;0; ...
0;0;0;0;0;Inf;Inf;Inf;Inf;0; ...
0;0;Inf;Inf;Inf;0; ...
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
