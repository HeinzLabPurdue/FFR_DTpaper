function x = p0001_calib
x = struct('General', {struct('program_name', {struct('name', {['CALrp2_PU.m ' ...
]} ...
,'date', {['16-Jan-2019 ' ...
]} ...
,'version', {['rp2 v3.0 ' ...
]} ...
)} ...
,'picture_number', {1 } ...
,'date', {['16-Jan-2019 ' ...
]} ...
,'time', {['12:01:30 ' ...
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
[94.1980857;94.46094213;94.5669295;94.57626337;94.88191872;94.61447122;94.61457958;94.74099553;94.75360194;94.79102579; ...
94.71790769;94.91394387;95.16235949;94.98228768;94.91644681;95.02300273;95.04669395;95.13981903;94.88841554;94.91127837; ...
95.29218773;95.22459897;95.18086223;95.28295511;95.30316311;95.1852522;95.27969785;95.21106244;95.17860581;95.1421996; ...
95.29234737;95.14024761;95.09904829;94.96727008;94.98136802;94.85042353;94.69201004;94.74565864;94.52575261;94.27717651; ...
94.20574089;94.20167785;94.0148002;93.93823128;93.84658416;93.63034394;93.60457127;93.35728027;93.12391336;93.12986692; ...
93.19944193;93.17621033;93.16454108;93.20015778;93.25735121;93.10215242;93.08642687;93.03672422;93.09374484;93.0760431; ...
92.8925116;92.96961255;92.87125651;92.86212368;92.67009672;92.79503563;92.77109988;92.81230118;92.81665141;92.84510037; ...
92.87599099;92.8174125;92.77821624;92.70148597;92.58119879;92.42319115;92.38446328;92.22348966;92.12517925;92.15297414; ...
92.05374772;91.94452967;91.92712919;91.74440449;91.67926711;91.53064079;91.41391799;91.31136945;91.19195378;91.21522334; ...
91.09097563;91.15799605;91.14200933;91.20279367;91.16989012;91.16794806;91.1194162;91.14588504;91.17753325;91.2686801; ...
91.37082207;91.55825626;91.71379121;91.93140165;92.16909918;92.29716656;92.50348499;92.64424859;92.8053493;92.96048307; ...
93.14237528;93.34946776;93.56355818;93.82206564;93.97640103;94.21706016;94.37156462;94.47944962;94.61715812;94.69476246; ...
94.84249492;94.97596199;95.12161528;95.2843842;95.59993104;95.9521496;96.28487652;96.59210069;96.85166915;97.06179809; ...
97.1013501;97.09950735;97.01100946;96.7974896;96.51265813;96.05501153;95.58258695;95.11693438;94.60482478;94.17281014; ...
93.74592266;93.50152127;93.12046691;92.88753699;92.62818939;92.52170886;92.93404072;93.52054004;94.16054921;94.64842165; ...
95.07755525;95.36993569;95.51195649;95.5189677;95.40964663;95.19876655;94.86008379;94.30327124;93.35796506;92.7568653; ...
93.92348461;95.63271754;96.42897239;96.77680457;97.17574445;97.69078472;98.62340698;99.61373175;100.3874213;100.6233975; ...
100.502851;100.2628609;100.1306025;100.391104;100.7940667;101.2031676;101.4568133;101.1510965;100.5455842;100.1028217; ...
99.76524826;99.52036141;99.34552461;98.92153951;98.17192883;97.69010855;97.16682323;97.13250266;97.22697896;97.2248068; ...
96.90437006;96.74192937;96.6961899;96.82365543;96.99029572;97.07419899;97.03954594;97.41234169;97.95533477;99.07293786; ...
99.12488988;98.60680284;97.98854665;97.91605897;98.51126088;99.01714629;99.68610831;100.0022483;100.2632958;100.7044219; ...
101.0959231;101.3773009;101.5668893;101.7672515;102.0586107;102.7744526;102.8210616;102.7740495;102.5672452;101.8084152; ...
100.4554167;99.37802872;98.24725193;96.23141605;94.59923565;94.0532738;92.87543968;91.54662827;91.61509448;91.74052877; ...
91.19623408;91.82757016;92.56904338;92.24620523;92.61748363;93.57769277;93.61949782;93.22201073;92.85791903;92.10233287; ...
91.84865931;91.53229806;91.70587535;91.27597057;91.09748569;90.30793945;90.89498771;90.5748758;89.04013411;86.62752475; ...
84.8446613;84.20995988;86.6294793;88.60750616;89.80463474;90.43002143;88.72540366;85.25555344;82.08760041;78.72148504; ...
74.15242546;70.96241305;68.18740482;64.3415572;58.83284212;56.77067241; ...
] ...
[-0.6501176595;-0.8599983611;-0.5146097646;-0.7233981748;0.5783931073;-0.6350752065;0.6510407502;-0.06621842441;-0.7867526031;-0.9618716258; ...
-0.2372226771;-0.9669351452;0.3005245022;0.9773209781;0.8274617255;-0.5214904823;0.7253115442;-0.03022633904;-0.1467815496;-0.200104395; ...
-0.973151697;0.2518911332;-0.5248064754;0.6969314126;-0.08271258494;0.4125881037;-0.3818970793;0.8226044508;0.02630001704;-0.7705955004; ...
-0.3572375618;0.8315625471;0.02022240412;-0.7910150402;-0.4479327972;0.7270099543;-0.09742212122;-0.9209550805;0.2566944058;-0.5641807029; ...
-0.3220116181;0.8445739572;0.01328730763;-0.81554254;-0.6476859512;0.5118420771;-0.3254271821;0.8408750336;0.01112939979;0.9028289981; ...
0.8836507868;0.8641373463;0.8442828168;0.8240812363;0.8035265384;0.7826125507;-0.4773340142;0.2190444249;0.8706059808;-0.5238178091; ...
-0.6378521358;-0.1046556968;0.3792903199;0.2111232769;0.6171619722;0.9704676518;0.7427154381;-0.9877975713;-0.7747576239;0.9320247777; ...
-0.9450743716;-0.8825808084;0.7523981088;0.7179963059;0.6189650534;0.1753123773;-0.02783489178;-0.517752816;-0.8297490018;0.631356478; ...
0.2055741445;-0.2954276817;-0.9298605823;0.4246168253;-0.2444097929;-0.9478696849;0.2531928688;-0.5114878504;0.710465047;-0.2755240833; ...
0.8806346197;0.02204323654;0.837009574;-0.09226312986;0.9622207655;-0.4345975682;0.5439836581;-0.4952892012;0.4472717699;0.7618219087; ...
-0.3785213576;0.4612025268;-0.719354858;-0.7173699433;0.01155281955;0.7182575437;0.4562926286;-0.9341068849;-0.3488101539;-0.8915102226; ...
-0.4103470677;0.04426728689;-0.7969930532;-0.4540153041;0.4758602135;0.701249509;0.8956191558;-0.506292162;-0.4379081244;0.3095713504; ...
0.2450637938;0.3588573302;0.3379387635;0.1787860317;-0.09776391716;-0.4739348425;-0.9981755036;0.3256747664;-0.5063123452;0.5018439319; ...
-0.6539685276;0.5654938514;-0.8515560596;-0.4441006572;-0.2166200042;-0.1736971629;0.6346476012;0.3808884985;-0.07030430675;0.4951320432; ...
-0.2747890464;0.1472354365;0.5416773931;-0.6855891523;-0.4455422499;-0.2362587857;0.04959703286;0.09273845823;-0.09991186325;-0.4130728889; ...
-0.9450743716;0.7237715615;-0.1508946354;0.7399947227;-0.6096559146;-0.2060835625;-0.1670093507;0.964494368;0.713472702;-0.4277156936; ...
0.411148289;-0.4114266984;0.4782175988;0.3242184121;-0.5865835029;-0.4644611808;-0.9835795321;-0.6973009343;0.6013287056;-0.3321040454; ...
0.9066643683;-0.4773209053;-0.2286873269;-0.4441006572;0.9974861493;0.0879125642;-0.1867211395;0.2664877959;-0.5878803834;-0.8191934406; ...
-0.6205634172;-0.5427747014;-0.9362581047;-0.4949044061;0.468647889;0.09096734129;-0.7055176635;0.1561910875;0.6012339145;0.6954168884; ...
-0.3522311494;-0.9721118934;0.9245526823;-0.6080042218;-0.9020716869;0.5366653574;-0.05566978355;0.70561796;-0.9302103989;-0.7577211689; ...
0.7038276298;-0.2582827174;-0.03186444314;0.2908647608;-0.2932917515;0.4834085041;0.7055892572;0.518912684;0.4810629645;-0.8831522119; ...
-0.2628600979;0.3501483517;0.5255752077;-0.1470021908;0.2821230704;0.2630283708;-0.7360980511;-0.146253205;-0.9394019172;0.1140470972; ...
-0.2692976955;0.2013569372;0.009806445436;-0.4949044061;-0.06597177552;-0.5811657775;-0.1165529904;0.3709538329;0.003525469931;-0.522915558; ...
0.4975209452;0.8950862459;-0.5652088014;0.5279620037;0.6841290247;0.9073330715;-0.1948442424;-0.9941576954;-0.7792490243;0.1017036481; ...
0.8548001687;0.306287962;-0.4824494611;0.2908647608;0.2402494912;-0.2748872438;0.6836953;0.08648544734;0.9873753097;-0.6568323633; ...
-0.8025453736;-0.3476049189;-0.8124484769;-0.8820131449;0.005599492708;-0.7791546068;0.1651481256;0.03768639069;0.2618047791;-0.5924911793; ...
-0.6095807154;-0.4194026339;0.03922578175;0.5152867817;0.8118873337;0.6910954921; ...
] ...
[0;0;0;Inf;Inf;0;0;0;0;0; ...
Inf;0;0;Inf;0;0;0;0;0;0; ...
0;0;Inf;0;0;0;0;0;Inf;0; ...
0;0;0;0;0;0;0;Inf;0;0; ...
0;0;0;0;Inf;0;0;0;0;0; ...
0;Inf;0;0;0;0;Inf;Inf;0;Inf; ...
0;Inf;0;0;0;Inf;Inf;Inf;0;0; ...
0;0;0;0;0;0;0;0;0;Inf; ...
0;0;0;0;0;0;0;0;Inf;0; ...
0;0;Inf;0;0;0;0;0;0;0; ...
0;Inf;0;0;Inf;0;Inf;Inf;Inf;0; ...
Inf;0;0;0;0;0;0;0;Inf;0; ...
0;0;0;Inf;0;0;0;0;0;0; ...
Inf;0;0;Inf;0;0;0;0;0;0; ...
0;0;0;0;0;Inf;0;Inf;0;0; ...
Inf;Inf;0;0;0;Inf;0;0;0;Inf; ...
0;Inf;0;0;0;0;0;0;0;0; ...
0;0;0;0;Inf;Inf;0;0;Inf;0; ...
0;0;Inf;0;Inf;0;0;0;0;0; ...
0;0;0;Inf;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;Inf; ...
0;0;0;0;0;0;0;0;Inf;0; ...
0;Inf;Inf;Inf;0;0;0;0;Inf;Inf; ...
0;Inf;Inf;Inf;0;0;Inf;0;Inf;Inf; ...
0;Inf;Inf;0;Inf;0;Inf;0;0;Inf; ...
0;Inf;Inf;Inf;0;Inf;0;0;0;0; ...
Inf;0;0;0;0;0; ...
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
