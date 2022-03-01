function x = p0001_calib
x = struct('General', {struct('program_name', {struct('name', {['CALrp2_PU.m ' ...
]} ...
,'date', {['19-Oct-2018 ' ...
]} ...
,'version', {['rp2 v3.0 ' ...
]} ...
)} ...
,'picture_number', {1 } ...
,'date', {['19-Oct-2018 ' ...
]} ...
,'time', {['15:33:59 ' ...
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
[92.10246676;92.39315173;92.26309461;91.94375302;92.20405738;92.61354624;92.42623918;92.56332925;92.21180473;92.69690337; ...
92.39357029;92.68016065;92.56751116;92.39488669;92.78484913;92.72321125;92.91515305;92.66947646;93.00711626;92.96941827; ...
92.94302419;92.79190598;92.93374262;92.91390407;92.84479736;92.96635678;93.0025551;92.75962866;92.92881538;92.65088003; ...
92.64978565;92.43448215;92.44001951;92.35454547;92.20860921;92.06022369;91.79384322;91.62552097;91.77649157;91.56929903; ...
91.291858;91.25059928;91.11549183;90.85565537;90.78790276;90.76843334;90.73390356;90.68278801;90.59717839;90.41068403; ...
90.47157628;90.60691127;90.59336354;90.53971914;90.73035373;90.73595583;90.76560578;90.91853118;90.73770898;90.88621974; ...
90.69792257;90.87283872;90.72767501;90.81033609;90.92320936;91.14567246;91.13029882;91.10335517;91.34011926;91.30015924; ...
91.28603592;91.21883666;91.08152992;91.05159053;90.98985502;90.80546356;90.66845362;90.45653279;90.25547798;90.1973031; ...
90.11688515;90.0590276;89.8315312;89.8272643;89.58028271;89.59978779;89.52931393;89.33541114;89.34687858;89.29336773; ...
89.26807172;89.33154216;89.42144036;89.40231281;89.5136649;89.68524393;89.73665483;89.74776692;89.76891055;89.92947435; ...
89.97124983;90.06897068;90.23293541;90.38685001;90.65655711;90.90880739;91.20470511;91.54495174;91.90721621;92.29889974; ...
92.81255289;93.19295744;93.70281244;94.19900438;94.69966586;95.23768414;95.68348437;96.12036012;96.57159408;96.98670867; ...
97.37579948;97.61216905;97.76434401;97.76786631;97.6538;97.37945119;96.93874319;96.49395907;96.10272835;95.82436405; ...
95.47732998;95.02835178;94.51970386;94.20577326;93.8290956;93.38072283;92.98883694;92.57186148;92.42868245;92.45076726; ...
92.44461264;92.41444868;92.28684541;92.43217292;92.58703354;92.97997951;93.1717186;93.50683698;93.65035229;93.57409121; ...
93.34748562;92.57639367;91.33877952;91.59363375;93.82098968;95.86675394;97.31270541;98.1364704;98.82996641;99.34709479; ...
99.41022622;99.21635398;98.8460208;98.32306528;98.07345455;98.26294609;98.84597191;99.81929366;100.6322695;100.9531478; ...
100.9477968;100.7341418;100.7656463;101.0104652;101.4317982;102.1353476;102.3258964;101.6775128;100.8604567;100.1563068; ...
99.92974065;100.0436722;100.0138787;99.90714765;99.55743614;99.24061269;99.64402502;100.2102582;100.6782655;100.9793212; ...
100.8169967;101.0882642;101.3129563;101.5729253;101.6361663;101.5525664;101.321713;101.1898273;101.4465629;101.3700478; ...
100.6608094;99.71583758;99.07219958;99.30025617;100.2229519;100.9816531;101.1661375;101.4544742;101.9584172;102.8953738; ...
103.7845359;104.5956333;104.9735665;105.1674391;105.276571;104.796936;103.4175648;101.8174774;100.2300279;98.25379334; ...
96.20024079;94.85360156;94.07232074;93.36395824;94.27340978;95.6687537;96.21838686;96.62266922;97.42057194;97.0004267; ...
96.14185832;96.9115804;97.77599745;97.37215463;98.37636074;98.64619778;98.1971952;98.20672774;97.5682682;96.950276; ...
96.60536027;96.04040071;95.99492837;94.77667574;93.39327394;91.02738608;89.97368998;89.38801416;88.56210335;88.72113217; ...
90.02783871;92.65132392;95.27297426;95.3684639;96.40070856;96.29532253;91.3921266;86.10359343;82.03448459;78.93587644; ...
76.78966709;74.05647778;70.44656734;65.58745022;60.4456416;53.80149386; ...
] ...
[0.411148289;0.7521445546;0.513551799;-0.2066369843;-0.4277171376;0.3412280514;-0.9113452187;0.8739842458;0.1349251361;-0.6076281287; ...
-0.7953972837;-0.1027977989;-0.2773340894;0.3893615963;0.2292249947;0.869815793;0.7253115442;-0.6603856016;0.5708699456;-0.200104395; ...
-0.973151697;0.927279108;-0.5248064754;-0.002274067582;-0.7941399617;0.4125881037;0.8815871117;0.8226044508;-0.7361889659;0.453587444; ...
-0.3572375618;0.0283863622;-0.7969930532;0.3774848279;-0.4479327972;0.7270099543;-0.09742212122;-0.9209550805;-0.6500629178;-0.5641807029; ...
-0.3220116181;-0.2002832091;0.985125744;0.173283333;-0.6354654615;0.5592355613;-0.2422458902;0.9604696932;0.1677739214;-0.6199291873; ...
0.5977771074;-0.1786774104;0.8442828168;0.8240812363;0.8035265384;0.7826125507;0.7613329929;0.739681475;0.7176514952;-0.6095271237; ...
0.01728871854;0.5968962021;0.5024516114;-0.9920548022;-0.5371285208;-0.6863548988;-0.3112213259;0.009761942971;0.2751983077;0.4836588888; ...
0.633681826;0.3301239314;0.3868290235;0.3809965698;0.3110340498;-0.1030417812;-0.276095211;-0.735393265;0.9837656453;0.4765708028; ...
0.08304224756;-0.4748563468;-0.9861837096;0.4022685713;-0.2199688136;-0.9004762007;0.4195554527;-0.3918931908;0.8671095685;0.1131599243; ...
-0.6539685276;0.293768544;-0.5401216933;0.6114119249;-0.6448323113;0.4349522291;-0.5013483134;0.02534784888;-0.9880312205;-0.01912384396; ...
0.931760351;-0.8372455743;0.02941933631;0.8762738985;-0.2970281666;-0.3862908234;0.3484191006;-0.9389878564;-0.2488982906;-0.6849737782; ...
-0.09532227713;0.4696767664;-0.2592693944;0.1979857514;-0.7558638008;-0.4121671249;0.3990985173;0.623146042;-0.81087883;-0.7128594894; ...
0.2450637938;0.3588573302;0.3379387635;0.1787860317;-0.1222048964;-0.568721811;0.8354619125;0.08648544734;-0.8196013883;0.5018439319; ...
-0.6539685276;0.02204323654;0.5255752077;-0.4441006572;-0.2166200042;0.6958526344;0.6346476012;-0.5778374014;-0.9409102876;-0.2858137095; ...
-0.9645073378;-0.4496607657;0.5416773931;0.9080546895;-0.754123236;-0.2362587857;-0.05827649518;0.09273845823;-0.09991186325;-0.1307288895; ...
-0.6300495811;0.7237715615;0.7078031643;0.3079778355;0.9268960569;-0.6182506874;-0.1530919048;0.09393257202;-0.8433475394;0.859424817; ...
-0.3948964935;-0.4114266984;0.4782175988;-0.9272838098;-0.5865835029;-0.4644611808;-0.9835795321;-0.6973009343;0.4558135838;0.4452639697; ...
-0.99301746;-0.4773209053;-0.2286873269;0.6676980283;0.9974861493;0.0879125642;-0.1867211395;-0.6415424062;-0.5878803834;-0.8191934406; ...
0.3476048211;0.3255088515;-0.8872258775;0.4381369492;0.468647889;0.4005126034;-0.7637941586;0.05857165783;0.8028203759;-0.7398951127; ...
0.9596033512;0.3426293689;0.9245526823;-0.6080042218;-0.6096559146;0.7219992144;-0.4036059308;0.70561796;0.8863595173;0.6679978099; ...
0.5458690324;0.7177146603;0.05526516167;0.2908647608;-0.2932917515;0.4834085041;-0.6107154427;-0.1783882503;-0.2910302436;-0.8831522119; ...
-0.2628600979;-0.6070368916;-0.8124484769;0.8208883886;-0.3022338761;0.2616097271;-0.7307047977;0.6792287969;0.4472717699;0.4951320432; ...
-0.4112683444;0.2013569372;0.05883867262;-0.1135349137;0.468647889;0.5091608096;0.0818351411;0.2733344032;0.003525469931;-0.522915558; ...
-0.1615865951;-0.4196550164;-0.6035785415;0.5279620037;0.976544797;0.4343303219;-0.08350467533;-0.07101126396;0.9837656453;-0.08572992932; ...
-0.552855091;-0.2102774728;-0.7906397126;-0.2788469856;-0.5418618461;-0.2748872438;0.7347812002;0.6918835787;0.9873753097;-0.6568323633; ...
-0.4012726868;-0.9520983777;0.03410027693;-0.8882013145;-0.2871507718;-0.7819918941;0.1759346326;-0.9340488163;0.05089185885;0.873433123; ...
-0.0563417221;-0.09389378239;0.09806445436;0.2678345786;-0.7524506651;0.8743116832; ...
] ...
[Inf;Inf;Inf;Inf;0;0;Inf;0;Inf;Inf; ...
Inf;0;Inf;0;0;Inf;Inf;0;0;0; ...
0;Inf;0;0;0;0;0;0;0;0; ...
0;0;0;0;0;0;Inf;0;Inf;0; ...
0;0;Inf;0;0;0;0;Inf;0;0; ...
Inf;0;Inf;0;0;0;0;0;0;0; ...
0;0;Inf;0;0;Inf;0;Inf;0;0; ...
0;0;Inf;Inf;0;0;Inf;0;0;0; ...
0;0;0;Inf;0;0;0;0;Inf;0; ...
0;0;0;0;0;Inf;0;Inf;Inf;0; ...
0;0;0;0;0;0;0;0;0;0; ...
Inf;0;0;0;0;0;0;Inf;0;0; ...
0;Inf;0;Inf;0;Inf;0;0;0;0; ...
0;0;Inf;Inf;0;0;0;0;Inf;Inf; ...
0;0;0;0;0;Inf;0;Inf;0;0; ...
Inf;0;Inf;0;0;0;0;0;0;0; ...
0;0;0;0;Inf;Inf;0;0;Inf;Inf; ...
0;Inf;0;Inf;0;0;Inf;0;0;0; ...
0;0;0;0;0;0;0;0;0;Inf; ...
0;0;Inf;0;0;0;0;0;0;Inf; ...
0;0;0;0;0;0;Inf;0;0;0; ...
0;0;0;0;Inf;0;0;0;Inf;0; ...
0;0;Inf;0;0;0;0;0;Inf;0; ...
0;0;0;0;Inf;Inf;Inf;0;0;0; ...
Inf;0;0;Inf;0;0;Inf;Inf;0;0; ...
0;Inf;0;0;0;0;Inf;0;0;0; ...
0;Inf;0;0;0;0; ...
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
