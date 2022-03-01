function x = p0004_calib
x = struct('General', {struct('program_name', {struct('name', {['CALrp2_PU.m ' ...
]} ...
,'date', {['25-Sep-2018 ' ...
]} ...
,'version', {['rp2 v3.0 ' ...
]} ...
)} ...
,'picture_number', {4 } ...
,'date', {['25-Sep-2018 ' ...
]} ...
,'time', {['11:35:14 ' ...
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
[94.74436215;94.92819815;94.94569078;95.21292855;95.07904391;95.19574523;95.3380085;95.30221605;95.28050438;95.45486293; ...
95.37635029;95.48990867;95.40821295;95.6422346;95.61551869;95.57213049;95.67416365;95.53507198;95.73430646;95.62711058; ...
95.63640525;95.53079984;95.41657865;95.55838816;95.42666486;95.28944393;95.29022917;95.17521258;95.1984108;95.10276167; ...
94.91880637;94.83990619;94.80328724;94.72802119;94.63545702;94.51696623;94.43589477;94.20709579;94.24634865;94.19193075; ...
94.1129567;93.9666424;93.81002994;93.81082151;93.71706416;93.58662005;93.52110539;93.41521871;93.41241893;93.45923396; ...
93.46335745;93.40000987;93.33290817;93.18876755;93.11515801;92.98363097;92.86467654;92.74225361;92.56383706;92.4429253; ...
92.30836656;92.25286469;92.08707967;92.01071597;91.92604592;91.78063168;91.62130469;91.57110859;91.55689594;91.5191975; ...
91.38340382;91.31806189;91.27360359;91.14886297;91.08030434;91.01176091;90.84889092;90.68201899;90.60319292;90.44044443; ...
90.32724924;90.32169537;90.1713774;90.09663693;89.93140846;90.05867764;90.01926992;90.01517786;90.03346918;90.06025441; ...
90.11809063;90.22213239;90.28268962;90.43779657;90.49337249;90.53747706;90.65114009;90.73893197;90.84007493;90.9319451; ...
91.0436787;91.20436164;91.4312953;91.62217734;91.82532418;92.07521712;92.38484862;92.57040739;92.84514442;92.96188926; ...
93.1201207;93.21226616;93.27257486;93.32693949;93.36247047;93.38120874;93.48047949;93.66012318;93.85868472;94.0781255; ...
94.43039627;94.71579264;94.89801452;95.08814161;95.29512246;95.42973355;95.45143817;95.4410872;95.44653968;95.90618294; ...
96.77491314;97.24195426;97.28450584;97.01120224;96.56383248;96.03310789;95.63536457;95.32452957;94.92934336;94.74157729; ...
94.64108605;94.99935983;95.34782133;95.69966064;95.88206851;96.14543081;96.22587821;96.01188327;95.8831589;95.5903522; ...
95.28964115;94.9162638;94.02986284;92.46271722;90.73653966;89.64401199;90.63896249;93.07582992;95.32924169;97.32462946; ...
98.80403192;99.81555693;100.4360298;100.8711481;101.3117471;101.9485637;102.7868356;103.5534164;104.2905175;104.5425605; ...
104.4601309;104.330662;104.271193;104.4542879;104.6931371;104.892424;104.8660883;104.5956675;103.9479548;103.5737206; ...
103.343855;103.1238659;102.8308222;102.2692627;101.2658598;100.3054862;99.7344522;99.38531499;99.26292929;98.853656; ...
98.43961508;97.86466796;98.00541765;98.29710133;98.49726676;98.24158285;97.79465241;97.50715081;98.01437465;99.12261961; ...
99.21881546;98.73504632;98.46706558;98.55796408;99.1839838;99.68649298;100.3645459;100.8052939;101.2613665;101.7670866; ...
102.4686437;102.9939605;103.3386912;103.8055353;104.0565183;105.0390656;105.1038715;105.1989722;105.3930426;104.9212318; ...
104.4025757;104.3205984;104.0032683;102.874234;102.0658212;101.7241657;100.465502;98.7552275;97.99819976;97.40030495; ...
96.02868026;95.73659082;95.49096561;94.49730551;94.53979562;94.85508387;93.97721113;93.99837592;94.06085609;93.82501984; ...
93.92109809;93.53950309;93.64378381;93.13227277;92.88986726;92.60514318;92.2534074;91.899308;90.81374029;89.62735419; ...
88.72627162;87.74718355;87.30899686;87.28381433;86.87591299;84.36545203;80.07409927;76.74441882;76.29002055;72.483543; ...
68.27087771;67.09680905;65.79881346;63.48926201;59.1712522;55.89176567; ...
] ...
[0.8805153147;-0.2927126117;-0.5146097646;-0.7233981748;-0.4277171376;-0.6350752065;0.6510407502;-0.5961170893;0.6740862665;-0.05904262774; ...
-0.7953972837;-0.1027977989;-0.855192681;0.3893615963;-0.3690117361;0.2611220684;-0.5133554629;0.7094551359;-0.7114785591;0.4951320432; ...
-0.3007221242;0.9011151838;0.1008064275;-0.7014795478;-0.2169947151;0.9648622872;0.1450713027;-0.6761753063;-0.2611669317;0.9019533329; ...
0.06400624059;-0.7747898227;0.3857914894;-0.4540153041;0.7060327046;-0.1338129664;-0.9732919616;-0.7033146315;0.4431797586;-0.4093950277; ...
0.7392543303;-0.1105688765;-0.9585511287;0.195631587;-0.6476859512;0.5118420771;-0.3254271821;-0.2189222962;0.9216777392;0.9028289981; ...
0.8836507868;0.8641373463;0.8442828168;0.8240812363;-0.3929469232;0.347837652;-0.9546680285;-0.3015926251;-0.411742524;0.171418629; ...
0.7070070099;-0.8062075957;-0.9950967773;-0.585698644;-0.8056927812;-0.4772581647;-0.2033477979;-0.4865773284;-0.2997355898;-0.6196093334; ...
-0.5238305693;-0.4889331783;-0.882032806;-0.9450039579;0.6189650534;0.4536665357;-0.02783489178;-0.300112367;-0.8297490018;0.631356478; ...
0.2055741445;-0.3851420143;-0.9298605823;0.4246168253;-0.2321893032;-0.9004762007;0.3363741607;-0.3918931908;0.8671095685;-0.08118207948; ...
-0.886666954;0.293768544;-0.8515560596;0.2595743975;-0.6448323113;0.0001773304563;-0.9786823276;0.02534784888;0.4472717699;-0.6286509676; ...
0.2766194967;-0.8372455743;-0.719354858;0.07945197758;0.8572623265;0.7182575437;-0.5976441354;0.06345262931;-0.3488101539;0.2117579996; ...
-0.4103470677;0.04426728689;-0.7969930532;-0.4540153041;-0.1400017936;0.701249509;0.8956191558;-0.94157306;-0.81087883;-0.403288139; ...
0.2450637938;0.1794286651;0.225292509;0.1340895238;-0.09776391716;-0.4739348425;-0.9981755036;0.5648640855;-0.193023302;-0.720788053; ...
0.2768251779;-0.8910555337;-0.2286873269;0.2595743975;0.5692738423;0.6958526344;-0.4106843703;-0.5778374014;-0.9409102876;-0.2858137095; ...
0.3457743708;-0.4496607657;0.03922578175;0.9080546895;-0.4455422499;-0.2362587857;0.04959703286;0.004880971486;-0.09991186325;-0.4130728889; ...
-0.9450743716;0.298362082;-0.6886182942;-0.9960242754;-0.6096559146;0.6804998036;-0.1530919048;-0.776629224;0.3457736306;0.2402821163; ...
-0.07897929871;-0.7702840285;0.2529250898;0.9833231744;0.729069077;-0.2748872438;-0.3254271821;-0.6973009343;-0.6265780863;0.4452639697; ...
-0.99301746;-0.4773209053;-0.2286873269;-0.4441006572;0.9974861493;0.08720324237;-0.1867211395;0.2664877959;0.6133263129;-0.8191934406; ...
-0.6205634172;0.3875848086;-0.9852903318;-0.6220275702;-0.3993309716;0.6367713891;-0.7055176635;0.0683336008;0.6029966494;-0.8706240022; ...
0.3295537701;0.05577621313;-0.7544731769;-0.6080042218;-0.6096559146;0.7219992144;-0.2922663637;0.70561796;0.951296936;-0.7577211689; ...
0.7038276298;-0.2582827174;0.05526516167;0.2908647608;-0.2932917515;0.4834085041;0.7055892572;0.518912684;0.4810629645;0.9052793944; ...
-0.06222375444;-0.3476049189;0.5255752077;0.9647964947;0.2821230704;-0.8677764928;-0.7360980511;-0.146253205;-0.1406086135;0.1140470972; ...
-0.6169025166;0.1392809801;0.05883867262;-0.4949044061;-0.2375914979;-0.03536172975;0.1983881315;0.1952388594;0.007050939863;0.3477084442; ...
-0.4620824061;0.3868966557;0.7186877619;-0.9440759926;-0.7310394307;0.8049980361;0.997943896;-0.07101126396;0.9837656453;-0.08572992932; ...
-0.552855091;-0.2102774728;-0.8777693174;-0.212139683;0.08944331376;-0.2748872438;0.6763973143;0.1621361834;-0.5568111065;-0.8462735745; ...
-0.4634964413;0.1763458923;-0.6419470923;0.382975901;0.7910644957;0.3488129695;-0.9120326837;-0.9528920116;-0.3515215338;-0.01947182723; ...
-0.7515513644;-0.2180456966;0.03922578175;0.2610404534;0.468647889;0.510442318; ...
] ...
[0;0;0;0;0;0;0;0;0;0; ...
0;0;0;Inf;0;0;0;0;0;Inf; ...
0;Inf;Inf;0;Inf;0;0;0;Inf;0; ...
Inf;0;Inf;0;0;0;0;0;0;0; ...
Inf;0;0;0;0;0;0;0;Inf;0; ...
0;0;Inf;0;Inf;0;Inf;0;Inf;Inf; ...
Inf;Inf;Inf;Inf;0;Inf;0;0;0;0; ...
Inf;Inf;Inf;0;0;0;0;0;0;0; ...
0;0;Inf;0;0;0;0;0;0;0; ...
0;0;0;0;0;Inf;0;0;Inf;Inf; ...
0;0;0;0;0;0;Inf;Inf;0;Inf; ...
0;Inf;0;0;0;0;Inf;0;0;Inf; ...
0;0;0;0;0;0;0;0;Inf;Inf; ...
Inf;0;0;0;0;0;0;0;0;Inf; ...
0;0;Inf;0;0;Inf;0;Inf;0;0; ...
Inf;0;0;0;0;0;0;0;0;0; ...
Inf;Inf;Inf;0;0;0;0;0;0;0; ...
0;0;Inf;Inf;Inf;0;Inf;0;0;Inf; ...
0;0;0;Inf;0;Inf;0;0;Inf;0; ...
Inf;0;0;0;Inf;0;0;0;0;0; ...
0;0;0;Inf;0;0;0;0;Inf;0; ...
0;0;0;0;0;0;Inf;0;0;0; ...
Inf;0;0;0;0;0;0;Inf;Inf;Inf; ...
0;0;0;0;Inf;0;0;0;0;0; ...
0;0;0;0;Inf;0;0;Inf;0;Inf; ...
0;Inf;0;0;0;0;0;0;0;Inf; ...
0;0;Inf;0;0;Inf; ...
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
