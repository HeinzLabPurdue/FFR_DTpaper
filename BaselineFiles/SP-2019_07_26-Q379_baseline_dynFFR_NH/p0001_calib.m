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
,'time', {['19:50:32 ' ...
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
[98.58799755;98.60137753;98.66735805;98.64761526;98.71709211;98.87415425;98.7512006;98.78670632;98.91187416;98.89933157; ...
98.96206075;98.99186259;98.80752855;98.78033647;98.86399994;98.92297438;98.83960398;98.68060712;98.64576912;98.70177743; ...
98.4912923;98.41247035;98.34165504;98.30429889;98.22676788;98.11074092;97.87370222;97.83697177;97.66520693;97.54237647; ...
97.47065158;97.42641762;97.27661396;97.25917624;97.19785904;97.11542909;97.20321887;97.10784108;97.09837111;97.02659805; ...
97.00963868;96.95950789;96.96663818;96.88509926;96.7136781;96.65368983;96.51295954;96.36850915;96.31132019;96.14669302; ...
96.02529542;95.88818373;95.68154505;95.53345241;95.38356057;95.30293936;95.24876275;95.17702465;95.10815868;95.0450141; ...
95.01152685;95.03150644;95.07614581;95.09640239;95.09644167;95.23335086;95.39647156;95.53515082;95.6784261;95.87681511; ...
96.12321357;96.36397082;96.57276809;96.68405155;96.82239712;97.00020797;97.08007524;97.16062245;97.13263022;97.14696505; ...
97.18407719;97.21544923;97.26283889;97.30438823;97.41096886;97.5054605;97.68613002;97.83715906;97.99659204;98.17157581; ...
98.35982387;98.56049634;98.72289658;98.85249549;99.01461665;99.18765176;99.36172253;99.55759564;99.75787258;99.96867804; ...
100.1727378;100.3894379;100.6012554;100.8232113;101.0616034;101.262169;101.5220526;101.7578141;102.0119601;102.2279803; ...
102.484444;102.734;102.9274344;103.2085299;103.4590293;103.7614338;104.0612382;104.3319744;104.6334437;104.8437056; ...
105.1390633;105.3095609;105.5439573;105.7737123;105.9988533;106.2330339;106.4471422;106.658516;106.7830774;106.8926065; ...
106.9335001;106.8985773;106.7635424;106.6023601;106.4477807;106.2085702;105.9943034;105.7911587;105.380378;105.0442351; ...
104.6833089;104.7821538;105.0708192;105.5430041;105.7974656;106.0198829;105.7322989;104.9939453;104.1859631;104.3444549; ...
105.5346635;107.1317702;108.2663154;109.0311645;109.5907672;110.0567859;110.5000514;110.9188125;111.5625175;112.1727451; ...
112.6599444;112.8534077;112.7770035;112.7688933;112.9060907;113.0451934;113.5468389;114.2101699;114.7138253;114.8095323; ...
114.8278216;114.7032711;114.757946;115.0240343;115.3503048;115.7957585;116.0660671;115.9039838;115.664179;115.3975399; ...
115.4008919;115.6212834;115.8401085;115.7646743;115.3937771;115.1913985;115.0555137;114.9829454;115.1992626;115.184655; ...
115.1017781;114.8933676;115.0292868;115.5587563;115.9729456;115.8379639;115.5740618;115.5031977;115.4410878;115.7050815; ...
115.8084935;115.3627928;115.0480945;114.6962544;114.5241526;114.467335;114.2071893;113.5366608;113.0811377;112.6769669; ...
112.4516653;112.0779935;111.3822932;110.7307037;110.7874157;110.5554699;109.795925;109.0738997;108.7079326;107.950832; ...
106.5589473;105.6947014;104.7266008;103.1893499;101.76667;101.7194335;100.4732348;98.98574181;98.75104171;98.50072565; ...
97.32209429;97.12772712;97.18846772;96.32825476;96.5023853;97.16231572;97.05036151;96.48259522;96.73685686;96.2087931; ...
95.46698;93.70465126;91.64474559;90.16507942;88.27572838;85.24392896;80.77894807;79.69183894;80.04803298;79.00773675; ...
75.8525289;73.86108438;81.15678883;85.32361072;85.80515508;86.09472666;85.98493663;79.99242723;78.09542758;81.13335341; ...
76.36907576;63.73217929;76.88751849;76.46528113;71.62751104;74.3572277; ...
] ...
[-0.6501176595;0.6624302221;-0.5146097646;-0.7233981748;0.07533798485;-0.6350752065;0.6510407502;-0.5961170893;0.6740862665;-0.05904262774; ...
-0.7953972837;-0.1027977989;-0.2773340894;0.9773209781;-0.3690117361;0.869815793;0.1059780406;-0.6603856016;0.5708699456;-0.8524861759; ...
0.3630630894;-0.4234968415;0.7879999761;-0.002274067582;-0.7941399617;0.4125881037;-0.3818970793;0.8226044508;0.02630001704;-0.7705955004; ...
-0.7784813641;-0.3652612679;0.8374378615;0.8719852237;0.09017069746;-0.6905212834;0.5301873999;-0.247415754;0.9769663766;0.2036407844; ...
-0.5670754119;0.6651452921;-0.09935894686;-0.8602390479;0.382865273;-0.3696742123;-0.1590645982;-0.9199356473;0.3244184429;-0.4255871836; ...
0.8304755338;0.09304789706;-0.6374152183;0.6395549766;-0.07555846902;0.7826125507;0.7613329929;0.739681475;0.7176514952;0.6952364382; ...
0.6724295728;0.6492240505;0.6256129028;0.6015890396;-0.8457095069;-0.8954516329;-0.946063236;-0.4963392714;-0.5749338974;-0.2065364445; ...
0.1062190118;-0.03176184939;0.1934145117;0.358998153;0.1555170249;0.2268332679;0.2343428734;-0.04123595899;-0.1351464717;-0.2973575731; ...
-0.6521491339;-0.9234280095;0.6758775269;0.2905273015;-0.1588663654;-0.6161152952;0.8354619125;0.3256747664;-0.3496678236;-0.9151300568; ...
0.2768251779;-0.3476049189;0.7056157721;0.01876203445;0.9269553817;0.1740518238;0.9306537294;0.1078960491;0.7060598083;-0.190542473; ...
0.8971829139;-0.7325898775;0.2757419193;0.4699177403;-0.6056091527;-0.5953875575;0.2405455725;-0.9438688279;0.801057641;-0.4784373337; ...
-0.9378098819;-0.3176184939;-0.990407565;-0.4760137209;0.01241218501;-0.9688754419;-0.5939427597;0.187865144;0.4431797586;0.9775691602; ...
-0.8931499789;-0.7965687126;0.1126462545;0.04469650793;-0.04888195858;-0.2843609055;-0.6654503357;0.8040534046;0.1202657411;-0.720788053; ...
0.2768251779;-0.3476049189;0.3941814058;0.9632494523;-0.6448323113;-0.4345975682;0.5439836581;0.4634366987;0.1884837316;0.9332405378; ...
0.82985849;0.9534430322;-0.4632258296;0.5016985313;0.9372957779;-0.4453555198;-0.05827649518;0.08785748674;0.3032610597;0.3581515544; ...
0.01452897958;-0.4254094795;0.9245526823;0.04399683364;0.9268960569;-0.4329168303;-0.1530919048;-0.776629224;-0.7784101207;0.9531416057; ...
0.411148289;-0.4114266984;0.4782175988;0.3242184121;-0.5865835029;-0.4644611808;-0.9835795321;-0.6973009343;0.4558135838;0.4452639697; ...
-0.99301746;-0.4773209053;-0.2286873269;-0.4441006572;0.9974861493;0.0879125642;-0.1867211395;-0.6415424062;-0.5878803834;-0.8191934406; ...
-0.4785927683;0.3255088515;-0.9362581047;-0.4949044061;0.468647889;0.09096734129;-0.7055176635;0.1561910875;0.6029966494;-0.8706240022; ...
0.3295537701;-0.850818959;-0.1125248952;-0.34402322;-0.4386236584;0.7634986252;-0.06958722944;0.446741552;0.4431797586;0.4805642325; ...
-0.3159171948;-0.5645706793;0.5058501796;-0.7818514294;0.08944331376;0.1042606302;0.007297985733;-0.6973009343;0.7468438274;0.6789595458; ...
-0.9549230427;0.2619754056;0.6960765923;0.00618816956;-0.8614523155;0.5218008106;0.3626557721;-0.3113496054;0.5430220061;-0.762169892; ...
0.1419706489;0.3255088515;0.06864511805;-0.615233445;0.3696902257;-0.4725175713;0.6707998142;0.3904777189;-0.3961219831;0.2169795547; ...
-0.1034706768;-0.8065516721;0.6419482816;-0.5120591053;-0.1462078861;-0.8243342498;-0.6680374026;-0.7901132627;0.5941411332;0.7657080283; ...
0.9732691168;-0.2582827174;0.2444613652;0.7151441268;0.6760268167;0.6919297644;0.0437879144;0.4324272367;-0.5568111065;-0.6568323633; ...
0.429202847;0.005086865678;-0.8124484769;-0.4348184029;-0.4580926199;0.09004052967;-0.01078650695;-0.9905784023;0.7539349264;0.4951320432; ...
0.4259119468;-0.5738126799;0.1372902361;-0.983014687;-0.7194647773;0.8743116832; ...
] ...
[0;Inf;Inf;0;0;0;0;Inf;0;0; ...
0;0;0;Inf;0;0;0;0;Inf;0; ...
Inf;0;Inf;Inf;0;Inf;Inf;Inf;Inf;Inf; ...
0;0;0;0;0;0;0;0;Inf;0; ...
0;Inf;0;Inf;0;Inf;Inf;0;0;0; ...
0;0;0;0;0;0;0;Inf;0;0; ...
0;Inf;0;0;0;Inf;0;0;0;0; ...
Inf;Inf;0;Inf;0;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;Inf;0;Inf;Inf;0;Inf;0;0; ...
Inf;0;0;0;0;0;0;0;0;Inf; ...
0;0;0;Inf;Inf;0;0;0;0;0; ...
0;0;0;0;0;0;0;0;0;0; ...
0;0;Inf;0;0;0;0;0;0;0; ...
Inf;0;Inf;0;Inf;0;0;0;Inf;0; ...
Inf;0;Inf;0;Inf;0;0;0;0;0; ...
Inf;0;0;0;0;0;0;0;0;Inf; ...
0;0;0;Inf;0;0;0;Inf;Inf;Inf; ...
Inf;0;Inf;Inf;0;0;Inf;0;Inf;0; ...
0;0;Inf;Inf;0;Inf;0;0;0;0; ...
0;0;0;Inf;0;0;0;Inf;0;0; ...
0;0;Inf;0;0;0;0;0;0;Inf; ...
0;0;0;Inf;0;0;0;Inf;Inf;0; ...
0;Inf;0;Inf;0;Inf;0;0;Inf;0; ...
Inf;Inf;0;0;0;0;Inf;Inf;0;0; ...
0;0;0;0;0;0;0;Inf;0;0; ...
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
