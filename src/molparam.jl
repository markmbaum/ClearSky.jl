const TMIN = 25.000000
const TMAX = 1000.000000
const MOLPARAM = [
  #1, molecule number
  [1,
    #2, molecule formula
    "H2O",
    #3, molecule name
    "Water",
    #4, global isotopologue numbers
    Int64[1, 2, 3, 4, 5, 6, 129],
    #5, isotopologue formulae
    String["H216O", "H218O", "H217O", "HD16O", "HD18O", "HD17O", "D216O"],
    #6, AFGL code
    Int64[161, 181, 171, 162, 182, 172, 262],
    #7, abundance fractions
    Float64[0.997317, 0.002, 0.000371884, 0.000310693, 6.23003e-07, 1.15853e-07, 2.4197e-08],
    #8, molecular masses (kg/mole)
    Float64[0.018010565, 0.020014811, 0.01901478, 0.019016739999999997, 0.021020985, 0.020020956000000003, 0.020022915000000002],
    #9, Qref
    Float64[174.58, 176.05, 1052.14, 864.74, 875.57, 5226.79, 1027.8],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[7, 7, 7, 8, 8, 8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0028, 0.0028, 0.0029, 0.0094, 0.0094, 0.0096, 0.0091],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[2.8854493216841015, 3.4632621262928915, 0.6001640690631624, 0.016178811904860996, 0.015211044164822699, -0.004364096635145624, 0.0012395758566148274],
      Float64[2.8871644108577477, 3.465944358793127, 0.6012775643496776, 0.01626203386609267, 0.015073634800670662, -0.004421502003156756, 0.0012358242435368538],
      Float64[2.8843598005839284, 3.4613544216657774, 0.5988802292683817, 0.015592997489868013, 0.015071636079156958, -0.0043773144167319105, 0.0012344781178038982],
      Float64[2.988927843477164, 3.630916093351383, 0.6832681181788238, 0.03771653748442871, 0.01698698700299332, -0.0046881118088947715, 0.0017425144860104983, -0.00048799452244198605],
      Float64[2.995264646453848, 3.641366946322687, 0.6888505901076387, 0.03935511770905783, 0.01702953291306411, -0.004802882155796924, 0.0017375357666057514, -0.0004781891512090551],
      Float64[3.001991615535067, 3.651228173736914, 0.6917712772263419, 0.03877308927826554, 0.016728580596165613, -0.004750977994398704, 0.0017126280124818902, -0.00047877803465772625],
      Float64[3.1511970042488793, 3.8905643121864615, 0.8073563192325093, 0.06660577751477327, 0.01805683243327309, -0.004508003223199252, 0.0020143284208121565, -0.0006707506150274156]
    ]
  ],
  #1, molecule number
  [2,
    #2, molecule formula
    "CO2",
    #3, molecule name
    "Carbon Dioxide",
    #4, global isotopologue numbers
    Int64[7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122],
    #5, isotopologue formulae
    String["12C16O2", "13C16O2", "16O12C18O", "16O12C17O", "16O13C18O", "16O13C17O", "12C18O2", "17O12C18O", "12C17O2", "13C18O2", "18O13C17O", "13C17O2"],
    #6, AFGL code
    Int64[626, 636, 628, 627, 638, 637, 828, 827, 727, 838, 837, 737],
    #7, abundance fractions
    Float64[0.984204, 0.011057, 0.003947, 0.000733989, 4.43446e-05, 8.24623e-06, 3.95734e-06, 1.4717999999999998e-06, 1.36847e-07, 4.4459999999999996e-08, 1.65354e-08, 1.5375000000000001e-09],
    #8, molecular masses (kg/mole)
    Float64[0.04398983, 0.044993185, 0.045994076, 0.044994044999999996, 0.046997431, 0.0459974, 0.047998322, 0.046998291000000005, 0.045998262, 0.049001675, 0.048001646, 0.047001618237800004],
    #9, Qref
    Float64[286.09, 576.64, 607.81, 3542.61, 1225.46, 7141.32, 323.42, 3766.58, 10971.57, 652.24, 7595.04, 22120.47],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true, true, true, true, true, true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[7, 6, 7, 7, 6, 6, 7, 7, 7, 6, 6, 6],
    #12, maximum relative errors of interpolation
    Float64[0.0083, 0.0097, 0.0084, 0.0084, 0.0097, 0.0097, 0.0085, 0.0085, 0.0084, 0.0097, 0.0097, 0.0097],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[3.629868764336615, 4.6480808871061505, 1.3583273380857648, 0.2691914939179858, 0.012300344743387518, 0.004326678712311664, -0.0004876212718173771],
      Float64[3.7179211688304883, 4.786969447800756, 1.420400066599882, 0.2814003436457875, 0.012916486545213601, 0.005068159441358943],
      Float64[3.679005578420316, 4.726613084391979, 1.3953922253638875, 0.27781440229728993, 0.01314456082627761, 0.004419817566995832, -0.0004970214933885941],
      Float64[3.655363377445139, 4.688828648648713, 1.3775475553307868, 0.2736558449760584, 0.012739727465199616, 0.0043752016281623325, -0.0004930857409695122],
      Float64[3.770167807416816, 4.870520321377293, 1.45990358750975, 0.2906789686290149, 0.013857386918618885, 0.005159901644174525],
      Float64[3.7449516621478636, 4.830204291937045, 1.4408492290337949, 0.28620850769525374, 0.013406177679559405, 0.005114771601030554],
      Float64[3.7311033311460196, 4.809810916666548, 1.4345733694792095, 0.28689967884368855, 0.014046000533340042, 0.004524547647267359, -0.0005104693856422907],
      Float64[3.705985669094004, 4.769708282307593, 1.4156946951706022, 0.2825250989881978, 0.013613061747109967, 0.0044736845305030455, -0.000503988378637151],
      Float64[3.6816588053606285, 4.730840055680833, 1.3973606497267912, 0.2782602087069206, 0.013191466220107065, 0.004426317951875092, -0.000498417536415848],
      Float64[3.8253535330410573, 4.958741134066428, 1.5015645373059177, 0.30044270685108937, 0.014856805465608858, 0.005263210781243188],
      Float64[3.798742219093048, 4.916205675937364, 1.4814840206725617, 0.29574048323158914, 0.014376502287099413, 0.0052127748195452735],
      Float64[3.772887856184446, 4.874864915118713, 1.461940449659132, 0.29114887205345086, 0.013907974624242314, 0.005165919678550068]
    ]
  ],
  #1, molecule number
  [3,
    #2, molecule formula
    "O3",
    #3, molecule name
    "Ozone",
    #4, global isotopologue numbers
    Int64[16, 17, 18, 19, 20],
    #5, isotopologue formulae
    String["16O3", "16O16O18O", "16O18O16O", "16O16O17O", "16O17O16O"],
    #6, AFGL code
    Int64[666, 668, 686, 667, 676],
    #7, abundance fractions
    Float64[0.992901, 0.003982, 0.001991, 0.000740475, 0.00037023700000000004],
    #8, molecular masses (kg/mole)
    Float64[0.047984744999999995, 0.049988990999999997, 0.049988990999999997, 0.04898896, 0.04898896],
    #9, Qref
    Float64[3483.71, 7465.68, 3647.08, 43330.85, 21404.96],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8, 8, 8, 8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0081, 0.0085, 0.0082, 0.0083, 0.0081],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[5.198813326073512, 7.183627282058914, 2.4304058215844266, 0.47946442974793274, 0.04864953304330084, -0.003877687472476156, 0.0035095252917969333, -0.001342048768093613],
      Float64[5.289447393725946, 7.328403504739264, 2.499786028390479, 0.4961019241043351, 0.05020868215804548, -0.0036478605934218778, 0.003467093548097568, -0.0013615638674566405],
      Float64[5.32163909580045, 7.380784814348417, 2.527176975966621, 0.5040656085549876, 0.05090642322592535, -0.0036226978332624276, 0.0035401534819499147, -0.0013961748984979547],
      Float64[5.24568623397049, 7.258502442396174, 2.46630032008053, 0.488079931032065, 0.049452989319160166, -0.0037593580768103658, 0.003494383612559509, -0.0013477264983363974],
      Float64[5.262065893189635, 7.285163964911438, 2.4802630772444894, 0.4921466852101934, 0.0498022530511696, -0.0037485193119807087, 0.003527212690246679, -0.0013719459179242222]
    ]
  ],
  #1, molecule number
  [4,
    #2, molecule formula
    "N2O",
    #3, molecule name
    "Nitrogen oxide",
    #4, global isotopologue numbers
    Int64[21, 22, 23, 24, 25],
    #5, isotopologue formulae
    String["14N216O", "14N15N16O", "15N14N16O", "14N218O", "14N217O"],
    #6, AFGL code
    Int64[446, 456, 546, 448, 447],
    #7, abundance fractions
    Float64[0.990333, 0.003641, 0.003641, 0.001986, 0.00036928000000000004],
    #8, molecular masses (kg/mole)
    Float64[0.044001062, 0.044998095999999994, 0.044998095999999994, 0.046005308, 0.045005277999999996],
    #9, Qref
    Float64[4984.9, 3362.01, 3458.58, 5314.74, 30971.79],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[6, 6, 6, 6, 6],
    #12, maximum relative errors of interpolation
    Float64[0.007, 0.006, 0.0067, 0.0069, 0.0069],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.063103273197024, 5.333914221331409, 1.6639844097596002, 0.33067143798925497, 0.018935088814124867, 0.005959907290501931],
      Float64[4.157765791743768, 5.4838751764969, 1.7297047842704245, 0.34323622510507956, 0.020258240973265628, 0.005979204910197211],
      Float64[4.108724068330206, 5.406344295333827, 1.6961335382408997, 0.33687644580232573, 0.019294311268007645, 0.005866268561468502],
      Float64[4.132106894555242, 5.44362298720071, 1.7136070490521895, 0.34075869422226673, 0.019555787367178824, 0.005909433476340098],
      Float64[4.103738130527359, 5.3982340580606385, 1.6921417741556322, 0.3357097168910115, 0.01903526376545912, 0.005863689357000013]
    ]
  ],
  #1, molecule number
  [5,
    #2, molecule formula
    "CO",
    #3, molecule name
    "Carbon Monoxide",
    #4, global isotopologue numbers
    Int64[26, 27, 28, 29, 30, 31],
    #5, isotopologue formulae
    String["12C16O", "13C16O", "12C18O", "12C17O", "13C18O", "13C17O"],
    #6, AFGL code
    Int64[26, 36, 28, 27, 38, 37],
    #7, abundance fractions
    Float64[0.986544, 0.011084, 0.001978, 0.000367867, 2.2225000000000005e-05, 4.13292e-06],
    #8, molecular masses (kg/mole)
    Float64[0.027994915, 0.028998270000000003, 0.029999161, 0.02899913, 0.031002516, 0.030002485],
    #9, Qref
    Float64[107.42, 224.69, 112.77, 661.17, 236.44, 1384.66],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4, 4, 4, 4, 4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.0066, 0.0066, 0.0066, 0.0066, 0.0067, 0.0067],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.7730339354052334, 1.7139901664342414, 0.04076391415723165, 0.012492844178307502],
      Float64[1.7765934677858093, 1.7196928711145105, 0.04337971749168412, 0.013094903152256995],
      Float64[1.7769442249020646, 1.7202283022733258, 0.043607118540340384, 0.013146398444950247],
      Float64[1.7750240019453913, 1.7171872823767103, 0.04223776885485675, 0.012834504758217532],
      Float64[1.7808256851317825, 1.7264850221557417, 0.04648621217831724, 0.013784369656134091],
      Float64[1.7787537718136903, 1.7231790283878674, 0.04498907480569617, 0.013455416168919024]
    ]
  ],
  #1, molecule number
  [6,
    #2, molecule formula
    "CH4",
    #3, molecule name
    "Methane",
    #4, global isotopologue numbers
    Int64[32, 33, 34, 35],
    #5, isotopologue formulae
    String["12CH4", "13CH4", "12CH3D", "13CH3D"],
    #6, AFGL code
    Int64[211, 311, 212, 312],
    #7, abundance fractions
    Float64[0.988274, 0.011103, 0.0006157510000000001, 6.917849999999999e-06],
    #8, molecular masses (kg/mole)
    Float64[0.016031300000000002, 0.017034655, 0.017037475, 0.01804083],
    #9, Qref
    Float64[590.48, 1180.82, 4794.73, 9599.16],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[11, 11, 8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0073, 0.0073, 0.0069, 0.0069],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.559320057536849, 6.246154790519078, 2.15851344284554, 0.563745455962757, 0.12015959690385385, 0.0061053536914517535, 0.004893088743789864, -0.0004958586065243687, 6.363692773234675e-05, -0.0001658264468598958, 0.0001301035427090369],
      Float64[4.5276900460578515, 6.192333563178976, 2.125650980876807, 0.549549228804292, 0.11577307616779402, 0.005066519270059899, 0.004703210444850381, -0.0005103723347237743, 6.382180973449892e-05, -0.00016759358231226428, 0.0001297215051426548],
      Float64[5.03072658646721, 7.023825144119229, 2.5828817683596705, 0.7081888710343947, 0.1496415451034128, 0.01136677561558308, 0.0055690189660779765, -0.0005441941308146982],
      Float64[5.038521129225728, 7.036277974255681, 2.5888898889443, 0.7096428906298792, 0.14975799107681567, 0.011399201669279182, 0.005580246938073178, -0.0005477396012294784]
    ]
  ],
  #1, molecule number
  [7,
    #2, molecule formula
    "O2",
    #3, molecule name
    "Oxygen",
    #4, global isotopologue numbers
    Int64[36, 37, 38],
    #5, isotopologue formulae
    String["16O2", "16O18O", "16O17O"],
    #6, AFGL code
    Int64[66, 68, 67],
    #7, abundance fractions
    Float64[0.995262, 0.003991, 0.000742235],
    #8, molecular masses (kg/mole)
    Float64[0.031989830000000004, 0.033994076, 0.032994045],
    #9, Qref
    Float64[215.73, 455.23, 2658.12],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4, 4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.0041, 0.0046, 0.0049],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.8169911476436258, 1.7787349638684251, 0.06038173512650354, 0.010527289874679694],
      Float64[1.8579285870104456, 1.8479587879352657, 0.09662133281030434, 0.021907776219117842],
      Float64[1.8534855385966251, 1.841048066429084, 0.09392781106403698, 0.021584982116504392]
    ]
  ],
  #1, molecule number
  [8,
    #2, molecule formula
    "NO",
    #3, molecule name
    "Nitric Oxide",
    #4, global isotopologue numbers
    Int64[39, 40, 41],
    #5, isotopologue formulae
    String["14N16O", "15N16O", "14N18O"],
    #6, AFGL code
    Int64[46, 56, 48],
    #7, abundance fractions
    Float64[0.993974, 0.003654, 0.001993],
    #8, molecular masses (kg/mole)
    Float64[0.029997989, 0.030995023, 0.032002234000000004],
    #9, Qref
    Float64[1142.13, 789.26, 1204.44],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[10, 10, 10],
    #12, maximum relative errors of interpolation
    Float64[0.0083, 0.0084, 0.0085],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[2.022879002635259, 2.1237886286875303, 0.12399689023022563, -0.009259037224789404, 0.01487593099804976, -0.007768913072230708, 0.002867728523709386, -0.0007316837794678621, -7.443930904828101e-05, 0.00013029548728265498],
      Float64[2.0273058689602084, 2.1309113128029358, 0.12723729874870546, -0.008638605216508043, 0.014726863586520296, -0.0077965256468788685, 0.0028961073794491693, -0.000737512056292195, -6.570282357958806e-05, 0.00013225984839040607],
      Float64[2.0294899939509414, 2.1344200568742244, 0.12882581216417885, -0.008331484385140185, 0.014653941163746952, -0.007814956390896736, 0.002906477403806661, -0.0007396933931758479, -6.357731128186433e-05, 0.00012902023942729102]
    ]
  ],
  #1, molecule number
  [9,
    #2, molecule formula
    "SO2",
    #3, molecule name
    "Sulfur Dioxide",
    #4, global isotopologue numbers
    Int64[42, 43],
    #5, isotopologue formulae
    String["32S16O2", "34S16O2"],
    #6, AFGL code
    Int64[626, 646],
    #7, abundance fractions
    Float64[0.945678, 0.04195],
    #8, molecular masses (kg/mole)
    Float64[0.063961901, 0.065957695],
    #9, Qref
    Float64[6340.3, 6368.98],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.008, 0.008],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[5.322899596644206, 7.361793665479842, 2.4685172475289034, 0.4593724914539494, 0.04693984072073576, -0.0022777800327913322, 0.0016853901083817568, -0.0009605272724346747, 0.0003893627846269787],
      Float64[5.322098471801372, 7.360509652089675, 2.4679029817288893, 0.4592255365312494, 0.04692738529504803, -0.0022790239425536374, 0.0016852856097759883, -0.0009602665910768415, 0.00038778403377648374]
    ]
  ],
  #1, molecule number
  [10,
    #2, molecule formula
    "NO2",
    #3, molecule name
    "Nitrogen Dioxide",
    #4, global isotopologue numbers
    Int64[44, 130],
    #5, isotopologue formulae
    String["14N16O2", "15N16O2"],
    #6, AFGL code
    Int64[646, 656],
    #7, abundance fractions
    Float64[0.991616, 0.003646],
    #8, molecular masses (kg/mole)
    Float64[0.045992904, 0.046989938],
    #9, Qref
    Float64[13577.48, 9324.7],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0098, 0.0098],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.28158486704412, 5.705703174793279, 1.6937108567035881, 0.28623161362735894, 0.03210370018398624, -0.004833093686365356, 0.002727672554545535, -0.0009386770019688129],
      Float64[4.28158486704412, 5.705703174793279, 1.6937108567035881, 0.28623161362735894, 0.03210370018398624, -0.004833093686365356, 0.002727672554545535, -0.0009386770019688129]
    ]
  ],
  #1, molecule number
  [11,
    #2, molecule formula
    "NH3",
    #3, molecule name
    "Ammonia",
    #4, global isotopologue numbers
    Int64[45, 46],
    #5, isotopologue formulae
    String["14NH3", "15NH3"],
    #6, AFGL code
    Int64[4111, 5111],
    #7, abundance fractions
    Float64[0.995872, 0.003661],
    #8, molecular masses (kg/mole)
    Float64[0.017026549, 0.018023583],
    #9, Qref
    Float64[1725.22, 1153.3],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0075, 0.0075],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[3.8555960403301373, 5.040288257810079, 1.4109697012367264, 0.24722432325106913, 0.04050656557366161, -0.0038690829501497603, 0.002775394410359233, -0.0008768219769410557, 0.00020945927498106087],
      Float64[3.7706796997985936, 4.891536114323181, 1.3123307039988636, 0.20002602773432243, 0.02604944595386982, -0.005835809118104995, 0.002916794286104585, -0.0008693703104207806, 0.00018585921410796402]
    ]
  ],
  #1, molecule number
  [12,
    #2, molecule formula
    "HNO3",
    #3, molecule name
    "Nitric Acid",
    #4, global isotopologue numbers
    Int64[47, 117],
    #5, isotopologue formulae
    String["H14N16O3", "H15N16O3"],
    #6, AFGL code
    Int64[146, 156],
    #7, abundance fractions
    Float64[0.98911, 0.003636],
    #8, molecular masses (kg/mole)
    Float64[0.062995644, 0.06399268],
    #9, Qref
    Float64[214000.0, 143000.0],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0072, 0.0073],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[23.230966874867196, 37.713046763492684, 20.728891814718033, 8.121744804012192, 2.305601327065478, 0.48245116851004877, 0.07614122727266803, 0.0068653086813306174, 0.0014066597448483265],
      Float64[23.63526946297278, 38.40170483557111, 21.15217371744953, 8.308610842642612, 2.3655429409129027, 0.49660656393444746, 0.07850397352590655, 0.0071504774322548315, 0.0014405894745461723]
    ]
  ],
  #1, molecule number
  [13,
    #2, molecule formula
    "OH",
    #3, molecule name
    "Hydroxyl",
    #4, global isotopologue numbers
    Int64[48, 49, 50],
    #5, isotopologue formulae
    String["16OH", "18OH", "16OD"],
    #6, AFGL code
    Int64[61, 81, 62],
    #7, abundance fractions
    Float64[0.997473, 0.002, 0.000155371],
    #8, molecular masses (kg/mole)
    Float64[0.01700274, 0.019006986, 0.018008914999999997],
    #9, Qref
    Float64[80.35, 80.88, 209.32],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0068, 0.0069, 0.0064],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.8348152488068417, 1.7541974761207322, 0.05921231834937474, -0.025587632394101223, 0.016697790226026132, -0.008941690431595317, 0.0052568726754337915, -0.003495533506610027, 0.0014773841302516688],
      Float64[1.8505560852233072, 1.7811349343619796, 0.07407620302602635, -0.020226248256200174, 0.017225839806468524, -0.009253742271296772, 0.005191375517307202, -0.003471462133196712, 0.0014838792162081837],
      Float64[1.8980931808868446, 1.8878476706620586, 0.06346999293280998, -0.025382216475793522, 0.01515875114784657, -0.007217914788709562, 0.003960092995426656, -0.0015407268617251596]
    ]
  ],
  #1, molecule number
  [14,
    #2, molecule formula
    "HF",
    #3, molecule name
    "Hydrogen Fluoride",
    #4, global isotopologue numbers
    Int64[51, 110],
    #5, isotopologue formulae
    String["H19F", "D19F"],
    #6, AFGL code
    Int64[19, 29],
    #7, abundance fractions
    Float64[0.999844, 0.000155741],
    #8, molecular masses (kg/mole)
    Float64[0.020006229, 0.021012404],
    #9, Qref
    Float64[41.47, 115.91],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8, 3],
    #12, maximum relative errors of interpolation
    Float64[0.0099, 0.0072],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.7159194314292705, 1.603109095729781, 0.007421783834501869, 7.083322701672393e-05, 0.00164524541490359, -0.0009050350115482955, 0.0008683741416758904, -0.0004347685534003633],
      Float64[1.7382395497979517, 1.6521889828315073, 0.016506244352682442]
    ]
  ],
  #1, molecule number
  [15,
    #2, molecule formula
    "HCl",
    #3, molecule name
    "Hydrogen Chloride",
    #4, global isotopologue numbers
    Int64[52, 53, 107, 108],
    #5, isotopologue formulae
    String["H35Cl", "H37Cl", "D35Cl", "D37Cl"],
    #6, AFGL code
    Int64[15, 17, 25, 27],
    #7, abundance fractions
    Float64[0.757587, 0.242257, 0.000118005, 3.7735000000000004e-05],
    #8, molecular masses (kg/mole)
    Float64[0.035976678, 0.037973729, 0.036982852999999996, 0.038979903999999996],
    #9, Qref
    Float64[160.65, 160.89, 462.78, 464.13],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[3, 3, 4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.0076, 0.0077, 0.0075, 0.0075],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.7389932647539617, 1.6540553376906317, 0.016823417474485014],
      Float64[1.7390283919416012, 1.6541608552427125, 0.016862589474273326],
      Float64[1.7769014767992009, 1.7170608001937298, 0.046270829696325975, 0.013326424837602602],
      Float64[1.7771378593252287, 1.7174573475033619, 0.04644272154650967, 0.013365234531843573]
    ]
  ],
  #1, molecule number
  [16,
    #2, molecule formula
    "HBr",
    #3, molecule name
    "Hydrogen Bromide",
    #4, global isotopologue numbers
    Int64[54, 55, 111, 112],
    #5, isotopologue formulae
    String["H79Br", "H81Br", "D79Br", "D81Br"],
    #6, AFGL code
    Int64[19, 11, 29, 21],
    #7, abundance fractions
    Float64[0.506781, 0.493063, 7.893840000000001e-05, 7.68016e-05],
    #8, molecular masses (kg/mole)
    Float64[0.07992616, 0.081924115, 0.08093233600000001, 0.082930289],
    #9, Qref
    Float64[200.17, 200.23, 586.4, 586.76],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4, 4, 4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.0075, 0.0075, 0.0072, 0.0072],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.7488750133890794, 1.669328183117771, 0.0263672556822101, 0.008005907904859702],
      Float64[1.748900710055193, 1.669366126576109, 0.026379118142379514, 0.008009591348277537],
      Float64[1.8028059428968322, 1.7589216134192316, 0.0639923517825675, 0.016990238558940145],
      Float64[1.8028645166164645, 1.7590220344851282, 0.06403776029402719, 0.01699860427688697]
    ]
  ],
  #1, molecule number
  [17,
    #2, molecule formula
    "HI",
    #3, molecule name
    "Hydrogen Iodide",
    #4, global isotopologue numbers
    Int64[56, 113],
    #5, isotopologue formulae
    String["H127I", "D127I"],
    #6, AFGL code
    Int64[17, 27],
    #7, abundance fractions
    Float64[0.999844, 0.000155741],
    #8, molecular masses (kg/mole)
    Float64[0.127912297, 0.128918472],
    #9, Qref
    Float64[388.99, 1147.06],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.0077, 0.0054],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.766635713209643, 1.699786979483175, 0.03930603336481905, 0.011626090261548741],
      Float64[1.8405398208009391, 1.8185964784146365, 0.08822013072731567, 0.020907340042985556]
    ]
  ],
  #1, molecule number
  [18,
    #2, molecule formula
    "ClO",
    #3, molecule name
    "Chlorine Monoxide",
    #4, global isotopologue numbers
    Int64[57, 58],
    #5, isotopologue formulae
    String["35Cl16O", "37Cl16O"],
    #6, AFGL code
    Int64[56, 76],
    #7, abundance fractions
    Float64[0.755908, 0.24172],
    #8, molecular masses (kg/mole)
    Float64[0.050963768, 0.052960819],
    #9, Qref
    Float64[3274.61, 3332.29],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[7, 7],
    #12, maximum relative errors of interpolation
    Float64[0.0096, 0.0095],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[2.726090979162171, 3.1747524473309494, 0.53340481837599, 0.00862471634302627, -0.0033188447116862716, 0.002818989455709797, -0.0014107531874335184],
      Float64[2.730927657002032, 3.182125195917262, 0.536151715425191, 0.008740020580801774, -0.0034008147472519936, 0.0028494736999678714, -0.001422681090598843]
    ]
  ],
  #1, molecule number
  [19,
    #2, molecule formula
    "OCS",
    #3, molecule name
    "Carbonyl Sulfide",
    #4, global isotopologue numbers
    Int64[59, 60, 61, 62, 63, 135],
    #5, isotopologue formulae
    String["16O12C32S", "16O12C34S", "16O13C32S", "16O12C33S", "18O12C32S", "16O13C34S"],
    #6, AFGL code
    Int64[622, 624, 632, 623, 822, 634],
    #7, abundance fractions
    Float64[0.937395, 0.041583, 0.010531, 0.007399, 0.00188, 0.00046750800000000005],
    #8, molecular masses (kg/mole)
    Float64[0.059966986, 0.06196278, 0.060970341, 0.060966371000000005, 0.061971231, 0.06296613599999999],
    #9, Qref
    Float64[1221.01, 1253.48, 2484.15, 4950.11, 1313.78, 2546.53],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[6, 6, 6, 6, 6, 6],
    #12, maximum relative errors of interpolation
    Float64[0.0057, 0.0058, 0.0054, 0.0057, 0.0051, 0.0057],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[5.126342127660176, 7.034492494351397, 2.46587731887598, 0.5193274408013411, 0.03999760328325479, 0.007938846511621734],
      Float64[5.156948279839978, 7.083476707015535, 2.4892336653688045, 0.5249638067698701, 0.040620760059749725, 0.008015867648682118],
      Float64[5.272543748922432, 7.269032398550013, 2.575970190347018, 0.5458344409651564, 0.043964679171464384, 0.008125632759047364],
      Float64[5.142002462968687, 7.059554966519485, 2.4778239245689733, 0.5222067488374783, 0.04031309017091687, 0.007976767222211124],
      Float64[5.2433970110477635, 7.222282031530504, 2.554774663298425, 0.5409890003043852, 0.043020671404242705, 0.008132056962323376],
      Float64[5.126342127660176, 7.034492494351397, 2.46587731887598, 0.5193274408013411, 0.03999760328325479, 0.007938846511621734]
    ]
  ],
  #1, molecule number
  [20,
    #2, molecule formula
    "H2CO",
    #3, molecule name
    "Formaldehyde",
    #4, global isotopologue numbers
    Int64[64, 65, 66],
    #5, isotopologue formulae
    String["H212C16O", "H213C16O", "H212C18O"],
    #6, AFGL code
    Int64[126, 136, 128],
    #7, abundance fractions
    Float64[0.986237, 0.01108, 0.001978],
    #8, molecular masses (kg/mole)
    Float64[0.030010565, 0.03101392, 0.032014811000000004],
    #9, Qref
    Float64[2844.53, 5837.69, 2986.44],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[7, 8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.009, 0.0066, 0.0066],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.099296685534555, 5.46254008647207, 1.6762575149353685, 0.3607112319829267, 0.0658749047199505, -0.0032135863867808943, 0.0031701326455504386],
      Float64[4.096562472963226, 5.45838322129773, 1.6745262086094586, 0.36054806548097645, 0.06622828358352943, -0.0027924541403665515, 0.0030753561215511077, -0.0005451817250504222],
      Float64[4.0966583506132235, 5.458540573571016, 1.6745828181967304, 0.3605577126652277, 0.06623155701519713, -0.0027929745447514065, 0.0030755155592398103, -0.0005453868311245154]
    ]
  ],
  #1, molecule number
  [21,
    #2, molecule formula
    "HOCl",
    #3, molecule name
    "Hypochlorous Acid",
    #4, global isotopologue numbers
    Int64[67, 68],
    #5, isotopologue formulae
    String["H16O35Cl", "H16O37Cl"],
    #6, AFGL code
    Int64[165, 167],
    #7, abundance fractions
    Float64[0.75579, 0.241683],
    #8, molecular masses (kg/mole)
    Float64[0.051971592999999996, 0.053968643999999996],
    #9, Qref
    Float64[19274.79, 19616.2],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0041, 0.0041],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.143619456416505, 5.464534302985132, 1.535937117235669, 0.21748079463269954, 0.01850070265914905, -0.003603611521264982, 0.0025436888861178897, -0.0012872390879090645, 0.0004250820383528975],
      Float64[4.143652887452801, 5.46458889118496, 1.5359558687538701, 0.21748297167936315, 0.018501348337276458, -0.0036037772770924903, 0.0025437650379011023, -0.0012874828201296928, 0.0004249718859172802]
    ]
  ],
  #1, molecule number
  [22,
    #2, molecule formula
    "N2",
    #3, molecule name
    "Nitrogen",
    #4, global isotopologue numbers
    Int64[69, 118],
    #5, isotopologue formulae
    String["14N2", "14N15N"],
    #6, AFGL code
    Int64[44, 45],
    #7, abundance fractions
    Float64[0.992687, 0.007478],
    #8, molecular masses (kg/mole)
    Float64[0.028006147999999998, 0.029003182],
    #9, Qref
    Float64[467.1, 644.1],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.006, 0.0061],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.7615954213555896, 1.6956422885022302, 0.03176548637294833, 0.01029121610063847],
      Float64[1.763674187848655, 1.6990430426037548, 0.03341123366974245, 0.010717926189910779]
    ]
  ],
  #1, molecule number
  [23,
    #2, molecule formula
    "HCN",
    #3, molecule name
    "Hydrogen Cyanide",
    #4, global isotopologue numbers
    Int64[70, 71, 72],
    #5, isotopologue formulae
    String["H12C14N", "H13C14N", "H12C15N"],
    #6, AFGL code
    Int64[124, 134, 125],
    #7, abundance fractions
    Float64[0.985114, 0.011068, 0.003622],
    #8, molecular masses (kg/mole)
    Float64[0.027010898999999998, 0.028014254000000002, 0.028007933000000002],
    #9, Qref
    Float64[892.2, 1830.97, 615.28],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[7, 5, 5],
    #12, maximum relative errors of interpolation
    Float64[0.0081, 0.01, 0.0079],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[3.1500573196671033, 3.871952656895569, 0.9740837519171931, 0.1661471760780874, -0.0004121953054499657, 0.0038889475688232977, -0.00040715469175870805],
      Float64[3.168146071291571, 3.8993492460929553, 0.9876728830526877, 0.17458143818914662, -0.0006551859591950038],
      Float64[3.3184960101262417, 4.140109825563163, 1.1038407893004953, 0.19960512535349317, -0.0016887611352265353]
    ]
  ],
  #1, molecule number
  [24,
    #2, molecule formula
    "CH3Cl",
    #3, molecule name
    "Methyl Chloride",
    #4, global isotopologue numbers
    Int64[73, 74],
    #5, isotopologue formulae
    String["12CH335Cl", "12CH337Cl"],
    #6, AFGL code
    Int64[215, 217],
    #7, abundance fractions
    Float64[0.748937, 0.239491],
    #8, molecular masses (kg/mole)
    Float64[0.049992328, 0.051989379],
    #9, Qref
    Float64[57916.12, 58833.9],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0083, 0.0083],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[7.29969230201383, 10.76299105271942, 4.6188215208386625, 1.4061580874699755, 0.3052046132618746, 0.04120328062335891, 0.009648711581299096, -0.0008605824437495357],
      Float64[7.299832398848403, 10.763224928473894, 4.61894630404787, 1.4062027724380906, 0.3052168648070337, 0.041205502783890936, 0.00964906200847285, -0.0008605280785930956]
    ]
  ],
  #1, molecule number
  [25,
    #2, molecule formula
    "H2O2",
    #3, molecule name
    "Hydrogen Peroxide",
    #4, global isotopologue numbers
    Int64[75],
    #5, isotopologue formulae
    String["H216O2"],
    #6, AFGL code
    Int64[1661],
    #7, abundance fractions
    Float64[0.994952],
    #8, molecular masses (kg/mole)
    Float64[0.03400548],
    #9, Qref
    Float64[9847.99],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[11],
    #12, maximum relative errors of interpolation
    Float64[0.0067],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[6.659148751736031, 9.483437660502407, 3.4240505334530247, 0.6560605220688904, 0.07013734708443575, 0.0001766124815127057, 0.000513028749733202, -0.00029800025544997056, 7.57280500629065e-05, -0.00016363556385172727, 0.00010136968662237678]
    ]
  ],
  #1, molecule number
  [26,
    #2, molecule formula
    "C2H2",
    #3, molecule name
    "Acetylene",
    #4, global isotopologue numbers
    Int64[76, 77, 105],
    #5, isotopologue formulae
    String["12C2H2", "H12C13CH", "H12C12CD"],
    #6, AFGL code
    Int64[1221, 1231, 1222],
    #7, abundance fractions
    Float64[0.977599, 0.021966, 0.00030455],
    #8, molecular masses (kg/mole)
    Float64[0.02601565, 0.027019005, 0.027021825],
    #9, Qref
    Float64[412.45, 1656.18, 1581.84],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8, 8, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0064, 0.0067, 0.005],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[6.606246082570006, 9.5021885391759, 3.8798354279414227, 1.0421227361830696, 0.15435826120434, 0.024812414831032723, 0.0013956143916372201, -0.0005373524520955186],
      Float64[6.691310653700157, 9.642511905249071, 3.9554948181251057, 1.0672351744197028, 0.1592037312638897, 0.02524828278728946, 0.0014209977344309874, -0.0005472408750750089],
      Float64[8.312818087515264, 12.296135611243283, 5.334539240450019, 1.5041555942106484, 0.25556636740202876, 0.040809216146493243, -5.996404029051661e-05, -0.0008823053560504945, 0.001325829564092551]
    ]
  ],
  #1, molecule number
  [27,
    #2, molecule formula
    "C2H6",
    #3, molecule name
    "Ethane",
    #4, global isotopologue numbers
    Int64[78, 106],
    #5, isotopologue formulae
    String["12C2H6", "12CH313CH3"],
    #6, AFGL code
    Int64[1221, 1231],
    #7, abundance fractions
    Float64[0.97699, 0.021953],
    #8, molecular masses (kg/mole)
    Float64[0.03004695, 0.031050305],
    #9, Qref
    Float64[70882.52, 36191.8],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0084, 0.0084],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[29.306413760745293, 48.87646045145203, 29.277666381786332, 13.421126486275961, 4.868450706061395, 1.4265486625599113, 0.350571677937209, 0.07303068189280282, 0.012110399026012075],
      Float64[29.31971200400752, 48.899765252850294, 29.29331784478448, 13.429280048978747, 4.87184061373862, 1.4276942043686383, 0.3508915525695748, 0.07310824356884638, 0.01212522620861023]
    ]
  ],
  #1, molecule number
  [28,
    #2, molecule formula
    "PH3",
    #3, molecule name
    "Phosphine",
    #4, global isotopologue numbers
    Int64[79],
    #5, isotopologue formulae
    String["31PH3"],
    #6, AFGL code
    Int64[1111],
    #7, abundance fractions
    Float64[0.999533],
    #8, molecular masses (kg/mole)
    Float64[0.033997238000000006],
    #9, Qref
    Float64[3249.44],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8],
    #12, maximum relative errors of interpolation
    Float64[0.0061],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.77754526457992, 6.561487284314362, 2.241225465013582, 0.5263044273306788, 0.09227919777852593, 0.002846217167035497, 0.004031049657860259, -0.0011050831525177987]
    ]
  ],
  #1, molecule number
  [29,
    #2, molecule formula
    "COF2",
    #3, molecule name
    "Carbonyl Fluoride",
    #4, global isotopologue numbers
    Int64[80, 119],
    #5, isotopologue formulae
    String["12C16O19F2", "13C16O19F2"],
    #6, AFGL code
    Int64[269, 369],
    #7, abundance fractions
    Float64[0.986544, 0.011083],
    #8, molecular masses (kg/mole)
    Float64[0.065991722, 0.066995083],
    #9, Qref
    Float64[70028.43, 140000.0],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0037, 0.0037],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[11.629931174695397, 17.881640555690623, 8.45377717908021, 2.6794696985587407, 0.5645965476906625, 0.07777414142631045, 0.00947097051224155, -0.001480514566241098, 0.0007799374263299796],
      Float64[11.635249499599023, 17.889802343475637, 8.457586722083246, 2.6806408754837325, 0.5648398483725394, 0.07780705463459148, 0.009474014430588262, -0.0014807057725274575, 0.0007800005503160179]
    ]
  ],
  #1, molecule number
  [30,
    #2, molecule formula
    "SF6",
    #3, molecule name
    "Sulfur Hexafluoride",
    #4, global isotopologue numbers
    Int64[126],
    #5, isotopologue formulae
    String["32S19F6"],
    #6, AFGL code
    Int64[29],
    #7, abundance fractions
    Float64[0.95018],
    #8, molecular masses (kg/mole)
    Float64[0.145962492],
    #9, Qref
    Float64[1620000.0],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[14],
    #12, maximum relative errors of interpolation
    Float64[0.00024],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1870.7117266249759, 3396.573083307492, 2545.0085260475485, 1579.432604203472, 815.1045425222118, 351.0061028885553, 126.3795862864786, 38.05037645536903, 9.559459105261725, 1.9936675698906514, 0.3424948763082424, 0.04771635428971002, 0.005364060141101408, 0.00045015125453154236]
    ]
  ],
  #1, molecule number
  [31,
    #2, molecule formula
    "H2S",
    #3, molecule name
    "Hydrogen Sulfide",
    #4, global isotopologue numbers
    Int64[81, 82, 83],
    #5, isotopologue formulae
    String["H232S", "H234S", "H233S"],
    #6, AFGL code
    Int64[121, 141, 131],
    #7, abundance fractions
    Float64[0.949884, 0.042137, 0.007498],
    #8, molecular masses (kg/mole)
    Float64[0.033987721, 0.035983514999999994, 0.034987105],
    #9, Qref
    Float64[505.79, 504.35, 2014.94],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8, 8, 8],
    #12, maximum relative errors of interpolation
    Float64[0.0056, 0.0083, 0.0083],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[3.1610034569623626, 3.9079200993763954, 0.8188618250430909, 0.07182197256747361, 0.018994732241263388, -0.004588514291410765, 0.00193425807974279, -0.0006258386386391075],
      Float64[3.136231167171715, 3.868597973183339, 0.801384127880046, 0.06811779813907153, 0.018316506534915036, -0.00478598833219644, 0.0020113628824416046, -0.0006519501360460518],
      Float64[3.136195916650622, 3.868549974556976, 0.8013728165558528, 0.06811697766351898, 0.0183160074539575, -0.004786003572910781, 0.002011194206114096, -0.0006519716060844973]
    ]
  ],
  #1, molecule number
  [32,
    #2, molecule formula
    "HCOOH",
    #3, molecule name
    "Formic Acid",
    #4, global isotopologue numbers
    Int64[84],
    #5, isotopologue formulae
    String["H12C16O16OH"],
    #6, AFGL code
    Int64[126],
    #7, abundance fractions
    Float64[0.983898],
    #8, molecular masses (kg/mole)
    Float64[0.04600548],
    #9, Qref
    Float64[39132.76],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9],
    #12, maximum relative errors of interpolation
    Float64[0.0038],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[11.020149393914501, 16.9755902316461, 8.167574402236214, 2.7594738274122346, 0.6634809411963838, 0.11285962933776061, 0.01800084865083562, -5.6527863588229366e-05, 0.000655575724366031]
    ]
  ],
  #1, molecule number
  [33,
    #2, molecule formula
    "HO2",
    #3, molecule name
    "Hydroperoxyl",
    #4, global isotopologue numbers
    Int64[85],
    #5, isotopologue formulae
    String["H16O2"],
    #6, AFGL code
    Int64[166],
    #7, abundance fractions
    Float64[0.995107],
    #8, molecular masses (kg/mole)
    Float64[0.032997655],
    #9, Qref
    Float64[4300.39],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8],
    #12, maximum relative errors of interpolation
    Float64[0.0074],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[3.4589368258406075, 4.387817669014019, 1.053710491228722, 0.12771732192601018, 0.018001614852028273, -0.006237091281306625, 0.002782788224954099, -0.000725510875379801]
    ]
  ],
  #1, molecule number
  [34,
    #2, molecule formula
    "O",
    #3, molecule name
    "Oxygen Atom",
    #4, global isotopologue numbers
    Int64[86],
    #5, isotopologue formulae
    String["16O"],
    #6, AFGL code
    Int64[6],
    #7, abundance fractions
    Float64[0.997628],
    #8, molecular masses (kg/mole)
    Float64[0.015994915000000002],
    #9, Qref
    Float64[6.72],
    #10, has interpolating chebyshev expansion?
    Bool[false],
    #11, lengths of interpolating chebyshev expansion
    Int64[0],
    #12, maximum relative errors of interpolation
    Float64[0],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[]
    ]
  ],
  #1, molecule number
  [35,
    #2, molecule formula
    "ClONO2",
    #3, molecule name
    "Chlorine Nitrate",
    #4, global isotopologue numbers
    Int64[127, 128],
    #5, isotopologue formulae
    String["35Cl16O14N16O2", "37Cl16O14N16O2"],
    #6, AFGL code
    Int64[5646, 7646],
    #7, abundance fractions
    Float64[0.74957, 0.239694],
    #8, molecular masses (kg/mole)
    Float64[0.096956672, 0.098953723],
    #9, Qref
    Float64[4790000.0, 4910000.0],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[10, 10],
    #12, maximum relative errors of interpolation
    Float64[0.0016, 0.0015],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[112.41039753461149, 192.5969283401371, 121.68354064743774, 57.30328926867914, 20.284701556330944, 5.410076429365737, 1.0808835468204576, 0.15853736954761644, 0.016535623193792363, 0.001150132727673281],
      Float64[112.45459127185921, 192.67265280226152, 121.7313921905935, 57.32582857162725, 20.292682333281764, 5.412206013024527, 1.081309309033366, 0.15859843363539816, 0.016538546940485805, 0.0011482640850822969]
    ]
  ],
  #1, molecule number
  [36,
    #2, molecule formula
    "NO+",
    #3, molecule name
    "Nitric Oxide Cation",
    #4, global isotopologue numbers
    Int64[87],
    #5, isotopologue formulae
    String["14N16O+"],
    #6, AFGL code
    Int64[46],
    #7, abundance fractions
    Float64[0.993974],
    #8, molecular masses (kg/mole)
    Float64[0.029997989],
    #9, Qref
    Float64[311.69],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4],
    #12, maximum relative errors of interpolation
    Float64[0.0059],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.7605530293617953, 1.6939645379456973, 0.030855228843472766, 0.010020286719835495]
    ]
  ],
  #1, molecule number
  [37,
    #2, molecule formula
    "HOBr",
    #3, molecule name
    "Hypobromous Acid",
    #4, global isotopologue numbers
    Int64[88, 89],
    #5, isotopologue formulae
    String["H16O79Br", "H16O81Br"],
    #6, AFGL code
    Int64[169, 161],
    #7, abundance fractions
    Float64[0.505579, 0.491894],
    #8, molecular masses (kg/mole)
    Float64[0.095921076, 0.097919027],
    #9, Qref
    Float64[28339.38, 28237.98],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0051, 0.0052],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[4.483055511505265, 5.995051967678824, 1.7647444347848933, 0.2533586767055279, 0.016936791901325243, -0.002838029482718052, 0.0023154774596215733, -0.0013133740502127011, 0.00048033005031822285],
      Float64[4.486274449888522, 6.000089813847239, 1.7669127031543173, 0.25370966952649165, 0.01695092315480551, -0.0028181141360982265, 0.0023235097509621827, -0.0013007580536887886, 0.00048488773683352804]
    ]
  ],
  #1, molecule number
  [38,
    #2, molecule formula
    "C2H4",
    #3, molecule name
    "Ethylene",
    #4, global isotopologue numbers
    Int64[90, 91],
    #5, isotopologue formulae
    String["12C2H4", "12CH213CH2"],
    #6, AFGL code
    Int64[221, 231],
    #7, abundance fractions
    Float64[0.977294, 0.021959],
    #8, molecular masses (kg/mole)
    Float64[0.028031300000000002, 0.029034655],
    #9, Qref
    Float64[11041.54, 45196.89],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[7, 7],
    #12, maximum relative errors of interpolation
    Float64[0.0091, 0.0091],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[10.306665304688702, 15.910097139939962, 7.822086138506613, 2.834536775621054, 0.7613961848840939, 0.1521884153819452, 0.030223301156906263],
      Float64[10.306853180349535, 15.910413449245775, 7.822252700108585, 2.83459834274206, 0.7614149949091343, 0.15219218454683414, 0.030224516959518628]
    ]
  ],
  #1, molecule number
  [39,
    #2, molecule formula
    "CH3OH",
    #3, molecule name
    "Methanol",
    #4, global isotopologue numbers
    Int64[92],
    #5, isotopologue formulae
    String["12CH316OH"],
    #6, AFGL code
    Int64[2161],
    #7, abundance fractions
    Float64[0.98593],
    #8, molecular masses (kg/mole)
    Float64[0.032026215000000004],
    #9, Qref
    Float64[70569.92],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[11],
    #12, maximum relative errors of interpolation
    Float64[0.0063],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[13.749331968563945, 21.61345904772461, 10.955577051534103, 3.977402542310325, 1.0933160614782178, 0.2226208170647183, 0.038152444149720924, 0.00540115477077876, 0.0001595952470339057, 1.367409072514647e-05, 8.94526436184151e-05]
    ]
  ],
  #1, molecule number
  [40,
    #2, molecule formula
    "CH3Br",
    #3, molecule name
    "Methyl Bromide",
    #4, global isotopologue numbers
    Int64[93, 94],
    #5, isotopologue formulae
    String["12CH379Br", "12CH381Br"],
    #6, AFGL code
    Int64[219, 211],
    #7, abundance fractions
    Float64[0.500995, 0.487433],
    #8, molecular masses (kg/mole)
    Float64[0.093941811, 0.095939764],
    #9, Qref
    Float64[83051.98, 83395.21],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9, 9],
    #12, maximum relative errors of interpolation
    Float64[0.0054, 0.0054],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[8.322478463681351, 12.445792571437991, 5.528514219140895, 1.7163686770516318, 0.37806637006999333, 0.05547356269277692, 0.010458543694655376, -0.0007773050759798394, 0.0004671678943388713],
      Float64[8.329461595021453, 12.457238194953582, 5.534590915534836, 1.7183801520233106, 0.3785478856826856, 0.05556730686168132, 0.010462922226282423, -0.0007745337034998911, 0.0004672506165874779]
    ]
  ],
  #1, molecule number
  [41,
    #2, molecule formula
    "CH3CN",
    #3, molecule name
    "Acetonitrile",
    #4, global isotopologue numbers
    Int64[95],
    #5, isotopologue formulae
    String["12CH312C14N"],
    #6, AFGL code
    Int64[2124],
    #7, abundance fractions
    Float64[0.973866],
    #8, molecular masses (kg/mole)
    Float64[0.041026549],
    #9, Qref
    Float64[88672.19],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[11],
    #12, maximum relative errors of interpolation
    Float64[0.0054],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[23.087630326330327, 37.50816558997893, 20.677521610415493, 8.221774398903325, 2.4554758737437163, 0.5635995654848657, 0.10208943244482249, 0.014556918325425272, 0.001771181894064, -0.00017207723102501403, 0.0001938628200505832]
    ]
  ],
  #1, molecule number
  [42,
    #2, molecule formula
    "CF4",
    #3, molecule name
    "PFC-14",
    #4, global isotopologue numbers
    Int64[96],
    #5, isotopologue formulae
    String["12C19F4"],
    #6, AFGL code
    Int64[29],
    #7, abundance fractions
    Float64[0.98889],
    #8, molecular masses (kg/mole)
    Float64[0.087993616],
    #9, Qref
    Float64[121000.0],
    #10, has interpolating chebyshev expansion?
    Bool[false],
    #11, lengths of interpolating chebyshev expansion
    Int64[0],
    #12, maximum relative errors of interpolation
    Float64[0],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[]
    ]
  ],
  #1, molecule number
  [43,
    #2, molecule formula
    "C4H2",
    #3, molecule name
    "Diacetylene",
    #4, global isotopologue numbers
    Int64[116],
    #5, isotopologue formulae
    String["12C4H2"],
    #6, AFGL code
    Int64[2211],
    #7, abundance fractions
    Float64[0.955998],
    #8, molecular masses (kg/mole)
    Float64[0.05001565],
    #9, Qref
    Float64[9818.97],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[10],
    #12, maximum relative errors of interpolation
    Float64[0.0044],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[139.66835957434776, 241.85305800353672, 157.95302972883783, 78.74096208726291, 30.283942023964958, 9.074668301483646, 2.130903025224642, 0.3944094168400372, 0.05767338336487329, 0.006565583269182045]
    ]
  ],
  #1, molecule number
  [44,
    #2, molecule formula
    "HC3N",
    #3, molecule name
    "Cyanoacetylene",
    #4, global isotopologue numbers
    Int64[109],
    #5, isotopologue formulae
    String["H12C314N"],
    #6, AFGL code
    Int64[1224],
    #7, abundance fractions
    Float64[0.963346],
    #8, molecular masses (kg/mole)
    Float64[0.051010899000000005],
    #9, Qref
    Float64[24786.84],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9],
    #12, maximum relative errors of interpolation
    Float64[0.0062],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[54.67987143122882, 91.66794189149792, 54.56630682138059, 23.48211282466563, 7.394320212397229, 1.7249800510739384, 0.2979283599543656, 0.03885721387024432, 0.0032754500365292927]
    ]
  ],
  #1, molecule number
  [45,
    #2, molecule formula
    "H2",
    #3, molecule name
    "Hydrogen",
    #4, global isotopologue numbers
    Int64[103, 115],
    #5, isotopologue formulae
    String["H2", "HD"],
    #6, AFGL code
    Int64[11, 12],
    #7, abundance fractions
    Float64[0.999688, 0.00031143200000000005],
    #8, molecular masses (kg/mole)
    Float64[0.00201565, 0.003021825],
    #9, Qref
    Float64[7.67, 29.87],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[11, 11],
    #12, maximum relative errors of interpolation
    Float64[0.0072, 0.0085],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.664734139290654, 1.541274894552553, -0.009725177039028837, 0.010711847479083847, -0.0017437301893213019, -0.0036796434870607795, 0.006507341594540717, -0.006895927202643417, 0.006015970973214247, -0.004996438617030119, 0.002286683531610123],
      Float64[1.702133003802882, 1.5472084650368458, 0.02093893129594614, -0.005247622987785912, 0.0065559250200964096, -0.004561928075496269, 0.0036220249507184165, -0.002945429432234235, 0.002444808920724935, -0.002149668002996563, 0.0010267991465589433]
    ]
  ],
  #1, molecule number
  [46,
    #2, molecule formula
    "CS",
    #3, molecule name
    "Carbon Monosulfide",
    #4, global isotopologue numbers
    Int64[97, 98, 99, 100],
    #5, isotopologue formulae
    String["12C32S", "12C34S", "13C32S", "12C33S"],
    #6, AFGL code
    Int64[22, 24, 32, 23],
    #7, abundance fractions
    Float64[0.939624, 0.041682, 0.010556, 0.007417],
    #8, molecular masses (kg/mole)
    Float64[0.043971036, 0.045966786999999995, 0.044974368, 0.044970399],
    #9, Qref
    Float64[253.62, 257.77, 537.5, 1022.97],
    #10, has interpolating chebyshev expansion?
    Bool[true, true, true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[4, 4, 4, 4],
    #12, maximum relative errors of interpolation
    Float64[0.0053, 0.0055, 0.0058, 0.0054],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[1.928548290410206, 1.9559277563753463, 0.13775785658135606, 0.025167997114126095],
      Float64[1.9320348490293633, 1.9612931183817734, 0.13973653631415997, 0.025297795999263812],
      Float64[1.941312753398263, 1.975502364252016, 0.1449170140435975, 0.025619198538681925],
      Float64[1.9303466528392785, 1.9586867922618854, 0.13877109254914913, 0.02523509205534813]
    ]
  ],
  #1, molecule number
  [47,
    #2, molecule formula
    "SO3",
    #3, molecule name
    "Sulfur trioxide",
    #4, global isotopologue numbers
    Int64[114],
    #5, isotopologue formulae
    String["32S16O3"],
    #6, AFGL code
    Int64[26],
    #7, abundance fractions
    Float64[0.9434],
    #8, molecular masses (kg/mole)
    Float64[0.07995682],
    #9, Qref
    Float64[7783.3],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[9],
    #12, maximum relative errors of interpolation
    Float64[0.0083],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[15.308429278661556, 24.00276207581289, 11.897538539905902, 3.9481327675633446, 0.8795595441612623, 0.12893108012856036, 0.011262184567776501, -0.0008464360123436876, 0.0011600801108988534]
    ]
  ],
  #1, molecule number
  [48,
    #2, molecule formula
    "C2N2",
    #3, molecule name
    "Cyanogen",
    #4, global isotopologue numbers
    Int64[123],
    #5, isotopologue formulae
    String["12C214N2"],
    #6, AFGL code
    Int64[4224],
    #7, abundance fractions
    Float64[0.970752],
    #8, molecular masses (kg/mole)
    Float64[0.052006148],
    #9, Qref
    Float64[15582.44],
    #10, has interpolating chebyshev expansion?
    Bool[true],
    #11, lengths of interpolating chebyshev expansion
    Int64[8],
    #12, maximum relative errors of interpolation
    Float64[0.0078],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[23.898246264067748, 38.30762066370932, 20.017340299203973, 6.9592760018001645, 1.616441488972797, 0.2564194297139437, 0.025346935562266384, 0.0020212521367994896]
    ]
  ],
  #1, molecule number
  [49,
    #2, molecule formula
    "COCl2",
    #3, molecule name
    "Phosgene",
    #4, global isotopologue numbers
    Int64[124, 125],
    #5, isotopologue formulae
    String["12C16O35Cl2", "12C16O35Cl37Cl"],
    #6, AFGL code
    Int64[2655, 2657],
    #7, abundance fractions
    Float64[0.566392, 0.362235],
    #8, molecular masses (kg/mole)
    Float64[0.09793261997960001, 0.0999296698896],
    #9, Qref
    Float64[1480000.0, 3040000.0],
    #10, has interpolating chebyshev expansion?
    Bool[true, true],
    #11, lengths of interpolating chebyshev expansion
    Int64[10, 10],
    #12, maximum relative errors of interpolation
    Float64[0.0067, 0.0067],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[28.334677111812944, 45.91333430423453, 24.71182133380345, 9.013419514341322, 2.223400830715287, 0.3666455762516466, 0.03834885535248765, 0.0015768709148049867, 0.0005194428496066747, -0.00021849483690535484],
      Float64[28.38289115870957, 45.99285456416629, 24.75674368149981, 9.030958287676107, 2.2280852109188753, 0.367484622231376, 0.03844429301198539, 0.0015817522400456913, 0.0005200295439730477, -0.00021895935986688327]
    ]
  ],
  [50], #no molecule has been assigned this number
  [51], #no molecule has been assigned this number
  [52], #no molecule has been assigned this number
  #1, molecule number
  [53,
    #2, molecule formula
    "CS2",
    #3, molecule name
    "Carbon disulfide",
    #4, global isotopologue numbers
    Int64[131, 132, 133, 134],
    #5, isotopologue formulae
    String["12C32S2", "32S12C34S", "32S12C33S", "13C32S2"],
    #6, AFGL code
    Int64[222, 224, 223, 232],
    #7, abundance fractions
    Float64[0.892811, 0.07926, 0.014094, 0.01031],
    #8, molecular masses (kg/mole)
    Float64[0.07594414000000001, 0.07793994000000001, 0.076943256, 0.076947495],
    #9, Qref
    Float64[1352.6, 2798.0, 1107.0, 2739.7],
    #10, has interpolating chebyshev expansion?
    Bool[false, false, false, false],
    #11, lengths of interpolating chebyshev expansion
    Int64[0, 0, 0, 0],
    #12, maximum relative errors of interpolation
    Float64[0, 0, 0, 0],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[],
      Float64[],
      Float64[],
      Float64[]
    ]
  ],
  [54], #no molecule has been assigned this number
  #1, molecule number
  [55,
    #2, molecule formula
    "NF3",
    #3, molecule name
    "Nitrogen trifluoride",
    #4, global isotopologue numbers
    Int64[136],
    #5, isotopologue formulae
    String["14N19F3"],
    #6, AFGL code
    Int64[4999],
    #7, abundance fractions
    Float64[0.996337],
    #8, molecular masses (kg/mole)
    Float64[0.070998284],
    #9, Qref
    Float64[346000.0],
    #10, has interpolating chebyshev expansion?
    Bool[false],
    #11, lengths of interpolating chebyshev expansion
    Int64[0],
    #12, maximum relative errors of interpolation
    Float64[0],
    #13, chebyshev expansion coefficients
    Vector{Float64}[
      Float64[]
    ]
  ]
]
