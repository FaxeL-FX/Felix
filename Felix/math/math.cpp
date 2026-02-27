#include "math.h"
namespace math {
	//	number
	const double
		pi = 3.1415926535897932384626433832795028841971693993751058209749445923,
		inf = std::numeric_limits<double>::infinity();
	const complex
		i(0, 1);
	const infsim
		infinity(acch + 1, 1.0),
		zero(acch - 1, 1.0),
		lnInf = -ln_sum(infinity);

	bool sign(double x) {
		unsigned long long i = *(unsigned long long*) & x;
		return i >> 63;
	}
	int E(double x) {
		unsigned long long i = *(unsigned long long*) & x;
		return (int)((i >> 52) % 2048) - 1023;
	}
	unsigned long long M(double x) {
		unsigned long long i = *(unsigned long long*) & x;
		return i % ((unsigned long long)1 << 52);
	}
	double fract(double x) {
		return (double)M(x) / ((unsigned long long)1 << 52);
	}

	complex getFctIntegralConstant() {
		const int n = 256;
		complex res, a = complex(0, pi / n), b = 8 / 256;
		for (int k = 0; k < n; k++) {
			complex t = exp(k * b - (n - k) * a) + 1;
			res = res + t * exp(-t) / ln(t) * (exp((k + 1) * b - (n - k - 1) * a) - exp(k * b - (n - k) * a));
		}
		return res;
	}

	//	complex
	complex::complex(infsim x) {
		this->R = x.getNum(acch).R;
		this->i = x.getNum(acch).i;
	}

	complex operator+(complex x, complex y) { return complex(x.R + y.R, x.i + y.i); }
	complex operator-(complex x) { return complex(-x.R, -x.i); }
	complex operator-(complex x, complex y) { return complex(x.R - y.R, x.i - y.i); }
	complex operator*(complex x, complex y) { return complex(x.R * y.R - x.i * y.i, x.R * y.i + x.i * y.R); }
	complex operator/(complex x, complex y) {
		if (x == y) return 1;
		if (x == 0) return 0;
		if (0 == y) {
			if (x.i == 0) return complex(x.R / 0.0);
			if (x.R == 0) return complex(0, x.i / 0.0);
			return complex(x.R / 0.0, x.i / 0.0);
		}
		if (y.R == inf)		return complex(0, 0);
		if (y.R == -inf)	return complex(-0, 0);
		if (y.i == inf)		return complex(0, 0);
		if (y.i == -inf)	return complex(0, -0);
		long double denominator = 1.0 / (y.R * y.R + y.i * y.i);
		return complex((x.R * y.R + x.i * y.i) * denominator, (y.R * x.i - x.R * y.i) * denominator);
	}
	complex operator%(complex x, complex y) { return x - y * floor(x / y); }
	complex operator<<(complex x, int n) {
		if (n < 0) {
			for (int i = 0; i < -n / 64; i++) x = x / ((unsigned long long)1 << 63) * 0.5;
			return x / ((unsigned long long)1 << -n % 64);
		}
		for (int i = 0; i < n / 64; i++) x = x * ((unsigned long long)1 << 63) * 2;
		return x * ((unsigned long long)1 << n % 64);
	}																						//	need to change for "flex_float" instead "long double"
	complex operator>>(complex x, int n) { return x << -n; }
	bool operator==(complex x, complex y) {
		if (x.R == y.R && x.i == y.i) return true;
		return false;
	}
	bool operator!=(complex x, complex y) { return !(x == y); }

	//	functions
	long double rand(int seed, std::vector<long double> v) {
		long double res = cos(seed);
		for (int i = 1; i <= 16; i++)
			for (auto n : v)
				res = res + sin(n) + res - 16 * floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}
	long double rand(int seed, std::vector<complex> v) {
		long double res = cos(seed);
		for (int i = 1; i <= 16; i++)
			for (auto n : v)
				res = res + sin(((complex)n).R) + res - 16 * floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}
	long double rand(int seed, complex z) {
		std::vector<long double> v = { z.R, z.i };
		return rand(seed, v);
	}
	long double rand(int n) {
		std::vector<long double> v = { (long double)n };
		return rand(n, v);
	}

	//	long double
	long double floor(long double x) { return std::floor(x); }
	long double ceil(long double x) { return std::ceil(x); }
	long double sign(long double x) {
		if (x < 0) return -1;
		return 1;
	}
	long double round(long double x) {
		return (long long)(x + 0.5);
	}

	long double exp(long double x) {
		long double res = x / ((unsigned long long)1 << 32);
		if (res < 0) res = -res;
		for (int i = 0; i < 32; i++) res = 2 * res + res * res;
		if (x < 0) return 1 / (res + 1);
		return res + 1;
	}
	long double log(long double x) {
		long double res = 0, part = 1;
		for (int k = 1; k < 32; k++) {
			part = part * (1 - x);
			res = res + (1 - part) / (k * ((unsigned long long)1 << k));
		}
		return res;
	}
	long double ln(long double x) { return log(fract(x)) + 0.6931471805599453 * E(x); }
	long double pow(long double x, int n) {
		long double res = 1;
		if (n > 0) for (;0 < n;n--) {
			res *= x;
		}
		if (n < 0) for (; n < 0; n++) {
			res /= x;
		}
		return res;
	}

	long double sqrt(long double x) { return exp(0.5 * ln(x)); }
	long double inv_sqrt(long double x) { return exp(-0.5 * ln(x)); }

	long double cos(long double x) {
		long double res = x / ((unsigned long long)1 << 16);
		res = 2 - res * res;
		for (int i = 0; i < 16; i++) res = res * res - 2;
		return res * 0.5;
	}
	long double sin(long double x) {
		long double res = x, sqr;
		for (int i = 1; i <= 32; i++) {
			sqr = x / (pi * i);
			res *= 1 - sqr * sqr;
		}
		return res * exp(-sqr * sqr * 32);
	}
	long double arccos(long double x) {
		if (x == 1) return 0;
		if (x < 0) return pi - arccos(-x);
		long double res = sqrt((x + 1) * 2);
		for (int i = 1; i < 8; i++) res = sqrt(2 + res);
		return sqrt(2 - res) * (1 << 8);
	}

	long double factorial(long double x) {
		long double res = 1.0;
		for (long double k = 2.0; k <= x; k++) res *= k;
		for (long double k = -1.0; x < k; k--) res *= k;
		if (x < -1.0) return 1.0 / res;
		return res;
	}
	long double Binom(long double n, long double k) {
		if (k == floor(k)) {
			long double res = 1;
			if (0 <= k) {
				for (int l = 0; l < k; l++) {
					res *= (n - l) / (l + 1);
				}
			}
			else return 0;
			return res;
		}
		return factorial(n) / (factorial(k) * factorial(n - k));
	}
	std::vector<long double> zbfValues = {};
	long double zbf(int x) {
		if (x == -1)	return -1;
		if (x < -1)		return 0;
		if (x < zbfValues.size()) return zbfValues[x];

		long double res = 0;
		for (int k = 1; k <= x + 1; k++) {
			res = res + zbf(x - k) * (2 * (k % 2) - 1) / factorial(k + 1);
		}
		zbfValues.push_back(res);
		return res;
	}
	std::vector<long double> zetaZeros = {
		14.134725141734693790457251983562470270784257115699243175685567460149,
		21.022039638771554992628479593896902777334340524902781754629520403587,
		25.010857580145688763213790992562821818659549672557996672496542006745,
		30.424876125859513210311897530584091320181560023715440180962146036993,
		32.935061587739189690662368964074903488812715603517039009280003440784,
		37.586178158825671257217763480705332821405597350830793218333001113622,
		40.918719012147495187398126914633254395726165962777279536161303667253,
		43.327073280914999519496122165406805782645668371836871446878893685521,
		48.005150881167159727942472749427516041686844001144425117775312519814,
		49.773832477672302181916784678563724057723178299676662100781955750433,
		52.970321477714460644147296608880990063825017888821224779900748140317,
		56.446247697063394804367759476706127552782264471716631845450969843958,
		59.347044002602353079653648674992219031098772806466669698122451754746,
		60.831778524609809844259901824524003802910090451219178257101348824808,
		65.112544048081606660875054253183705029348149295166722405966501086675,
		67.079810529494173714478828896522216770107144951745558874196669551694,
		69.546401711173979252926857526554738443012474209602510157324539999663,
		72.067157674481907582522107969826168390480906621456697086683306151488,
		75.704690699083933168326916762030345922811903530697400301647775301574,
		77.144840068874805372682664856304637015796032449234461041765231453151,
		79.337375020249367922763592877116228190613246743120030878438720497101,
		82.910380854086030183164837494770609497508880593782149146571306283235,
		84.735492980517050105735311206827741417106627934240818702735529689045,
		87.425274613125229406531667850919213252171886401269028186455557938439,
		88.809111207634465423682348079509378395444893409818675042199871618814,
		92.491899270558484296259725241810684878721794027730646175096750489181,
		94.651344040519886966597925815208153937728027015654852019592474274513,
		95.870634228245309758741029219246781695256461224987998420529281651651,
		98.831194218193692233324420138622327820658039063428196102819321727565,
		101.31785100573139122878544794029230890633286638430089479992831871523,
		103.72553804047833941639840810869528083448117306949576451988516579403,
		105.44662305232609449367083241411180899728275392853513848056944711418,
		107.16861118427640751512335196308619121347670788140476527926471042155,
		111.02953554316967452465645030994435041534596839007305684619079476550,
		111.87465917699263708561207871677059496031174987338587381661941961969,
		114.32022091545271276589093727619107980991765772382989228772843104130,
		116.22668032085755438216080431206475512732985123238322028386264231147,
		118.79078286597621732297913970269982434730621059280938278419371651419,
		121.37012500242064591894553297049992272300131063172874230257513263573,
		122.94682929355258820081746033077001649621438987386351721195003491528,
		124.25681855434576718473200796612992444157353877469356114035507691395,
		127.51668387959649512427932376690607626808830988155498248279977930068,
		129.57870419995605098576803390617997360864095326465943103047083999886,
		131.08768853093265672356637246150134905920354750297504538313992440777,
		133.49773720299758645013049204264060766497417494390467501510225885516,
		134.75650975337387133132606415716973617839606861364716441697609317354,
		138.11604205453344320019155519028244785983527462414623568534482856865,
		139.73620895212138895045004652338246084679005256538260308137013541090,
		141.12370740402112376194035381847535509030066087974762003210466509596,
		143.11184580762063273940512386891392996623310243035463254859852295728,
		146.00098248676551854740250759642468242897574123309580363697688496658,
		147.42276534255960204952118501043150616877277525047683060101046081273,
		150.05352042078488035143246723695937062303732155952820044842911127506,
		150.92525761224146676185252467830562760242677047299671770031135495336,
		153.02469381119889619825654425518544650859043490414550667519976756379,
		156.11290929423786756975018931016919474653530850094292080385607815839,
		157.59759181759405988753050315849876573072389951914173353824961760978,
		158.84998817142049872417499477554027141433508304942696625772418341154,
		161.18896413759602751943734412936955436491579032747546657918809379411,
		163.03070968718198724331103900068799489696446141647768311520959169590,
		165.53706918790041883003891935487479732836725174506860447895315460558,
		167.18443997817451344095775624621037873646076924261676736110699343540,
		169.09451541556882148950587118143183479666764858044162508738214912188,
		169.91197647941169896669984359582179228839443712534137301854144160780,
		173.41153651959155295984611864934559525415606606342011793368228539153,
		174.75419152336572581337876245586691793875571762057166344561154743789,
		176.44143429771041888889264105786093352811849710880971534761261578625,
		178.37740777609997728583093541418442618313236146127250370148904080374,
		179.91648402025699613934003661205123745368760755301840654130067065381,
		182.20707848436646191540703722698779869079745777823990876663006454018,
		184.87446784838750880096064661723425841335102291195066777317864468070,
		185.59878367770747146652770426839264661293471764951328308891979623038,
		187.22892258350185199164154058613124301681073460399031915146420316373,
		189.41615865601693708485228909984532449135710302319335435541994217710,
		192.02665636071378654728363142558343010583992029797709691628912343221,
		193.07972660384570404740220579437605460402061581054886013850435583088,
		195.26539667952923532146318781486225092690505245228692406011097663218,
		196.87648184095831694862226391469620773574602869194221548282317318163,
		198.01530967625191242491991870220886715506269543857099672153480159423,
		201.26475194370378873301613342754817322240286363918673408063271979951,
		202.49359451414053427768666063786431582102024489942005390906915428511,
		204.18967180310455433071643838631368513653452922874190735095968021739,
		205.39469720216328602521237939069309092372291477204840700213409541714,
		207.90625888780620986150196790775364426865940376888399985865752750992,
		209.57650971685625985283564428988675217539078318132616246897745334620,
		211.69086259536530756390748673071929425339403098293564373621001482077,
		213.34791935971266619063912202107260882189718327663306905985370458536,
		214.54704478349142322294420107259069104559988805308307640008161991904,
		216.16953850826370026586956335449812857545371427416411097637615056594,
		219.06759634902137898567725659043724124514918292701135137355787499323,
		220.71491883931400336911559263390633965676114507766196570161193204082,
		221.43070555469333873209747511927607795022233107731990937941995151378,
		224.00700025460433521172887552850489535608598994959552976295036068233,
		224.98332466958228750378252368052865677209005448558742698847775254720,
		227.42144427967929131046143616065963996396914832197662836489382008238,
		229.33741330552534810776008330605574008275234138781851753263649248435,
		231.25018870049916477380618677001037260670849584312337140680603034414,
		231.98723525318024860377166853919786220541983399456249648472682389683,
		233.69340417890830064070449473256978817953722775456583636301480873894,
		236.52422966581620580247550795566297868952949521218912370091896098781,
	};
	long double zetaZero(int n) {
		if (0 <= n && n < zetaZeros.size()) return zetaZeros[n];
		return 0;
	}

	//	complex
	complex floor(complex x) { return complex(floor(x.R), 0); }
	complex ceil(complex x) { return complex(ceil(x.R), 0); }
	long double abs(complex x) {
		if (x.i == 0) return sign(x.R) * x.R;
		if (x.R == 0) return sign(x.i) * x.i;
		return sqrt(x.R * x.R + x.i * x.i);
	}
	long double inv_abs(complex x) { return inv_sqrt(x.R * x.R + x.i * x.i); }
	complex normalize(complex x) {
		if (x == 0) return 1;
		if (x.i == 0)
			if (x.R < 0)	return -1;
			else			return 1;
		if (x.R == 0)
			if (x.i < 0)	return -i;
			else			return i;
		return x * inv_sqrt(x.R * x.R + x.i * x.i);
	}
	complex mul_i(complex x) { return complex(-x.i, x.R); }
	long double arg(complex x) {
		long double cosine;
		if (x.i == 0)
			if (x.R < 0)	cosine = -1;
			else			cosine = 1;
		else				cosine = normalize(x).R;
		if (x.i < 0) return -arccos(cosine);
		return arccos(cosine);
	}
	bool exist(complex x) {
		if (x.R == x.R && x.i == x.i) return true;
		return false;
	}

	complex exp(complex x) {
		if (x.R < -709) return 0;
		complex res = x >> 64;
		if (res.R < 0) res = -res;
		for (int i = 0; i < 64; i++) res = (res << 1) + res * res;
		if (x.R < 0) return 1 / (res + 1);
		return res + 1;
	}
	complex ln(complex x) { return complex(ln(abs(x)), arg(x)); }
	complex pow(complex x, complex y) {
		if (abs(y) == inf) return exp(y * ln(x));
		complex res = 1;
		for (; 0 < y.R;) {
			res = res * x;
			y.R = y.R - 1;
		}
		for (; y.R < 0;) {
			res = res / x;
			y.R = y.R + 1;
		}
		return res * exp(y * ln(x));
	}
	complex sqrt(complex x) { return exp(0.5 * ln(x)); }
	complex inv_sqrt(complex x) { return exp(-0.5 * ln(x)); }

	complex Re(complex x) { return x.R; }
	complex Im(complex x) { return x.i; }
	complex grid(complex x) {
		complex res = x - math::floor(x + complex(0.5, 0.5));
		if (res.R < 0) res.R = -res.R;
		if (res.i < 0) res.i = -res.i;
		if (res.R > res.i) return (res.i - res.i * res.i) * 4;
		return (res.R - res.R * res.R) * 4;
	}

	complex conjugate(complex x) { return complex(x.R, -x.i); }

	complex cosh(complex x) { return (exp(x) + exp(-x)) * 0.5; }
	complex sinh(complex x) { return (exp(x) - exp(-x)) * 0.5; }
	complex coth(complex x) {
		complex positive = exp(x), negative = exp(-x);
		return (positive + negative) / (positive - negative);
	}
	complex tanh(complex x) {
		complex positive = exp(x), negative = exp(-x);
		return (positive - negative) / (positive + negative);
	}
	complex arccosh(complex x) { return ln(x + sqrt(x * x - 1)); }
	complex arcsinh(complex x) { return ln(x + sqrt(x * x + 1)); }
	complex arccoth(complex x) { return arctanh(1 / x); }
	complex arctanh(complex x) { return -0.5 * ln((1 - x) / (1 + x)); }

	complex cos(complex x) { return        cosh(mul_i(x)); }
	complex sin(complex x) { return -mul_i(sinh(mul_i(x))); }
	complex cot(complex x) { return  mul_i(coth(mul_i(x))); }
	complex tan(complex x) { return -mul_i(tanh(mul_i(x))); }
	complex arccos(complex x) { return -mul_i(arccosh(x)); }
	complex arcsin(complex x) { return -mul_i(arcsinh(mul_i(x))); }
	complex arccot(complex x) { return  mul_i(arccoth(mul_i(x))); }
	complex arctan(complex x) { return -mul_i(arctanh(mul_i(x))); }

	complex sin1(complex x) {
		x.R = fmod(x.R + 1.0, 2.0) - 1.0;
		/**/ if (x.R > 0.5) return sin1(1 - x);
		else if (x.R < -0.5) return sin1(-1 - x);
		complex res = x;
		for (int i = 1; i < 64; i++) {
			res = res * (1 - x * x / (i * i));
		}
		return res;
	}
	complex cos1(complex x) {
		x.R = fmod(x.R + 1.0, 2.0) - 1.0;
		/**/ if (x.R > 0.5) return -cos1(1 - x);
		else if (x.R < -0.5) return -cos1(-1 - x);
		complex res = 1;
		for (double i = 0.5; i < 64; i++) {
			res = res * (1 - x * x / (i * i));
		}
		return res;
	}
	complex Binom(complex n, complex k) {
		if (k.i == 0 && k.R == floor(k.R)) {
			complex res = 1;
			if (0 <= k.R) {
				for (int l = 0; l < k.R; l++) {
					res = res * (n - l) / (l + 1);
				}
			}
			else return 0;
			return res;
		}
		complex d = n - k;
		if (d.i == 0 && d.R == floor(d.R) && k.R > n.R) return 0;
		return fct(n) / (fct(k) * fct(n - k));
	}

	complex fct(complex x) {
		if (x.R < -0.5) return -1 / (sin1(x) * fct(-1 - x));
		float n = 2048.5;
		complex res = exp(-x + (x + 0.5) * ln(x + n) - 0.5 * ln(n));
		for (int i = 1; i < n; i++) res = res * i / (x + i) / n * (x + n);
		return res;
	}
	complex inv_fct(complex x) {
		if (x.R < -0.5) return -sin1(x) * fct(-1 - x);
		float n = 2048.5;
		complex res = exp(x - (x + 0.5) * ln(x + n) + 0.5 * ln(n));
		for (int i = 1; i < n; i++) res = res * (x + i) / i * n / (x + n);
		return res;
	}
	complex gamma(complex x) {
		complex res = exp(-0.5772156649 * x) / x;
		for (int k = 1; k != 2048; k++) res = res * exp(x / k) / (1 + x / k);
		return res;
	}
	complex fctIntegral(complex x, complex N) {
		int n = 2048;
		double invN = 1 / (double)n;
		complex res, t = ln(complex(-1.0)), dt, w = 0.0;
		dt = (ln(n) - t) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			w = exp(t) + 1.0;
			res = res + exp(-w + x * ln(w) - N * ln(ln(w))) * (exp(t + 0.5 * dt) - exp(t - 0.5 * dt));
		}
		return res;
	}
	complex Harmonic(complex x) {
		if (x.R < -0.5) return Harmonic(-x - 1) - cos1(x) / sin1(x);
		int n = 2048;
		complex res = ln(1 + x / n);
		for (int k = 1; k <= n; k++)
			res = res + x / (k * (x + k));
		return res;
	}
	complex zeta(complex x) {
		if (x.R > 3) return -fct(x) * exp(-x * 1.83787706640934548356) * sin1(x * 0.5) * zeta(-x - 1);
		if (x.R < -170) return 1;
		int n = 64;
		long double pow2 = 1;
		complex res1 = 0;

		for (int i = 0; i < n; i++) {
			complex res2 = 0;
			for (int j = 0; j <= i; j++) {
				res2 = res2 + (1 - 2 * (j % 2)) * factorial(i) / factorial(i - j) * pow(j + 1, x) / factorial(j);
			}
			pow2 *= 0.5;
			res1 = res1 + res2 * pow2;
		}

		return res1 / (1 - 2 * pow(2, x));
	}
	complex zetaByFct(complex x) {
		if (x.i == 0 && x.R == floor(x.R)) return zbf(x.R);
		return zeta(x) / fct(x);
	}

	complex USumN(complex x, complex n) {
		if (x.R < 1) return USumN(x + 1, n) - pow(x + 1, n);
		complex res;
		if (0 <= n.R) {
			res = pow(x, n + 1) / (n + 1) + pow(x, n) * 0.5;
			if (n.i == 0 && floor(n.R) == n.R) {
				for (int k = 1; k <= n.R;) {
					res = res - zetaByFct(k) * fct(n) / fct(n - k) * pow(x, n - k);
					k += 2;
				}
				return res;
			}
			else {
				for (int k = 1; k <= n.R;) {
					res = res - zetaByFct(k) * fct(n) / fct(n - k) * pow(x, n - k);
					k += 2;
				}
				return res;
			}
		}
		else {
			int m = 64;
			if (n == -1)	res = ln(x + m + 0.5);
			else			res = pow(x + m + 0.5, n + 1) / (n + 1);
			for (int k = 1; k <= m; k++) res = res - pow(x + k, n);
			return res;
		}
	}

	// infsim
	infsim operator+(infsim x, infsim y) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) + y.getNum(i));
		return x;
	}
	infsim operator-(infsim x) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, -x.getNum(i));
		return x;
	}
	infsim operator-(infsim x, infsim y) { return x + -y; }
	infsim operator*(infsim x, infsim y) {
		infsim res;
		for (int i = 0; i < accuracy; i++)
			for (int j = acch - i; i + j - acch < accuracy; j++)
				res.setNum(i + j - acch, res.getNum(i + j - acch) + x.getNum(i) * y.getNum(j));
		return res;
	}
	infsim inv(infsim x) {
		int maxPower = accuracy - 1;
		for (; 0 <= maxPower && x.getNum(maxPower).isZero(); maxPower--) {}
		int invPower = acch - maxPower;
		complex denominator = 1.0 / (x.getNum(maxPower));
		x =  mul(x * infsim(acch + invPower, 1.0), denominator);
		infsim numerator(acch + invPower, 1.0), res;
		for (int i = acch + invPower; 0 <= i; i--) {
			infsim mul(i, numerator.getNum(i));
			res = res + mul;
			numerator = numerator - mul * x;
		}
		return mul(res, denominator);
	}
	infsim operator/(infsim x, infsim y) { return x * inv(y); }
	infsim mul(infsim x, complex y) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) * y);
		return x;
	}
	infsim div(infsim x, complex y) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) / y);
		return x;
	}
	infsim operator<<(infsim x, int n) {
		if (n < 0)	for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) >> -n);
		else		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) << n);
		return x;
	}
	infsim operator>>(infsim x, int n) { return x << -n; }


	infsim Re(infsim x) {
		infsim res;
		for (int i = 0; i < accuracy; i++) res.setNum(i, x.getNum(i).R);
		return res;
	}
	infsim Im(infsim x) {
		infsim res;
		for (int i = 0; i < accuracy; i++) res.setNum(i, x.getNum(i).i);
		return res;
	}
	complex grid(infsim x) {
		complex res = x.getNum(acch) - math::floor(x.getNum(acch) + complex(0.5, 0.5));
		if (res.R < 0) res.R = -res.R;
		if (res.i < 0) res.i = -res.i;
		if (res.R > res.i) return res.i;
		return res.R;
	}
	infsim mul_i(infsim x) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) * math::i);
		return x;
	}
	infsim floor(infsim x) {
		x.Pol[acch] = floor(x.Pol[acch]);
		return x;
	}
	infsim ceil(infsim x) {
		x.Pol[acch] = ceil(x.Pol[acch]);
		return x;
	}

	infsim exp(infsim x) {
		infsim res = x >> 64;
		if (res.getNum(acch).R < 0) res = -res;
		for (int i = 0; i < 64; i++) res = (res << 1) + res * res;
		if (x.getNum(acch).R < 0) return 1 / (res + 1);
		return res + 1;
	}
	infsim ln_sum(infsim x) {
		complex c = x.getNum(acch);
		infsim res, s = 1, t = 1 - infsim(c) / x;
		for (int k = 1; k < 32; k++) {
			s = s * t;
			res = res + div(s, k);
		}
		return res + infsim(ln(c));
	}
	infsim ln_integral(infsim x) {
		int n = 128;
		double invN = 1 / (double)n;
		infsim res, t = 1.0, th = i * abs(x.getNum(acch)), dt;
		if (x.getNum(acch).i < 0.0) th = -th;
		if (x.getNum(acch).R < 0.0) {
			dt = (th - 1.0) * invN;
			for (int i = 0; i < n; i++) {
				t = t + dt;
				res = res + inv(t) * dt;
			}
		}
		else th = 1.0;
		dt = (x - th) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			res = res + inv(t) * dt;
		}
		return res;
	}
	infsim ln(infsim x) {
		int maxPower = accuracy - 1;
		for (; 0 <= maxPower && x.getNum(maxPower).isZero(); maxPower--) {}
		int maxPowerVal = maxPower - acch;
		x = infsim(acch - maxPowerVal, 1.0) * x;
		//return mul(lnInf, maxPowerVal) + ln_integral(x);
		return mul(lnInf, maxPowerVal) + ln_sum(x);
	}
	infsim pow(infsim x, infsim y) {
		if (y.getNum(acch).R >= 1) return pow(x, y - 1) * x;
		if (y.getNum(acch).R < 0) return pow(x, y + 1) / x;
		return exp(y * ln(x));
	}
	infsim sqrt(infsim x) { return exp(0.5 * ln(x)); }
	infsim inv_sqrt(infsim x) { return exp(-0.5 * ln(x)); }

	infsim cosh(infsim x) { return (exp(x) + exp(-x)) * 0.5; }
	infsim sinh(infsim x) { return (exp(x) - exp(-x)) * 0.5; }
	infsim coth(infsim x) {
		infsim positive = exp(x), negative = exp(-x);
		return (positive + negative) / (positive - negative);
	}
	infsim tanh(infsim x) {
		infsim positive = exp(x), negative = exp(-x);
		return (positive - negative) / (positive + negative);
	}
	infsim arccosh(infsim x) { return ln(x + sqrt(x * x - 1)); }
	infsim arcsinh(infsim x) { return ln(x + sqrt(x * x + 1)); }
	infsim arccoth(infsim x) { return arctanh(1 / x); }
	infsim arctanh(infsim x) { return -0.5 * ln((1 - x) / (1 + x)); }

	infsim cos(infsim x) { return        cosh(mul_i(x)); }
	infsim sin(infsim x) { return -mul_i(sinh(mul_i(x))); }
	infsim cot(infsim x) { return  mul_i(coth(mul_i(x))); }
	infsim tan(infsim x) { return -mul_i(tanh(mul_i(x))); }
	infsim arccos(infsim x) { return -mul_i(arccosh(x)); }
	infsim arcsin(infsim x) { return -mul_i(arcsinh(mul_i(x))); }
	infsim arccot(infsim x) { return  mul_i(arccoth(mul_i(x))); }
	infsim arctan(infsim x) { return -mul_i(arctanh(mul_i(x))); }

	infsim sin1(infsim x) {
		x.setNum(acch, fmod(x.getNum(acch).R + 1.0, 2.0) - 1.0);
		/**/ if (x.getNum(acch).R > 0.5) return sin1(1 - x);
		else if (x.getNum(acch).R < -0.5) return sin1(-1 - x);
		infsim res = x;
		for (int i = 1; i < 64; i++) {
			res = res * (1 - div(x * x, i * i));
		}
		return res;
	}
	infsim cos1(infsim x) {
		x.setNum(acch, fmod(x.getNum(acch).R + 1.0, 2.0) - 1.0);
		/**/ if (x.getNum(acch).R > 0.5) return -cos1(1 - x);
		else if (x.getNum(acch).R < -0.5) return -cos1(-1 - x);
		infsim res = 1;
		for (double i = 0.5; i < 64; i++) {
			res = res * (1 - div(x * x, i * i));
		}
		return res;
	}
	infsim Binom(infsim n, infsim k) {
		return fct(n) / (fct(k) * fct(n - k));
	}

	infsim fct(infsim x) {
		if (x.getNum(acch).R <= -1) return fct(x + 1) / (x + 1);
		if (x.getNum(acch).R > 0) return fct(x - 1) * x;
		float n = 512;
		infsim res = exp(-x + (x + 0.5) * ln(x + n) - 0.5 * ln(n));
		for (int i = 1; i < n; i++) res = res * i / n * (x + n) / (x + i);
		return res;
	}
	infsim inv_fct(infsim x) {
		if (x.getNum(acch).R < -0.5) return -sin1(x) * fct(-1 - x);
		float n = 2048.5;
		infsim res = exp(x - (x + 0.5) * ln(x + n) + 0.5 * ln(n));
		for (int i = 1; i < n; i++) res = res * (x + i) / i * n / (x + n);
		return res;
	}
	infsim gamma(infsim x) {
		infsim res = exp(mul(x, -0.5772156649)) / x;
		for (int k = 1; k != 256; k++) res = res * exp(div(x, k)) / (1 + div(x, k));
		return res;
	}
	infsim fctIntegral(infsim x, infsim y) {
		int n = 32;
		double invN = 1 / (double)n;
		infsim res, t = ln(infsim(-1.0)), dt, w = 0.0;
		dt = (infsim(256) - t) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			w = exp(t) + 1.0;
			res = res + exp(-w + x * ln(w) - y * ln(ln(w))) * (exp(t) - exp(t - dt));
		}
		return res;
	}
	infsim Harmonic(infsim x) {
		if (x.getNum(acch).R < -0.5) return Harmonic(-x - 1) - cos1(x) / sin1(x);
		int n = 64;
		infsim res = ln(1 + x / n);
		for (int k = 1; k <= n; k++)
			res = res + x / (k * (x + k));
		return res;
	}
	infsim zeta(infsim x) {
		if (x.getNum(acch).R > 0) return -pow(2 * pi, -x) * sin1(x * 0.5) * fct(x) * zeta(-x - 1);
		infsim res1 = 0;
		for (int i = 0; i < 12; i++) {
			infsim res2 = 0;
			for (int j = 0; j <= i; j++) {
				res2 = res2 + (1 - 2 * (j % 2)) * factorial(i) * pow(j + 1, x) / (factorial(j) * factorial(i - j));
			}
			res1 = res1 + mul(res2, pow(complex(2.0), -(i + 1)));
		}
		return res1 / (1 - pow(2, x + 1));
	}
}