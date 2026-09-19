#include "math.h"
namespace math {
	//	number
	const real
		pi = (long double)3.1415926535897932384626433832795028841971693993751058209749445923,
		inf = std::numeric_limits<long double>::infinity();
	const complex
		i(0, 1);

	//	complex
	complex operator+(complex x, complex y) { return complex(x.R + y.R, x.i + y.i); }
	complex operator-(complex x) { return complex(-x.R, -x.i); }
	complex operator-(complex x, complex y) { return complex(x.R - y.R, x.i - y.i); }
	complex operator*(complex x, complex y) { return complex(x.R * y.R - x.i * y.i, x.R * y.i + x.i * y.R); }
	complex operator/(complex x, complex y) {
		if (abs(x) == abs(y)) {
			x = normalize(x);
			y = normalize(y);
		}
		if (x == (real)0) return (real)0;
		if ((real)0 == y) {
			if (x.i == (real)0) return complex(x.R / (real)0.0);
			if (x.R == (real)0) return complex(0, x.i / (real)0.0);
			return complex(x.R / (real)0.0, x.i / (real)0.0);
		}
		if (y.R == inf)		return complex(0, 0);
		if (y.R == -inf)	return complex(-0, 0);
		if (y.i == inf)		return complex(0, 0);
		if (y.i == -inf)	return complex(0, -0);
		real denominator = (real)1.0 / (y.R * y.R + y.i * y.i);
		return complex((x.R * y.R + x.i * y.i) * denominator, (y.R * x.i - x.R * y.i) * denominator);
	}
	complex operator%(complex x, complex y) { return x - y * floor(x / y); }
	bool operator==(complex x, complex y) { return x.R == y.R && x.i == y.i; }
	bool operator!=(complex x, complex y) { return !(x == y); }

	//	FUNCTIONS

	//	real
	real rand(int seed) {
		float res = std::cos((float)seed);
		for (int i = 1; i <= 16; i++)
			res = res + std::sin(res) + res - 16 * std::floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}

	real factorial(real x) {
		real res = 1.0;
		for (long long k = 2.0; (real)k <= x; k++) res = res * (real)k;
		for (long long k = -1.0; x < (real)k; k--) res = res * (real)k;
		if (x < (real)-1.0) return (real)1.0 / res;
		return res;
	}
	real Binom(real n, real k) {
		if (k == (real)std::floor(k)) {
			real res = 1;
			if ((real)0 <= k) {
				for (int l = 0; (real)l < k; l++) {
					res = res * (n - (real)l) / (real)(l + 1);
				}
			}
			else return 0;
			return res;
		}
		return factorial(n) / (factorial(k) * factorial(n - k));
	}
	std::vector<real> zbfValues = { (real)(-0.5) };
	std::mutex zbfMutex;
	real zbf(int x) {
		if (x == -1)	return -1;
		if (x < -1)		return 0;
		{
			std::lock_guard<std::mutex> l{ zbfMutex };
			if (x < zbfValues.size()) return zbfValues[x];
		}

		real res = 0;
		if (x % 2) for (int k = 1; k <= x + 1; k++) {
			res = res + zbf(x - k) * (2 * (k % 2) - 1) / factorial(k + 1);
		}
		{
			std::lock_guard<std::mutex> l{ zbfMutex };
			if (x >= zbfValues.size()) zbfValues.resize(x + 1);
			zbfValues[x] = res;
		}
		return res;
	}
	std::vector<real> zetaZeros = {
		(real)14.134725141734693790457251983562470270784257115699243175685567460149,
		(real)21.022039638771554992628479593896902777334340524902781754629520403587,
		(real)25.010857580145688763213790992562821818659549672557996672496542006745,
		(real)30.424876125859513210311897530584091320181560023715440180962146036993,
		(real)32.935061587739189690662368964074903488812715603517039009280003440784,
		(real)37.586178158825671257217763480705332821405597350830793218333001113622,
		(real)40.918719012147495187398126914633254395726165962777279536161303667253,
		(real)43.327073280914999519496122165406805782645668371836871446878893685521,
		(real)48.005150881167159727942472749427516041686844001144425117775312519814,
		(real)49.773832477672302181916784678563724057723178299676662100781955750433,
		(real)52.970321477714460644147296608880990063825017888821224779900748140317,
		(real)56.446247697063394804367759476706127552782264471716631845450969843958,
		(real)59.347044002602353079653648674992219031098772806466669698122451754746,
		(real)60.831778524609809844259901824524003802910090451219178257101348824808,
		(real)65.112544048081606660875054253183705029348149295166722405966501086675,
		(real)67.079810529494173714478828896522216770107144951745558874196669551694,
		(real)69.546401711173979252926857526554738443012474209602510157324539999663,
		(real)72.067157674481907582522107969826168390480906621456697086683306151488,
		(real)75.704690699083933168326916762030345922811903530697400301647775301574,
		(real)77.144840068874805372682664856304637015796032449234461041765231453151,
		(real)79.337375020249367922763592877116228190613246743120030878438720497101,
		(real)82.910380854086030183164837494770609497508880593782149146571306283235,
		(real)84.735492980517050105735311206827741417106627934240818702735529689045,
		(real)87.425274613125229406531667850919213252171886401269028186455557938439,
		(real)88.809111207634465423682348079509378395444893409818675042199871618814,
		(real)92.491899270558484296259725241810684878721794027730646175096750489181,
		(real)94.651344040519886966597925815208153937728027015654852019592474274513,
		(real)95.870634228245309758741029219246781695256461224987998420529281651651,
		(real)98.831194218193692233324420138622327820658039063428196102819321727565,
		(real)101.31785100573139122878544794029230890633286638430089479992831871523,
		(real)103.72553804047833941639840810869528083448117306949576451988516579403,
		(real)105.44662305232609449367083241411180899728275392853513848056944711418,
		(real)107.16861118427640751512335196308619121347670788140476527926471042155,
		(real)111.02953554316967452465645030994435041534596839007305684619079476550,
		(real)111.87465917699263708561207871677059496031174987338587381661941961969,
		(real)114.32022091545271276589093727619107980991765772382989228772843104130,
		(real)116.22668032085755438216080431206475512732985123238322028386264231147,
		(real)118.79078286597621732297913970269982434730621059280938278419371651419,
		(real)121.37012500242064591894553297049992272300131063172874230257513263573,
		(real)122.94682929355258820081746033077001649621438987386351721195003491528,
		(real)124.25681855434576718473200796612992444157353877469356114035507691395,
		(real)127.51668387959649512427932376690607626808830988155498248279977930068,
		(real)129.57870419995605098576803390617997360864095326465943103047083999886,
		(real)131.08768853093265672356637246150134905920354750297504538313992440777,
		(real)133.49773720299758645013049204264060766497417494390467501510225885516,
		(real)134.75650975337387133132606415716973617839606861364716441697609317354,
		(real)138.11604205453344320019155519028244785983527462414623568534482856865,
		(real)139.73620895212138895045004652338246084679005256538260308137013541090,
		(real)141.12370740402112376194035381847535509030066087974762003210466509596,
		(real)143.11184580762063273940512386891392996623310243035463254859852295728,
		(real)146.00098248676551854740250759642468242897574123309580363697688496658,
		(real)147.42276534255960204952118501043150616877277525047683060101046081273,
		(real)150.05352042078488035143246723695937062303732155952820044842911127506,
		(real)150.92525761224146676185252467830562760242677047299671770031135495336,
		(real)153.02469381119889619825654425518544650859043490414550667519976756379,
		(real)156.11290929423786756975018931016919474653530850094292080385607815839,
		(real)157.59759181759405988753050315849876573072389951914173353824961760978,
		(real)158.84998817142049872417499477554027141433508304942696625772418341154,
		(real)161.18896413759602751943734412936955436491579032747546657918809379411,
		(real)163.03070968718198724331103900068799489696446141647768311520959169590,
		(real)165.53706918790041883003891935487479732836725174506860447895315460558,
		(real)167.18443997817451344095775624621037873646076924261676736110699343540,
		(real)169.09451541556882148950587118143183479666764858044162508738214912188,
		(real)169.91197647941169896669984359582179228839443712534137301854144160780,
		(real)173.41153651959155295984611864934559525415606606342011793368228539153,
		(real)174.75419152336572581337876245586691793875571762057166344561154743789,
		(real)176.44143429771041888889264105786093352811849710880971534761261578625,
		(real)178.37740777609997728583093541418442618313236146127250370148904080374,
		(real)179.91648402025699613934003661205123745368760755301840654130067065381,
		(real)182.20707848436646191540703722698779869079745777823990876663006454018,
		(real)184.87446784838750880096064661723425841335102291195066777317864468070,
		(real)185.59878367770747146652770426839264661293471764951328308891979623038,
		(real)187.22892258350185199164154058613124301681073460399031915146420316373,
		(real)189.41615865601693708485228909984532449135710302319335435541994217710,
		(real)192.02665636071378654728363142558343010583992029797709691628912343221,
		(real)193.07972660384570404740220579437605460402061581054886013850435583088,
		(real)195.26539667952923532146318781486225092690505245228692406011097663218,
		(real)196.87648184095831694862226391469620773574602869194221548282317318163,
		(real)198.01530967625191242491991870220886715506269543857099672153480159423,
		(real)201.26475194370378873301613342754817322240286363918673408063271979951,
		(real)202.49359451414053427768666063786431582102024489942005390906915428511,
		(real)204.18967180310455433071643838631368513653452922874190735095968021739,
		(real)205.39469720216328602521237939069309092372291477204840700213409541714,
		(real)207.90625888780620986150196790775364426865940376888399985865752750992,
		(real)209.57650971685625985283564428988675217539078318132616246897745334620,
		(real)211.69086259536530756390748673071929425339403098293564373621001482077,
		(real)213.34791935971266619063912202107260882189718327663306905985370458536,
		(real)214.54704478349142322294420107259069104559988805308307640008161991904,
		(real)216.16953850826370026586956335449812857545371427416411097637615056594,
		(real)219.06759634902137898567725659043724124514918292701135137355787499323,
		(real)220.71491883931400336911559263390633965676114507766196570161193204082,
		(real)221.43070555469333873209747511927607795022233107731990937941995151378,
		(real)224.00700025460433521172887552850489535608598994959552976295036068233,
		(real)224.98332466958228750378252368052865677209005448558742698847775254720,
		(real)227.42144427967929131046143616065963996396914832197662836489382008238,
		(real)229.33741330552534810776008330605574008275234138781851753263649248435,
		(real)231.25018870049916477380618677001037260670849584312337140680603034414,
		(real)231.98723525318024860377166853919786220541983399456249648472682389683,
		(real)233.69340417890830064070449473256978817953722775456583636301480873894,
		(real)236.52422966581620580247550795566297868952949521218912370091896098781,
	};
	real zetaZero(int n) {
		if (0 <= n && n < zetaZeros.size()) return zetaZeros[n];
		return 0;
	}

	//	complex
	complex floor(complex x) { return complex(std::floor(x.R), 0); }
	complex ceil(complex x) { return complex(std::ceil(x.R), 0); }
	real abs(complex x) {
		if (x.i == (real)0) return std::abs(x.R);
		if (x.R == (real)0) return std::abs(x.i);
		return (real)std::sqrt((long double)(x.R * x.R + x.i * x.i));
	}
	complex normalize(complex x) {
		if (x == (real)0) return (real)1;
		if (x.i == 0)
			if (x.R < 0)	return (real)-1;
			else			return (real)1;
		if (x.R == 0)
			if (x.i < 0)	return -i;
			else			return i;
		if (abs(x) == inf) {
			if (x.R == inf) {
				if (x.i < 0)	return ((real)1 - i) * (real)std::sqrt(0.5);
				else			return ((real)1 + i) * (real)std::sqrt(0.5);
			}
			else {
				if (x.i < 0)	return ((real)-1 - i) * (real)std::sqrt(0.5);
				else			return ((real)-1 + i) * (real)std::sqrt(0.5);
			}
		}
		return x * ((real)1.0 / (real)std::sqrt(x.R * x.R + x.i * x.i));
	}
	complex mul_i(complex x) { return complex(-x.i, x.R); }
	real arg(complex x) {
		real cosine = normalize(x).R;
		if (x.i < 0) return (real)-std::acos(cosine);
		return (real)std::acos(cosine);
	}
	bool exist(complex x) {
		if (x.R == x.R && x.i == x.i) return true;
		return false;
	}

	complex exp(complex x) {
		complex res{ (real)std::cos(x.i), (real)std::sin(x.i) };
		res.R = res.R * std::exp(x.R);
		res.i = res.i * std::exp(x.R);
		return res;
	}
	complex ln(complex x) {
		//if (abs(x) == inf) return complex(710, arg(x));
		//if (abs(x) == 0.0) return complex(-710, arg(x));
		return complex(std::log(abs(x)), arg(x));
	}
	complex pow(complex x, complex y) {
		if (abs(y) == inf) return exp(y * ln(x));
		if (x == (real)0) {
			if (y.R < 0) return inf;
			if (0 < y.R) return (real)0;
		}

		complex res = (real)1;
		if (y.R != y.R - 1) for (; 0 < y.R;) {
			res = res * x;
			y.R = y.R - 1;
		}
		if (y.R != y.R + 1) for (; y.R < 0;) {
			res = res / x;
			y.R = y.R + 1;
		}
		if (y == (real)0) return res;

		return res * exp(y * ln(x));
	}
	complex sqrt(complex x) {
		if (x == (real)0) return (real)0;
		return exp((real)0.5 * ln(x));
	}
	complex inv_sqrt(complex x) { return exp((real)-0.5 * ln(x)); }

	real Re(complex x) { return x.R; }
	real Im(complex x) { return x.i; }
	real grid(complex x) {
		complex res = x - math::floor(x + complex(0.5, 0.5));
		if (res.R < 0) res.R = -res.R;
		if (res.i < 0) res.i = -res.i;
		if (res.R > res.i) return (res.i - res.i * res.i) * 4;
		return (res.R - res.R * res.R) * 4;
	}
	real axis(complex x) {
		if (x.R < 0) x.R = -x.R;
		if (x.i < 0) x.i = -x.i;
		return std::min(x.R, x.i);
	}

	complex sin1(complex x) {
		x.R = fmod(x.R + 1.0, 2.0) - 1.0;
		/**/ if (x.R >  0.5) return sin1((real)1 - x);
		else if (x.R < -0.5) return sin1((real)-1 - x);
		return sin(pi * x) / pi;
	}
	complex cos1(complex x) {
		x.R = fmod((long double)x.R + 1.0, 2.0) - 1.0;
		/**/ if (x.R >  0.5) return -cos1((real)1 - x);
		else if (x.R < -0.5) return -cos1((real)-1 - x);
		return cos(pi * x);
	}
	complex Binom(complex n, complex k) {
		if (k.i == 0 && k.R == std::floor(k.R)) {
			complex res = (real)1;
			if (0 <= k.R) {
				for (int l = 0; l < k.R; l++) {
					res = res * (n - (real)l) / (real)(l + 1);
				}
			}
			else return (real)0;
			return res;
		}
		complex d = n - k;
		if (d.i == 0 && d.R == std::floor(d.R) && k.R > n.R) return (real)0;
		return fct(n) / (fct(k) * fct(n - k));
	}

	complex fct(complex x) {
		if (x.R == inf) return inf;
		if (x.R == -inf) return (real)0;
		if (abs(x.i) == inf) return (real)0;
		if (x.R < -0.5) return (real)-1 / (sin1(x) * fct((real)-1 - x));
		complex 
			res = (real)1;
		if (x != x - (real)1 && !(std::isnan(res.R) || std::isnan(res.i))) {
			for (; 1 <= x.R;) {
				res = res * x;
				x = x - (real)1;
			}
		}
		else {
			return sqrt((real)2 * pi) * exp((x + (real)0.5) * (ln(x + (real)0.5) - (real)1.0));
		}
		if (x == (real)0) return res;
		real n = 128.5;
		res = res * exp(-x + (x + (real)0.5) * ln(x + n) - (real)(0.5 * std::log(n)));
		for (long long i = 1; i < n; i++) res = res * (real)i / (x + (real)i) / n * (x + n);
		return res;
	}
	complex gamma(complex x) {
		return fct(x - (real)1);

		complex res = exp((real)-0.5772156649 * x) / x;
		for (int k = 1; k != 2048; k++) res = res * exp(x / (real)k) / ((real)1 + x / (real)k);
		return res;
	}
	complex fctIntegral(complex x, complex N) {
		int n = 1024;
		real invN = 1 / (double)n;
		complex res, t = ln(complex(-1.0)), dt, w = (real)0.0;
		dt = ((real)std::log(n) - t) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			w = exp(t) + (real)1.0;
			res = res + exp(-w + x * ln(w) - N * ln(ln(w))) * (exp(t + (real)0.5 * dt) - exp(t - (real)0.5 * dt));
		}
		return res;
	}
	complex Harmonic(complex x) {
		if (x.R < -0.5) return Harmonic(-x - (real)1) - cos1(x) / sin1(x);
		int n = 2048;
		complex res = ln((real)1 + x / (real)n);
		for (int k = 1; k <= n; k++)
			res = res + x / ((real)k * (x + (real)k));
		return res;
	}
	complex zeta(complex x) {
		if (0 < x.R) return -fct(x) * pow((real)2 * pi, -x) * sin1(x * (real)0.5) * zeta(-x - (real)1);
		int n = 128;
		complex res = (real)0;
		for (int k = 1; k <= n; k++) res = res + pow((real)k, x);
		return res - pow((real)n + (real)0.5, x + (real)1) / (x + (real)1);
	}
	complex zetaByFct(complex x) {
		if (x.i == 0 && x.R == std::floor(x.R)) return zbf(x.R);
		return zeta(x) / fct(x);
	}

	complex USumN(complex x, complex n) {
		complex res;
		if (x != x + (real)1) for (; x.R < 1 + std::abs(n.i);) {
			x = x + (real)1;
			res = res - pow(x, n);
		}
		if (0 <= n.R) {
			res = res + pow(x, n + (real)1) / (n + (real)1) + pow(x, n) * (real)0.5;
			for (int k = 1; k < n.R + 8; k++) {
				complex multiplier = zbf(k);
				for (int l = 0; l < k; l++) multiplier = multiplier * (n - (real)l);
				res = res - multiplier * pow(x, n - (real)k);
			}
			return res;
		}
		else {
			int m = 128;
			if (n == (real)-1)	res = res + ln(x + (real)m + (real)0.5);
			else				res = res + pow(x + (real)m + (real)0.5, n + (real)1) / (n + (real)1);
			for (int k = 1; k <= m; k++) res = res - pow(x + (real)k, n);
			return res;
		}
	}

	complex Integral(complex(*f)(complex), complex a, complex b) {
		const long long n = 16 * (1 + abs(b - a));
		complex
			res = (real)0,
			dx = (b - a) / (real)n,
			x = a + dx * (real)0.5;
		for (int k = 0; k < n; k++) {
			res = res + f(x);
			x = x + dx;
		}
		return res * dx;
	}
	complex Integral(complex(*f)(complex), complex(*path)(real)) {
		const long long n = 512;
		const real dt = 1.0 / n;
		complex res = (real)0;
		for (real k = dt * (real)0.5; k < 1.0;) {
			res = res + f(path(k)) * (path(k + dt * (real)0.5) - path(k - dt * (real)0.5));
			k = k + dt;
		}
		return res;
	}
	complex Derivative(complex(*f)(complex), complex x) {
		const real dx = 1.0 / 4096;
		return (f(x + dx) - f(x - dx)) / ((real)2 * dx);
	}
}