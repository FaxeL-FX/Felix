//	v1.10.2

#include <iostream>
#include "include.h"
#include "object/object.h"
#include "bmp/BMPWriter.h"
#include <thread>
#include <Windows.h>

const unsigned int prc_count = std::thread::hardware_concurrency();

struct Color {
	float R, G, B, A;

	Color() {
		this->R = 0;
		this->G = 0;
		this->B = 0;
		this->A = 1;
	}
	Color(float R, float G, float B, float A) {
		this->R = R;
		this->G = G;
		this->B = B;
		this->A = A;
	}
	Color(std::vector<float> rgb) {
		this->R = rgb[0];
		this->G = rgb[1];
		this->B = rgb[2];
		this->A = 1;
	}
	Color(float W) {
		this->R = W;
		this->G = W;
		this->B = W;
		this->A = 1;
	}
};
Color pen(Color a, Color b) {
	a.R = b.A * b.R + (1 - b.A) * a.R;
	a.G = b.A * b.G + (1 - b.A) * a.G;
	a.B = b.A * b.B + (1 - b.A) * a.B;
	a.A = 1;
	return a;
}
Color penAdd(Color a, Color b) {
	a.R += b.R * b.A;
	a.G += b.G * b.A;
	a.B += b.B * b.A;
	if (a.R > 1) a.R = 1;
	if (a.G > 1) a.G = 1;
	if (a.B > 1) a.B = 1;
	if (a.R < 0) a.R = 0;
	if (a.G < 0) a.G = 0;
	if (a.B < 0) a.B = 0;
	return a;
}
float absMul(float c, float absolute) { return c * absolute / (absolute + 1); }
Color toCol(math::complex x) {
	Color c;
	float angle = math::arg(x), absolute = math::abs(x);
	c.A = 1;
	c.R = (math::cos(angle) + 1) * 0.5;
	c.G = (math::cos(angle + 2.0943951) + 1) * 0.5;
	c.B = (math::cos(angle - 2.0943951) + 1) * 0.5;
	c.R = absMul(c.R, absolute);
	c.G = absMul(c.G, absolute);
	c.B = absMul(c.B, absolute);
	return c;
}
std::vector<float> toVecF(Color c) { return { c.R, c.G, c.B }; }
Color fncColor(Color col, long double im) {
	if (im == im)	col.A = 1 / (1 + im * im);
	else			col.A = 0;
	return col;
}
Color RGBColor(double col) {
	return Color(
		(2 + cos(col * 1.8)) / 3,
		(2 + cos(col * 1.8 + 2.095)) / 3,
		(2 + cos(col * 1.8 - 2.095)) / 3,
		1);
}

math::complex plot_center = 0;
long double plot_radius = 8;

bool run_command(std::string c) {
	std::vector<std::string> args;
	for (;;) {
		while (c[0] == ' ') c = c.substr(1);
		int index = c.find(' ');
		if (-1 < index) {
			args.push_back(c.substr(0, index));
			c = c.substr(index + 1);
		}
		else {
			args.push_back(c);
			break;
		}
	}

	if (args[0] == "def") {
		if (args.size() < 3) return false;
		Function f;

		int bracket_index = args[1].find('(');
		if (-1 < bracket_index) {
			f.name = args[1].substr(0, bracket_index);
			std::string arguments = args[1].substr(bracket_index + 1);
			while (arguments.length() > 0) {
				int index1 = arguments.find(','),
					index2 = arguments.find(')');
				Variable var;
				/**/ if (-1 < index1) var.name = arguments.substr(0, index1);
				else if (-1 < index2) var.name = arguments.substr(0, index2);
				else /*------------*/ var.name = arguments;

				var = Variable(var.name);

				f.args.push_back(var);
				arguments = arguments.substr(f.args[f.args.size() - 1].name.length() + 1);
			}
		}
		else f.name = args[1];

		f.id = 0;
		for (char c : f.name)
			f.id = (f.id << 1) + (int)c;

		std::string expr = "";
		for (int i = 2 + (args[2] == "="); i < args.size(); i++) expr += args[i] + ' ';
		f.objects = parse_expr(expr);

		for (int i = 0; i < functions_list.size(); i++)
			if (functions_list[i].name == f.name) {
				functions_list[i] = f;
				return true;
			}
		functions_list.push_back(f);
		return true;
	}
	if (args[0] == "del") {
		if (args.size() < 2) return false;
		for (int i = 0; i < functions_list.size(); i++)
			if (functions_list[i].name == args[1])
				functions_list.erase(functions_list.begin() + i);
		return true;
	}
	if (args[0] == "list") {
		for (auto f : functions_list) {
			std::cout << " -> " << f.name;
			if (f.args.size() > 0) {
				std::cout << "(";
				for (auto a : f.args) std::cout << a.name << ",";
				std::cout << "\b)";
			}
			std::cout << " [ID:" << f.id << "]\n";
		}
		return true;
	}
	if (args[0] == "clear") {
		functions_list.clear();
		return true;
	}
	if (args[0] == "scale") {
		if (args.size() < 3) return false;
		plot_center = value(parse_expr(args[1]), {});
		plot_radius = math::abs(value(parse_expr(args[2]), {}));
		return true;
	}
	if (args[0] == "print") {
		bool grid = true, cmplx = false, eq = false;
		if (args.size() < 2) return false;
		for (int i = 3; i < args.size(); i++) {
			/**/ if (args[i] == "noGrid") grid = false;
			else if (args[i] == "c") cmplx = true;
			else if (args[i] == "eq") eq = true;
		}
		int resolution;
		if (args.size() < 3)	resolution = 400;
		else					resolution = std::stoi(args[2]);
		if (resolution < 0) resolution = -resolution;
		if (1600 < resolution) resolution = 1600;
		int prc = resolution + (prc_count < resolution) * (prc_count - resolution);
		std::vector<std::vector<float>> img;
		std::vector<std::vector<float>>* im = &img;
		for (int i = 0; i < resolution * resolution; i++) img.push_back({ 0, 0, 0 });

		std::vector<Function> target_fncs;
		args[1] += ';';
		for (; 0 < args[1].size();) {
			if (args[1][0] == ';' || args[1][0] == ',') break;
			bool no_func = true;
			std::string token = parse_token(args[1]);
			for (auto fnc : functions_list)
				if (fnc.name == token) {
					target_fncs.push_back(fnc);
					no_func = false;
					break;
				}
			if (no_func) {
				target_fncs.push_back(Function(-1, token,
					{
						Variable("x")
					},
					{
						Object(nameToType(token), { 1 }),
						Object("x")
					}
				));
			}
			args[1] = args[1].substr(args[1].find_first_of(";,") + 1);
		}

		double colorNum = 0;
		Color col(1);
		for (auto f : target_fncs) {
			if (target_fncs.size() > 1) {
				col = RGBColor(colorNum);
			}

			int progress = 0, * pr = &progress, * resol = &resolution;
			bool go = true, * g = &go;
			std::thread counter = std::thread([](int* progress, int* resolution, bool* go) {
				int max_steps = (*resolution) * (*resolution);
				std::cout << " render [" << std::to_string(100.0).substr(1) << "%]\b\b";
				while (*go) {
					float percent = 100 + 100 * (float)(*progress) / max_steps;
					std::cout << "\b\b\b\b\b\b\b\b\b" << std::to_string(percent).substr(1);
				}
				std::cout << "\b\b\b\b\b\b\b\b\b\bcomplete     \n";
			}, pr, resol, g);

			switch (f.args.size()) {
			case(0): {
				Color res = toCol(f.return_value());
				res.R = res.R * res.A;
				res.G = res.G * res.A;
				res.B = res.B * res.A;
				for (int i = 0; i < resolution * resolution; i++) {
					img[i][0] = res.R;
					img[i][1] = res.G;
					img[i][2] = res.B;
				}
				break;
			}
			case(1): {
				if (cmplx) {
					std::vector<std::thread> thr(prc);
					for (int i = 0; i < prc; i++) {
						thr[i] = std::thread([](std::vector<std::vector<float>>* img, int begin, int end, int resolution, Function f, int* progress) {
							long double
								startX = plot_center.R - plot_radius,
								startY = plot_center.i - plot_radius,
								step = 2 * plot_radius / resolution;
							math::complex res;
							for (int iX = begin; iX < end && iX < resolution; iX++) {
								for (int iY = 0; iY < resolution; iY++) {
									f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
									res = f.return_value();
									(*img)[iY * resolution + iX] = toVecF(penAdd((*img)[iY * resolution + iX], toCol(res)));
									(*progress)++;
								}
							}
							}, im, i * resolution / prc, (i + 1) * resolution / prc, resolution, f, pr);
					}
					for (int i = 0; i < prc; i++) thr[i].join();
				}
				else if (eq) {
					std::vector<std::thread> thr(prc);
					for (int i = 0; i < prc; i++) {
						thr[i] = std::thread([](std::vector<std::vector<float>>* img, int begin, int end, int resolution, Function f, int* progress) {
							long double
								startX = plot_center.R - plot_radius,
								startY = plot_center.i - plot_radius,
								step = 2 * plot_radius / resolution;
							math::complex res;
							for (int iX = begin; iX < end && iX < resolution; iX++) {
								for (int iY = 0; iY < resolution; iY++) {
									f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
									res = 1 / (1 + math::abs(f.return_value()));
									(*img)[iY * resolution + iX] = toVecF(penAdd((*img)[iY * resolution + iX], Color(res.R)));
									(*progress)++;
								}
							}
							}, im, i * resolution / prc, (i + 1) * resolution / prc, resolution, f, pr);
					}
					for (int i = 0; i < prc; i++) thr[i].join();
				}
				else {
					std::vector<std::thread> thr(prc);
					for (int i = 0; i < prc; i++) {
						thr[i] = std::thread([](std::vector<std::vector<float>>* img, int begin, int end, int resolution, Function f, int* progress, Color color) {
							int iY, iYp;
							long double start = plot_center.R - plot_radius, step = 2 * plot_radius / resolution;
							math::complex res;
							f.args[0].value = start + step * (begin - 1);
							res = f.return_value();
							iYp = ((res.R - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
							if (resolution <= iYp) iYp = resolution - 1;
							if (iYp < 0) iYp = 0;
							for (int iX = begin; iX < end && iX < resolution; iX++) {
								f.args[0].value = start + step * iX;
								res = f.return_value();
								iY = ((res.R - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
								if (resolution <= iY) iY = resolution - 1;
								if (iY < 0) iY = 0;
								int up, down;
								if (iYp > iY) {
									up = iYp;
									down = iY;
								}
								else {
									up = iY;
									down = iYp;
								}
								if (down < up) down++;
								for (int j = down; j <= up && j < resolution; j++) {
									(*img)[j * resolution + iX] = toVecF(penAdd((*img)[j * resolution + iX], fncColor(color, res.i)));
								}
								iYp = iY;
								(*progress) += resolution;
							}
							}, im, i * resolution / prc, (i + 1) * resolution / prc, resolution, f, pr, col);
					}
					for (int i = 0; i < prc; i++) thr[i].join();
				}
				break;
			}
			case(2): {
				if (cmplx) {
					std::vector<std::thread> thr(prc);
					for (int i = 0; i < prc; i++) {
						thr[i] = std::thread([](std::vector<std::vector<float>>* img, int begin, int end, int resolution, Function f, int* progress) {
							long double
								startX = plot_center.R - plot_radius,
								startY = plot_center.i - plot_radius,
								step = 2 * plot_radius / resolution;
							math::complex res;
							for (int iX = begin; iX < end && iX < resolution; iX++) {
								for (int iY = 0; iY < resolution; iY++) {
									f.args[0].value = startX + step * iX;
									f.args[1].value = startY + step * iY;
									res = f.return_value();
									(*img)[iY * resolution + iX] = toVecF(penAdd((*img)[iY * resolution + iX], toCol(res)));
									(*progress)++;
								}
							}
							}, im, i * resolution / prc, (i + 1) * resolution / prc, resolution, f, pr);
					}
					for (int i = 0; i < prc; i++) thr[i].join();
				}
				else {
					std::vector<std::thread> thr(prc);
					for (int i = 0; i < prc; i++) {
						thr[i] = std::thread([](std::vector<std::vector<float>>* img, int begin, int end, int resolution, Function f, int* progress) {
							long double
								startX = plot_center.R - plot_radius,
								startY = plot_center.i - plot_radius,
								step = 2 * plot_radius / resolution;
							math::complex res;
							for (int iX = begin; iX < end && iX < resolution; iX++) {
								for (int iY = 0; iY < resolution; iY++) {
									f.args[0].value = startX + step * iX;
									f.args[1].value = startY + step * iY;
									res = f.return_value();
									Color c(1, 1, 1, 1 / (1 + math::abs(res)));
									(*img)[iY * resolution + iX] = toVecF(penAdd((*img)[iY * resolution + iX], c));
									(*progress)++;
								}
							}
							}, im, i * resolution / prc, (i + 1) * resolution / prc, resolution, f, pr);
					}
					for (int i = 0; i < prc; i++) thr[i].join();
				}
				break;
			}
			}
			go = false;
			counter.join();
			colorNum++;
		}

		if (grid) {
			for (int i = 1; i < 2 * plot_radius; i++) {
				long double x = math::floor(plot_center.R - plot_radius) + i;
				int iX = ((x - plot_center.R) / plot_radius + 1) * 0.5 * resolution;
				for (int j = 0; j < resolution; j++) {
					if (x == 0) img[j * resolution + iX] = toVecF(penAdd(img[j * resolution + iX], Color(0.5)));
					else		img[j * resolution + iX] = toVecF(penAdd(img[j * resolution + iX], Color(0.25)));
				}
			}
			for (int i = 1; i < 2 * plot_radius; i++) {
				long double y = math::floor(plot_center.i - plot_radius) + i;
				int iY = ((y - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
				for (int j = 0; j < resolution; j++) {
					if (y == 0) img[iY * resolution + j] = toVecF(penAdd(img[iY * resolution + j], Color(0.5)));
					else		img[iY * resolution + j] = toVecF(penAdd(img[iY * resolution + j], Color(0.25)));
				}
			}
		}

		BMPWriter::write_image(img, resolution, (char*)"graph.bmp");
		return true;
	}
	if (args[0] == "help") {
		std::string response = "";
		if (args.size() == 1) {
			response += " Choose a Category\n";
			response += "   >help operators\n";
			response += "   >help commands\n";
			response += "   >help functions\n";
			response += "   >help syntax <function name>\n";
			std::cout << response;
			return true;
		}
		if (args.size() >= 2) {
			if (args[1] == "operators") {
				response += " => x+y -> summation\n";
				response += " => x-y -> subtraction\n";
				response += " =>  -y -> opposite\n";
				response += " => x*y -> multiplication\n";
				response += " => x/y -> division\n";
				response += " => x^y -> power\n";
				response += " => x%y -> modulus (remainder of the division)\n";
				response += " => x!  -> factorial\n";
				std::cout << response;
				return true;
			}
			if (args[1] == "commands") {
				response += " Define\n";
				response += "  defining a constant\n";
				response += "   >def n <expression>\n";
				response += "   >def n = <expression>\n";
				response += "  defining a function\n";
				response += "   >def f(x,...) <expression>\n";
				response += "   >def f(x,...) = <expression>\n";

				response += "\n Delete\n";
				response += "  >del <function_name>   -> deletes the function\n";

				response += "\n List\n";
				response += "  >list    -> shows functions & constants list\n";

				response += "\n Clear\n";
				response += "  >clear   -> clears the list\n";

				response += "\n Scale\n";
				response += "  >scale <center> <radius>   -> defines position & scale of plot\n";
				response += "     center is a complex number (\"1-2i\" translates to (1;-2))\n";

				response += "\n Print\n";
				response += "  >print <functions:\"f;g;h\"|\"f,g,h\"> <*resolution> <*modificators>\n";
				response += "     resolution is one integer number (400 means image with size 400x400)\n";
				response += "     modificators:\n";
				response += "     =>     eq -> equation mode (f(z) = 0)\n";
				response += "     =>      c -> complex spectrum mode (abs(z) - brightness, arg(z) - hue)\n";
				response += "     => noGrid -> disables the grid\n";

				response += "\n Help\n";
				response += "  >help <*category>   :D\n";
				std::cout << response;
				return true;
			}
			if (args[1] == "functions") {
				response += " How To Call Functions\n";
				response += " =>           n -> n = const\n";
				response += " =>        f(x) -> f(x, ...)\n";
				response += " =>     f[x](y) -> f(x, y)\n";
				response += " => f{x;y;z}[w] -> f(x, y, z, w)\n";

				response += "\n Power Functions\n";
				response += " =>      exp(x) -> exponent\n";
				response += " =>       ln(x) -> natural logarithm\n";
				response += " =>     sqrt(x) -> square root\n";
				response += " => inv_sqrt(x) -> inverse square root\n";
				response += " =>  root[y](x) -> y-th root of x\n";
				response += " =>   log[y](x) -> logarithm of x to base y\n";

				response += "\n Trigonometric Functions\n";
				response += " =>     cos(x)       cot(x)\n";
				response += " =>    cosh(x)      coth(x)\n";
				response += " =>  arccos(x)    arccot(x)\n";
				response += " => arccosh(x)   arccoth(x)\n";
				response += " =>     sin(x)       tan(x)\n";
				response += " =>    sinh(x)      tanh(x)\n";
				response += " =>  arcsin(x)    arctan(x)\n";
				response += " => arcsinh(x)   arctanh(x)\n";
				response += " =>    sin1(x)\n";
				response += " =>    cos1(x)\n";

				response += "\n FOR-Like Functions\n";
				response += " =>    S{t;begin;end}[f(t)] -> Sum\n";
				response += " =>           FD{t;x}[f(t)] -> Forward Difference\n";
				response += " =>         FD{t;x;n}[f(t)] -> n-th Forward Difference\n";
				response += " =>           BD{t;x}[f(t)] -> Backward Difference\n";
				response += " =>         BD{t;x;n}[f(t)] -> n-th Backward Difference\n";
				response += " =>    P{t;begin;end}[f(t)] -> Product\n";
				response += " =>          R{t;x;n}[f(t)] -> Return\n";
				response += " =>          I{t;a;b}[f(t)] -> Integral\n";
				response += " =>            D{t;x}[f(t)] -> Derivative\n";
				response += " =>          D{t;x;n}[f(t)] -> n-th Derivative\n";
				response += " =>       Iexp{t;a;b}[f(t)] -> Integral along exp\n";
				response += " => Poly{t;f(t);x}[a;b;...] -> Polynomial\n";

				response += "\n Value Functions\n";
				response += " =>     abs(x) -> absolute value\n";
				response += " => inv_abs(x) -> inverse absolute value\n";
				response += " =>     arg(x) -> argument of x\n";
				response += " =>    sign(x) -> normalized number\n";
				response += " =>      Re(x) -> real part\n";
				response += " =>      Im(x) -> imaginary part\n";
				response += " =>   floor(x)\n";
				response += " =>    ceil(x)\n";
				response += " =>   round(x)\n";

				response += "\n Logic Functions\n";
				response += " =>   exist(x) -> if number exist returns 0 else 1\n";
				response += " =>    grid(x) -> if number on the grid returns 0 else (0 < ...)\n";

				response += "\n Other Functions\n";
				response += " =>   gamma(x) -> gamma function\n";
				response += " =>  fctI(x,n) -> n-th integral of factorial\n";
				response += " => rand | rand(...)\n";

				response += "\n Experimental Functions\n";
				response += " =>   zbf(x) -> zeta(x)/x!\n";

				std::cout << response;
				return true;
			}
			if (args.size() >= 3) {
				if (args[1] == "syntax") {
					ObjType func = nameToType(args[2]);
					switch (func) {
					case(_exp): {
						response += " exp(x)\n";
						response += "   Экспонента\n";
						response += "   -> e^x\n";
						response += " Пример:\n";
						response += "   exp(1) = e = 2.718282\n";
						response += "   exp(0) = 1\n";
						std::cout << response;
						return true;
					}
					case(_ln): {
						response += " ln(x)\n";
						response += "   Натуральный логарифм x\n";
						response += "   -> log[e](x)\n";
						std::cout << response;
						return true;
					}
					case(_sqrt): {
						response += " sqrt(x)\n";
						response += "   Квадратный корень из x\n";
						response += "   -> root[2](x)\n";
						response += " Пример:\n";
						response += "   sqrt(4) = 2\n";
						response += "   sqrt(9) = 3\n";
						std::cout << response;
						return true;
					}
					case(_inv_sqrt): {
						response += " inv_sqrt(x)\n";
						response += "   Обратный квадратный корень из x\n";
						response += "   -> 1/sqrt(x)\n";
						std::cout << response;
						return true;
					}
					case(_root): {
						response += " root[y](x)\n";
						response += "   Корень степени y из x\n";
						response += "   -> x^(1/y)\n";
						response += " Пример:\n";
						response += "   root[3](8) = 2\n";
						std::cout << response;
						return true;
					}
					case(_log): {
						response += " log[y](x)\n";
						response += "   Логарифм x по основанию y\n";
						response += "   -> ln(x)/ln(y)\n";
						response += " Пример:\n";
						response += "   log[2](8) = 3\n";
						std::cout << response;
						return true;
					}

					case(_sin1): {
						response += " sin1(x)\n";
						response += "   -> sin(πx)/π\n";
						std::cout << response;
						return true;
					}
					case(_cos1): {
						response += " cos1(x)\n";
						response += "   -> cos(πx)\n";
						response += " Пример:\n";
						response += "   можно использовать вместо (-1)^n\n";
						response += "   cos1(0) = 1\n";
						response += "   cos1(1) = -1\n";
						response += "   cos1(2) = 1\n";
						response += "   cos1(3) = -1\n";
						std::cout << response;
						return true;
					}

					case(_Sum): {
						response += " S{k;a;b}[f(k)]\n";
						response += "   Сумма по k от a до b\n";
						response += "      k -> переменная\n";
						response += "      a -> нижний предел\n";
						response += "      b -> верхний предел\n";
						response += "   f(k) -> функция от переменной k\n";
						response += " Пример:\n";
						response += "   S{k;1;5}[k] = 1 + 2 + 3 + 4 + 5 = 15\n";
						response += "   S{k;1;4}[k^2] = 1^2 + 2^2 + 3^2 + 4^2 = 30\n";
						std::cout << response;
						return true;
					}
					case(_Product): {
						response += " P{k;a;b}[f(k)]\n";
						response += "   Произведение по k от a до b\n";
						response += "      k -> переменная\n";
						response += "      a -> нижний предел\n";
						response += "      b -> верхний предел\n";
						response += "   f(k) -> функция от переменной k\n";
						response += " Пример:\n";
						response += "   P{k;1;4}[k] = 1 * 2 * 3 * 4 = 24\n";
						std::cout << response;
						return true;
					}
					case(_Return): {
						response += " R{t;a;n}[f(t)]\n";
						response += "   Итерирует функцию f(t)\n";
						response += "      t -> переменная\n";
						response += "      a -> начальное значение переменной t\n";
						response += "      n -> количество итераций\n";
						response += "   f(t) -> итерируемая функция\n";
						response += "   R{t;a;n}[f(t)] = f(f(...f(a)...))\n";
						response += " Пример:\n";
						response += "   R{t;13;0}[t-5] = 13\n";
						response += "   R{t;13;3}[t-5] = (((13-5)-5)-5) = -2\n";
						std::cout << response;
						return true;
					}
					case(_Integral): {
						response += " I{t;a;b}[f(t)]\n";
						response += "   Интеграл функции f(t)\n";
						response += "      t -> переменная\n";
						response += "      a -> нижний предел\n";
						response += "      b -> верхний предел\n";
						response += "   f(t) -> функция от переменной t\n";
						response += " Пример:\n";
						response += "   I{t;0;1}[t] = 0.5\n";
						std::cout << response;
						return true;
					}
					case(_IntegralAlongExp): {
						response += " Iexp{t;a;b}[f(t)]\n";
						response += "   Криволинейный интеграл функции f(t)\n";
						response += "      t -> переменная\n";
						response += "      a -> нижний предел\n";
						response += "      b -> верхний предел\n";
						response += "   f(t) -> функция от переменной t\n";
						response += " Пример:\n";
						response += "   Iexp{t;1;-1}[1/t] = πi\n";
						std::cout << response;
						return true;
					}
					case(_Derivative): {
						response += " D{t;x}[f(t)]\n";
						response += "   Производная функции f(t) в точке x\n";
						response += " D{t;x;n}[f(t)]\n";
						response += "   n-ная производная функции f(t) в точке x\n";
						response += "      t -> переменная\n";
						response += "      x -> \n";
						response += "      n -> порядок производной\n";
						response += "   f(t) -> функция от переменной t\n";
						response += " Пример:\n";
						response += "   D{t;2}[t^2] = 2t|t=2 = 4\n";
						response += "   D{t;1;2}[t^3] = 6t|t=1 = 6\n";
						std::cout << response;
						return true;
					}
					case(_ForwardDifference): {
						response += " FD{t;x}[f(t)] = f(t+1) - f(t)\n";
						response += "   Правая разность функции f(t) в точке x\n";
						response += " FD{t;x;n}[f(t)]\n";
						response += "   n-ная правая разность функции f(t) в точке x\n";
						response += "      t -> переменная\n";
						response += "      x -> \n";
						response += "      n -> порядок разности\n";
						response += "   f(t) -> функция от переменной t\n";
						response += " Пример:\n";
						response += "   FD{t;2}[t^2] = (t+1)^2-t^2 = 3^2-2^2 = 5\n";
						std::cout << response;
						return true;
					}
					case(_BackwardDifference): {
						response += " BD{t;x}[f(t)] = f(t) - f(t-1)\n";
						response += "   Левая разность функции f(t) в точке x\n";
						response += " BD{t;x;n}[f(t)]\n";
						response += "   n-ная левая разность функции f(t) в точке x\n";
						response += "      t -> переменная\n";
						response += "      x -> \n";
						response += "      n -> порядок разности\n";
						response += "   f(t) -> функция от переменной t\n";
						response += " Пример:\n";
						response += "   BD{t;2}[t^2] = t^2-(t-1)^2 = 2^2-1^2 = 3\n";
						std::cout << response;
						return true;
					}
					case(_Polynomial): {
						response += " Poly{t;f(t);x}[...]\n";
						response += "    Полиноминальная аппроксимация функции f(t) в точке x\n";
						response += "      t -> переменная\n";
						response += "   f(t) -> функция от переменной t\n";
						response += "      x -> \n";
						response += "    ... -> список значений переменной t\n";
						response += " Пример:\n";
						response += "   >def f(x) = Poly{t;cos(t);x}[1;2;3;4;5;6]\n";
						response += "   -> f(x) = ax^5 + bx^4 + cx^3 + ...\n";
						std::cout << response;
						return true;
					}

					default: {
						response += "  unknown function\n";
						std::cout << response;
						return true;
					}
					}
				}
			}
		}
		return true;
	}

	return false;
}

int main() {
	SetConsoleOutputCP(CP_UTF8);
	std::string expression;
	for (;;) {
		std::getline(std::cin, expression);
		if (expression[0] == '>') {
			if (run_command(expression.substr(1)))	std::cout << " ----> done";
			else									std::cout << " --x-> unknown command";
		}
		else {
			math::number answer = value(parse_expr(expression), {});
			std::cout << " -> " << answer.toString();
		}
		std::cout << "\n\n";
	}
}