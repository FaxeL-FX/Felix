﻿//	v1.6.2

#include <iostream>
#include "include.h"
#include "object/object.h"
#include "bmp/BMPWriter.h"
#include <thread>

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
	a.R = a.R * a.A + b.R * b.A;
	a.G = a.G * a.A + b.G * b.A;
	a.B = a.B * a.A + b.B * b.A;
	if (a.R > 1) a.R = 1;
	if (a.G > 1) a.G = 1;
	if (a.B > 1) a.B = 1;
	return a;
}
Color toCol(math::complex_linear x) {
	Color c;
	float angle = math::arg(x), absolute = math::abs(x);
	c.A = 1;
	c.R = (math::cos(angle) + 1) * 0.5;
	c.G = (math::cos(angle + 2.0943951) + 1) * 0.5;
	c.B = (math::cos(angle - 2.0943951) + 1) * 0.5;
	c.R = 1 - 0.5 / (0.5 + c.R * c.R * absolute);
	c.G = 1 - 0.5 / (0.5 + c.G * c.G * absolute);
	c.B = 1 - 0.5 / (0.5 + c.B * c.B * absolute);
	return c;
}
Color toCol(math::complex_exponential x) {
	Color c;
	c.A = 1;
	c.R = (math::cos(x.a) + 1) * 0.5;
	c.G = (math::cos(x.a + 2.0943951) + 1) * 0.5;
	c.B = (math::cos(x.a - 2.0943951) + 1) * 0.5;
	c.R = 1 - 0.5 / (0.5 + c.R * c.R * x.r);
	c.G = 1 - 0.5 / (0.5 + c.G * c.G * x.r);
	c.B = 1 - 0.5 / (0.5 + c.B * c.B * x.r);
	return c;
}
std::vector<float> toVecF(Color c) {
	return { c.R, c.G, c.B };
}

math::complex_linear plot_center = 0;
math::number plot_radius = 4;

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
				/**/ if (-1 < index1) f.args.push_back(Variable(f.args.size(), arguments.substr(0, index1), 0));
				else if (-1 < index2) f.args.push_back(Variable(f.args.size(), arguments.substr(0, index2), 0));
				else /*------------*/ f.args.push_back(Variable(f.args.size(), arguments, 0));
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
		for (int i = 0; i < f.objects.size(); i++) {
			for (auto arg : f.args)
				if (f.objects[i].name == arg.name) {
					f.objects[i].argID = arg.id;
					break;
				}
		}

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
		bool monochrome = false, grid = true, cmplx = false, eq = false;
		if (args.size() < 3) return false;
		for (int i = 3; i < args.size(); i++) {
			/**/ if (args[i] == "m") monochrome = true;
			else if (args[i] == "noGrid") grid = false;
			else if (args[i] == "c") cmplx = true;
			else if (args[i] == "eq") eq = true;
		}
		Function f;
		for (auto fnc : functions_list)
			if (fnc.name == args[1]) {
				f = fnc;
				break;
			}

		int resolution = std::stoi(args[2]);
		int prc = resolution + (prc_count < resolution) * (prc_count - resolution);
		std::vector<std::vector<float>> img;
		std::vector<std::vector<float>>* im = &img;
		for (int i = 0; i < resolution * resolution; i++) {
			img.push_back({ 0, 0, 0 });
		}
		if (grid) {
			for (int i = 1; i < 2 * plot_radius; i++) {
				math::number x = math::floor(plot_center.R - plot_radius) + i;
				int iX = ((x - plot_center.R) / plot_radius + 1) * 0.5 * resolution;
				for (int j = 0; j < resolution; j++) {
					if (x == 0) img[j * resolution + iX] = { 0.5,  0.5,  0.5  };
					else		img[j * resolution + iX] = { 0.25, 0.25, 0.25 };
				}
			}
			for (int i = 1; i < 2 * plot_radius; i++) {
				math::number y = math::floor(plot_center.i - plot_radius) + i;
				int iY = ((y - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
				for (int j = 0; j < resolution; j++) {
					if (y == 0) img[iY * resolution + j] = { 0.5,  0.5,  0.5 };
					else		img[iY * resolution + j] = { 0.25, 0.25, 0.25 };
				}
			}
		}

		int progress = 0, * pr = &progress, * resol = &resolution;
		bool go = true, *g = &go;
		std::thread counter = std::thread([](int* progress, int* resolution, bool* go) {
			int max_steps = (*resolution) * (*resolution);
			std::cout << " render [" << std::to_string(100.0).substr(1) << "%]\b\b";
			while(*go) {
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
						math::number
							startX = plot_center.R - plot_radius,
							startY = plot_center.i - plot_radius,
							step = 2 * plot_radius / resolution;
						math::complex res;
						for (int iX = begin; iX < end && iX < resolution; iX++) {
							for (int iY = 0; iY < resolution; iY++) {
								f.args[0].value = math::complex_linear(startX + step * iX, startY + step * iY);
								res = f.return_value();
								(*img)[iY * resolution + iX] = toVecF(penAdd((*img)[iY * resolution + iX], toCol(res)));
								(*progress)++;
							}
						}
					}, im, i* resolution / prc, (i + 1)* resolution / prc, resolution, f, pr);
				}
				for (int i = 0; i < prc; i++) thr[i].join();
			}
			else if (eq) {
				std::vector<std::thread> thr(prc);
				for (int i = 0; i < prc; i++) {
					thr[i] = std::thread([](std::vector<std::vector<float>>* img, int begin, int end, int resolution, Function f, int* progress) {
						math::number
							startX = plot_center.R - plot_radius,
							startY = plot_center.i - plot_radius,
							step = 2 * plot_radius / resolution;
						math::complex_linear res;
						for (int iX = begin; iX < end && iX < resolution; iX++) {
							for (int iY = 0; iY < resolution; iY++) {
								f.args[0].value = math::complex_linear(startX + step * iX, startY + step * iY);
								res = 1 / (1 + math::abs(f.return_value()));
								(*img)[iY * resolution + iX] = toVecF(penAdd((*img)[iY * resolution + iX], Color(res.R, res.R, res.R, 1)));
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
						int iY, iYp;
						math::number start = plot_center.R - plot_radius, step = 2 * plot_radius / resolution;
						math::complex_linear res;
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
								(*img)[j * resolution + iX] = toVecF(penAdd((*img)[j * resolution + iX], Color(1 / (1 + res.i * res.i))));
							}
							iYp = iY;
							(*progress) += resolution;
						}
					}, im, i * resolution / prc, (i + 1) * resolution / prc, resolution, f, pr);
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
						math::number
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
						math::number
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
		default: return false;
		}
		go = false;
		counter.join();
		BMPWriter::write_image(img, resolution, (char*)"graph.bmp");
		return true;
	}
	if (args[0] == "help") {
		std::string response = "";
		if (args.size() == 1) {
			response += " Choose a Category\n";
			response += "   >help operators\n";
			response += "   >help functions\n";
			response += "   >help commands\n";
			std::cout << response;
			return true;
		}
		if (args.size() == 2) {
			if (args[1] == "operators") {
				response += " => x+y -> sum\n";
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
				response += " =>     cos(x)       cot(x) = ctg(x)\n";
				response += " =>    cosh(x)      coth(x) = ctgh(x)\n";
				response += " =>  arccos(x)    arccot(x) = arcctg(x)\n";
				response += " => arccosh(x)   arccoth(x) = arcctgh(x)\n";
				response += " =>     sin(x)       tan(x) = tg(x)\n";
				response += " =>    sinh(x)      tanh(x) = tgh(x)\n";
				response += " =>  arcsin(x)    arctan(x) = arctg(x)\n";
				response += " => arcsinh(x)   arctanh(x) = arctgh(x)\n";
				response += "\n Other Functions\n";
				response += " =>   gamma(x) -> gamma function\n";
				response += " =>    fctI(x) -> integral of factorial\n";
				response += " =>     abs(x) -> absolute value\n";
				response += " => inv_abs(x) -> inverse absolute value\n";
				response += " =>     arg(x) -> argument of x\n";
				response += " =>      Re(x) -> real part\n";
				response += " =>      Im(x) -> imaginary part\n";
				response += " =>   floor(x)\n";
				response += " =>   round(x)\n";
				response += " =>   exist(x) -> if number exist returns 0 else 1\n";
				response += " =>    grid(x) -> if number on the grid returns 0 else (0 < ...)\n";
				response += " => S{var;begin;end}[f(var)] = Sum{var;begin;end}[f(var)]\n";
				response += " => P{var;begin;end}[f(var)] = Product{var;begin;end}[f(var)]\n";
				response += " =>       R{var;x;n}[f(var)] = Return{var;x;n}[f(var)]\n";
				response += " =>         D{var;x}[f(var)] = Derivative{var;x}[f(var)]\n";
				response += " =>       I{var;a;b}[f(var)] = Integral{var;a;b}[f(var)]\n";
				response += " =>    Iexp{var;a;b}[f(var)] = IntegralAlongExp{var;a;b}[f(var)]\n";
				response += " => rand | rand(x) | rand(x, ...)\n";
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
				response += "     center is a complex number (1-2i translates to (1;-2))\n";
				response += "\n Print\n";
				response += "  >print <function_name> <resolution> <*modificators>\n";
				response += "     resolution is one integer number (400 means image with size 400x400)\n";
				response += "     modificators:\n";
				response += "     =>     eq -> equation mode (f(z) = 0)\n";
				response += "     =>      c -> complex spectrum mode (abs(z) - brightness, arg(z) - hue)\n";
				response += "     => noGrid -> disables the grid\n";
				std::cout << response;
				return true;
			}
		}
		return true;
	}

	return false;
}

int main() {
	std::string expression;
	for (;;) {
		std::getline(std::cin, expression);
		if (expression[0] == '>') {
			if (run_command(expression.substr(1)))	std::cout << " ----> done";
			else									std::cout << " --x-> unknown command";
		}
		else {
			math::complex answer = value(parse_expr(expression), {});
			std::cout << " -> " << math::toString((math::complex_linear)answer) << " = " << math::toString((math::complex_exponential)answer);
		}
		std::cout << "\n\n";
	}
}