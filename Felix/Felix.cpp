//	v1.5.0

#include <iostream>
#include "include.h"
#include "object/object.h"
#include "bmp/BMPWriter.h"

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
};
Color pen(Color a, Color b) {
	a.R = b.A * b.R + (1 - b.A) * a.R;
	a.G = b.A * b.G + (1 - b.A) * a.G;
	a.B = b.A * b.B + (1 - b.A) * a.B;
	a.A = 1;
	return a;
}
Color penAdd(Color a, Color b) {
	a.R = a.R + b.R;
	a.G = a.G + b.G;
	a.B = a.B + b.B;
	a.A = 1;
	if (a.R > 1) a.R = 1;
	if (a.G > 1) a.G = 1;
	if (a.B > 1) a.B = 1;
	return a;
}
Color toCol(math::complex x) {
	Color c;
	float angle = math::ln(x).i, absolute = 1 - 1 / (1 + math::abs(x));
	c.A = absolute;
	absolute *= absolute;
	c.R = 1 + (math::cos(angle) - 0) * (1 - absolute);
	c.G = 1 + (math::cos(angle + 2.0943951) - 0) * (1 - absolute);
	c.B = 1 + (math::cos(angle - 2.0943951) - 0) * (1 - absolute);
	return c;
}
std::vector<float> toVecF(Color c) {
	return { c.R, c.G, c.B };
}

math::complex plot_center = 0;
math::number plot_radius = 4;

bool run_command(std::string c) {
	std::vector<std::string> args;
	for (;;) {
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

	if (args[0] == "add") {
		if (args.size() < 3) return false;
		Function f;

		int bracket_index = args[1].find('(');
		if (-1 < bracket_index) {
			f.name = args[1].substr(0, bracket_index);
			std::string arguments = args[1].substr(bracket_index + 1);
			while (arguments.length() > 0) {
				int index1 = arguments.find(','),
					index2 = arguments.find(')');
				/**/ if (-1 < index1) f.args.push_back(Variable(arguments.substr(0, index1), 0));
				else if (-1 < index2) f.args.push_back(Variable(arguments.substr(0, index2), 0));
				else /*------------*/ f.args.push_back(Variable(arguments, 0));
				arguments = arguments.substr(f.args[f.args.size() - 1].name.length() + 1);
			}
		}
		else f.name = args[1];

		std::string expr = "";
		for (int i = 2; i < args.size(); i++) expr += args[i];
		f.objects = parse_expr(expr);

		functions_list.push_back(f);
		return true;
	}
	if (args[0] == "remove") {
		if (args.size() < 2) return false;
		if (args[1] == "all") {
			functions_list.clear();
			return true;
		}
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
			std::cout << "\n";
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
				math::number
					startX = plot_center.R - plot_radius,
					startY = plot_center.i - plot_radius,
					step = 2 * plot_radius / resolution;
				math::complex res;
				for (int iX = 0; iX < resolution; iX++) {
					for (int iY = 0; iY < resolution; iY++) {
						f.args[0].value.R = startX + step * iX;
						f.args[0].value.i = startY + step * iY;
						res = f.return_value();
						img[iY * resolution + iX] = toVecF(pen(img[iY * resolution + iX], toCol(res)));
					}
				}
			}
			else if (eq) {
				math::number
					startX = plot_center.R - plot_radius,
					startY = plot_center.i - plot_radius,
					step = 2 * plot_radius / resolution;
				math::complex res;
				for (int iX = 0; iX < resolution; iX++) {
					for (int iY = 0; iY < resolution; iY++) {
						f.args[0].value.R = startX + step * iX;
						f.args[0].value.i = startY + step * iY;
						res = 1 / (1 + math::abs(f.return_value()));
						img[iY * resolution + iX] = toVecF(penAdd(img[iY * resolution + iX], Color(res.R, res.R, res.R, 1)));
					}
				}
			}
			else {
				int iY, iYp;
				math::number start = plot_center.R - plot_radius, step = 2 * plot_radius / resolution;
				math::complex res;
				for (int iX = 0; iX < resolution; iX++) {
					f.args[0].value = start + step * iX;
					res = f.return_value();
					iY = ((res.R - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
					if (resolution <= iY) iY = resolution - 1;
					if (iY < 0) iY = 0;
					if (iX == 0) iYp = iY;
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
						img[j * resolution + iX] = toVecF(pen(img[j * resolution + iX], Color(1, 1, 1, 0.75)));
					}
					iYp = iY;
				}
			}
			break;
		}
		case(2): {
			math::number
				startX = plot_center.R - plot_radius,
				startY = plot_center.i - plot_radius,
				step = 2 * plot_radius / resolution;
			math::complex res;
			if (cmplx) {
				for (int iX = 0; iX < resolution; iX++) {
					for (int iY = 0; iY < resolution; iY++) {
						f.args[0].value = startX + step * iX;
						f.args[1].value = startY + step * iY;
						res = f.return_value();
						img[iY * resolution + iX] = toVecF(pen(img[iY * resolution + iX], toCol(res)));
					}
				}
			}
			else {
				for (int iX = 0; iX < resolution; iX++) {
					for (int iY = 0; iY < resolution; iY++) {
						f.args[0].value = startX + step * iX;
						f.args[1].value = startY + step * iY;
						res = f.return_value();
						Color c(1,1,1, 1 / (1 + math::abs(res)));
						img[iY * resolution + iX] = toVecF(penAdd(img[iY * resolution + iX], c));
					}
				}
			}
			break;
		}
		default: return false;
		}
		BMPWriter::write_image(img, resolution, (char*)"graph.bmp");
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
			std::cout << " -> " << math::toString(value(parse_expr(expression), {}));
		}
		std::cout << "\n\n";
	}
}