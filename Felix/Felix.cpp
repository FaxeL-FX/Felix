//	v1.15.4

#include "Felix.h"
#include "math/math.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <ratio>
#include <thread>
#include <chrono>

#include <print>
#include <vector>

const unsigned int prc_count = std::thread::hardware_concurrency();

namespace
{

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
	Color(std::array<float, 3> rgb) {
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

Color RGBColor(double col) {
	return Color(
		(2 + cos(col * 1.8)) / 3,
		(2 + cos(col * 1.8 + 2.095)) / 3,
		(2 + cos(col * 1.8 - 2.095)) / 3,
		1);
}

struct RenderStripTask
{
	std::function<void(RenderStripTask &)> render_function;

	void do_work()
	{
		if (done) {
			return;
		}
		render_function(*this);
		done_function();
		return;
	}

	bool done = false;
	std::vector<Function> target_fncs;
	size_t target_fnc_idx = 0;

	int resolution;
	Color color;
	math::number parameter;

	size_t strip_width;
	size_t strip_index;

	std::vector<Function> functions_list;
	math::complex plot_center;
	long double plot_radius = 8;
	bool gif = false;
	std::vector<std::array<float, 3>> strip;

	float progress = 0;

	void done_function()
	{
		target_fnc_idx++;
		progress = static_cast<float>(target_fnc_idx) / target_fncs.size();
		if (target_fnc_idx == target_fncs.size()) {
			done = true;
			return;
		}
		color = RGBColor(target_fnc_idx);
	}
};

struct ThreadedRenderStripTask
{
	ThreadedRenderStripTask(std::reference_wrapper<RenderStripTask> task) : _task{task} {}

	bool do_work()
	{
		Guard in_work_guard{_mutex};
		if (!in_work_guard.own) return false;
		_task.get().do_work();
		std::atomic_ref{_progress}.store(_task.get().progress);
		if (_task.get().done) std::atomic_ref{_done}.store(true);
		return true;
	}

	float progress()
	{
		return std::atomic_ref{_progress}.load();
	}

	bool done()
	{
		return std::atomic_ref{_done}.load();
	}

private:
	struct Guard : std::atomic_ref<bool>
	{
		Guard(bool &mutex) : std::atomic_ref<bool>{mutex}
		{
			own = !exchange(true);
		}
		~Guard()
		{
			if(!own) return;
			store(false);
		}
		bool own;
	};

	std::reference_wrapper<RenderStripTask> _task;
	bool _mutex = false;
	bool _done = false;
	float _progress = 0;
};

struct RenderFrameTask
{
	struct PostProcessInfo
	{
		math::complex parameter;
		size_t resolution;
		size_t strip_width;
		math::complex plot_center;
		long double plot_radius = 8;
	};

	RenderFrameTask(std::vector<RenderStripTask> &&strip_tasks, PostProcessInfo i) :
		_i{i},
		_strip_tasks{std::move(strip_tasks)}
	{
		_img.resize(_i.resolution * _i.resolution);
		for (auto &task : _strip_tasks) _strip_task_wraps.emplace_back(std::ref(task));
	}

	bool do_work()
	{
		if (done()) return false;
		if (!std::atomic_ref{_strips_done}.load()) return do_work_strips();

		if (std::atomic_ref{_postprocess_in_progress}.exchange(true)) return false;
		postprocess();

		std::atomic_ref{_postprocess_done}.store(true);
		return true;
	}

	float progress()
	{
		float total = 0;
		for (auto &task : _strip_task_wraps) total += task.progress();
		total /= _strip_task_wraps.size();
		total *= 0.8;
		total += std::atomic_ref{_postprocess_progress}.load() * 0.2;
		return total;
	}

	bool done()
	{
	    return std::atomic_ref{_postprocess_done}.load();
	}

	std::vector<std::function<void(std::vector<std::array<float, 3>> &img, RenderFrameTask::PostProcessInfo &I)>> postprocess_functions;

	auto &image()
	{
		return _img;
	}

private:
	PostProcessInfo _i;
	std::vector<RenderStripTask> _strip_tasks;

	std::vector<std::array<float, 3>> _img;
	std::vector<ThreadedRenderStripTask> _strip_task_wraps;
	bool _strips_done = false;
	bool _postprocess_done = false;
	bool do_work_strips()
	{
		bool all_done = true;
		for (auto &task : _strip_task_wraps) {
			if (task.done()) continue;
			all_done = false;
			if (task.do_work()) return true;
		}
		if (all_done) std::atomic_ref{_strips_done}.store(true);
		return false;
	}

	float _postprocess_progress = 0;
	bool _postprocess_in_progress = false;
	void postprocess()
	{
		std::atomic_ref pp{_postprocess_progress};
		merge_strips();
		pp = 0.333333333333333;
		external_postproc();
		pp = 0.666666666666666;
		std::atomic_ref{_postprocess_done}.store(true);
		pp = 1;
	}

	void merge_strips()
	{
		for (size_t y = 0; y != _i.resolution; ++y) {
			for (size_t x = 0; x != _i.resolution; ++x) {
				auto &image_pixel = _img[y * _i.resolution + x];
				size_t strip_index = x / _i.strip_width;
				auto &strip = _strip_tasks.at(strip_index).strip;
				auto &strip_pixel = strip[y * _i.strip_width + (x % _i.strip_width)];
				image_pixel = strip_pixel;
			}
		}
	}

	size_t _postproc_done = 0;
	void external_postproc()
	{
		std::atomic_ref pp_done{_postproc_done};
		for (auto &f : postprocess_functions) { f(_img, _i); pp_done++;};
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
float absMul(float c, math::real absolute) {
	float t;
	if (absolute > 1.0) t = 1.0 - 0.5 / std::pow((double)absolute, 0.25);
	else				t = 0.5 * std::pow((double)absolute, 0.25);
	return (1.0 - 2.0 * c) * t * t + 2.0 * c * t;
}
Color toCol(math::complex x) {
	Color c;
	float angle = math::arg(x);
	math::real absolute = math::abs(x);
	c.A = 1;
	c.R = (std::cos(angle) + 1) * 0.5;
	c.G = (std::cos(angle + 2.0943951) + 1) * 0.5;
	c.B = (std::cos(angle - 2.0943951) + 1) * 0.5;
	c.R = absMul(c.R, absolute);
	c.G = absMul(c.G, absolute);
	c.B = absMul(c.B, absolute);
	return c;
}
std::array<float, 3> toVecF(Color c) { return { c.R, c.G, c.B }; }
Color fncColor(Color col, long double im) {
	if (im == im)	col.A = 1 / (1 + im * im);
	else			col.A = 0;
	return col;
}
Color inverse(Color c) {
	c.R = 1.0 - c.R;
	c.G = 1.0 - c.G;
	c.B = 1.0 - c.B;
	return c;
}
} // namespace


struct Felix::PrintTask::Impl
{
	std::vector<RenderFrameTask> frame_tasks;
	uint32_t gif_delay;
	uint32_t resolution;
	bool all_done = false;
};

bool Felix::PrintTask::do_work()
{
	if (std::atomic_ref{I->all_done}.load()) return false;
	bool all_done = true;
	for (auto &task : I->frame_tasks){
		if (task.done()) continue;
		all_done = false;
		if (task.do_work()) return true;
	}
	if (all_done) std::atomic_ref{I->all_done}.store(true);
	return false;
}

float Felix::PrintTask::progress()
{
	float total = 0;
	for (auto &task : I->frame_tasks) total += task.progress();
	total /= I->frame_tasks.size();
	return total;
}

bool Felix::PrintTask::done()
{
	return std::atomic_ref{I->all_done}.load();
}

Felix::ImageList Felix::PrintTask::images()
{
	ImageList images;
	images.gif_delay = I->gif_delay;
	for (auto &task : I->frame_tasks) images.emplace_back(task.image(), I->resolution);
	return images;
}

Felix::PrintTask::~PrintTask() = default;

Felix::CommandResult Felix::run_command(std::string c) {
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
		if (args.size() < 3) return ErrorMessage{"Bad arguments count"};
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
		f.expression = expr;

		if (f.args.size() == 0) {
			auto objs = parse_expr(expr);
			std::vector<Variable> vars = {};
			f.objects = { Object(ObjType::_Const, {}, value(objs, vars, functions_list)) };
		}
		else f.objects = parse_expr(expr);

		for (int i = 0; i < functions_list.size(); i++)
			if (functions_list[i].name == f.name) {
				functions_list[i] = f;
				return std::monostate{};
			}
		functions_list.push_back(f);
		return std::monostate{};
	}
	if (args[0] == "del") {
		if (args.size() < 2) return ErrorMessage{"Bad arguments count"};
		for (int i = 0; i < functions_list.size(); i++)
			if (functions_list[i].name == args[1])
				functions_list.erase(functions_list.begin() + i);
		return std::monostate{};
	}
	if (args[0] == "list") {
		std::string out;

		for (auto f : functions_list) {
			out += " -> " + f.name;
			if (f.args.size() > 0) {
				out += "(";
				for (auto a : f.args) out += a.name += ",";
				out.erase(out.end() - 1, out.end());
				out += ") = " + f.toString() + "\n";
			}
			else out += " = " + f.return_value(functions_list).toString() + "\n";
		}
		return out;
	}
	if (args[0] == "clear") {
		functions_list.clear();
		return std::monostate{};
	}
	if (args[0] == "scale") {
		if (args.size() < 3) return ErrorMessage{"Bad arguments count"};

		auto objs1 = parse_expr(args[1]);
		auto objs2 = parse_expr(args[2]);
		std::vector<Variable> vars = {};

		plot_center = value(objs1, vars, functions_list);
		plot_radius = math::abs(value(objs2, vars, functions_list));
		return std::monostate{};
	}
	if (args[0] == "print") {
		bool grid = true, cmplx = false, eq = false, gif = false, one_thread = false, debug = false, param = false;
		if (args.size() < 2) return ErrorMessage{"Bad arguments count"};

		int resolution = 400;	// bmp/gif
		int frames = 64;		// gif
		uint32_t delay = 5;		// gif
		math::number start_value = (math::real)0, end_value = (math::real)1;	// parameter

		// render options flags
		for (int i = 2; i < args.size(); i++) {
			/**/ if (args[i] == "no_grid") grid = false;
			else if (args[i] == "c") cmplx = true;
			else if (args[i] == "eq") eq = true;
			else if (args[i] == "gif") gif = true;
			else if (args[i] == "debug") debug = true;
			else if (args[i] == "param") param = true;
			else if (args[i] == "one_thread") one_thread = true;

			else if (args[i].substr(0, 4) == "res=") {	// substr(length)
				resolution = std::stoi(args[i].substr(4));
				if (resolution < 0)	resolution = -resolution;
			}
			else if (args[i].substr(0, 7) == "frames=") {
				frames = std::stoi(args[i].substr(7));
				if (frames < 0)	frames = -frames;
				if (frames < 1)	frames = 1;
			}
			else if (args[i].substr(0, 6) == "delay=") {
				delay = std::stoi(args[i].substr(6));
			}
			else if (args[i].substr(0, 8) == "p_start=") {
				auto objs = parse_expr(args[i].substr(8));
				std::vector<Variable> vars = {};
				start_value = value(objs, vars, functions_list);
			}
			else if (args[i].substr(0, 6) == "p_end=") {
				auto objs = parse_expr(args[i].substr(6));
				std::vector<Variable> vars = {};
				end_value = value(objs, vars, functions_list);
			}
		}

		if (!gif) frames = 1;

		int prc = prc_count;
		if (prc > resolution) prc = resolution;
		if (one_thread) prc = 1;

		std::vector<Function> target_fncs;
		if (args[1] == "all") {
			for (auto f : functions_list) {
				if (f.args.size() > 0) target_fncs.push_back(f);
			}
		}
		else {
			args[1] += ';';
			for (; 0 < args[1].size();) {
				if (args[1][0] == ';' || args[1][0] == ',') break;
				bool no_func = true;
				std::string token = parse_token(args[1]);
				for (auto fnc : functions_list)
					if (fnc.name == token) {
						if (fnc.args.size() == 0) fnc.args.push_back(Variable());
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
		}

		// render functions
		auto render_graph = [](RenderStripTask &C) {
				auto &resolution = C.resolution;
				auto &f = C.target_fncs.at(C.target_fnc_idx);
				auto &color = C.color;
				auto &parameter = C.parameter;

				auto &functions_list = C.functions_list;
				auto &plot_center = C.plot_center;
				auto &plot_radius = C.plot_radius;
				auto &gif = C.gif;

				long long iY, iYp;
				long double start = plot_center.R - plot_radius, step = 2 * plot_radius / resolution;
				math::complex res;
				if (gif && f.args.size() == 2) {
					f.args[1].value = parameter;
				}

				auto compute_point = [&](int64_t offset)
				{
					f.args[0].value = (math::real)(start + step * offset);
					return f.return_value(functions_list);
				};

				auto remap = [&](math::complex P)
				{
					return ((long double)(P.R - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
				};

				std::optional<math::complex> prev;
				for (size_t x = 0; x < C.strip_width; x++) {
					int64_t point_index = C.strip_width * C.strip_index + x;
					if (point_index >= resolution) {
						continue;
					}
					if (!prev.has_value()) prev = compute_point(point_index - 1);
					res = compute_point(point_index);
					iYp = remap(prev.value());
					iY = remap(res);
					int64_t down = std::min(iY, iYp);
					down = std::max<int64_t>(0, down);
					int64_t up = std::max(iY, iYp);
					up = std::max<int64_t>(0, up);
					if (down < up) down++;
					for (size_t y = down; y <= up && y < resolution; y++) {
						auto &pixel = C.strip[y * C.strip_width + x];
						pixel = toVecF(penAdd(pixel, fncColor(color, res.i)));
					}
					prev = res;
				}
			};
		auto render_eq = [](RenderStripTask &C) {
				auto &resolution = C.resolution;
				auto &f = C.target_fncs.at(C.target_fnc_idx);
				auto &parameter = C.parameter;

				auto &functions_list = C.functions_list;
				auto &plot_center = C.plot_center;
				auto &plot_radius = C.plot_radius;
				auto &gif = C.gif;

				long double
					startX = plot_center.R - plot_radius,
					startY = plot_center.i - plot_radius,
					step = 2 * plot_radius / resolution;
				math::complex res;
				for (int iY = 0; iY < resolution; iY++) {
					for (size_t x = 0; x < C.strip_width; x++) {
						size_t iX = C.strip_width * C.strip_index + x;
						auto &pixel = C.strip[iY * C.strip_width + x];
						switch (f.args.size()) {
						case(1):
							f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
							break;
						case(2):
							if (gif) {
								f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
								f.args[1].value = parameter;
							}
							else {
								f.args[0].value = math::complex(startX + step * iX);
								f.args[1].value = math::complex(startY + step * iY);
							}
							break;
						case(3):
							if (gif) {
								f.args[0].value = math::complex(startX + step * iX);
								f.args[1].value = math::complex(startY + step * iY);
								f.args[2].value = parameter;
							}
							break;
						}
						res = math::complex(1 / (1 + math::abs(f.return_value(functions_list))));
						pixel = toVecF(penAdd(pixel, Color(res.R)));
					}
				}
			};
		auto render_complex = [](RenderStripTask &C) {
				std::vector<std::array<float, 3>> strip;
				auto &resolution = C.resolution;
				auto &f = C.target_fncs.at(C.target_fnc_idx);
				auto &parameter = C.parameter;

				auto &functions_list = C.functions_list;
				auto &plot_center = C.plot_center;
				auto &plot_radius = C.plot_radius;
				auto &gif = C.gif;

				long double
					startX = plot_center.R - plot_radius,
					startY = plot_center.i - plot_radius,
					step = 2 * plot_radius / resolution;
				math::complex res;
				for (int iY = 0; iY < resolution; iY++) {
					for (size_t x = 0; x < C.strip_width; x++) {
						size_t iX = C.strip_width * C.strip_index + x;
						auto &pixel = C.strip[iY * C.strip_width + x];
						switch (f.args.size()) {
						case(1):
							f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
							break;
						case(2):
							if (gif) {
								f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
								f.args[1].value = parameter;
							}
							else {
								f.args[0].value = math::complex(startX + step * iX);
								f.args[1].value = math::complex(startY + step * iY);
							}
							break;
						case(3):
							if (gif) {
								f.args[0].value = math::complex(startX + step * iX);
								f.args[1].value = math::complex(startY + step * iY);
								f.args[2].value = parameter;
							}
							break;
						}
						res = f.return_value(functions_list);
						pixel = toVecF(penAdd(pixel, toCol(res)));
					}
				}
			};
		auto render_debug = [](RenderStripTask &C) {
				std::vector<std::array<float, 3>> strip;
				auto &resolution = C.resolution;
				auto &f = C.target_fncs.at(C.target_fnc_idx);
				auto &parameter = C.parameter;

				auto &functions_list = C.functions_list;
				auto &plot_center = C.plot_center;
				auto &plot_radius = C.plot_radius;
				auto &gif = C.gif;

				long double
					startX = plot_center.R - plot_radius,
					startY = plot_center.i - plot_radius,
					step = 2 * plot_radius / resolution;
				math::complex res;
				for (int iY = 0; iY < resolution; iY++) {
					for (size_t x = 0; x < C.strip_width; x++) {
						size_t iX = C.strip_width * C.strip_index + x;
						auto &pixel = C.strip[iY * C.strip_width + x];
						switch (f.args.size()) {
						case(1):
							f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
							break;
						case(2):
							if (gif) {
								f.args[0].value = math::complex(startX + step * iX, startY + step * iY);
								f.args[1].value = parameter;
							}
							else {
								f.args[0].value = math::complex(startX + step * iX);
								f.args[1].value = math::complex(startY + step * iY);
							}
							break;
						case(3):
							if (gif) {
								f.args[0].value = math::complex(startX + step * iX);
								f.args[1].value = math::complex(startY + step * iY);
								f.args[2].value = parameter;
							}
							break;
						}
						auto start = std::chrono::high_resolution_clock::now();
						f.return_value(functions_list);
						auto end = std::chrono::high_resolution_clock::now();
						float dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 100.0;
						dur = dur / (dur + 1);
						pixel = toVecF(penAdd(pixel, Color(dur)));
					}
				}
			};

		std::function<void(RenderStripTask &)>
			render_function = render_graph;
		/**/ if (eq)	render_function = render_eq;
		else if (cmplx) render_function = render_complex;
		else if (debug) render_function = render_debug;

		size_t n_strips = prc;
		size_t strip_width = resolution / prc;
		if (strip_width * prc < resolution) n_strips++;

		Felix::PrintTask print_task;
		print_task.I = std::make_unique<Felix::PrintTask::Impl>();
		auto &I = *print_task.I;
		I.gif_delay = delay;
		I.resolution = resolution;
		std::vector<RenderFrameTask> &frame_tasks = I.frame_tasks;

		for (int frame = 0; frame < frames; frame++) {
			math::number p = start_value;
			if (frames > 1) p = start_value * (math::real)(1 - (double)frame / (frames - 1)) + end_value * (math::real)((double)frame / (frames - 1));

			std::vector<RenderStripTask> strip_tasks;
			for (size_t strip_index = 0; strip_index != n_strips; ++strip_index) {
				RenderStripTask &C = strip_tasks.emplace_back();
				C.render_function = render_function;

				C.strip_width = strip_width;
				C.strip_index = strip_index;
				C.strip.resize(strip_width * resolution);
				std::fill(C.strip.begin(), C.strip.end(), std::array<float,3>{0, 0, 0});

				C.functions_list = functions_list;

				C.plot_center = plot_center;
				C.plot_radius = plot_radius;

				C.gif = gif;
				C.resolution = resolution;
				C.target_fncs = target_fncs;
				C.parameter = p;
				C.target_fnc_idx = 0;
				C.done = false;
				C.color = target_fncs.size() == 1 ? Color{1} : RGBColor(0);
			}

			RenderFrameTask::PostProcessInfo I;
			I.parameter = p;
			I.resolution = resolution;
			I.strip_width = strip_width;
			I.plot_center = plot_center;
			I.plot_radius = plot_radius;

			RenderFrameTask frame_task{std::move(strip_tasks), I};

			if (grid) frame_task.postprocess_functions.push_back([](std::vector<std::array<float, 3>> &img, RenderFrameTask::PostProcessInfo &I) {
				auto &plot_radius = I.plot_radius;
				auto &plot_center = I.plot_center;
				auto &resolution = I.resolution;
				for (int i = 1; i <= 2 * plot_radius; i++) {
					long double x = std::floor(plot_center.R - plot_radius) + i;
					int iX = ((x - plot_center.R) / plot_radius + 1) * 0.5 * resolution;
					for (int j = 0; j < resolution; j++) {
						int index = j * resolution + iX;
						if (index < resolution * resolution) {
							Color gridLine = Color(1.0);
							gridLine.A = 0.125;
							if (x == 0) gridLine.A = 0.25;
							img[index] = toVecF(penAdd(img[index], gridLine));
						}
					}
				}
				for (int i = 1; i <= 2 * plot_radius; i++) {
					long double y = std::floor(plot_center.i - plot_radius) + i;
					int iY = ((y - plot_center.i) / plot_radius + 1) * 0.5 * resolution;
					for (int j = 0; j < resolution; j++) {
						int index = iY * resolution + j;
						if (index < resolution * resolution) {
							Color gridLine = Color(1.0);
							gridLine.A = 0.125;
							if (y == 0) gridLine.A = 0.25;
							img[index] = toVecF(penAdd(img[index], gridLine));
						}
					}
				}
			});
			if (param) frame_task.postprocess_functions.push_back([](std::vector<std::array<float, 3>> &img, RenderFrameTask::PostProcessInfo &I) {
				auto &plot_radius = I.plot_radius;
				auto &plot_center = I.plot_center;
				auto &resolution = I.resolution;
				auto &p = I.parameter;
				int
					scale = 2,
					pos_x = 1 * scale,
					pos_y = resolution - 1 - 1 * scale;
				ImageString s{ p.toString(8) };

				// rectangle
				for (int iX = pos_x - scale; iX < pos_x + scale * 4 * s.str.size(); iX++) {
					for (int iY = pos_y + scale; iY >= pos_y - 6 * scale + 1; iY--) {
						if (iX < 0 || resolution <= iX) continue;
						if (iY < 0 || resolution <= iY) continue;
						img.at(iX + resolution * iY) = toVecF(0.0);
					}
				}
				// number
				for (int i = 0; i < s.str.size(); i++) {
					for (int iX = 0; iX < s.str.at(i).c.size() * scale; iX++) {
						for (int iY = 0; iY < s.str.at(i).c.at(iX / scale).size() * scale; iY++) {
							if (s.str.at(i).c.at(iX / scale).at(iY / scale)) {
								int
									index_x = scale * 4 * i + iX + pos_x,
									index_y = pos_y - iY;
								if (index_x < 0 || resolution <= index_x) continue;
								if (index_y < 0 || resolution <= index_y) continue;
								img.at(index_x + resolution * index_y) = toVecF(1.0);
							}
						}
					}
				}
			});
			frame_tasks.push_back(std::move(frame_task));
		}

		auto worker = [&print_task](){ while(!print_task.done()) print_task.do_work();};
		std::vector<std::thread> thread_pool;
		for (size_t i = 0; i != prc; ++i) thread_pool.emplace_back(worker);

		std::cout << " render [" << std::to_string(100.0).substr(1) << "%]\b\b";
		while (!print_task.done()){
			float percent = 100 + 100 * print_task.progress();
			std::cout << "\b\b\b\b\b\b\b\b\b" << std::to_string(percent).substr(1);
			std::this_thread::sleep_for(std::chrono::milliseconds{10});
		}
		for(auto &T : thread_pool) T.join();
		std::cout << "\b\b\b\b\b\b\b\b\b\bcomplete     \n";

		return print_task.images();
	}
	if (args[0] == "help") {
		std::string response = "";
		if (args.size() == 1) {
			response += " Choose a Category\n";
			response += "   >help operators\n";
			response += "   >help commands\n";
			response += "   >help functions\n";
			response += "   >help syntax <function name>\n";
			return response;
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
				return response;
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
				return response;
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
				response += " =>    norm(x) | normalize(x)\n";

				response += "\n Logic Functions\n";
				response += " =>   exist(x) -> if number exist returns 0 else 1\n";
				response += " =>    grid(x) -> if number on the grid returns 0 else (0 < ...)\n";
				response += " => if{x}[a;b] -> if x is 0 then returns a else returns b\n";

				response += "\n Other Functions\n";
				response += " =>   gamma(x) -> gamma function\n";
				response += " =>  fctI(x,n) -> n-th integral of factorial\n";
				response += " => rand | rand(...)\n";

				response += "\n Experimental Functions\n";
				response += " => Poly{t;f(t);x}[a;b;...] -> Polynomial\n";
				response += " =>                  zbf(x) -> zeta(x)/x!\n";
				response += " =>              USumN(x,n) -> Undefined sum of x^n\n";

				return response;
			}
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
					return response;
				}
				case(_ln): {
					response += " ln(x)\n";
					response += "   Натуральный логарифм x\n";
					response += "   -> log[e](x)\n";
					return response;
				}
				case(_sqrt): {
					response += " sqrt(x)\n";
					response += "   Квадратный корень из x\n";
					response += "   -> root[2](x)\n";
					response += " Пример:\n";
					response += "   sqrt(4) = 2\n";
					response += "   sqrt(9) = 3\n";
					return response;
				}
				case(_inv_sqrt): {
					response += " inv_sqrt(x)\n";
					response += "   Обратный квадратный корень из x\n";
					response += "   -> 1/sqrt(x)\n";
					return response;
				}
				case(_root): {
					response += " root[y](x)\n";
					response += "   Корень степени y из x\n";
					response += "   -> x^(1/y)\n";
					response += " Пример:\n";
					response += "   root[3](8) = 2\n";
					return response;
				}
				case(_log): {
					response += " log[y](x)\n";
					response += "   Логарифм x по основанию y\n";
					response += "   -> ln(x)/ln(y)\n";
					response += " Пример:\n";
					response += "   log[2](8) = 3\n";
					return response;
				}

				case(_sin1): {
					response += " sin1(x)\n";
					response += "   -> sin(πx)/π\n";
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
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
					return response;
				}
				case(_If): {
					response += " if{x}[a;b]\n";
					response += "    Функция ветвления\n";
					response += "      x -> условие (логическое выражение)\n";
					response += "      a -> значение если условие истинно\n";
					response += "      b -> значение если условие ложно\n";
					response += " Пример:\n";
					response += "   >def f(x) = if{x >= 0}[cos(x);sin(x)]   - кусочно-заданная функция\n";
					response += "   >def f(x) = if{x = 0}[1;sin(x)/x]       - устранение разрыва у функции\n";
					return response;
				}

				default: {
					response += "  unknown function\n";
					return response;
				}
				}
			}
		}
		return std::monostate{};
	}

	return ErrorMessage{std::format("Unknown command [{}]", args[0])};
}

std::string Felix::evaluate(std::string expression)
{
	auto objs = parse_expr(expression);
	std::vector<Variable> vars = {};
	math::number answer = value(objs, vars, functions_list);
	return answer.toString();
}
