#pragma once

#include "../include.h"

struct ImageSymbol {
	std::vector<std::vector<bool>> c;

	ImageSymbol(char c);
};

struct ImageString {
	std::vector<ImageSymbol> str;

	ImageString() {}
	ImageString(std::string text);
};