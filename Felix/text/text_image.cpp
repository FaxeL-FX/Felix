#include "text_image.h"

ImageSymbol::ImageSymbol(char c) {
	switch (c) {
	case'0':
		this->c = {
			{1,1,1,1,1},
			{1,0,0,0,1},
			{1,1,1,1,1}
		};
		break;
	case'1':
		this->c = {
			{1,0,0,0,1},
			{1,1,1,1,1},
			{0,0,0,0,1}
		};
		break;
	case'2':
		this->c = {
			{1,0,1,1,1},
			{1,0,1,0,1},
			{1,1,1,0,1}
		};
		break;
	case'3':
		this->c = {
			{1,0,1,0,1},
			{1,0,1,0,1},
			{1,1,1,1,1}
		};
		break;
	case'4':
		this->c = {
			{1,1,1,0,0},
			{0,0,1,0,0},
			{1,1,1,1,1}
		};
		break;
	case'5':
		this->c = {
			{1,1,1,0,1},
			{1,0,1,0,1},
			{1,0,1,1,1}
		};
		break;
	case'6':
		this->c = {
			{1,1,1,1,1},
			{1,0,1,0,1},
			{1,0,1,1,1}
		};
		break;
	case'7':
		this->c = {
			{1,0,0,0,0},
			{1,0,0,0,0},
			{1,1,1,1,1}
		};
		break;
	case'8':
		this->c = {
			{1,1,1,1,1},
			{1,0,1,0,1},
			{1,1,1,1,1}
		};
		break;
	case'9':
		this->c = {
			{1,1,1,0,1},
			{1,0,1,0,1},
			{1,1,1,1,1}
		};
		break;
	case'.':
		this->c = {
			{0,0,0,0,0},
			{0,0,0,0,1},
			{0,0,0,0,0}
		};
		break;
	case'+':
		this->c = {
			{0,0,1,0,0},
			{0,1,1,1,0},
			{0,0,1,0,0}
		};
		break;
	case'-':
		this->c = {
			{0,0,1,0,0},
			{0,0,1,0,0},
			{0,0,1,0,0}
		};
		break;
	case'i':
		this->c = {
			{0,0,1,0,0},
			{1,0,1,1,1},
			{1,0,0,0,1}
		};
		break;
	case' ':
		this->c = {
			{0,0,0,0,0},
			{0,0,0,0,0},
			{0,0,0,0,0}
		};
		break;
	default:
		this->c = {
			{1,1,1,1,1},
			{1,1,1,1,1},
			{1,1,1,1,1}
		};
	}
}

ImageString::ImageString(std::string text) {
	for (auto c : text) this->str.push_back(c);
}