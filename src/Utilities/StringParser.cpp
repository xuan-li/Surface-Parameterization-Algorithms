#include "StringParser.h"

void Split(std::string delimiter, std::string s, std::vector<std::string> &result)
{
	result.clear();
	size_t pos = 0;
	std::string token;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		if (token.length() > 0)
			result.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	if (s.length() > 0)
	result.push_back(s);
}
