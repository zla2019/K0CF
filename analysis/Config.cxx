#ifndef _CONFIG_
#define _CONFIG_
#include "Config.h"
#include <iostream>
#include <stdlib.h>

Config::Config(std::ifstream* ifConfig)
{
	mIfConfig = ifConfig;
	if(!mIfConfig || !mIfConfig->is_open()) {
		std::cout << "ERROR: Config file load failed" << std::endl;
		exit(1);
	}
	std::cout << "CONFIG: Config file load successful" << std::endl;
	loadConfig();
	printConfig();
}

Config::Config(std::string config)
{
	mIfConfig = new std::ifstream(config);
	if(!mIfConfig || !mIfConfig->is_open()) {
		std::cout << "ERROR: Config file load failed" << std::endl;
		exit(1);
	}
	std::cout << "CONFIG: Config file load successful" << std::endl;
	loadConfig();
	printConfig();
}

Config::~Config()
{
}

void Config::loadConfig()
{
	std::cout << "CONFIG: start load cuts..." << std::endl;
	while(!mIfConfig->eof()) {
		std::string type, name, slower, supper;//slower: string lower, supper: string upper
		*mIfConfig >> type;
		if(type == "Switch:") {
			*mIfConfig >> name;
			*mIfConfig >> slower;
			bool _switch;
			if(slower == "0" || slower == "false") {
				_switch = false;
			} else if(slower == "1" || slower == "true") {
				_switch = true;
			}
			setSwitch(name, _switch);
		} else if(type == "Cut:") {
			*mIfConfig >> name;
			*mIfConfig >> slower >> supper;
			float lower = std::stof(slower), upper = std::stof(supper);
			setCuts(name, lower, upper);
		} else if(type == "THS:") {
			*mIfConfig >> name;
			*mIfConfig >> slower;
			*mIfConfig >> supper;
			float ths = std::stof(slower);
			bool pattern;
			if(supper == "0" || supper == "false") {
				pattern = false;
			} else if(supper == "1" || supper == "true") {
				pattern = true;
			}
			setTHS(name, ths, pattern);
		} else if(type == "SET:") {
			*mIfConfig >> name;
			*mIfConfig >> slower;
			mSetList[name] = slower;
		} else if(type == "") {
			std::cout << "CONFIG: end with empty char \"" << type << "\"" << std::endl;
		} else {
			std::cout << "ERROR: Ilegal type \"" << type << "\"" << std::endl;
			exit(1);
		}
	}
}

void Config::setCuts(std::string cutName, float cutLower, float cutUpper)
{
	mCutList[cutName].first = cutLower;
	mCutList[cutName].second = cutUpper;
}

void Config::setSwitch(std::string switchName, bool _switch)
{
	mSwitchList[switchName] = _switch;
}

void Config::setTHS(std::string thsName, float ths, bool pattern)
{
	mThsList[thsName].first = ths;
	mThsList[thsName].second = pattern;
}

void Config::printConfig()
{
	for(auto cut:mCutList) {
		std::cout << "Cut: " << cut.first << " " << cut.second.first << "~" << cut.second.second << std::endl;
	}
	for(auto _switch:mSwitchList) {
		std::cout << "Switch: " << _switch.first << " " << _switch.second << std::endl;
	}
	for(auto ths:mThsList) {
		std::cout << "THS: " << ths.first << " ths: " << ths.second.first << " pattern: " << ths.second.second << std::endl;
	}
	for(auto set:mSetList) {
		std::cout << "SET: " << set.first << " " << set.second << std::endl;
	}
}

#endif
