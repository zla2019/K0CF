#include <unordered_map>
#include <fstream>

class Config
{
public:
	Config() {};
	Config(std::ifstream* ifConfig);
	Config(std::string config);
	~Config();
	void printConfig();

	//TODO: optimize the interface to use those cuts
	std::unordered_map<std::string, std::pair<float, float>> mCutList;
        std::unordered_map<std::string, bool> mSwitchList;
	std::unordered_map<std::string, std::pair<float, bool>> mThsList;
	std::unordered_map<std::string, std::string> mSetList;
private:
	std::ifstream* mIfConfig;

	void loadConfig();
	void setCuts(std::string cutName, float cutLower, float cutUpper);
	void setSwitch(std::string switchName, bool _switch);
	void setTHS(std::string thsName, float ths, bool pattern);
};
