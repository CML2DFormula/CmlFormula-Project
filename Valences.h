#include <map>
#include <string>
#include <vector>

// Singlton class
class Valences
{
	//element => vector of valencies in increasing order
	std::map<std::string, std::vector<int> > valencies;

	static Valences* instance;

	Valences(); // private constructor
	static const Valences& getInstance();

public:

	static int getCeilValence(const std::string& element, int lowerValence);
};
