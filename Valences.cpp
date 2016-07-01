#include "Valences.h"
#include <algorithm>
#include <sstream>

using std::map;
using std::pair;
using std::make_pair;
using std::vector;
using std::string;

Valences* Valences::instance = nullptr;

Valences::Valences()
{
	vector<int> v1;
	v1.push_back(1);

	valencies.insert(make_pair("Li",v1));
	valencies.insert(make_pair("Na",v1));
	valencies.insert(make_pair("K",v1));
	valencies.insert(make_pair("Rb",v1));
	valencies.insert(make_pair("Cs",v1));
	valencies.insert(make_pair("Fr",v1));
	valencies.insert(make_pair("F",v1));

	vector<int> v2;
	v2.push_back(2);

	valencies.insert(make_pair("Be",v2));
	valencies.insert(make_pair("Mg",v2));
	valencies.insert(make_pair("Ca",v2));
	valencies.insert(make_pair("Sr",v2));
	valencies.insert(make_pair("Ba",v2));
	valencies.insert(make_pair("Ra",v2));

	vector<int> v3;
	v3.push_back(3);

	valencies.insert(make_pair("B",v3));
	valencies.insert(make_pair("Al",v3));
	valencies.insert(make_pair("Ga",v3));
	valencies.insert(make_pair("In",v3));

	vector<int> v13;
	v13.push_back(1);
	v13.push_back(3);

	valencies.insert(make_pair("Tl",v13));

	vector<int> v4;
	v4.push_back(4);

	valencies.insert(make_pair("C",v4));
	valencies.insert(make_pair("Si",v4));
	valencies.insert(make_pair("Ge",v4));

	vector<int> v24;
	v24.push_back(2);
	v24.push_back(4);

	valencies.insert(make_pair("Pb",v24));
	valencies.insert(make_pair("Sn",v24));
	valencies.insert(make_pair("O",v24));

	vector<int> v35;
	v35.push_back(3);
	v35.push_back(5);

	valencies.insert(make_pair("N",v35));
	valencies.insert(make_pair("P",v35));
	valencies.insert(make_pair("As",v35));
	valencies.insert(make_pair("Sb",v35));
	valencies.insert(make_pair("Bi",v35));
	
	vector<int> v246;
	v246.push_back(2);
	v246.push_back(4);
	v246.push_back(6);
	
	valencies.insert(make_pair("S",v246));
	valencies.insert(make_pair("Se",v246));
	valencies.insert(make_pair("Te",v246));
	valencies.insert(make_pair("Po",v246));
	
	vector<int> v1357;
	v1357.push_back(1);
	v1357.push_back(3);
	v1357.push_back(5);
	v1357.push_back(7);
	
	valencies.insert(make_pair("Cl",v1357));
	valencies.insert(make_pair("Br",v1357));
	valencies.insert(make_pair("I",v1357));
	valencies.insert(make_pair("At",v1357));
}

const Valences& Valences::getInstance()
{
	if (Valences::instance == nullptr)
	{
		Valences::instance = new Valences();
	}

	return *Valences::instance;
}

int Valences::getCeilValence(const string& element, int lowerValence)
{
	const Valences& inst = Valences::getInstance();

	auto it = inst.valencies.find(element);
		
	if (it == inst.valencies.end())
	{
		return lowerValence;
	}
	else
	{
		auto& valencies_vector = it->second;
			
		auto pos = std::find_if(valencies_vector.begin(), 
								valencies_vector.end(),
								[lowerValence](int x){ return x >= lowerValence; }
								);
		if (pos == valencies_vector.end())
		{
			std::ostringstream oss;
			oss << "Not found valence >= " << lowerValence << " for element " << element << ".";
			throw std::logic_error(oss.str());
		}
		return *pos;
	}
}

