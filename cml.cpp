#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <tuple>
#include <algorithm>
#include <set>
#include <queue>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include "cml.h"
#include "Valences.h"

using boost::property_tree::ptree;
using std::map;
using std::set;
using std::queue;
using std::string;
using std::pair;
using std::vector;
using std::ostringstream;
using std::logic_error;

class Atom
{
	string id;
	string element;
	int valence; // 0 is it must be deduced from bonds number and standard valences
	int bondsCnt;   // to accumulate number of explicit bonds
	int charge;
public:
	Atom(const string& id, const string& element, int valence = 0, int charge = 0) : id(id), element(element), valence(valence), bondsCnt(0), charge(charge)
	{
	}

	void addBonds(int order)
	{
		bondsCnt += order;

		if (valence != 0 && bondsCnt > valence)
		{
			ostringstream oss;
			oss << "Atom with id '" << id << "' has more bounds than valence.";
			throw logic_error(oss.str());
		}
	}

	// if valence == 0 then deduce it from bonds and standard valences for that element
	void correct_valence()
	{
		if (this->valence == 0)
			this->valence = Valences::getCeilValence(this->element, this->bondsCnt);
	}

	int getImplicitHydrogenCnt() const {
		return valence - bondsCnt + charge;
	}

	string getElement() const { return element; }
};

class Atoms
{
	map<string, Atom> atoms; // atom id => Atom

	// atom id => list of bound atom ids
	map<string, vector<string> > graph;

	void addAtom(const string& id, const string& element, int valence = 0, int charge = 0)
	{
		if (atoms.find(id) != atoms.end())
		{
			ostringstream oss;
			oss << "Atom with id '" << id << "' already exists.";
			throw logic_error(oss.str());
		}

		atoms.insert(make_pair(id, Atom(id, element, valence, charge)));
		graph.insert(make_pair(id, vector<string>()));
	}

	void addBond(const string& id1, const string& id2, int order)
	{
		graph[id1].push_back(id2);
		graph[id2].push_back(id1);
		Atom& a1 = atoms.find(id1)->second;
		Atom& a2 = atoms.find(id2)->second;
		a1.addBonds(order);
		a2.addBonds(order);
	}

	vector<vector<string>> get_connected_components() const
	{
		vector<vector<string> > connected_components;
		set<string> visited;
		
		for (auto it = graph.cbegin(); it != graph.cend(); ++it)
		{
			const string& currentAtomID = it->first;
			
			if (visited.find(currentAtomID) == visited.end())
			{
				vector<string> component;
				queue<string> q;
				q.push(currentAtomID);
				visited.insert(currentAtomID);
				component.push_back(currentAtomID);

				while (!q.empty())
				{
					string neighborID = q.front();
					q.pop();
					const vector<string>& adjList = graph.find(neighborID)->second;

					for (auto adjIt = adjList.cbegin(); adjIt != adjList.cend(); ++adjIt)
					{
						if (visited.find(*adjIt) == visited.end())
						{
							visited.insert(*adjIt);
							q.push(*adjIt);
							component.push_back(*adjIt);
						}
					}
				}

				connected_components.push_back(component);
			}
		}

		return connected_components;
	}

	static string getFormula(const map<string, int>& molecule)
	{
		std::ostringstream oss;
		std::map<string, int>::const_iterator tC, tH;
		tH = molecule.end();
		tC = molecule.find("C");
		if(tC != molecule.end())
		{
			const string& element = tC->first;
			int cnt = tC->second;
			
			oss << element;

			if (cnt > 1)
				oss << cnt;
			
			tH = molecule.find("H");
			
			const string& element1 = tH->first;
			cnt = tH->second;
			
			oss << element1;

			if (cnt > 1)
				oss << cnt;
		}
		for (std::map<string, int>::const_iterator el = molecule.begin(); el != molecule.end(); ++el)
		{
			if(el == tC || el == tH)
				continue;
			const string& element = el->first;
			int cnt = el->second;

			oss << element;

			if (cnt > 1)
				oss << cnt;
		}

		return oss.str();
	}

	// vector of molecules, where molecule = { element => element count }
	vector<map<string, int>> getMolecules() const
	{
		vector<map<string, int>> molecules;

		vector<vector<string>> connected_components = this->get_connected_components();
		for (std::vector<vector<string>>::const_iterator component = connected_components.begin(); component != connected_components.end(); ++component)
		{
			map<string, int> molecule;
			int implicitHydrogenCnt = 0;
			for (std::vector<string>::const_iterator id = (*component).begin(); id != (*component).end(); ++id)
			{
				const Atom& a = atoms.find(*id)->second;
				molecule[a.getElement()]++;
				implicitHydrogenCnt += a.getImplicitHydrogenCnt();
			}

			if (implicitHydrogenCnt)
				molecule["H"] += implicitHydrogenCnt;

			molecules.push_back(molecule);
		}

		return molecules;
	}


public:
	void readFromCML(const string& filename)
	{
		ptree pt;
		try
		{
			read_xml(filename, pt);
		}
		catch(std::exception &e)
		{
			throw string("Error: no file exists");
		}
		std::ifstream chem_atom ("table.txt");
		if(!chem_atom.is_open())
			throw string("Error: source file problem");
		vector<string> atom_tab;
		while(true)
		{
			if(chem_atom.eof())
				break;
			string str;
			getline(chem_atom, str, '\n');
			atom_tab.push_back(str);
		}
		BOOST_FOREACH(ptree::value_type &i, pt.get_child("cml.molecule.atomArray"))
		{
			std::string name = i.first;
			ptree& sub_pt = i.second;
			int valence = 0;
			int charge = 0;

			if (name == "atom")
			{
				string id = sub_pt.get<string>("<xmlattr>.id");
				string elementType = sub_pt.get<string>("<xmlattr>.elementType");
				
				if(find(atom_tab.begin(), atom_tab.end(), elementType) == atom_tab.end())
				{
					throw string("Error: non chemical atom in file");
				}

				boost::optional<int> mrvValence = sub_pt.get_optional<int>("<xmlattr>.mrvValence");
				
				if (mrvValence)
				{
					valence = mrvValence.value();
				}

				boost::optional<int> formalcharge = sub_pt.get_optional<int>("<xmlattr>.formalCharge");
				
				if (formalcharge)
				{
					charge = formalcharge.value();
				}

				boost::optional<int> hydrogenCount = sub_pt.get_optional<int>("<xmlattr>.hydrogenCount");
				
				if (hydrogenCount)
				{
					valence -= hydrogenCount.value();
				}

				this->addAtom(id, elementType, valence, charge);
			}
		}

		BOOST_FOREACH(ptree::value_type &i, pt.get_child("cml.molecule.bondArray"))
		{
			std::string name = i.first;
			ptree& sub_pt = i.second;

			if (name == "bond")
			{
				string atomRefs2 = sub_pt.get<string>("<xmlattr>.atomRefs2");
				int order = sub_pt.get<int>("<xmlattr>.order");

				vector<string> ids;
				boost::split(ids, atomRefs2, boost::is_any_of(" \t\n"));
				if (ids.size() != 2)
				{
					ostringstream oss;
					oss << "Wrong bond '" << atomRefs2 << "'.";
					throw logic_error(oss.str());
				}
				this->addBond(ids.at(0), ids.at(1), order);
			}
		}

	}

	void correct_valences()
	{
		for (std::map<string, Atom>::iterator id_atom = atoms.begin(); id_atom != atoms.end(); ++id_atom)
		{
			id_atom->second.correct_valence();
		}
	}
	
	string getFormula()
	{
		vector<map<string, int>> molecules = this->getMolecules();

		vector<string> molecules_str;
		for (std::vector<map<string, int>>::const_iterator molecule = molecules.begin(); molecule != molecules.end(); ++molecule)
		{
			molecules_str.push_back(getFormula(*molecule));
		}

		string formula = "";
		for (unsigned i = 0; i < molecules_str.size(); ++i)
		{
			if (i > 0)
				formula += ".";
			formula += molecules_str[i];
		}
		return formula;
	}
};

string getFormulaFromCML(const string& xmlfile)
{
	Atoms atoms;
	atoms.readFromCML(xmlfile);
	atoms.correct_valences();
	return atoms.getFormula();
}
