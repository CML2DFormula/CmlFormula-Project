#include "cml.h"
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

int main(int argc, char** argv)
{
	if(argc > 1)
	{
		std::string filename = argv[1];
		
		try
		{
			string formula = getFormulaFromCML(filename);
			cout << "formula = " << formula << endl;
		}
		catch(std::exception& e)
		{
			cout << "Error: undefined" << endl;
		}
		catch(string& s)
		{
			cout << s << endl;
		}
	}
	else
		cout << "Error: No filename" << endl;
	return 0;
}