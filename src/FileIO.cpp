#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

void saveStates(vector<unsigned long int> &indexTable, string indexTableName)
{
    unsigned long int numStates_m = indexTable.size();
    ofstream wi(indexTableName, ios::out | ios::binary);
    if(!wi)
    {
        cout << "Cannot open file!" << endl;
    }
    wi.write((char *) &numStates_m, sizeof(numStates_m));
    for (int i = 0; i < numStates_m; i++)
        wi.write((char *) &indexTable[i], sizeof(indexTable[0]));
    wi.close();
    if (!wi.good()) 
    {
        cout << "Error occured at index writing time!" << endl;
    }
}


void readStates(vector<unsigned long int> &indexRead, string indexTableName)
{
    
    unsigned long int numStates_m;
    unsigned long int item;
    string item_str;
    ifstream rf(indexTableName, ios::out | ios::binary);
    if(!rf)
    {    cout << "Cannot open file!" << endl;
    }
    rf.read((char *) &numStates_m, sizeof(item));
    cout << "output from file "  << endl;
    for (int i = 0; i < numStates_m ; i++)
    {
        rf.read((char *) &item, sizeof(item));
        indexRead.push_back(item);
    }
    rf.close();
    if(!rf.good()) {
      cout << "Error occurred at reading time!" << endl;    
    }
//    cout << "numStates_m in save" << numStates_m << endl;
    
//    for (int i = 0; i < numStates_m; i++)
//    {
//        cout << indexRead[i] << endl;
        
//    }

    cout << "end reading ..." << endl;
}


