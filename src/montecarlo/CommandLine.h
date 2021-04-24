//
// Created by grady on 4/24/21.
//

#ifndef ATOM_COMMANDLINE_H
#define ATOM_COMMANDLINE_H

#include<iostream>
#include<sstream>

using namespace std;

class CommandLine {

public:
    static long parseNumPoints(int argc, char **argv, long defaultNumPoints) {
        if(argc != 2) {
            cout << "Using " << defaultNumPoints << " evaluation points since no number specified on the command line." << endl;
            return defaultNumPoints;
        } else {
            stringstream sstr(argv[1]);
            double realN;
            sstr >> realN;
            if (sstr.fail() || realN < 1) {
                cout << "Using " << defaultNumPoints << " evaluation points since I could not parse \"" << argv[1] << "\" as an integer." << endl;
                return defaultNumPoints;
            } else {
                return realN;
            }
        }
    }
};

#endif //ATOM_COMMANDLINE_H
