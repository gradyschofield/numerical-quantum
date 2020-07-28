//
// Created by grady on 7/25/20.
//

#ifndef ATOM_ENERGY_H
#define ATOM_ENERGY_H

class Energy {
public:
    double kineticEnergy;
    double ionElectronEnergy;
    double electronElectronEnergy;
    double exchangeCorrelationEnergy;
    double ionIonEnergy;
    double totalEnergy;

    Energy(double kineticEnergy, double ionElectronEnergy, double electronElectronEnergy, double exchangeCorrelationEnergy, double ionIonEnergy)
        : kineticEnergy(kineticEnergy), ionElectronEnergy(ionElectronEnergy), electronElectronEnergy(electronElectronEnergy),
            exchangeCorrelationEnergy(exchangeCorrelationEnergy), ionIonEnergy(ionIonEnergy)
    {
        totalEnergy = kineticEnergy + ionElectronEnergy + electronElectronEnergy + exchangeCorrelationEnergy + ionIonEnergy;
    }
};

#endif //ATOM_ENERGY_H
