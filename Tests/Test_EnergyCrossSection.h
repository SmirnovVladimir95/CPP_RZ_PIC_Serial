#ifndef CPP_RZ_PIC_TEST_ENERGYCROSSSECTION_H
#define CPP_RZ_PIC_TEST_ENERGYCROSSSECTION_H

#include "../ElementaryProcesses/EnergyCrossSection.h"
#include <iostream>
#include <cassert>
#define EV 1.6021766208e-19


void test_EnergyCrossSection() {
    cout << "test_EnergyCrossSection:" << endl;

    //Constructor from data in code
    vector<array<scalar, 2>> energy_cross_section = {{0, 2}, {1, 3}};
    EnergyCrossSection cross_section(energy_cross_section);
    assert(not cross_section.empty());
    cout << cross_section.get_cross_section(0.5*EV) << endl;

    //Constructor default
    EnergyCrossSection* new_cross_section;
    //assert(new_cross_section->empty());
    new_cross_section = &cross_section;
    assert(not new_cross_section->empty());
    cout << new_cross_section->get_cross_section(0.75*EV) << endl;
    EnergyCrossSection another_new_cross_section;
    assert(another_new_cross_section.empty());

    //Constructor from file
    EnergyCrossSection new_new_cross_section("../ElementaryProcesses/CrossSectionData/e-He_ionization.txt");
    assert(not new_new_cross_section.empty());
    cout << new_new_cross_section.get_cross_section(30.5*EV) << endl;

    //Constructor from file for He+ - He elastic
    EnergyCrossSection helium_elastic("../ElementaryProcesses/CrossSectionData/He+-He_elastic.txt");
    assert(not helium_elastic.empty());
    cout << helium_elastic.get_cross_section(5.25*EV) << endl;

    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_ENERGYCROSSSECTION_H
