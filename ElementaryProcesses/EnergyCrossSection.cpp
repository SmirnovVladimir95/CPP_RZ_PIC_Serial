#include "EnergyCrossSection.h"
#include <fstream>
#include "../Tools/Helpers.h"
#include <iostream>
#define EV 1.6021766208e-19


scalar EnergyCrossSection::get_cross_section(scalar energy) {
    energy /= EV; // convert Joule to ElectronVolt
    int idx = binary_search(energy);
    if (idx == energy_cross_section.size() - 1) {
        cout << "Warning: cross-section interpolation is out of bounds (cross-section is set to max element in energy_cross-section)!!!" << endl;
        return energy_cross_section[idx][1];
    }
    if (idx == -1) {
        //cout << "Warning: cross-section interpolation is out of bounds (cross-section is set to min element in energy_cross-section)!!!" << endl;
        return energy_cross_section[0][1];
    }
    return linear_interp(idx, idx + 1, energy);
}

EnergyCrossSection::EnergyCrossSection(const vector<array<scalar, 2>>& energy_cross_section) : energy_cross_section(energy_cross_section) {}

EnergyCrossSection::EnergyCrossSection(const string& filename) {
    ifstream input(filename);
    scalar energy, cross_section;
    if (input) {
        while (input) {
            input >> energy;
            input >> cross_section;
            energy_cross_section.push_back({energy, cross_section});
        }
    } else {
        cout << "Warning: cross-section data file not found!!!" << endl;
    }
}

int EnergyCrossSection::binary_search(scalar energy) {
    int idx_1 = 0;
    int idx_2 = energy_cross_section.size() - 1;
    if (energy >= energy_cross_section[idx_2][0]) {
        return idx_2;
    }
    if (energy <= energy_cross_section[idx_1][0]) {
        return -1;
    }
    int idx;
    while (idx_2 - idx_1 > 1) {
        idx = (idx_1 + idx_2) / 2;
        if (energy > energy_cross_section[idx][0]) {
            idx_1 = idx;
        }
        else {
            idx_2 = idx;
        }
    }
    return idx_1;
}

scalar EnergyCrossSection::linear_interp(int ptr_1, int ptr_2, scalar energy) {
    scalar delta_cross = energy_cross_section[ptr_2][1] - energy_cross_section[ptr_1][1];
    scalar delta_energy = energy_cross_section[ptr_2][0] - energy_cross_section[ptr_1][0];
    return energy_cross_section[ptr_1][1] + (energy - energy_cross_section[ptr_1][0]) * delta_cross / delta_energy;
}

bool EnergyCrossSection::empty() {
    return energy_cross_section.empty();
}

array<scalar, 2> EnergyCrossSection::get_energy_cross_section(int idx) {
    array<scalar, 2> energy_cross_section_item = energy_cross_section[idx];
    energy_cross_section_item[0] *= EV; // convert from ElectronVolt to Joule
    return energy_cross_section_item;
}

EnergyCrossSection::EnergyCrossSection() = default;
