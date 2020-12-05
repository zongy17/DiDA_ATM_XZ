//
// Created by lenovo on 2020/12/3.
//

#include "State.h"

void State::Print_coordinates() {
    x_full.raw_print("x_full");
    x_half.raw_print("x_half");
    A_eta_half.raw_print("A_half");
    B_eta_half.raw_print("B_half");
    eta_half.raw_print("eta_half");
    A_eta_full.raw_print("A_full");
    B_eta_full.raw_print("B_full");
    eta_full.raw_print("eta_full");
}

void State::Print_Console_All_Prognostic(const char* extra_txt) {
    printf(extra_txt);
    phs.print("phs");
    pt.print("pt");
    u.print("u");
    w.print("w");
    geo_potential.print("gz");
}

void State::Write_dat(const Col<type_f> &col_val, const char* filename) {
    col_val.save(filename, raw_ascii);
}
void State::Write_dat(const Row<type_f> &row_val, const char *filename) {
    row_val.save(filename, raw_ascii);
}
void State::Write_dat(const Mat<type_f> &mat_val, const char *filename) {
    mat_val.save(filename, raw_ascii);
}

void State::Save_all(){
    phs.save("phs_last.txt", raw_ascii);
    u.save("u_last.txt", raw_ascii);
    w.save("w_last.txt", raw_ascii);
    pt.save("pt_last.txt", raw_ascii);
    geo_potential.save("geo_last.txt", raw_ascii);
}
