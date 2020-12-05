//
// Created by lenovo on 2020/12/3.
//

#ifndef XZ_V11_INTEGRATE_H
#define XZ_V11_INTEGRATE_H

#include "State.h"
void reset_state_4(State* & curr, State* & next, State* inter_1, State* inter_2, State* inter_3);

void reset_state_3(State* & curr, State* & next, State* inter_1, State* inter_2);

void integrate_per_substep(State* curr, State* last_sub, State* next_sub, type_f dt);
void advance_phs(State* curr, State* last_sub, State* next_sub, type_f dt);
void advance_pt(State* curr, State* last_sub, State* next_sub, type_f dt);
void advance_u(State* curr, State* last_sub, State* next_sub, type_f dt);
//void advance_w(State* curr, State* ET, State* next, type_f dt);
//void advance_geo_potential(State* curr, State* ET, State* next, type_f dt);
void advance_w_gz_implicit(State* curr, State* last_sub, State* next_sub, type_f dt);

#endif //XZ_V11_INTEGRATE_H
