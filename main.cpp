#include "State.h"
#include "integrate.h"

void Parse_case_info(const char* case_control_file, int & method, type_f & x_lo, type_f & x_hi,\
                    int & NLEV_full, int & NX_full, int & nstep, \
                    int & nstep_to_write, type_f & tot_delta_t, type_f & u_ref, type_f w_ref, type_f & zh,\
                    type_f & tau_0, type_f & kexi){
    FILE * fp = fopen(case_control_file, "r");
    char option_name[100], format[10];
    if (sizeof(type_f) == sizeof(double)) sprintf(format, "%s", "%lf");
    else sprintf(format, "%f");
    while (fscanf(fp, "%s", option_name) != EOF){
        if (strcmp(option_name, "method:") == 0)                        fscanf(fp, "%d", &method);
        else if (strcmp(option_name, "x_lo:") == 0)                     fscanf(fp, format, &x_lo);
        else if (strcmp(option_name, "x_hi:") == 0)                     fscanf(fp, format, &x_hi);
        else if (strcmp(option_name, "NX_full:") == 0)                  fscanf(fp, "%d", &NX_full);
        else if (strcmp(option_name, "NLEV_full:") == 0)                fscanf(fp, "%d", &NLEV_full);
        else if (strcmp(option_name, "nstep:") == 0)                    fscanf(fp, "%d", &nstep);
        else if (strcmp(option_name, "nstep_to_write:") == 0)           fscanf(fp, "%d", &nstep_to_write);
        else if (strcmp(option_name, "delta_t:") == 0)                  fscanf(fp, format, &tot_delta_t);
        else if (strcmp(option_name, "u_ref:") == 0)                    fscanf(fp, format, &u_ref);
        else if (strcmp(option_name, "w_ref:") == 0)                    fscanf(fp, format, &w_ref);
        else if (strcmp(option_name, "Rayleigh_layer_height:") == 0)    fscanf(fp, format, &zh);
        else if (strcmp(option_name, "tau_0:") == 0)                    fscanf(fp, format, &tau_0);
        else if (strcmp(option_name, "kexi:") == 0)                     fscanf(fp, format, &kexi);
        else {printf("an unrecognized param!\n"); exit(1);}
    }
    fclose(fp);
}


bool Check_step_size(State* const state, type_f delta_t){
    if (!state->geo_potential_vtx_valid) printf("Check_step_size: geo_potential_vtx unusable!\n");
    int NX_half = state->NX_half, NLEV_full = state->NLEV_full;
    type_f min_delta_z = Datum<type_f>::inf;
    int min_i = -1, min_k = -1;
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_half; i++){
            type_f local_delta_z = (state->geo_potential_vtx(k,i) - state->geo_potential_vtx(k+1,i)) / gravity;
            if (local_delta_z < min_delta_z){ min_delta_z = local_delta_z; min_i = i; min_k = k; }
        }
    }
    printf("k:%d, i:%d, min_delta_z: %lf, delta_t: %lf\n", min_k, min_i, min_delta_z, delta_t);
    if (min_delta_z < delta_t * sound_velocity ) return false;
    else return true;
}


void write_err_to_log(State* const curr, State* const next, FILE* fp){
    Mat<type_f> w_err = abs(next->w - curr->w);
    Mat<type_f> u_err = abs(next->u - curr->u);
    int max_u_err_idx = u_err.index_max(), max_w_err_idx = w_err.index_max();
    //fprintf(fp, "max_u_err: %lf, max_w_err: %lf, max_dlayer_wdt_err: %lf\n", max(max(u_err)), max(max(w_err)),\
                                                                            max(max(dlayer_wdt_err)));
    fprintf(fp, "max_u_err: %lf,  max_w_err: %lf\n", u_err.max(), w_err.max());
}

int main()
{
    type_f x_lo, x_hi, tot_delta_t, u_ref, w_ref, zh, tau_0, kexi;
    int vert_init_method, NX_full, NLEV_full, nstep, nstep_to_write;

    Parse_case_info("case-control.txt", vert_init_method, x_lo, x_hi, NLEV_full, NX_full, nstep, \
        nstep_to_write, tot_delta_t,u_ref, w_ref, zh, tau_0, kexi);

    State* curr_state = new State(NLEV_full, NX_full);
    State* in_1_state = new State(NLEV_full, NX_full);
#ifdef mod_RK4
    State* in_2_state = new State(NLEV_full, NX_full);
    State* in_3_state = new State(NLEV_full, NX_full);
    printf("Time Method: mod RK-4\n");
#endif
#ifdef mod_RK3
    State* in_2_state = new State(NLEV_full, NX_full);
    printf("Time Method: mod RK-3\n");
#endif
    State* next_state = new State(NLEV_full, NX_full);

    type_f curr_t = 0.0;

    if (vert_init_method == 1){
        curr_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt", "A_full.txt", "B_full.txt");
        in_1_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt", "A_full.txt", "B_full.txt");
#ifdef mod_RK4
        in_2_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "topo.txt", "A_half.txt", "B_half.txt", "A_full.txt", "B_full.txt");
        in_3_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "topo.txt", "A_half.txt", "B_half.txt", "A_full.txt", "B_full.txt");
#endif
#ifdef mod_RK3
        in_2_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt", "A_full.txt", "B_full.txt");
#endif
        next_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt", "A_full.txt", "B_full.txt");
    } else {
        curr_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt");
        in_1_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt");
#ifdef mod_RK4
        in_2_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "topo.txt", "A_half.txt", "B_half.txt");
        in_3_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "topo.txt", "A_half.txt", "B_half.txt");
#endif
#ifdef mod_RK3
        in_2_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt");
#endif
        next_state->Pre_Process(x_lo, x_hi, u_ref, w_ref, zh, tau_0, kexi, "zs_half.txt", "zs_full.txt", "A_half.txt", "B_half.txt");
    }

    printf("All States finished Pre_Process.\n");

    curr_state->initiate("u0.txt","w0.txt","pt0.txt","geo0.txt","phs0.txt");
    curr_state->diagnose();
    curr_state->tend_all_explicit();


    printf("Initial State set.\n");

    //记录residual
    FILE* err_log = fopen("err_log.txt", "w+");
    //LOOP BEGIN HERE
    char filename[100];
    sprintf(filename, "eta_full.dat\0");
    State::Write_dat(curr_state->eta_full, filename);
    sprintf(filename, "eta_half.dat\0");
    State::Write_dat(curr_state->eta_half, filename);
    sprintf(filename, "x_full.dat\0");
    State::Write_dat(curr_state->x_full, filename);
    sprintf(filename, "x_half.dat\0");
    State::Write_dat(curr_state->x_half, filename);

    for (int step = 0; step < nstep; step++)
    {
        printf("\nstep: %d, curr_t: %.2f\n", step, curr_t);
        Check_step_size(curr_state, tot_delta_t);

        //写出和输出
#ifdef OUTPUT
        if (step % nstep_to_write == 0 && step != 0) {
            printf("writing to .nc file...\n");
            sprintf(filename, "ph_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->ph_cell, filename);
            sprintf(filename, "u_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->u, filename);
            sprintf(filename, "u_ns_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->u_ns, filename);
            sprintf(filename, "w_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->w, filename);
            sprintf(filename, "pt_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->pt, filename);
            sprintf(filename, "gz_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->geo_potential, filename);
            sprintf(filename, "K_hori_adv%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->K_hori_adv, filename);
            sprintf(filename, "u_vert_adv%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->u_vert_adv, filename);
            sprintf(filename, "hpgf_1_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->hpgf_1, filename);
            sprintf(filename, "hpgf_2_%.2fs.dat\0", curr_t);
            State::Write_dat(curr_state->hpgf_2, filename);
            printf("finish writing.\n");
            //输出到控制台观察
            //curr_state->u.print("curr->u");
            //curr_state->w.print("curr->w");
        }
#endif

        //假定curr已经完成了diagnose和tend_all_explicit
#ifdef mod_RK4
        integrate_per_substep(curr_state, curr_state, in_1_state, tot_delta_t / 4);
        integrate_per_substep(curr_state, in_1_state, in_2_state, tot_delta_t / 3);
        integrate_per_substep(curr_state, in_2_state, in_3_state, tot_delta_t / 2);
        integrate_per_substep(curr_state, in_3_state, next_state, tot_delta_t);
#endif
#ifdef mod_RK3
        integrate_per_substep(curr_state, curr_state, in_1_state, tot_delta_t / 3.0);
        integrate_per_substep(curr_state, in_1_state, in_2_state, tot_delta_t / 2.0);
        integrate_per_substep(curr_state, in_2_state, next_state, tot_delta_t);
#endif


#ifdef DEBUG
        next_state->Print_Console_All_Prognostic("\nnext_state\n");
#endif
        //damping
        //next_state->Dh_damping(tot_delta_t);//用n+1时刻的水平速度场对自己做水平散度耗散
        //next_state->Rayleigh_damping_u(tot_delta_t);//用n+1时刻的速度场对自己做(隐式)Rayleigh耗散

        //write error to log
        write_err_to_log(curr_state, next_state, err_log);

        //switch State pointers
#ifdef mod_RK4
        reset_state_4(curr_state, next_state, in_1_state, in_2_state, in_3_state);
#endif
#ifdef mod_RK3
        reset_state_3(curr_state, next_state, in_1_state, in_2_state);
#endif
        //time advance
        curr_t += tot_delta_t;
    }

    fclose(err_log);
    printf("Computing finished.\n");
    curr_state->Save_all();
#ifdef OUTPUT
    printf("writing to .nc file...\n");
    sprintf(filename, "ph_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->ph_cell, filename);
    sprintf(filename, "u_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->u, filename);
    sprintf(filename, "u_ns_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->u_ns, filename);
    sprintf(filename, "w_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->w, filename);
    sprintf(filename, "pt_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->pt, filename);
    sprintf(filename, "gz_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->geo_potential, filename);
    sprintf(filename, "K_hori_adv%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->K_hori_adv, filename);
    sprintf(filename, "u_vert_adv%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->u_vert_adv, filename);
    sprintf(filename, "hpgf_1_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->hpgf_1, filename);
    sprintf(filename, "hpgf_2_%.2fs.dat\0", curr_t);
    State::Write_dat(curr_state->hpgf_2, filename);
    printf("finish writing.\n");
    FILE * fp = fopen("netcdf-control.txt", "w+");
    fprintf(fp, "NX_full: %d\nNLEV_full: %d\n", NX_full, NLEV_full);
    fprintf(fp, "n_time: %d\nt_per_output: %f\n", nstep/nstep_to_write + 1, tot_delta_t*nstep_to_write);
    fclose(fp);
    //输出到控制台观察
    curr_state->u.print("curr->u");
    curr_state->w.print("curr->w");
#endif
#ifdef modmod_RK4
    delete curr_state; delete in_1_state; delete in_2_state; delete in_3_state; delete next_state;
#endif
#ifdef mod_RK3
    delete curr_state; delete in_1_state; delete in_2_state; delete next_state;
#endif
    return 0;
}