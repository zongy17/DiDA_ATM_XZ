//
// Created by lenovo on 2020/12/3.
//

#include "State.h"

void State::calc_h_div_from_u_at_cell() {
    if (h_div_valid) {printf("calc_h_div_from_u: h_div already done.\n"); return;}
    if (!u_valid) {printf("calc_h_div_from_u: u unusable!\n"); exit(1);}
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++)
            h_div(k,i) = (u(k, i+1) - u(k, i)) / (x_half(i+1) - x_half(i));
    h_div_valid = true;
}

void State::Dh_damping(type_f dt) {
    if (!u_valid) {printf("damping: u unusable!\n"); exit(1);}
    //if (!last->h_div_valid) last->calc_h_div_from_u_at_cell();//先要计算水平散度，才能继续用于下一步的耗散
    if (!this->h_div_valid) this->calc_h_div_from_u_at_cell();

    //static type_f beta = 0.5;

    for (int k = 0; k < NLEV_full; k++){
        //middle part
        for (int i = 1; i < NX_half-1; i++){
            // 搞错了，注意水平散度耗散系数μ=C*Δx/Δt
            // u(k,i) = u(k,i) + Coeff(k) * dt * ((1.-beta)*(this->h_div(k,i)-this->h_div(k,i-1))\
                                              +    beta *(last->h_div(k,i)-last->h_div(k,i-1)))/(x_full(i)-x_full(i-1));
            type_f delta_x = x_full(i)-x_full(i-1);
            type_f miu = Coeff(k) * delta_x / dt;
            u(k,i) = u(k,i) + miu * dt * (this->h_div(k,i) - this->h_div(k,i-1)) / delta_x;
        }
#ifdef PERIODIC
        //i == 0 and NX_half-1: 周期性边界条件
        type_f delta_x = x_full(0) - x_half(0) + x_half(NX_half-1) - x_full(NX_full-1);
        // 搞错了，注意水平散度耗散系数μ=C*Δx/Δt
        // u(k,0) = u(k,0) + Coeff(k) * dt * ((1.-beta)*(this->h_div(k,0)-this->h_div(k,NX_full-1))\
                                              +    beta *(last->h_div(k,0)-last->h_div(k,NX_full-1)))/ delta_x;
        type_f miu = Coeff(k) * delta_x / dt;
        u(k,0) =  u(k,0) + miu * dt * (this->h_div(k,0) - this->h_div(k,NX_full-1)) / delta_x;
        u(k,NX_half-1) = u(k,0);
#endif
#ifdef UPSTREAM
        //上游来流的边界条件，不用对最左和最右的进行耗散，它们的值本来就是直接给定的
#endif
    }
    //damping过后应该还需要通过更新后的u更新插值的u_cell,u_ns,u_vtx等
    u_valid = true;
    interp_from_we_edges_to_cell(u, u_cell); u_cell_valid = true;
    interp_from_cell_to_ns_edges(u_cell, u_ns); u_ns_valid= true;
    interp_from_we_edges_to_vtxs(u, u_vtx);    u_vtx_valid= true;
    /*
    //以及更新底部free-slip条件的w
    for (int i = 0; i < NX_full; i++){
        int idx = idx_half(NLEV_half-1, i);
        w[idx] = u_ns[idx] * dzsdx_full[i];//bottom
    }
    */
}

type_f calc_Rayleigh_damping_coeff(type_f z, type_f zh, type_f z_top){
    /* z    : height at current position
     * zh   : height above which Rayleigh damping plays a role
     * z_top: height of the model top in this column
     */
    if (z > zh){
        type_f _sin = sin(pi/2.0 * (z-zh) / (z_top-zh));
        return _sin * _sin;
    }
    else return 0.0;
}

void State::Rayleigh_damping_u(type_f dt) {
    if (!this->u_valid) {printf("Rayleigh_damping: this->u unusable!\n"); exit(1);}
    if (!this->geo_potential_vtx_valid) {printf("Rayleigh_damping: this->geo_potential_vtx unusable!\n"); exit(1);}

    //对水平速度u做耗散
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_half; i++){
            type_f gz_this = 0.5 * (this->geo_potential_vtx(k+1,i) + this->geo_potential_vtx(k,i));
            type_f fz_this = calc_Rayleigh_damping_coeff(gz_this/gravity, this->zh, this->ztop_half(i));
            this->u(k,i) = this->u(k,i) - fz_this / tau_0 * (this->u(k,i)-this->u_ref) * dt ;
        }
    }
    u_valid = true;
    this->interp_from_we_edges_to_cell(this->u, this->u_cell); this->u_cell_valid = true;
    this->interp_from_cell_to_ns_edges(this->u_cell, this->u_ns); this->u_ns_valid= true;
    this->interp_from_we_edges_to_vtxs(this->u, this->u_vtx);    this->u_vtx_valid= true;
}

void State::Rayleigh_damping_w(type_f dt) {
    if (!this->w_valid) {printf("Rayleigh_damping: this->w unusable!\n"); exit(1);}
    if (!this->geo_potential_valid) {printf("Rayleigh_damping: this->geo_potential unusable!\n"); exit(1);}

    //对垂直速度w做耗散
    for (int k = 0; k < NLEV_half; k++)
        for (int i = 0; i < NX_full; i++){
            type_f fz_this = calc_Rayleigh_damping_coeff(this->geo_potential(k,i)/gravity, this->zh, this->geo_potential(0,i)/gravity);
            this->w(k,i) = this->w(k,i) - fz_this / tau_0 * (this->w(k,i)-this->w_ref) * dt;
        }
    this->w_valid = true;
    this->interp_from_ns_edges_to_cell(this->w, this->w_cell); this->w_cell_valid = true;
    this->interp_from_ns_edges_to_vtxs(this->w, this->w_vtx);  this->w_vtx_valid  = true;
}

void State::implicit_Rayleigh_damping_w(type_f dt){
    if (!this->w_valid) {printf("Rayleigh_damping: this->w unusable!\n"); exit(1);}
    if (!this->geo_potential_valid) {printf("Rayleigh_damping: this->geo_potential unusable!\n"); exit(1);}

    //对垂直速度w做耗散
    for (int k = 0; k < NLEV_half; k++)
        for (int i = 0; i < NX_full; i++){
            type_f fz_this = calc_Rayleigh_damping_coeff(this->geo_potential(k,i)/gravity, this->zh, this->geo_potential(0,i)/gravity);
            type_f miu = fz_this / tau_0;
            this->w(k,i) = this->w(k,i) / (1.0 + dt * miu);
        }
    this->w_valid = true;
    this->interp_from_ns_edges_to_cell(this->w, this->w_cell); this->w_cell_valid = true;
    this->interp_from_ns_edges_to_vtxs(this->w, this->w_vtx);  this->w_vtx_valid  = true;
}

void State::Damp_Pressure(State* last, type_f beta_d){
    if (!last->p_valid) {printf("Damp_Pressure: last->p unusable!\n"); exit(1);}
    if (!this->p_valid) {printf("Damp_Pressure: this->p unusable!\n"); exit(1);}
    this->p = this->p + beta_d * (this->p - last->p);// damping

    //if (!last->p_we_valid) {printf("Damp_Pressure: last->p_we unusable!\n"); exit(1);}
    //if (!this->p_we_valid) {printf("Damp_Pressure: this->p_we unusable!\n"); exit(1);}
    //this->p_we = this->p_we + beta_d * (this->p_we - last->p_we);

    if (!last->p_ns_valid) {printf("Damp_Pressure: last->p_ns unusable!\n"); exit(1);}
    if (!this->p_ns_valid) {printf("Damp_Pressure: this->p_ns unusable!\n"); exit(1);}
    this->p_ns = this->p_ns + beta_d * (this->p_ns - last->p_ns);

    //if (!last->p_vtx_valid) {printf("Damp_Pressure: last->p_vtx unusable!\n"); exit(1);}
    //if (!this->p_vtx_valid) {printf("Damp_Pressure: this->p_vtx unusable!\n"); exit(1);}
    //this->p_vtx = this->p_vtx + beta_d * (this->p_vtx - last->p_vtx);
}