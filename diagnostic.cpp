//
// Created by lenovo on 2020/12/3.
//

#include "State.h"

//diagnose诊断计算当前时刻state的非预报变量
void State::diagnose() {
    tend_dphsdt();//诊断计算当前时刻的地面静力气压倾向
    calc_m_detadt_interp_to_cell_vtx();//诊断计算垂直坐标速度和m的乘积: mη'，并插值
    calc_K_from_u_at_cell();//诊断计算水平动能
    calc_rho_from_dphidph_at_cell_interp_to_we();//诊断计算密度
    calc_p_from_rho_at_cell_interp_to_we_ns_vtx();//诊断计算气压，并插值
}

void State::calc_ph_at_cell_ns_we_vtx() {
    if (ph_ns_valid&&ph_cell_valid&&ph_we_valid&&ph_vtx_valid) {printf("calc_ph_at_cell_ns_we_vtx: already done.\n"); return;}
    if (!phs_valid) {printf("calc_ph_layer_ph_interp_to_ns_we_vtx: phs unusable!\n"); exit(1);}
    // calc ph on half levels(i.e., north/south interface) from definitions
    for (int k = 0; k < NLEV_half; k++){//from top to bottom
        for (int i = 0; i < NX_full; i++)//from left to right
            ph_ns(k,i) = A_eta_half(k) * ph0 + B_eta_half(k) * phs(i);
    }
    ph_ns_valid = true;
    // calc ph on full levels(i.e., cell centroid) from definitions
    for (int k = 0; k < NLEV_full; k++) {//from top to bottom
        for (int i = 0; i < NX_full; i++)//from left to right
            ph_cell(k,i) = A_eta_full(k) * ph0 + B_eta_full(k) * phs(i);
    }
    ph_cell_valid = true;
    //interp to west/east interface from centroid
    interp_from_cell_to_we_edges(ph_cell, ph_we); ph_we_valid = true;
    //interp to vtx from north/south interface
    interp_from_ns_edges_to_vtxs(ph_ns, ph_vtx); ph_vtx_valid = true;
}

void State::calc_layer_ph_at_cell_ns_we_vtx() {
    if (layer_ph_cell_valid&&layer_ph_ns_valid&&layer_ph_we_valid&&layer_ph_vtx_valid) {printf("calc_layer_ph_at_cell_ns_we: already done.\n"); return;}
    //calc layer_ph at cell centroid
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++)
            layer_ph_cell(k,i) = ph_ns(k+1,i) - ph_ns(k,i);
    layer_ph_cell_valid = true;
    //calc layer_ph at n/s interface
    for (int k = 1; k < NLEV_half-1; k++)//middle part
        for (int i = 0; i < NX_full; i++)
            layer_ph_ns(k,i) = ph_cell(k,i) - ph_cell(k-1,i);
    //NOTE: instead of extrapolation, use half layer (layer_ph_ns at bottom and top refer to only half layer)
    for (int i = 0; i < NX_full; i++){
        layer_ph_ns(NLEV_half-1,i) = ph_ns(NLEV_half-1,i) - ph_cell(NLEV_full-1,i);//bottom
        layer_ph_ns(0,i) = ph_cell(0,i) - ph_ns(0,i);//top
    }
    layer_ph_ns_valid = true;
    //calc layer_ph at w/e interface
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_half; i++)
            layer_ph_we(k,i) = ph_vtx(k+1,i) - ph_vtx(k,i);
    layer_ph_we_valid = true;
    //calc layer_ph at vtx
    for (int k = 1; k < NLEV_half-1; k++)        //middle part
        for (int i = 0; i < NX_half; i++)
            layer_ph_vtx(k,i) = ph_we(k,i) - ph_we(k-1,i);
    //NOTE: instead of extrapolation, use half layer (layer_ph_ns at bottom and top refer to only half layer)
    for (int i = 0; i < NX_half; i++){
        layer_ph_vtx(NLEV_half-1,i) = ph_vtx(NLEV_half-1,i) - ph_we(NLEV_full-1,i);//bottom
        layer_ph_vtx(0,i) = ph_we(0,i) - ph_vtx(0,i);//top
    }
    layer_ph_vtx_valid = true;
}

void State::tend_dphsdt() {
    if (dphsdt_valid) {printf("calc_dphsdt: dphsdt already done.\n"); return;}
    if (!layer_ph_we_valid) {printf("calc_dphsdt: layer_ph_we unusable!\n"); exit(1);}
    if (!u_valid) {printf("calc_dphsdt: u unusable!\n"); exit(1);}
    //version-2
    for (int i = 0; i < NX_full; i++){
        type_f sum = 0.0;
        type_f delta_x = x_half(i+1) - x_half(i);
        for (int k = 0; k < NLEV_full; k++)
            sum += u(k,i+1) * layer_ph_we(k,i+1) - u(k,i) * layer_ph_we(k,i);
        dphsdt(i) = - sum / delta_x;
    }
    dphsdt_valid = true;
}

void State::calc_m_detadt_interp_to_cell_vtx() {
    if (m_detadt_valid&&m_detadt_cell_valid) {printf("calc_m_detadt_interp_to_cell_vtx: m_detadt already done.\n"); return;}
    if (!layer_ph_we_valid) {printf("calc_m_detadt_interp_to_cell_vtx: layer_ph_we unusable!\n"); exit(1);}
    if (!u_valid) {printf("calc_m_detadt_interp_to_cell_vtx: u unusable!\n"); exit(1);}
    if (!dphsdt_valid) {printf("calc_m_detadt_interp_to_cell_vtx: dphsdt unusable!\n"); exit(1);}

    for (int i = 0; i < NX_full; i++){
        m_detadt(0,i) = 0.0;//top
        type_f delta_x = x_half(i+1) - x_half(i);
        type_f part_sum = 0.0;//从顶部开始累计求和
        for (int k = 1; k < NLEV_half; k++){
            //因为k指的是half level的值，紧挨着上面的full level的指标应该是k-1
            part_sum += u(k-1,i+1) * layer_ph_we(k-1,i+1) - u(k-1,i) * layer_ph_we(k-1,i);
            m_detadt(k,i) = - B_eta_half(k) * dphsdt(i) - part_sum / delta_x;
        }
    }
    m_detadt_valid = true;
    interp_from_ns_edges_to_cell(m_detadt, m_detadt_cell); m_detadt_cell_valid = true;
    //interp_from_cell_to_we_edges(m_detadt_cell, m_detadt_we); m_detadt_we_valid= true;
    interp_from_ns_edges_to_vtxs(m_detadt, m_detadt_vtx);  m_detadt_vtx_valid  = true;
}

void State::calc_K_from_u_at_cell() {
    if (K_valid) {printf("calc_K_from_u: K already done.\n"); return;}
    if (!u_valid) {printf("calc_K_from_u: u unusable!\n"); exit(1);}
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++)
            K(k,i) = ( 0.5*u(k,i+1)*u(k,i+1) + 0.5*u(k,i)*u(k,i) ) / 2.0;
    K_valid = true;
}

void State::calc_rho_from_dphidph_at_cell_interp_to_we() {
    if (rho_valid) {printf("calc_rho_from_dphidph_at_cell: already done.\n"); return;}
    if (!layer_ph_cell_valid) {printf("calc_rho_from_dphidph_at_cell: layer_ph_cell unusable!\n"); exit(1);}
    if (!geo_potential_valid) {printf("calc_rho_from_dphidph_at_cell: geo_potential unusable!\n"); exit(1);}
    //calc rho at cell centroid
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++)
            rho(k,i) = - layer_ph_cell(k,i) / (geo_potential(k+1,i) - geo_potential(k,i));
    rho_valid = true;
    /*
    if (!layer_ph_we_valid) {printf("calc_rho_from_dphidph_at_cell: layer_ph_we unusable!\n"); exit(1);}
    if (!geo_potential_vtx_valid) {printf("calc_rho_from_dphidph_at_cell: geo_potential_vtx unusable!\n"); exit(1);}
    //calc rho at w/e interface
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_half; i++)
            rho_we(k,i) = - layer_ph_we(k,i) / (geo_potential_vtx(k+1,i) - geo_potential_vtx(k,i));
    rho_we_valid = true;
    */
    interp_from_cell_to_we_edges(rho, rho_we);
    rho_we_valid = true;
}

void State::calc_p_from_rho_at_cell_interp_to_we_ns_vtx() {
    if (p_valid&&p_we_valid&&p_vtx_valid&&p_ns_valid) {printf("calc_p_from_rho_interp_to_ns_vtx: already done.\n"); return;}

    //calc p at cell centroid
    if (!pt_valid) {printf("calc_p_from_rho_interp_to_ns_vtx: pt unusable!\n"); exit(1);}
    if (!rho_valid){printf("calc_p_from_rho_interp_to_ns_vtx: rho unusable!\n"); exit(1);}
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++){
            type_f base = Rd * pt(k,i) * rho(k,i) / p0;
            p(k,i) = p0 * pow(base, gamma);
        }
    p_valid = true;

    interp_from_cell_to_we_edges(p, p_we);  p_we_valid = true;
/*
    //calc p at n/s interface
    if (!pt_ns_valid){printf("calc_p_from_rho_interp_to_ns_vtx: pt_ns unusable!\n"); exit(1);}
    if (!geo_potential_cell_valid){printf("calc_p_from_rho_interp_to_ns_vtx: geo_potential_cell unusable!\n"); exit(1);}
    for (int k = 1; k < NLEV_half-1; k++){
        for (int i = 0; i < NX_full; i++){
            type_f rho_ns = - layer_ph_ns(k,i) / (geo_potential_cell(k,i) - geo_potential_cell(k-1,i));
            type_f base   = Rd * pt_ns(k,i) * rho_ns / p0;
            p_ns(k,i) = p0 * pow(base, gamma);
        }
    }
    int k;
    //NOTE: instead of extrapolation, use part-diff here instead of center-diff
    if (!geo_potential_valid) {printf("calc_p_from_rho_interp_to_ns_vtx: geo_potential unusable!\n"); exit(1);}
    if (!layer_ph_ns_valid) {printf("calc_p_from_rho_interp_to_ns_vtx: layer_ph_ns unusable!\n"); exit(1);}
    k = 0;//top
    for (int i = 0; i < NX_full; i++){
        type_f rho_ns = - layer_ph_ns(k,i) / (geo_potential_cell(k,i) - geo_potential(k,i));
        type_f base   = Rd * pt_ns(k,i) * rho_ns / p0;
        p_ns(k,i) = p0 * pow(base, gamma);
    }
    k = NLEV_half-1;//bottom
    for (int i = 0; i < NX_full; i++){
        type_f rho_ns = - layer_ph_ns(k,i) / (geo_potential(k,i) - geo_potential_cell(k-1,i));
        type_f base   = Rd * pt_ns(k,i) * rho_ns / p0;
        p_ns(k,i) = p0 * pow(base, gamma);
    }
    p_ns_valid = true;

    //calc p_vtx, same as p_ns
    interp_from_ns_edges_to_vtxs(p_ns, p_vtx); p_vtx_valid = true;
*/
    interp_from_we_edges_to_vtxs(p_we, p_vtx);  p_vtx_valid= true;
    interp_from_cell_to_ns_edges(p, p_ns);      p_ns_valid = true;
}

void State::calc_p_at_cell_linear_interp_to_ns_we_vtx(State* last){
    //calc p at cell
    if (!last->p_valid) {printf("calc_p_linear: last->p unusable!\n"); exit(1);}
    if (!this->pt_valid) {printf("calc_p_linear: this->pt unusable!\n"); exit(1);}
    if (!this->layer_ph_cell_valid) {printf("calc_p_linear: this->layer_ph_cell unusable!\n"); exit(1);}
    if (!last->pt_valid) {printf("calc_p_linear: last->pt unusable!\n"); exit(1);}
    if (!last->layer_ph_cell_valid) {printf("calc_p_linear: last->layer_ph_cell unusable!\n"); exit(1);}
    if (!this->geo_potential_valid) {printf("calc_p_linear: this->geo_potential unusable!\n"); exit(1);}
    if (!last->geo_potential_valid) {printf("calc_p_linear: last->geo_potential unusable!\n"); exit(1);}
    this->p = last->p % (1.0 + gamma * (this->layer_ph_cell % this->pt) / (last->layer_ph_cell % last->pt));
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_full; i++){
            this->p(k,i) += - gamma * last->p(k,i) * (this->geo_potential(k+1,i) - this->geo_potential(k,i))\
                                                    /(last->geo_potential(k+1,i) - last->geo_potential(k,i));
        }
    }
    p_valid = true;
    interp_from_cell_to_ns_edges(p, p_ns);  p_ns_valid = true;
    interp_from_cell_to_we_edges(p, p_we);  p_we_valid = true;
    interp_from_we_edges_to_vtxs(p_we, p_vtx);  p_vtx_valid= true;
    /*
    //calc p_we
    if (!last->p_we_valid) {printf("calc_p_linear: last->p_we unusable!\n"); exit(1);}
    if (!this->pt_we_valid) {printf("calc_p_linear: this->pt_we unusable!\n"); exit(1);}
    if (!this->layer_ph_we_valid) {printf("calc_p_linear: this->layer_ph_we unusable!\n"); exit(1);}
    if (!last->pt_we_valid) {printf("calc_p_linear: last->pt_we unusable!\n"); exit(1);}
    if (!last->layer_ph_we_valid) {printf("calc_p_linear: last->layer_ph_we unusable!\n"); exit(1);}
    if (!this->geo_potential_vtx_valid) {printf("calc_p_linear: this->geo_potential_vtx unusable!\n"); exit(1);}
    if (!last->geo_potential_vtx_valid) {printf("calc_p_linear: last->geo_potential_vtx unusable!\n"); exit(1);}
    this->p_we = last->p_we % (1.0 + gamma * (this->layer_ph_we % this->pt_we) / (last->layer_ph_we % last->pt_we));
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_half; i++){
            this->p_we(k,i) += - gamma * last->p_we(k,i) * (this->geo_potential_vtx(k+1,i) - this->geo_potential_vtx(k,i))\
                                                          /(last->geo_potential_vtx(k+1,i) - last->geo_potential_vtx(k,i));
        }
    }
    //if (!last_sub->p_we_valid) {printf("calc_p_linear: last_sub->p_we unusable!\n"); exit(1);}
    //this->p_we = this->p_we + beta_d * (this->p_we - last_sub->p_we);// damping
    this->p_we_valid = true;
    */
}