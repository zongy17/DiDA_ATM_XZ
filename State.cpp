//
// Created by lenovo on 2020/12/3.
//

#include "State.h"


#include "State.h"

State::State(int nlev_full, int nx_full) {
    NLEV_full = nlev_full;      NX_full = nx_full;
    NLEV_half = nlev_full+1;    NX_half = nx_full+1;
    x_full = Row<type_f>(NX_full, fill::zeros);
    x_half = Row<type_f>(NX_half, fill::zeros);
    zs_full= Row<type_f>(NX_full, fill::zeros);
    zs_half= Row<type_f>(NX_half, fill::zeros);
    dzsdx_full = Row<type_f>(NX_full, fill::zeros);
    ztop_half  = Row<type_f>(NX_half, fill::zeros);
    dztopdx_full = Row<type_f>(NX_full, fill::zeros);
    eta_full = Col<type_f>(NLEV_full, fill::zeros);
    eta_half = Col<type_f>(NLEV_half, fill::zeros);
    A_eta_full = Col<type_f>(NLEV_full, fill::zeros);
    A_eta_half = Col<type_f>(NLEV_half, fill::zeros);
    B_eta_full = Col<type_f>(NLEV_full, fill::zeros);
    B_eta_half = Col<type_f>(NLEV_half, fill::zeros);

    //horizontal velocity
    u                   = Mat<type_f>(NLEV_full, NX_half);
    u_cell              = Mat<type_f>(NLEV_full, NX_full);
    u_vtx               = Mat<type_f>(NLEV_half, NX_half);
    u_ns                = Mat<type_f>(NLEV_half, NX_full);
    //vertical velocity
    w                   = Mat<type_f>(NLEV_half, NX_full);
    w_cell              = Mat<type_f>(NLEV_full, NX_full);
    w_vtx               = Mat<type_f>(NLEV_half, NX_half);
    //kinetic energy
    K                   = Mat<type_f>(NLEV_full, NX_full);
    //potential temperature
    pt                  = Mat<type_f>(NLEV_full, NX_full);
    pt_ns               = Mat<type_f>(NLEV_half, NX_full);
    pt_we               = Mat<type_f>(NLEV_full, NX_half);
    //pt_vtx              = Mat<type_f>(NLEV_half, NX_half);
    //hydrostatic pressure
    ph_cell             = Mat<type_f>(NLEV_full, NX_full);
    ph_ns               = Mat<type_f>(NLEV_half, NX_full);
    //ph_we               = Mat<type_f>(NLEV_full, NX_half);
    //ph_vtx              = Mat<type_f>(NLEV_half, NX_half);
    layer_ph_vtx        = Mat<type_f>(NLEV_half, NX_half);
    layer_ph_cell       = Mat<type_f>(NLEV_full, NX_full);
    layer_ph_we         = Mat<type_f>(NLEV_full, NX_half);
    layer_ph_ns         = Mat<type_f>(NLEV_half, NX_full);
    layer_ph_vtx        = Mat<type_f>(NLEV_half, NX_half);
    phs                 = Row<type_f>(NX_full);
    //pressure
    p                   = Mat<type_f>(NLEV_full, NX_full);
    p_we                = Mat<type_f>(NLEV_full, NX_half);
    p_ns                = Mat<type_f>(NLEV_half, NX_full);
    p_vtx               = Mat<type_f>(NLEV_half, NX_half);
    //density
    rho                 = Mat<type_f>(NLEV_full, NX_full);
    rho_we              = Mat<type_f>(NLEV_full, NX_half);
    //vertical coordinate velocity
    m_detadt            = Mat<type_f>(NLEV_half, NX_full);
    m_detadt_vtx        = Mat<type_f>(NLEV_half, NX_half);
    m_detadt_cell       = Mat<type_f>(NLEV_full, NX_full);
    //m_detadt_we         = Mat<type_f>(NLEV_full, NX_half);
    //geo potential
    geo_potential       = Mat<type_f>(NLEV_half, NX_full);
    geo_potential_cell  = Mat<type_f>(NLEV_full, NX_full);
    geo_potential_vtx   = Mat<type_f>(NLEV_half, NX_half);
    //horizontal divergence
    h_div               = Mat<type_f>(NLEV_full, NX_full);
    Coeff               = Col<type_f>(NLEV_full);
    //tendencies for prognostic vars
    dphsdt              = Row<type_f>(NX_full);
    u_lhs               = Mat<type_f>(NLEV_full, NX_half);
    u_rhs               = Mat<type_f>(NLEV_full, NX_half);
    w_lhs               = Mat<type_f>(NLEV_half, NX_full);
    layer_pt_lhs        = Mat<type_f>(NLEV_full, NX_full);
    gz_lhs              = Mat<type_f>(NLEV_half, NX_full);

    hpgf_1              = Mat<type_f>(NLEV_full, NX_half);
    hpgf_2              = Mat<type_f>(NLEV_full, NX_half);
    u_vert_adv          = Mat<type_f>(NLEV_full, NX_half);
    K_hori_adv          = Mat<type_f>(NLEV_full, NX_half);

    Set_flags_false();
}

void State::Set_flags_false() {
    //valid signals
    dztopdx_full_valid = false;
    u_valid = false; u_cell_valid = false; u_vtx_valid = false; u_ns_valid = false;
    w_valid = false; w_cell_valid = false; w_vtx_valid = false;
    K_valid = false;
    pt_valid = false; pt_ns_valid = false; pt_we_valid = false;// pt_vtx_valid = false;
    ph_cell_valid = false; ph_ns_valid = false;// ph_we_valid = false; ph_vtx_valid = false;
    layer_ph_cell_valid = false; layer_ph_we_valid = false; layer_ph_ns_valid = false; layer_ph_vtx_valid = false;
    phs_valid = false; dphsdt_valid = false;
    p_valid = false; p_we_valid = false; p_vtx_valid = false; p_ns_valid = false;
    rho_valid = false; rho_we_valid = false;
    m_detadt_valid = false; m_detadt_cell_valid = false; m_detadt_vtx_valid = false;
    geo_potential_valid = false; geo_potential_cell_valid = false; geo_potential_vtx_valid = false;
    h_div_valid = false;
    dphsdt_valid = false; u_lhs_valid = false; u_rhs_valid = false; w_lhs_valid = false;
    layer_pt_lhs_valid = false; gz_lhs_valid = false;
}

//用方案二的坐标的初始化
void State::Pre_Process(type_f x_lo, type_f x_hi, type_f u_ref, type_f w_ref, type_f zh, type_f tau_0, type_f kexi,\
                        const char* zs_half_f, const char* zs_full_f, const char *A_half_f, const char *B_half_f) {
    /*  x_lo: x at left-most half, x_hi: x at right-most half
     *  topo: surface height on half
     *  A_half, B_half: params defined previously
     * */
    // prepare horizontal coordinates
    type_f x_delta = (x_hi - x_lo) / NX_full;
    for (int i = 0; i < NX_half; i++) x_half(i) = x_lo + x_delta * i;
    for (int i = 0; i < NX_full; i++) x_full(i) = (x_half(i) + x_half(i+1)) / 2.0;//_full always lays on the center!

    this->zh = zh;
    this->tau_0 = tau_0;
    this->kexi = kexi;
    this->u_ref = u_ref;
    this->w_ref = w_ref;

    //prepare topological data: zs_half, zs_full and dzs/dx_full
    zs_half.load(zs_half_f);
    zs_full.load(zs_full_f);
    for (int i = 0; i < NX_full; i++){
        dzsdx_full(i) = (zs_half(i+1) - zs_half(i)) / (x_half(i+1) - x_half(i));
    }

    //prepare vertical coeffs
    A_eta_half.load(A_half_f); B_eta_half.load(B_half_f);
    for (int k = 0; k < NLEV_half; k++)
        eta_half(k) = A_eta_half(k) + B_eta_half(k);
    for (int k = 0; k < NLEV_full; k++){
        A_eta_full(k) = (A_eta_half(k) + A_eta_half(k+1)) / 2.0;//_full always lays on the center!
        B_eta_full(k) = (B_eta_half(k) + B_eta_half(k+1)) / 2.0;
        eta_full(k)   = A_eta_full(k) + B_eta_full(k);
    }

    //determine the ph_top and z_top from relations of hydrostatic balance
    ph_top = eta_half(0) * p0;//assume p0 corresponds to z=0

    //prepare Coeff for horizontal divergence
    for (int k = 0; k < NLEV_full; k++){
        type_f p_ref = eta_full(k) * p0;
        Coeff(k) = 1.0/128.0 * max(1.0, 8.0*(1.0+tanh(log(ph_top/p_ref))) );
    }
}

//用方案一的坐标的初始化
void State::Pre_Process(type_f x_lo, type_f x_hi, type_f u_ref, type_f w_ref, type_f zh, type_f tau_0, type_f kexi,\
                        const char* zs_half_f, const char* zs_full_f, const char *A_half_f, const char *B_half_f, const char* A_full_f, const char* B_full_f) {
    /*  x_lo: x at left-most half, x_hi: x at right-most half
     *  topo: surface height on half
     *  A_half, B_half: params defined previously
     * */
    // prepare horizontal coordinates
    type_f x_delta = (x_hi - x_lo) / NX_full;
    for (int i = 0; i < NX_half; i++) x_half(i) = x_lo + x_delta * i;
    for (int i = 0; i < NX_full; i++) x_full(i) = (x_half(i) + x_half(i+1)) / 2.0;//_full always lays on the center!

    this->zh = zh;
    this->tau_0 = tau_0;
    this->kexi = kexi;
    this->u_ref = u_ref;
    this->w_ref = w_ref;

    //prepare topological data: zs_half, zs_full and dzs/dx_full
    zs_half.load(zs_half_f);
    zs_full.load(zs_full_f);
    for (int i = 0; i < NX_full; i++){
        dzsdx_full(i) = (zs_half(i+1) - zs_half(i)) / (x_half(i+1) - x_half(i));
    }

    //prepare vertical coeffs
    A_eta_half.load(A_half_f); B_eta_half.load(B_half_f);
    A_eta_full.load(A_full_f); B_eta_full.load(B_full_f);
    type_f delta_eta = 1. / NLEV_full;
    for (int k = 0; k < NLEV_half; k++)
        eta_half(k) = k * delta_eta;
    for (int k = 0; k < NLEV_full; k++)
        eta_full(k) = eta_half(k) + 0.5 * delta_eta;

    //determine the ph_top and z_top from relations of hydrostatic balance
    ph_top = eta_half(0) * p0;//assume p0 corresponds to z=0

    //prepare Coeff for horizontal divergence
    for (int k = 0; k < NLEV_full; k++){
        type_f p_ref = eta_full(k) * p0;
        Coeff(k) = 1.0/128.0 * max(1.0, 8.0*(1.0+tanh(log(ph_top/p_ref))) );
    }
}


void State::initiate(const char *u0, const char *w0, const char *pt0, const char *geo_potential0, const char *phs0) {
    phs.load(phs0); phs_valid = true;
    calc_ph_at_cell_ns_we_vtx();
    calc_layer_ph_at_cell_ns_we_vtx();

    u.load(u0); u_valid = true;
    //interp_from_we_edges_to_cell(u, u_cell); u_cell_valid = true;
    //interp_from_we_edges_to_vtxs(u, u_vtx);    u_vtx_valid= true;
    //interp_from_cell_to_ns_edges(u_cell, u_ns); u_ns_valid= true;

    w.load(w0); w_valid = true;
    //interp_from_ns_edges_to_cell(w, w_cell); w_cell_valid = true;
    //interp_from_ns_edges_to_vtxs_UPWIND(w, w_vtx);  w_vtx_valid  = true;

    pt.load(pt0); pt_valid = true;
    //interp_from_cell_to_ns_edges(pt, pt_ns); pt_ns_valid = true;
    //interp_from_cell_to_we_edges_UPWIND(pt, pt_we); pt_we_valid = true;

    geo_potential.load(geo_potential0); geo_potential_valid = true;
    //interp_from_ns_edges_to_cell(geo_potential, geo_potential_cell); geo_potential_cell_valid = true;
    //interp_from_ns_edges_to_vtxs_UPWIND(geo_potential, geo_potential_vtx);  geo_potential_vtx_valid  = true;

    interp_pt_w_gz_u();

    //顺带更新z_top_half
    ztop_half = geo_potential_vtx.row(0) / gravity;
    for (int i = 0; i < NX_full; i++)
        dztopdx_full(i) = (ztop_half(i+1) - ztop_half(i)) / (x_half(i+1) - x_half(i));
    dztopdx_full_valid = true;
}

