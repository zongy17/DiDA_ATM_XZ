//
// Created by lenovo on 2020/12/3.
//

#ifndef XZ_V11_STATE_H
#define XZ_V11_STATE_H

#include "params.h"

class State {
public:
    int NLEV_half, NLEV_full;
    int NX_half, NX_full;
    type_f ph_top;//height, hydrostatic pressure of model top
    type_f u_ref, w_ref, zh, tau_0;//zh: height of the Rayleigh damping layer, tau_0: τ_0 for Rayleight damping
    type_f kexi;

    //coordinate
    //suffix _full: refer to centroid, _half: refer to interface
    Row<type_f> x_full, x_half;
    Row<type_f> zs_half, zs_full;
    Row<type_f> ztop_half;
    Row<type_f> dzsdx_full; //dzs/dx:surface gradient
    Row<type_f> dztopdx_full;//dztop/dx:model tio gradient
    bool dztopdx_full_valid;
    Col<type_f> eta_full, eta_half;
    Col<type_f> A_eta_full, A_eta_half;
    Col<type_f> B_eta_full, B_eta_half;

    //suffix _we: on west/east edges, _ns: north/south edges, _vtx: on vertexes, _cell: on cells
    //_valid: imply if the interpolated vars are valid in current step
    //horizontal velocity
    Mat<type_f> u, u_vtx, u_cell, u_ns;
    bool u_valid, u_vtx_valid, u_cell_valid, u_ns_valid;
    //vertical velocity
    Mat<type_f> w, w_cell, w_vtx;
    bool w_valid, w_cell_valid, w_vtx_valid;
    //kinetic energy
    Mat<type_f> K;
    bool K_valid;
    //potential temperature
    Mat<type_f> pt, pt_ns, pt_we;
    bool pt_valid, pt_ns_valid, pt_we_valid;
    //hydrostatic pressure
    Mat<type_f> ph_cell, ph_ns;// ph_we, ph_vtx,
    Mat<type_f> layer_ph_cell, layer_ph_we, layer_ph_ns, layer_ph_vtx;
    bool ph_cell_valid, ph_ns_valid;// ph_we_valid, ph_vtx_valid,
    bool layer_ph_cell_valid, layer_ph_we_valid, layer_ph_ns_valid, layer_ph_vtx_valid;
    Row<type_f> phs;//surface hydrostatic pressure
    bool phs_valid;
    //pressure
    Mat<type_f> p, p_we, p_vtx, p_ns;
    bool p_valid, p_we_valid, p_vtx_valid, p_ns_valid;
    //density
    Mat<type_f> rho, rho_we;
    bool rho_valid, rho_we_valid;
    //vertical coordinate velocity
    Mat<type_f> m_detadt, m_detadt_cell, m_detadt_vtx;
    bool m_detadt_valid, m_detadt_cell_valid, m_detadt_vtx_valid;
    //geo potential (phi)
    Mat<type_f> geo_potential, geo_potential_cell, geo_potential_vtx;
    bool geo_potential_valid, geo_potential_cell_valid, geo_potential_vtx_valid;
    //horizontal divergence 水平散度
    Mat<type_f> h_div;
    Col<type_f> Coeff;//Coeff是算一次之后都不变的，在预处理中计算
    bool h_div_valid;
    //tendencies of prognostic vars
    Row<type_f> dphsdt;
    Mat<type_f> u_lhs, u_rhs, w_lhs, layer_pt_lhs, gz_lhs;//u_rhs实际上就是水平气压梯度力
    bool dphsdt_valid, u_lhs_valid, u_rhs_valid, w_lhs_valid, layer_pt_lhs_valid, gz_lhs_valid;

    Mat<type_f> hpgf_1, hpgf_2, u_vert_adv, K_hori_adv;

    State(int nlev_full, int nx_full);
    void Set_flags_false();//set all flags as false
    //Pre-Process functions
    void Pre_Process(type_f x_lo, type_f x_hi, type_f u_ref, type_f w_ref, type_f zh, type_f tau_0, type_f kexi,\
                    const char* zs_half_f, const char* zs_full_f, const char* A_half_f, const char* B_half_f);
    void Pre_Process(type_f x_lo, type_f x_hi, type_f u_ref, type_f w_ref, type_f zh, type_f tau_0, type_f kexi,\
                    const char* zs_half_f, const char* zs_full_f, const char* A_half_f, const char* B_half_f, const char* A_full_f, const char* B_full_f);
    //initiate function
    void initiate(const char* u0_f, const char* w0_f, const char* pt0_f,\
                  const char* geo_potential0_f, const char* phs0_f);//用文件输入
    //interp functions
    void interp_pt_w_gz_u();
    void interp_from_cell_to_we_edges_UPWIND(Mat<type_f> const & cell, Mat<type_f> & edges);
    void interp_from_ns_edges_to_vtxs_UPWIND(Mat<type_f> const & edges, Mat<type_f> & vtxs);
    void interp_from_cell_to_ns_edges(Mat<type_f> const & cell, Mat<type_f> & edges);
    void interp_from_cell_to_we_edges(Mat<type_f> const & cell, Mat<type_f> & edges);
    void interp_from_ns_edges_to_cell(Mat<type_f> const & edges, Mat<type_f> & cell);
    void interp_from_we_edges_to_cell(Mat<type_f> const & edges, Mat<type_f> & cell);
    void interp_from_ns_edges_to_vtxs(Mat<type_f> const & edges, Mat<type_f> & vtxs);
    void interp_from_we_edges_to_vtxs(Mat<type_f> const & edges, Mat<type_f> & vtxs);
    //void interp_from_ns_edges_to_vtxs_1D(Row<type_f> const & edges, Row<type_f> & vtxs);
    //diagnostic诊断计算，准备之后可以为下一时刻或中间时刻计算explicit tendency服务
    void diagnose();//状态方程用原始的非线性的
    void calc_ph_at_cell_ns_we_vtx();
    void calc_layer_ph_at_cell_ns_we_vtx();
    void tend_dphsdt();//计算静力气压的倾向
    void calc_m_detadt_interp_to_vtx_cell_UPWIND();//诊断计算垂直方向质量通量，mη'
    void calc_K_from_u_at_cell();
    void calc_rho_from_dphidph_at_cell_interp_to_we();//diagnose rho from δΦ/δΠ
    void calc_p_from_rho_at_cell_interp_to_we_ns_vtx();
    void calc_p_at_cell_linear_interp_to_ns_we_vtx(State* last);
    void calc_h_div_from_u_at_cell();//计算水平散度
    //tendencies(tend_dphsdt已经在诊断中计算过了)
    void tend_all_explicit();
    void tend_u_lhs();
    void tend_u_rhs();
    void tend_w_lhs();
    void tend_layer_pt_lhs();
    void tend_gz_lhs();
    //damping
    void Dh_damping(type_f dt);//水平散度耗散
    void Rayleigh_damping_u(type_f dt);
    void Rayleigh_damping_w(type_f dt);//文献Wong中的耗散方法
    void implicit_Rayleigh_damping_w(type_f dt);//Klemp 2008
    void Damp_Pressure(State* last, type_f beta_d);

    //Print functions for debug
    void Print_coordinates();
    void Print_Console_All_Prognostic(const char* extra_txt);
    static void Write_dat(Row<type_f> const & row_var, const char* filename);
    static void Write_dat(Col<type_f> const & col_val, const char* filename);
    static void Write_dat(Mat<type_f> const & mat_val, const char* filename);
    //Save for cont'd
    void Save_all();
};

#endif //XZ_V11_STATE_H
