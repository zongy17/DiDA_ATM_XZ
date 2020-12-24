//
// Created by lenovo on 2020/12/3.
//

#include "State.h"

void State::tend_layer_pt_lhs() {
    if (layer_pt_lhs_valid) {printf("tend_dlayer_pt_lhs: already done.\n"); return; }
    if (!layer_ph_we_valid) {printf("tend_dlayer_pt_lhs: layer_ph_we unusable!\n"); exit(1);}
    if (!u_valid) {printf("tend_dlayer_pt_lhs: u unusable!\n"); exit(1);}
    if (!pt_we_valid) {printf("tend_dlayer_pt_lhs: pt_we unusable!\n"); exit(1);}
    if (!m_detadt_valid) {printf("tend_dlayer_pt_lhs: m_detadt unusable!\n"); exit(1);}
    if (!pt_ns_valid) {printf("tend_dlayer_pt_lhs: pt_ns unusable!\n"); exit(1);}

    /*
    static type_f beta = 1.0;
    type_f beta_div2 = beta / 2.0;

    int l1, l2, cen, r1, r2;
    type_f pt_l, pt_r;
    Mat<type_f> sign_u = sign(u);//迎风标志
    */
    type_f hori_adv, vert_adv;
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_full; i++){
            /*
            cen = i;
            l1 = (cen-1)%NX_full; l2 = (cen-2)%NX_full;
            r1 = (cen+1)%NX_full; r2 = (cen+2)%NX_full;
            pt_l = 7.0/12.0*(pt(k,cen)+pt(k,l1)) - 1.0/12.0*(pt(k,r1)+pt(k,l2))\
                 + sign_u(k,i)*beta_div2*(pt(k,r1)-pt(k,l2) - 3.0*(pt(k,cen)-pt(k,l1)));
            pt_r = 7.0/12.0*(pt(k,r1)+pt(k,cen)) - 1.0/12.0*(pt(k,r2)+pt(k,l1))\
                 + sign_u(k,i+1)*beta_div2*(pt(k,r2)-pt(k,l1) - 3.0*(pt(k,r1)-pt(k,cen)));
            hori_adv =  (layer_ph_we(k,i+1)*u(k,i+1)*pt_r - layer_ph_we(k,i)*u(k,i)*pt_l )\
                        / (x_half(i+1) - x_half(i));
            */
            hori_adv = ( layer_ph_we(k,i+1)*u(k,i+1)*pt_we(k,i+1)\
                              - layer_ph_we(k,i)*u(k,i)*pt_we(k,i) ) / (x_half(i+1) - x_half(i));
            vert_adv = m_detadt(k+1,i)*pt_ns(k+1,i) - m_detadt(k,i)*pt_ns(k,i);
            layer_pt_lhs(k,i) = hori_adv + vert_adv;
        }
    }
    layer_pt_lhs_valid = true;
}

void State::tend_u_lhs() {
    if (u_lhs_valid) {printf("tend_u_lhs: already done.\n"); return; }
    if (!u_valid) {printf("tend_u_lhs: u unusable!\n"); exit(1);}
    if (!u_vtx_valid) {printf("tend_u_lhs: u_vtx unusable!\n"); exit(1);}
    if (!K_valid) {printf("tend_u_lhs: K unusable!\n"); exit(1);}
    if (!m_detadt_vtx_valid) {printf("tend_u_lhs: m_detadt_vtx unusable!\n"); exit(1);}
    if (!layer_ph_we_valid) {printf("tend_u_lhs: layer_ph_we unusable!\n"); exit(1);}

    type_f delta_x, hori_adv, vert_adv;

    for (int k = 0; k < NLEV_full; k++){
        for (int i = 1; i < NX_half-1; i++){//middle part
            delta_x = x_full(i) - x_full(i-1);
            hori_adv = (K(k,i) - K(k,i-1)) / delta_x;
            /*
            //张祎方法
            vert_adv = (m_detadt_vtx(k+1,i)*(u_vtx(k+1,i)-u(k,i))\
                       +m_detadt_vtx(k,i)*(u(k,i)-u_vtx(k,i))) / layer_ph_we(k,i);
            */
            //能量守恒形式
            vert_adv = 0.0;
            vert_adv += (k==0)           ? 0 : ( m_detadt_vtx(k,i)*(u(k,i)-u(k-1,i)) );
            vert_adv += (k==NLEV_full-1) ? 0 : ( m_detadt_vtx(k+1,i)*(u(k+1,i)-u(k,i)) );
            vert_adv /= (2.0*layer_ph_we(k,i));

            u_lhs(k,i) = hori_adv + vert_adv;
            u_vert_adv(k,i) = vert_adv;
            K_hori_adv(k,i) = hori_adv;
        }
#ifdef PERIODIC
        //i == 0 and i == NX_half-1
        int i = 0;
        delta_x = x_full(i) - x_half(i) + x_half(NX_half-1) - x_full(NX_full-1);
        hori_adv = (K(k,i) - K(k,NX_full-1)) / delta_x;
        /*
        //张祎方法
        vert_adv = (m_detadt_vtx(k+1,i)*(u_vtx(k+1,i)-u(k,i))\
                   +m_detadt_vtx(k,i)*(u(k,i)-u_vtx(k,i))) / layer_ph_we(k,i);
        */
        //能量守恒形式
        vert_adv = 0.0;
        vert_adv += (k==0)           ? 0 : ( m_detadt_vtx(k,i)*(u(k,i)-u(k-1,i)) );
        vert_adv += (k==NLEV_full-1) ? 0 : ( m_detadt_vtx(k+1,i)*(u(k+1,i)-u(k,i)) );
        vert_adv /= (2.0*layer_ph_we(k,i));

        u_lhs(k,i) = hori_adv + vert_adv;
        u_lhs(k,NX_half-1) = u_lhs(k,i);//周期性边界条件
        u_vert_adv(k,i) = vert_adv;
        u_vert_adv(k,NX_half-1) = vert_adv;
        K_hori_adv(k,i) = hori_adv;
        K_hori_adv(k,NX_half-1) = hori_adv;
#endif
#ifdef UPSTREAM
        u_lhs(k,0) = 0.0;//应该用不到这个倾向，因为直接赋为0
        u_lhs(k,NX_half-1) = 0.0;//应该用不到这个倾向，因为直接速度外推
#endif
    }
    u_lhs_valid = true;
}

void State::tend_u_rhs() {
    if (u_rhs_valid) {printf("tend_u_rhs: already done.\n"); return; }
    if (!p_vtx_valid) {printf("tend_u_rhs: p_vtx unusable!\n"); exit(1);}
    if (!layer_ph_we_valid) {printf("tend_u_rhs: layer_ph_we unusable!\n"); exit(1);}
    if (!ph_ns_valid) {printf("tend_u_rhs: ph_ns unusable!\n"); exit(1);}
    if (!geo_potential_valid) {printf("tend_u_rhs: geo_potential unusable!\n"); exit(1);}
    if (!p_ns_valid) {printf("tend_u_rhs: p_ns unusable!\n"); exit(1);}
    if (!rho_we_valid)  {printf("tend_u_rhs: rho_we unusable!\n"); exit(1);}

    for (int k = 0; k < NLEV_full; k++){
        for (int i = 1; i < NX_half-1; i++){//middle part
            type_f delta_x = x_full(i) - x_full(i-1);
            /*
            //原方法存在大项小差问题
            type_f rhs_1    = - (p(k,i) - p(k,i-1)) / delta_x / rho_we(k,i);
            type_f rhs_2    = - (p_vtx(k+1,i)-p_vtx(k,i)) / layer_ph_we(k,i)\
                            *(geo_potential_cell(k,i)-geo_potential_cell(k,i-1))/delta_x;
            dudt(k,i) = - hori_adv - vert_adv + rhs_1 + rhs_2;
            */
            //林仙剑有限体积方法
            type_f line_int_ph = (ph_ns(k+1,i)-ph_ns(k,i-1) + ph_ns(k+1,i-1)-ph_ns(k,i));
            type_f pgf_1 = ((ph_ns(k+1,i)-ph_ns(k,i-1))*(geo_potential(k+1,i-1)-geo_potential(k,i))\
                           +(ph_ns(k+1,i-1)-ph_ns(k,i))*(geo_potential(k,i-1)-geo_potential(k+1,i))\
                           ) * (p_vtx(k+1,i)-p_vtx(k,i)) / layer_ph_we(k,i) / delta_x / line_int_ph;
            type_f pgf_2 = ((ph_ns(k+1,i)-ph_ns(k,i-1))*(p_ns(k+1,i-1)-p_ns(k,i))\
                           +(ph_ns(k+1,i-1)-ph_ns(k,i))*(p_ns(k,i-1)-p_ns(k+1,i))\
                           ) / rho_we(k,i) / delta_x / line_int_ph;
            u_rhs(k,i) = pgf_1 + pgf_2;
            hpgf_1(k,i)= pgf_1;
            hpgf_2(k,i)= pgf_2;
        }
#ifdef PERIODIC
        //i == 0 and i == NX_half-1
        int i = 0;
        type_f delta_x = x_full(i) - x_half(i) + x_half(NX_half-1) - x_full(NX_full-1);
        /*
        //原方法存在大项小差问题
        type_f rhs_1 = - (p(k,i) - p(k,NX_full-1)) / delta_x / rho_we(k,i);
        type_f rhs_2 = - (p_vtx(k+1,i)-p_vtx(k,i)) / layer_ph_we(k,i)\
                       *(geo_potential_cell(k,i)-geo_potential_cell(k,NX_full-1))/delta_x;
        dudt(k,i) = - hori_adv - vert_adv + rhs_1 + rhs_2;
        */
        //林仙剑有限体积方法
        type_f pgf_1 = (p_vtx(k+1,i)-p_vtx(k,i))/layer_ph_we(k,i)\
                     * ((ph_ns(k+1,i)-ph_ns(k,NX_full-1))*(geo_potential(k+1,NX_full-1)-geo_potential(k,i))\
                       +(ph_ns(k+1,NX_full-1)-ph_ns(k,i))*(geo_potential(k,NX_full-1)-geo_potential(k+1,i))\
                       ) / delta_x / (ph_ns(k+1,i)-ph_ns(k,NX_full-1) + ph_ns(k+1,NX_full-1)-ph_ns(k,i));
        type_f pgf_2 = ((ph_ns(k+1,i)-ph_ns(k,NX_full-1))*(p_ns(k+1,NX_full-1)-p_ns(k,i))\
                       +(ph_ns(k+1,NX_full-1)-ph_ns(k,i))*(p_ns(k,NX_full-1)-p_ns(k+1,i))\
                       ) / rho_we(k,i) / delta_x / (ph_ns(k+1,i)-ph_ns(k,NX_full-1) + ph_ns(k+1,NX_full-1)-ph_ns(k,i));
        u_rhs(k,i) = pgf_1 + pgf_2;
        u_rhs(k,NX_half-1) = u_rhs(k,i);//周期性边界条件
        hpgf_1(k,i)= pgf_1;
        hpgf_2(k,i)= pgf_2;
#endif
#ifdef UPSTREAM
        u_rhs(k,0) = 0.0;//应该用不到这个倾向，因为直接赋为0
        u_rhs(k,NX_half-1) = 0.0;//应该用不到这个倾向，因为直接速度外推
#endif
    }
    u_rhs_valid = true;
}

void State::tend_w_lhs() {
    if (w_lhs_valid) {printf("tend_w_lhs: already done.\n"); return ;}
    if (!layer_ph_vtx_valid) {printf("tend_w_lhs: layer_ph_vtx unusable!\n"); exit(1);}
    if (!u_vtx_valid) {printf("tend_w_lhs: u_vtx unusable!\n"); exit(1);}
    if (!w_valid) {printf("tend_w_lhs: w unusable!\n"); exit(1);}
    if (!w_vtx_valid) {printf("tend_w_lhs: w_vtx unusable!\n"); exit(1);}
    if (!layer_ph_ns_valid) {printf("tend_w_lhs: layer_ph_ns unusable!\n"); exit(1);}
    if (!m_detadt_cell_valid) {printf("tend_w_lhs: m_detadt_cell unusable!\n"); exit(1);}
    if (!w_cell_valid) {printf("tend_w_lhs: w_cell unusable!\n"); exit(1);}

    type_f hori_adv, vert_adv, delta_x;
    type_f delta_mdetadt_w, delta_mdetadt, deta_up, deta_down;
    //k == 0: top, 这个倾向应该不会用来算
    int k = 0;//改中心差分为单侧差分，实验发现w_lhs在顶部和底部都可以不用算，而直接设为0
    for (int i = 0; i < NX_full; i++) {

        delta_x    = x_half(i+1) - x_half(i);
        hori_adv   = (layer_ph_vtx(k,i+1)*u_vtx(k,i+1)*(w_vtx(k,i+1)-w(k,i))\
                     +layer_ph_vtx(k,i)*u_vtx(k,i)*(w(k,i)-w_vtx(k,i)) )\
                     / delta_x / layer_ph_ns(k,i);
        //垂直对流项改为单侧差分
        vert_adv   = m_detadt_cell(k,i)*(w_cell(k,i)-w(k,i))/ layer_ph_ns(k,i);
        w_lhs(k,i) = hori_adv + vert_adv;

        //w_lhs(0, i) = 0.0;
    }
    //middle part
    for (int k = 1; k < NLEV_half-1; k++){
        for (int i = 0; i < NX_full; i++){
            delta_x    = x_half(i+1) - x_half(i);
            hori_adv   = (layer_ph_vtx(k,i+1)*u_vtx(k,i+1)*(w_vtx(k,i+1)-w(k,i))\
                         +layer_ph_vtx(k,i)*u_vtx(k,i)*(w(k,i)-w_vtx(k,i)) )\
                         / delta_x / layer_ph_ns(k,i);
            /*
            vert_adv   = (m_detadt_cell(k,i)*(w_cell(k,i)-w(k,i))\
                         +m_detadt_cell(k-1,i)*(w(k,i)-w_cell(k-1,i)) )\
                         / layer_ph_ns(k,i);
            */
            deta_up = eta_full(k)-eta_half(k);
            deta_down = eta_half(k) - eta_full(k-1);
            delta_mdetadt = deta_down/deta_up * (m_detadt_cell(k,i)-m_detadt(k,i))\
                          - deta_up/deta_down * (m_detadt_cell(k-1,i)-m_detadt(k,i));
            delta_mdetadt_w = deta_down/deta_up\
                  * (m_detadt_cell(k,i)*w_cell(k,i) - m_detadt(k,i)*w(k,i))\
                             - deta_up/deta_down\
                  * (m_detadt_cell(k-1,i)*w_cell(k-1,i) - m_detadt(k,i)*w(k,i));
            vert_adv = (delta_mdetadt_w - w(k,i) * delta_mdetadt) / layer_ph_ns(k,i);
            w_lhs(k,i) = hori_adv + vert_adv;
        }
    }
    //k == NLEV_half-1: bottom，这个倾向应该不会用来算
    k = NLEV_half - 1;
    for (int i = 0; i < NX_full; i++){

        //水平对流项不用变
        delta_x    = x_half(i+1) - x_half(i);
        hori_adv   = (layer_ph_vtx(k,i+1)*u_vtx(k,i+1)*(w_vtx(k,i+1)-w(k,i))\
                     +layer_ph_vtx(k,i)*u_vtx(k,i)*(w(k,i)-w_vtx(k,i)) )\
                     / delta_x / layer_ph_ns(k,i);
        //垂直对流项改为单侧差分
        vert_adv   = m_detadt_cell(k-1,i)*(w(k,i)-w_cell(k-1,i)) / layer_ph_ns(k,i);
        w_lhs(k,i) = hori_adv + vert_adv;

        //w_lhs(NLEV_half-1,i) = 0.0;
    }
    w_lhs_valid = true;
}

void State::tend_gz_lhs() {
    if (gz_lhs_valid) {printf("tend_gz_lhs: already done.\n"); return ;}
    if (!layer_ph_vtx_valid) {printf("tend_gz_lhs: layer_ph_vtx unusable!\n"); exit(1);}
    if (!u_vtx_valid) {printf("tend_gz_lhs: u_vtx unusable!\n"); exit(1);}
    if (!geo_potential_vtx_valid) {printf("tend_gz_lhs: geo_potential_vtx unusable!\n"); exit(1);}
    if (!geo_potential_valid) {printf("tend_gz_lhs: geo_potential unusable!\n"); exit(1);}
    if (!m_detadt_cell_valid) {printf("tend_gz_lhs: m_detadt_cell unusable!\n"); exit(1);}
    if (!geo_potential_cell_valid) {printf("tend_gz_lhs: geo_potential_cell unusable!\n"); exit(1);}
    if (!layer_ph_ns_valid) {printf("tend_gz_lhs: layer_ph_ns unusable!\n"); exit(1);}

    type_f hori_adv, vert_adv, delta_x;
    type_f delta_mdetadt_gz, delta_mdetadt, deta_up, deta_down;
    //k == 0: top，将中心差分改成单侧差分
    int k = 0;
    for (int i = 0; i < NX_full; i++){
        //水平对流项不用变
        delta_x    = x_half(i+1) - x_half(i);
        hori_adv   = (layer_ph_vtx(k,i+1)*u_vtx(k,i+1)*(geo_potential_vtx(k,i+1)-geo_potential(k,i))\
                     +layer_ph_vtx(k,i)*u_vtx(k,i)*(geo_potential(k,i)-geo_potential_vtx(k,i)) )\
                     / delta_x / layer_ph_ns(k,i);
        //垂直对流项改为单侧差分
        vert_adv   = m_detadt_cell(k,i)*(geo_potential_cell(k,i) - geo_potential(k,i))/layer_ph_ns(k,i);

        gz_lhs(k,i) = hori_adv + vert_adv;
    }
    //middle part
    for (int k = 1; k < NLEV_half-1; k++){
        for (int i = 0; i < NX_full; i++){
            delta_x = x_half(i+1) - x_half(i);
            hori_adv   = (layer_ph_vtx(k,i+1)*u_vtx(k,i+1)*(geo_potential_vtx(k,i+1)-geo_potential(k,i))\
                         +layer_ph_vtx(k,i)*u_vtx(k,i)*(geo_potential(k,i)-geo_potential_vtx(k,i)) )\
                         / delta_x / layer_ph_ns(k,i);
            /*
            vert_adv   = ( m_detadt_cell(k,i)*(geo_potential_cell(k,i)-geo_potential(k,i))\
                         + m_detadt_cell(k-1,i)*(geo_potential(k,i)-geo_potential_cell(k-1,i)) )\
                         / layer_ph_ns(k,i);
            */
            deta_up = eta_full(k)-eta_half(k);
            deta_down = eta_half(k) - eta_full(k-1);
            delta_mdetadt = deta_down/deta_up * (m_detadt_cell(k,i)-m_detadt(k,i))\
                          - deta_up/deta_down * (m_detadt_cell(k-1,i)-m_detadt(k,i));
            delta_mdetadt_gz = deta_down/deta_up\
                  * (m_detadt_cell(k,i)*geo_potential_cell(k,i) - m_detadt(k,i)*geo_potential(k,i))\
                             - deta_up/deta_down\
                  * (m_detadt_cell(k-1,i)*geo_potential_cell(k-1,i) - m_detadt(k,i)*geo_potential(k,i));
            vert_adv = (delta_mdetadt_gz - geo_potential(k,i) * delta_mdetadt) / layer_ph_ns(k,i);
            gz_lhs(k,i) = hori_adv + vert_adv;
        }
    }
    k = NLEV_half - 1;// bottom，地表的gz取决于地形，一直不变，倾向为0，似乎这个倾向应该不会用来算
    //但如果简单地设gz_lhs(NLEV_half-1,i) = 0.0;会导致隐式计算出现相差和值偏大的问题
    //这是因为在解三对角的时候，在算离底部次近的会用到gz_1的最底部的一项
    //对于显式计算，是不会有此问题的
    for (int i = 0; i < NX_full; i++){
        //水平对流项不用变
        delta_x = x_half(i+1) - x_half(i);
        hori_adv   = (layer_ph_vtx(k,i+1)*u_vtx(k,i+1)*(geo_potential_vtx(k,i+1)-geo_potential(k,i))\
                     +layer_ph_vtx(k,i)*u_vtx(k,i)*(geo_potential(k,i)-geo_potential_vtx(k,i)) )\
                     / delta_x / layer_ph_ns(k,i);
        //垂直对流项改为单侧差分
        vert_adv   = m_detadt_cell(k-1,i)*(geo_potential(k,i) - geo_potential_cell(k-1,i))/layer_ph_ns(k,i);
        gz_lhs(k,i) = hori_adv + vert_adv;
    }
    gz_lhs_valid = true;
}

void State::tend_all_explicit() {
    tend_layer_pt_lhs();
    tend_u_lhs();
    tend_w_lhs();
    tend_gz_lhs();
}

