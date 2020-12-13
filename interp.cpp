//
// Created by lenovo on 2020/12/3.
//

#include "State.h"

void State::interp_from_cell_to_ns_edges(const Mat<type_f> &cell, Mat<type_f> &edges) {
    interp2(x_full, eta_full, cell, x_full, eta_half, edges);
    /*
    //first: calc middle part
    for (int k = 1; k < NLEV_half-1; k++){
        for (int i = 0; i < NX_full; i++){
            edges(k,i) = ( cell(k,i)          * (eta_half(k)-eta_full(k-1)) \
                         + cell(k-1,i) * (eta_full(k)-eta_half(k)     ) \
                         ) / (eta_full(k) - eta_full(k-1));
        }
    }*/
    //second: extrapolate the bottom and top
    Col<type_f>::fixed<3> etas, vars, coeff;
    type_f x;
    for (int i = 0; i < NX_full; i++){
        //二次外插
        //printf("quadratic extrapolation...\n");
        //bottom
        etas = {eta_full(NLEV_full-1), eta_half(NLEV_half-2), eta_full(NLEV_full-2)};
        vars = {cell(NLEV_full-1,i), edges(NLEV_half-2,i), cell(NLEV_full-2,i)};
        polyfit(coeff, etas, vars, 2);
        x = eta_half(NLEV_half-1);
        edges(NLEV_half-1,i) = coeff(0)*x*x + coeff(1)*x + coeff(2);
        //top
        etas = {eta_full(0), eta_half(1), eta_full(1)};
        vars = {cell(0,i), edges(1,i), cell(1,i)};
        polyfit(coeff, etas, vars, 2);
        x = eta_half(0);
        edges(0,i) = coeff(0)*x*x + coeff(1)*x + coeff(2);
        /*
        //线性外插
        printf("linear extrapolation...\n");
        edges(0,i) = 2.0 * cell(0,i) - edges(1,i);
        edges(NLEV_half-1,i) = 2.0 * cell(NLEV_full-1,i) - edges(NLEV_half-2,i);
        */
    }
}

void State::interp_from_cell_to_we_edges(const Mat<type_f> &cell, Mat<type_f> &edges) {//主要是用来插值pt
    interp2(x_full, eta_full, cell, x_half, eta_full, edges);
#ifdef UPSTREAM
    Col<type_f>::fixed<3> xs, vars, coeff;
    type_f x;
#endif
    for (int k = 0; k < NLEV_full; k++){//处理边界
#ifdef PERIODIC
        /*
        //middle part
        for (int i = 1; i < NX_half-1; i++){
            edges(k,i) = ( cell(k,i)          * (x_half(i)-x_full(i-1)) \
                         + cell(k,i-1) * (x_full(i)-x_half(i))\
                         ) / (x_full(i) - x_full(i-1));
        }
         */
        // i == 0
        type_f right_dx = x_full(0) - x_half(0), left_dx = x_half(NX_half-1) - x_full(NX_full-1);
        type_f right_val = cell(k,0), left_val = cell(k,NX_full-1);
        edges(k,0) = (right_val*left_dx + left_val*right_dx) / (left_dx+right_dx);
        // and i == NX_half-1
        edges(k,NX_half-1) = edges(k,0);//周期性边界条件
#endif
#ifdef UPSTREAM
        //i == 0: 二次外插
        xs = {x_full(0), x_half(1), x_full(1)};
        vars = {cell(k,0), edges(k,1), cell(k,1)};
        polyfit(coeff, xs, vars, 2);
        x = x_half(0);
        edges(k,0) = coeff(0)*x*x + coeff(1)*x + coeff(2);
        //i == NX_half-1: 线性外推，即梯度为0
        edges(k,NX_half-1) = 2.0*cell(k,NX_full-1) - edges(k,NX_half-2);
#endif
    }
}

void State::interp_from_ns_edges_to_cell(Mat<type_f> const & edges, Mat<type_f> & cell) {
    interp2(x_full, eta_half, edges, x_full, eta_full, cell);
    /*
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++)
            cell(k,i) = (edges(k,i) + edges(k+1,i)) / 2.0;
    */
}

void State::interp_from_we_edges_to_cell(Mat<type_f> const & edges, Mat<type_f> & cell) {
    interp2(x_half, eta_full, edges, x_full, eta_full, cell);
    /*
    for (int k = 0; k < NLEV_full; k++)
        for (int i = 0; i < NX_full; i++)
            cell(k,i) = (edges(k,i) + edges(k,i+1)) / 2.0;
    */
}

void State::interp_from_ns_edges_to_vtxs(Mat<type_f> const & edges, Mat<type_f> & vtxs) {
    //edges: [NLEV_half * NX_full] , vtxs: [NLEV_half * NX_half]
    interp2(x_full, eta_half, edges, x_half, eta_half, vtxs);
#ifdef UPSTREAM
    Col<type_f>::fixed<3> xs, vars, coeff;
    type_f x;
#endif
    for (int k = 0; k < NLEV_half; k++){//处理边界条件
#ifdef PERIODIC
        /*
        //middle part
        for (int i = 1; i < NX_half-1; i++){
            vtxs(k,i) = ( edges(k,i)*(x_half(i) - x_full(i-1)) + edges(k,i-1)*(x_full(i) - x_half(i)) \
                        ) / (x_full(i) - x_full(i-1));
        }
         */
        // i == 0
        type_f right_dx = x_full[0] - x_half[0], left_dx = x_half[NX_half-1] - x_full[NX_full-1];
        type_f right_val = edges(k,0), left_val = edges(k,NX_full-1);
        vtxs(k,0) = (right_val*left_dx + left_val*right_dx) / (left_dx+right_dx);
        // i == NX_half-1
        vtxs(k,NX_half-1) = vtxs(k,0);//周期性边界条件
#endif
#ifdef UPSTREAM
        //i == 0: 二次外插
        xs = {x_full(0), x_half(1), x_full(1)};
        vars = {edges(k,0), vtxs(k,1), edges(k,1)};
        polyfit(coeff, xs, vars, 2);
        x = x_half(0);
        vtxs(k,0) = coeff(0)*x*x + coeff(1)*x + coeff(2);
        //i == NX_half-1: 线性外推，即梯度为0
        vtxs(k,NX_half-1) = 2.0*edges(k,NX_full-1) - vtxs(k,NX_half-2);
#endif
    }
}

void State::interp_from_we_edges_to_vtxs(Mat<type_f> const & edges, Mat<type_f> & vtxs) {
    //edges: [NLEV_full * NX_half] , vtxs: [NLEV_half * NX_half]
    interp2(x_half, eta_full, edges, x_half, eta_half, vtxs);
    /*
    //middle part
    for (int k = 1; k < NLEV_half-1; k++){
        for (int i = 0; i < NX_half; i++){
            vtxs(k,i) = ( edges(k,i)          * (eta_half(k)-eta_full(k-1))\
                        + edges(k-1,i) * (eta_full(k)-eta_half(k))\
                        ) / (eta_full(k) - eta_full(k-1));
        }
    }*/
    //extrapolate to k == 0 and k == NLEV_half-1
    Col<type_f>::fixed<3> etas, vars, coeff;
    type_f x;
    for (int i = 0; i < NX_half; i++){
        //二次外插
        //bottom
        etas = {eta_full(NLEV_full-1), eta_half(NLEV_half-2), eta_full(NLEV_full-2)};
        vars = {edges(NLEV_full-1,i), vtxs(NLEV_half-2,i), edges(NLEV_full-2,i)};
        polyfit(coeff, etas, vars, 2);
        x = eta_half(NLEV_half-1);
        vtxs(NLEV_half-1,i) = coeff(0)*x*x + coeff(1)*x + coeff(2);
        //top
        etas = {eta_full(0), eta_half(1), eta_full(1)};
        vars = {edges(0,i), vtxs(1,i), edges(1,i)};
        polyfit(coeff, etas, vars, 2);
        x = eta_half(0);
        vtxs(0,i) = coeff(0)*x*x + coeff(1)*x + coeff(2);
        /*
        //线性外插
        printf("linear extrapolation...\n");
        vtxs(0,i) = 2.0 * edges(0,i) - vtxs(1,i);
        vtxs(NLEV_half-1,i) = 2.0 * edges(NLEV_full-1,i) - vtxs(NLEV_half-2,i);
        */
    }
}

/*
void State::interp_pt_we_upwind() {
    static type_f beta = 1.0;
    type_f beta_div2 = beta / 2.0;
    int l1, l2, r1, r2;
    Mat<type_f> sign_u = sign(u);//迎风标志

    for (int k = 0; k < NLEV_full; k++) {
        for (int i = 0; i < NX_half; i++) {
            r1 = i;
            r2 = (r1+1)%NX_full;
            l1 = (r1-1)%NX_full;
            l2 = (l1-2)%NX_full;
            pt_we(k,i) = 7.0/12.0*(pt(k,r1)+pt(k,l1)) - 1.0/12.0*(pt(k,r2)+pt(k,l2))\
                       + sign_u(k,i)*beta_div2*(pt(k,r2)-pt(k,l2) - 3.0*(pt(k,r1)-pt(k,l1)));
        }
    }
    pt_we_valid = true;
}
*/

void State::interp_from_cell_to_we_edges_UPWIND(const Mat<type_f> &cell, Mat<type_f> &edges) {
    if (!u_valid) {printf("interp_from_cell_to_we_edges_UPWIND: u unusable!\n"); exit(1); }
    Mat<type_f> sign_u = sign(u);//迎风标志

    //O3/O4
    //static type_f beta = 0.25;//这个beta是对于位温的，可以取0.25试试
    static type_f beta = 1.0;
    type_f beta_div12 = beta / 12.0;
    int l1, l2, r1, r2;
    for (int k = 0; k < NLEV_full; k++) {
        for (int i = 0; i < NX_half; i++) {
            r1 =  i           %NX_full;
            r2 = (i+1+NX_full)%NX_full;
            l1 = (i-1+NX_full)%NX_full;
            l2 = (i-2+NX_full)%NX_full;
            edges(k,i) = 7.0/12.0*(cell(k,r1)+cell(k,l1))\
                       - 1.0/12.0*(cell(k,r2)+cell(k,l2))\
                       + sign_u(k,i)*beta_div12*(\
                                      cell(k,r2)-cell(k,l2)\
                               - 3.0*(cell(k,r1)-cell(k,l1))\
                       );
        }
    }
    /*
    //O5/O6
    static type_f beta = 1.0;
    type_f beta_div60 = beta / 60.0;
    int l1,l2,l3,r1,r2,r3;
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_half; i++){
            r1 =  i           %NX_full;
            r2 = (i+1+NX_full)%NX_full;
            r3 = (i+2+NX_full)%NX_full;
            l1 = (i-1+NX_full)%NX_full;
            l2 = (i-2+NX_full)%NX_full;
            l3 = (i-3+NX_full)%NX_full;
            edges(k,i) = 37.0/60.0*(cell(k,r1)+cell(k,l1))\
                       - 2.0/15.0* (cell(k,r2)+cell(k,l2))\
                       + 1.0/60.0* (cell(k,r3)+cell(k,l3))\
                       - sign_u(k,i)*beta_div60*(\
                                      cell(k,r3)-cell(k,l3)
                               - 5.0*(cell(k,r2)-cell(k,l2))\
                               +10.0*(cell(k,r1)-cell(k,l1))\
                       );
        }
    }
    */
}

void State::interp_from_ns_edges_to_vtxs_UPWIND(const Mat<type_f> &edges, Mat<type_f> &vtxs) {
    if (!u_vtx_valid) {printf("interp_from_ns_edges_to_vtxs_UPWIND: u_vtx unusable!\n"); exit(1); }
    Mat<type_f> sign_u_vtx = sign(u_vtx);//迎风标志

    //O3/O4
    static type_f beta = 1.0;//这个beta是对于w和Φ的，取1.0
    type_f beta_div12 = beta / 12.0;
    int l1, l2, r1, r2;
    for (int k = 0; k < NLEV_half; k++) {
        for (int i = 0; i < NX_half; i++) {
            r1 =  i           %NX_full;
            r2 = (i+1+NX_full)%NX_full;
            l1 = (i-1+NX_full)%NX_full;
            l2 = (i-2+NX_full)%NX_full;
            vtxs(k,i) = 7.0/12.0*(edges(k,r1)+edges(k,l1))
                      - 1.0/12.0*(edges(k,r2)+edges(k,l2))\
                      + sign_u_vtx(k,i)*beta_div12*(\
                                   edges(k,r2)-edges(k,l2)\
                            - 3.0*(edges(k,r1)-edges(k,l1))\
                      );
        }
    }
    /*
    //O5/O6
    static type_f beta = 1.0;
    type_f beta_div60 = beta / 60.0;
    int l1,l2,l3,r1,r2,r3;
    for (int k = 0; k < NLEV_full; k++){
        for (int i = 0; i < NX_half; i++){
            r1 =  i           %NX_full;
            r2 = (i+1+NX_full)%NX_full;
            r3 = (i+2+NX_full)%NX_full;
            l1 = (i-1+NX_full)%NX_full;
            l2 = (i-2+NX_full)%NX_full;
            l3 = (i-3+NX_full)%NX_full;
            vtxs(k,i)  = 37.0/60.0*(edges(k,r1)+edges(k,l1))\
                       - 2.0/15.0* (edges(k,r2)+edges(k,l2))\
                       + 1.0/60.0* (edges(k,r3)+edges(k,l3))\
                       - sign_u_vtx(k,i)*beta_div60*(\
                                      edges(k,r3)-edges(k,l3)
                               - 5.0*(edges(k,r2)-edges(k,l2))\
                               +10.0*(edges(k,r1)-edges(k,l1))\
                       );
        }
    }
    */
}

void State::interp_pt_w_gz_u() {
    //另一个预报变量ph因为需要在forward step中用，所以每次算出来phs之后立刻计算，不在此
    if (!u_valid) {printf("interp_pt_w_gz_u: u unusable!\n"); exit(1);}
    interp_from_we_edges_to_cell(u, u_cell);        u_cell_valid = true;
    interp_from_we_edges_to_vtxs(u, u_vtx);         u_vtx_valid= true;
    interp_from_cell_to_ns_edges(u_cell, u_ns);     u_ns_valid= true;

    if (!pt_valid) {printf("interp_pt_w_gz_u: pt unusable!\n"); exit(1);}
    interp_from_cell_to_we_edges_UPWIND(pt, pt_we); pt_we_valid = true;
    interp_from_cell_to_ns_edges(pt, pt_ns);        pt_ns_valid = true;

    if (!w_valid) {printf("interp_pt_w_gz_u: w unusable!\n"); exit(1);}
    interp_from_ns_edges_to_vtxs_UPWIND(w, w_vtx);  w_vtx_valid = true;
    interp_from_ns_edges_to_cell(w, w_cell);        w_cell_valid= true;

    if (!geo_potential_valid) {printf("interp_pt_w_gz_u: gz unusable!\n"); exit(1);}
    interp_from_ns_edges_to_vtxs_UPWIND(geo_potential, geo_potential_vtx);  geo_potential_vtx_valid = true;
    interp_from_ns_edges_to_cell(geo_potential, geo_potential_cell);        geo_potential_cell_valid= true;
}