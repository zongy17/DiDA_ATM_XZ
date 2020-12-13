//
// Created by lenovo on 2020/12/3.
//

#include "integrate.h"

void reset_state_4(State* & curr, State* & next, State* inter_1, State* inter_2, State* inter_3){
    //clear flags
    curr->Set_flags_false();
    inter_1->Set_flags_false();
    inter_2->Set_flags_false();
    inter_3->Set_flags_false();

    swap(curr, next);
}

void reset_state_3(State* & curr, State* & next, State* inter_1, State* inter_2){
    //clear flags
    curr->Set_flags_false();
    inter_1->Set_flags_false();
    inter_2->Set_flags_false();

    swap(curr, next);
}

void integrate_per_substep(State* curr, State* last_sub, State* next_sub, type_f dt){
    //假定curr已经完成了diagnose, last_sub已经完成了tend_all(它提供显式倾向)，因此可以进行推进计算
    //forward
    advance_phs(curr, last_sub, next_sub, dt);//从当前时刻的地面静力气压倾向推进计算phs，内含预报phs后计算ph和layer_ph
    advance_pt(curr, last_sub, next_sub, dt);//推进计算cell上的pt
    //backward
    advance_w_gz_implicit(curr, last_sub, next_sub, dt);//垂直隐式计算w和gz

    //得到w和gz之后立刻计算rho和p，这时候不是用线性状态方程
    next_sub->calc_rho_from_dphidph_at_cell_interp_to_we();//诊断计算密度
    //next->calc_p_from_rho_at_cell_we_interp_to_vtx();//诊断计算气压，并插值
    //注意可能有误：Zhang Yi中用的是n时刻的状态curr_state来线性化外推next_sub的气压，而不是上一中间状态last_sub
    next_sub->calc_p_at_cell_linear_interp_to_ns_we_vtx(curr);
    next_sub->Damp_Pressure(last_sub, 0.0);

    next_sub->tend_u_rhs();//得到p之后立刻计算水平气压梯度力以便计算u
    advance_u(curr, last_sub, next_sub, dt);

    next_sub->interp_pt_w_gz_u();

    //顺带更新z_top_half和dztop/dx_full
    if (!next_sub->geo_potential_vtx_valid) {printf("advance_w_gz_implicit: next_sub->geo_potential_vtx unusable!\n"); exit(1);}
    next_sub->ztop_half = next_sub->geo_potential_vtx.row(0) / gravity;
    for (int i = 0; i < next_sub->NX_full; i++)
        next_sub->dztopdx_full(i) = (next_sub->ztop_half(i+1)-next_sub->ztop_half(i))\
                                  / (next_sub->x_half(i+1) - next_sub->x_half(i));
    next_sub->dztopdx_full_valid = true;

    //一旦算出一个新的state的所有预报变量之后就立刻令其进行诊断量和倾向的计算
    //next_sub->interp_from_cell_to_we_edges_UPWIND(next_sub->rho, next_sub->rho_we); next_sub->rho_we_valid = true;
    //next_sub->interp_from_cell_to_we_edges_UPWIND(next_sub->p, next_sub->p_we); next_sub->p_we_valid = true;
    //next_sub->interp_from_ns_edges_to_vtxs_UPWIND(next_sub->p_ns, next_sub->p_vtx); next_sub->p_vtx_valid = true;
    next_sub->tend_dphsdt();//诊断计算当前时刻的地面静力气压倾向
    next_sub->calc_m_detadt_interp_to_vtx_cell_UPWIND();//诊断计算垂直坐标速度和m的乘积: mη'，并插值
    next_sub->calc_K_from_u_at_cell();//诊断计算水平动能

    //计算所有显式倾向，为下一小步的循环服务
    next_sub->tend_all_explicit();//计算除了dphsdt的其它倾向
}

//advance_phs仍然是显式的
void advance_phs(State* curr, State* last_sub, State* next_sub, type_f dt){
    if (!curr->phs_valid) {printf("integrate_phs: curr->phs unusable!\n"); exit(1);}
    if (!last_sub->dphsdt_valid) {printf("integrate_phs: last_sub->dphsdt unusable!\n"); exit(1);}

    next_sub->phs = curr->phs + last_sub->dphsdt * dt;

    next_sub->phs_valid = true;
    next_sub->calc_ph_at_cell_ns_we_vtx();//由phs计算各位置的ph，以便于计算后面的pt u w geo_potential
    next_sub->calc_layer_ph_at_cell_ns_we_vtx();
}

//advance_pt仍然是显式的:但注意应该是从last_sub也就是ET那里开始推进的
void advance_pt(State* curr, State* last_sub, State* next_sub, type_f dt) {
    if (!curr->pt_valid) {printf("advance_pt: curr->pt unusable!\n"); exit(1);}
    if (!last_sub->layer_pt_lhs_valid) {printf("advance_pt: last_sub->layer_pt_lhs unusable!\n"); exit(1);}
    if (!next_sub->layer_ph_cell_valid) {printf("advance_pt: next->layer_ph_cell unusable!\n"); exit(1);}

    next_sub->pt = (curr->pt % curr->layer_ph_cell - last_sub->layer_pt_lhs * dt) / next_sub->layer_ph_cell;

    next_sub->pt_valid = true;
    //next_sub->interp_from_cell_to_ns_edges(next_sub->pt, next_sub->pt_ns); next_sub->pt_ns_valid = true;
    //next_sub->interp_from_cell_to_we_edges(next_sub->pt, next_sub->pt_we); next_sub->pt_we_valid = true;
}

//垂直隐式计算w和gz
void advance_w_gz_implicit(State* curr, State* last_sub, State* next_sub, type_f dt) {
    static type_f Ca = 1.0 - curr->kexi;
    static type_f Cb = curr->kexi;//Ca + Cb = 1
    type_f CB = Cb * dt * gravity, CA = Ca * dt * gravity;
    int NX_full = curr->NX_full, NLEV_full = curr->NLEV_full;
    int NX_half = NX_full + 1,   NLEV_half = NLEV_full + 1;

    //准备gz_1
    if (!curr->geo_potential_valid) {printf("advance_w_gz_implicit: curr->gz unusable!\n"); exit(1); }
    if (!last_sub->gz_lhs_valid) {printf("advance_w_gz_implicit: last_sub->gz_lhs unusable!\n"); exit(1); }
    if (!last_sub->w_valid) {printf("advance_w_gz_implicit: last_sub->w unusable!\n"); exit(1); }
    Mat<type_f> gz_1 = curr->geo_potential - dt * last_sub->gz_lhs + CA * last_sub->w;

    //准备w_1
    if (!curr->w_valid) {printf("advance_w_gz_implicit: curr->w unusable!\n"); exit(1); }
    if (!last_sub->w_lhs_valid) {printf("advance_w_gz_implicit: last_sub->w_lhs unusable!\n"); exit(1); }
    if (!last_sub->p_valid) {printf("advance_w_gz_implicit: last_sub->p unusable!\n"); exit(1); }
    if (!last_sub->p_ns_valid) {printf("advance_w_gz_implicit: last_sub->p_ns unusable!\n"); exit(1); }
    if (!last_sub->layer_ph_ns_valid) {printf("advance_w_gz_implicit: last_sub->layer_ph_ns unusable!\n"); exit(1); }
    Mat<type_f> w_1 = curr->w - dt * last_sub->w_lhs - gravity * dt;
    int k = 0;//改中心差分为单侧差分
    for (int i = 0; i < NX_full; i++)
        w_1(k,i) += CA * (last_sub->p(k,i) - last_sub->p_ns(k,i)) / last_sub->layer_ph_ns(k,i);
    for (int k = 1; k < NLEV_half-1; k++)
        for (int i = 0; i < NX_full; i++){
            w_1(k,i) += CA * (last_sub->p(k,i) - last_sub->p(k-1,i)) / last_sub->layer_ph_ns(k,i);
        }
    k = NLEV_half - 1;//改中心差分为单侧差分
    for (int i = 0; i < NX_full; i++)
        w_1(k,i) += CA * (last_sub->p_ns(k,i) - last_sub->p(k-1,i)) / last_sub->layer_ph_ns(k,i);

    //通过线性化的状态方程近似地求next_sub的p
    //准备(δp/δph)_1(version-1) 或者 (δp)_1(version-2)
    if (!curr->p_valid) {printf("advance_w_gz_implicit: curr->p unusable!\n"); exit(1);}
    if (!next_sub->layer_ph_cell_valid) {printf("advance_w_gz_implicit: next_sub->layer_ph_cell unusable!\n"); exit(1);}
    if (!next_sub->pt_valid) {printf("advance_w_gz_implicit: next_sub->pt unusable!\n"); exit(1);}
    if (!curr->layer_ph_cell_valid) {printf("advance_w_gz_implicit: curr->layer_ph_cell unusable!\n"); exit(1);}
    if (!curr->pt_valid) {printf("advance_w_gz_implicit: curr->pt unusable!\n"); exit(1);}
    if (!next_sub->layer_ph_ns_valid) {printf("advance_w_gz_implicit: next_sub->layer_ph_ns unusable!\n"); exit(1);}

    //version-1: 没有乘静力气压项
    Mat<type_f> layerp_div_ph_1(NLEV_half, NX_full);
    for (int k = 1; k < NLEV_half-1; k++)
        for (int i = 0; i < NX_full; i++){
            type_f layer_p = curr->p(k,i) - curr->p(k-1,i);
            type_f p_pt_ratio = curr->p(k,i)*(next_sub->layer_ph_cell(k,i)*next_sub->pt(k,i))\
                                            /(curr->layer_ph_cell(k,i)*curr->pt(k,i))\
                              - curr->p(k-1,i)*(next_sub->layer_ph_cell(k-1,i)*next_sub->pt(k-1,i))\
                                            /(curr->layer_ph_cell(k-1,i)*curr->pt(k-1,i));
            layerp_div_ph_1(k,i) = (layer_p + gamma * p_pt_ratio) / next_sub->layer_ph_ns(k,i);
        }
    /*
    //version-2: Zhang Yi
    Mat<type_f> layer_p1(NLEV_half, NX_full);
    for (int k = 1; k < NLEV_half-1; k++)
        for (int i = 0; i < NX_full; i++){
            type_f layer_p = curr->p(k,i) - curr->p(k-1,i);
            type_f p_pt_ratio = curr->p(k,i)*(next_sub->layer_ph_cell(k,i)*next_sub->pt(k,i))\
                                            /(curr->layer_ph_cell(k,i)*curr->pt(k,i))\
                              - curr->p(k-1,i)*(next_sub->layer_ph_cell(k-1,i)*next_sub->pt(k-1,i))\
                                            /(curr->layer_ph_cell(k-1,i)*curr->pt(k-1,i));
            layer_p1(k,i) = layer_p + gamma * p_pt_ratio;
        }
    */

    //准备系数G和H
    if (!curr->p_valid) {printf("advance_w_gz_implicit: curr->p unusable!\n"); exit(1);}
    if (!next_sub->layer_ph_ns_valid) {printf("advance_w_gz_implicit: next_sub->layer_ph_ns unusable!\n"); exit(1);}
    if (!curr->geo_potential_valid) {printf("advance_w_gz_implicit: curr->geo_potential unusable!\n"); exit(1);}
    Mat<type_f> G(NLEV_half, NX_full), H(NLEV_half, NX_full);

    //version-1: 没有乘静力气压项
    for (int k = 1; k < NLEV_half-1; k++)
        for (int i = 0; i < NX_full; i++){
            G(k,i) = gamma * curr->p(k,i) / next_sub->layer_ph_ns(k,i)\
                                          /(curr->geo_potential(k+1,i) - curr->geo_potential(k,i));
            H(k,i) = gamma * curr->p(k-1,i) / next_sub->layer_ph_ns(k,i)\
                                          /(curr->geo_potential(k,i) - curr->geo_potential(k-1,i));
        }
    /*
    //version-2: Zhang Yi
    for (int k = 1; k < NLEV_half-1; k++)
        for (int i = 0; i < NX_full; i++){
            G(k,i) = gamma * curr->p(k,i) /(curr->geo_potential(k+1,i) - curr->geo_potential(k,i));
            H(k,i) = gamma * curr->p(k-1,i) /(curr->geo_potential(k,i) - curr->geo_potential(k-1,i));
        }
    */
    //求解w的垂直三对角方程：写成标准格式A(k)*w(k) = B(k)*w(k+1) + C(k)*w(k-1) + D(k)
    if (!last_sub->u_ns_valid) {printf("advance_w_gz_implicit: last_sub->u_ns unusable!\n"); exit(1);}
    if (!last_sub->dztopdx_full_valid) {printf("advance_w_gz_implicit: last_sub->dztop/dx unusable!\n"); exit(1);}
    Col<type_f> A(NLEV_half), B(NLEV_half), C(NLEV_half), D(NLEV_half), res(NLEV_half);
    for (int i = 0; i < NX_full; i++){

        //version-1: 没有乘静力气压项
        A = 1.0 - CB * CB * (G.col(i) + H.col(i));
        B =     - CB * CB *  G.col(i);
        C =     - CB * CB *  H.col(i);
        D = w_1.col(i) + CB * layerp_div_ph_1.col(i);
        for (int k = 1; k < NLEV_half-1; k++)
            D(k) = D(k) - CB * G(k,i) * (gz_1(k+1,i)-gz_1(k,i)) + CB * H(k,i) * (gz_1(k,i)-gz_1(k-1,i));
        //设置边界条件
        A(0) = 1.0; B(0) = 0.0; C(0) = 0.0;//top
        //D(0) = last_sub->u_ns(0,i) * last_sub->dztopdx_full(i);
        D(0) = 0.0;//顶部强行赋为0
        A(NLEV_half-1) = 1.0; B(NLEV_half-1) = 0.0; C(NLEV_half-1) = 0.0;//bottom
        D(NLEV_half-1) = last_sub->u_ns(NLEV_half-1,i) * last_sub->dzsdx_full(i);
        /*
        //version-2: Zhang Yi
        A = next_sub->layer_ph_ns.col(i) - CB * CB * (G.col(i) + H.col(i));
        B =                              - CB * CB *  G.col(i);
        C =                              - CB * CB *  H.col(i);
        D = w_1.col(i) % next_sub->layer_ph_ns.col(i) + CB * layer_p1.col(i);
        for (int k = 1; k < NLEV_half-1; k++)
            D(k) = D(k) - CB * G(k,i) * (gz_1(k+1,i)-gz_1(k,i)) + CB * H(k,i) * (gz_1(k,i)-gz_1(k-1,i));
        //设置边界条件
        A(0) = 1.0; B(0) = 0.0; C(0) = 0.0;//top
        D(0) = last_sub->u_ns(0,i) * last_sub->dztopdx_full(i);
        //D(0) = 0.0;
        A(NLEV_half-1) = 1.0; B(NLEV_half-1) = 0.0; C(NLEV_half-1) = 0.0;//bottom
        D(NLEV_half-1) = last_sub->u_ns(NLEV_half-1,i) * last_sub->dzsdx_full(i);
        */
        //调用Thomas算法求解w
        Thomas_solve(A, B, C, D, res);
        //结果拷贝回去
        next_sub->w.col(i) = res;
    }
    next_sub->w_valid = true;

    //有w之后求解gz
    next_sub->geo_potential = gz_1 + CB * next_sub->w;
    //gz的下边界条件（上边界在之前一并处理了）
    next_sub->geo_potential.row(NLEV_half-1) = gravity * next_sub->zs_full;
    next_sub->geo_potential_valid = true;
    //Damping: modiﬁes the state of w immediately after the vertically implicit solver
    next_sub->implicit_Rayleigh_damping_w(dt);//implicit Rayleigh damping内含将w插值到w_cell和w_vtx

    next_sub->geo_potential = gz_1 + CB * next_sub->w;//再重新算gz
    next_sub->geo_potential.row(NLEV_half-1) = gravity * next_sub->zs_full;//gz的下边界条件（上边界在之前一并处理了）
    next_sub->geo_potential_valid = true;
}

//advance_u要等这个状态的w和gz得到update了之后才能计算:注意可能也是从last_sub也就是ET那里开始推进的，而不用n时刻的curr
void advance_u(State* curr, State* last_sub, State* next_sub, type_f dt){
    if (!curr->u_valid) {printf("advance_u: curr->u unusable!\n"); exit(1);}
    if (!next_sub->u_rhs_valid) {printf("advance_u: next_sub->u_rhs unusable!\n"); exit(1);}
    if (!last_sub->u_lhs_valid) {printf("advance_u: last_sub->u_lhs unusable!\n"); exit(1);}

    //if (!last_sub->u_rhs_valid) {printf("advance_u: last_sub->u_rhs unusable!\n"); exit(1);}
    //static type_f Ca = 1.0 - curr->kexi;
    //static type_f Cb = curr->kexi;//Ca + Cb = 1

    //next_sub->u = curr->u + (Cb * next_sub->u_rhs + Ca * last_sub->u_rhs - last_sub->u_lhs) * dt;
    next_sub->u = curr->u + (next_sub->u_rhs - last_sub->u_lhs) * dt;

#ifdef UPSTREAM
    //treat the boundary
    int NLEV_full = curr->NLEV_full, NX_half = curr->NX_half;
    for (int k = 0; k < NLEV_full; k++){
        next->u(k,0) = next->u_ref;//上游来流速度给定
        //下游出口速度充分发展
        next->u(k, NX_half-1) = next->u(k,NX_half-2) - (next->u(k,NX_half-3) - next->u(k,NX_half-2))\
            * (next->x_half(NX_half-2) - next->x_half(NX_half-1)) / (next->x_half(NX_half-3) - next->x_half(NX_half-2));
    }
#endif
    next_sub->u_valid = true;
    //next_sub->interp_from_we_edges_to_cell(next_sub->u, next_sub->u_cell);      next_sub->u_cell_valid = true;
    //next_sub->interp_from_cell_to_ns_edges(next_sub->u_cell, next_sub->u_ns);   next_sub->u_ns_valid= true;
    //next_sub->interp_from_we_edges_to_vtxs(next_sub->u, next_sub->u_vtx);       next_sub->u_vtx_valid= true;
}

void Thomas_solve(Col<type_f> & A, Col<type_f> & B, Col<type_f> & C, Col<type_f> & D, Col<type_f> & X){
    // Solve: A(i)*X(i) = B(i)*X(i+1) + C(i)*X(i-1) + D(i)
    // Note: A, B, C and D to be modified in the function
    //  [  A(0)  -B(0)                                      ] [X(0)  ]   [D(0)  ]    默认C(0) = 0.0
    //  [ -C(1)   A(1) -B(1)                                ] [X(1)  ]   [D(1)  ]
    //  [        -C(2)  A(2) -B(2)                          ] [X(2)  ]   [D(2)  ]
    //  [            ...  ...  ...                          ] [...   ] = [...   ]
    //  [                 ...  ...  ...                     ] [...   ]   [...   ]
    //  [                                                   ] [      ]   [      ]
    //  [                            -C(n-2)  A(n-2) -B(n-2)] [X(n-2)]   [D(n-2)]
    //  [                                    -C(n-1)  A(n-1)] [X(n-1)]   [D(n-1)]    默认B(n-1) = 0.0
    int N = A.size();
    //消元
    B(0) = B(0) / A(0);
    D(0) = D(0) / A(0);
    for(int i = 1; i < N; i++){
        B(i) = B(i) / (A(i) - C(i)*B(i-1));
        D(i) = (D(i) + C(i)*D(i-1)) / (A(i) - C(i)*B(i-1));
    }
    //回代
    X(N-1) = D(N-1);
    for(int i = N-2; i >= 0; i--)
        X(i) = B(i)*X(i+1) + D(i);
}

/*
void advance_w(State* curr, State* ET, State* next, type_f dt){
    if (!curr->w_valid) {printf("advance_w: curr->w unusable!\n"); exit(1);}
    if (!ET->dwdt_valid) {printf("advance_w: ET->dwdt unsable!\n"); exit(1);}
    if (!next->u_ns_valid) {printf("advance_w: next->u_ns unusable!\n"); exit(1);}
    if (!next->dztopdx_full_valid){printf("advance_w: next->dztopdx_full unusable!\n"); exit(1);}

    next->w = curr->w + ET->dwdt * dt;
    //treat the bottom
    int bot_row = curr->NLEV_half - 1;
    next->w.row(bot_row) = next->u_ns.row(bot_row) % next->dzsdx_full;
    //treat the top
    next->w.row(0) = next->u_ns.row(0) % next->dztopdx_full;

    next->w_valid = true;
    next->interp_from_ns_edges_to_cell(next->w, next->w_cell); next->w_cell_valid = true;
    next->interp_from_ns_edges_to_vtxs(next->w, next->w_vtx);  next->w_vtx_valid  = true;
}

void advance_geo_potential(State* curr, State* ET, State* next, type_f dt){
    if (!curr->geo_potential_valid) {printf("advance_geo_potential: curr->geo_potential unusable!\n"); exit(1);}
    if (!ET->dgzdt_valid) {printf("advance_geo_potential: ET->dgzdt unusable!\n"); exit(1);}

    next->geo_potential = curr->geo_potential + ET->dgzdt * dt;
    //treat the bottom
    int bot_row = curr->NLEV_half - 1;
    next->geo_potential.row(bot_row) = gravity * next->zs_full;

    next->geo_potential_valid = true;
    next->interp_from_ns_edges_to_cell(next->geo_potential, next->geo_potential_cell); next->geo_potential_cell_valid = true;
    next->interp_from_ns_edges_to_vtxs(next->geo_potential, next->geo_potential_vtx);  next->geo_potential_vtx_valid  = true;

    //顺带更新z_top_half和dztop/dx_full
    next->ztop_half = next->geo_potential_vtx.row(0) / gravity;
    for (int i = 0; i < next->NX_full; i++)
        next->dztopdx_full(i) = (next->ztop_half(i+1) - next->ztop_half(i)) / (next->x_half(i+1) - next->x_half(i));
    next->dztopdx_full_valid = true;
}
*/
