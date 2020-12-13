# DiDA_ATM_XZ
2D (XZ) nonhydrostatic model 

Temporal: modified RK3

Spatial: upwind(O3)/cent-diff(O4) for pt, w, gz and m*deta/dt 's horizontal advection terms; cent-diff(O2) for all other terms

Algorithm: Horizontal Explicit and Vertical Implicit (HEVI)
