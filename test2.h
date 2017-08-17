template <class SYS, typename real_type, typename vect_type>
void calc_force_on_i_p2(SYS &sys, const int nbody, const int i, const real_type eps2) {
vect_type acc0(real_type(0));
vect_type acc1(real_type(0));
vect_type icoord0 = sys.template get_coord<0>(i);
vect_type icoord1 = sys.template get_coord<1>(i);
for(int j=0; j<nbody; j++){
vect_type dr0 = sys.template get_coord<0>(j) - icoord0;
vect_type dr1 = sys.template get_coord<1>(j) - icoord1;

real_type s0 = eps2 + (dr0*dr0);
if(eps2==s0) continue;
real_type s1 = (dr0*dr1);

#if 1
real_type rinv1 = rsqrt(s0);
real_type rinv2 = rinv1 * rinv1;
#else
real_type rinv2 = 1.0 / s0;
real_type rinv1 = sqrt(rinv2);
#endif
real_type rinv3 = rinv1 * rinv2;
rinv3 *= sys.get_mass(j);

real_type q1 = rinv2 * (-3.0*s1);

acc0 += rinv3 * (dr0);
acc1 += rinv3 * (dr1 +(q1)*dr0);
} // for(j) 
sys.template set_force<0>(i, acc0);
sys.template set_force<1>(i, acc1);
} // end of function 
