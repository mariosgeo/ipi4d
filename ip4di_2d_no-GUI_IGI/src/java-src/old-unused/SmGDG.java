/**
 * Short description...
 * Package 'interp': set of functions to apply the matrix-vector product G'DGx
 * for an input vector x, where G is the gradient operator and D the diffusion
 * tensor field contained in an input tensor array.
 * The matrix-vector product G'DG makes use of Dave Hale's routine
 * LocalDiffusionKernel in JTK library edu.mines.jtk.dsp.
 * The input tensor array must be defined previously using JTK.
 * Comments refer to the use of the package in the frame of image-guided inversion.
 * <p>
 * Function 'getEigenvalues' gets the eigenvalues of tensors
 * in all points of a 2D array (while Dave Hale's routines
 * only extract some at a given point).
 * <p>
 * Function 'getLaplacianI(tensor T,vector Gm, matrix J'J, scalar scale, vector s, scalar beta)'
 * computes the model update
 *    dm = Hessien^-1 * Gradient
 * by solving the linear system
 *    Hessien * dm = Gradient
 * ie (J'J + beta*Cm^-1) dm = Gm,
 * with the gradient Gm = J'*Cd*[G(X)-d)].
 * Vectors Gm and s are in array format of model size [m1,m2].
 * <p>
 * Function 'applyGDG(tensor T, float[][] x_in, float scale,float[][] meshs)'
 * computes the matrix-vector product y=G'DGx,
 *
 * ...
 *
 * @author Jieyi Zhou, Colorado School of Mines
 *         Francois Lavoue, Colorado School of Mines (applyGDG, applyIGDG)
 * @version October 20, 2015
 */

 package interp;  // FL: FIND MORE APPROPRIATE NAME FOR PACKAGE?

import interp.SmoothCovariance.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;
import dnp.*;


public class SmGDG{   //FL: FIND MORE APPROPRIATE NAME FOR CLASS?? (e.g. applyLocalDiffusionKernel?)


  // Get eigenvalues of structure tensors in all points of the model
  // (while Dave Hale's routine only get them locally at a given point)
  public static float[][][] getEigenvalues(EigenTensors2 t, float[][] au, float[][] av){
    t.getEigenvalues(au,av);          // get eigenvalues of structure tensor mesh.et2 defined in initcm2.m
    return new float[][][]{au,av};
  }


  //FL: FIND A MORE APPROPRIATE NAME FOR THIS FUNCTION?? (e.g. solveIGIModelUpdate?)
  //            interp.Sm.getLaplacianI(  mesh.et2,         Gm1,           J'J,  lagrn?, mesh.s=[]?,    lagrn?)
  public static float[][] getLaplacianI(Tensors2 t, float[][] x, float[][] jtj, float c,float[][] s,float beta){
    int n1 = x[0].length;                          // [n1,n2] = size(Gm1)
    int n2 = x.length;                             //         = [mesh.m1,mesh.m2]
    float[][] y = new float[n2][n1];               // new array of size (m1,m2)
    VecArrayFloat2 vp = new VecArrayFloat2(x);     //  "    "   "   "      "    (FL: x=Gm1 or empty??)
    VecArrayFloat2 vq = new VecArrayFloat2(y);     //  "    "   "   "      "
    Multiply m = new Multiply(x,t,jtj,c,s,beta);   // m = (J'J+beta*G'DG)*Gm1   %FL: should be only (J'J+beta*G'DG) ...??
    CgSolver cgs = new CgSolver(1e-4,10000);       // define new CG solver with params 1e-4,1e4 (FL: params=?)
    CgSolver.Info info = cgs.solve(m,vp,vq);       // Solve(A a, Vec b, Vec x), solve Ax=b, ie m*vq=vp, ie (J'J+beta*G'DG)*Gm1 * dm = Gm1
    return y;
  }


  // Apply (J'J+G'DG)
  private static class Multiply implements CgSolver.A {

    //               Gm1,   mesh.et2,           J'J,   scale,   mesh.s=[],      lagrn
    Multiply(float[][] x, Tensors2 t, float[][] jtj, float c, float[][] s, float beta)
    {
      n1 = x[0].length;
      n2 = x.length;
      m = jtj.length;    // m = length(J'J) = mesh.num_param (FL: NOT the same 'm' as in 'Multiply m = new Multiply(x,t,jtj,c,s,beta)', I guess... ??)
      _t = t;
      _jtj = jtj;
      _c = c;
      _s = s;
      _beta = beta;
    }

    public void apply(Vec x, Vec y) { 
      float[][] ax = ((VecArrayFloat2)x).getArray();   // ax = Gm1 = J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m
      float[][] ay = ((VecArrayFloat2)y).getArray();   // ay = output = G'DG*ax
      //copy(ax,ay);                                   // ay = ax   %init. y=x to apply (I+G'DG) instead of G'DG

      //-- http://dhale.github.io/jtk/api/index.html doc for Class LocalDiffusionKernel --
      // Local diffusion kernel for use in anisotropic diffusion filtering.
      //   This kernel is a filter that computes y += G'DGx = y+G'DG*x where G is a gradient operator, G' is its adjoint,
      // and D is a local diffusion tensor field that determines for each image sample the filter coefficients.
      //   Given y = 0, this kernel computes y = G'DGx. Given y = x, it computes y = (I+G'DG)x.
      LocalDiffusionKernel ldk = new LocalDiffusionKernel();
      //  apply(et2,c, s,ax,ay)   with   Tensors2 t, float c, float[][] s, float[][] ax, float[][] ay,   ay = output
      ldk.apply(_t,_c,_s,ax,ay);                   // ay = ay+G'DG*ax = G'DG * Gm1   or   (I+G'DG)*Gm1
      mul(_beta,ay,ay);                            // ay = beta * G'DG * Gm1
      for (int i2=0; i2<m; ++i2) {                 // for i2 = 1 : mesh.num_param
	int k2 = i2/n1;                            //     k2 = ii2   %i2 index in model of size [m1,m2]
	int k1 = (i2%n1!=0)?(i2%n1-1):(n1-1);      //     k1 = ii1   %i1   "   "    "   "   "      "
	for (int i1=0; i1<m; ++i1) {               //     for i1 = 1 : mesh.num_param
	  int j2 = i1/n1;                          //         j2 = jj2   %i2 index in model of size [m1,m2]
	  int j1 = (i1%n1!=0)?(i1%n1-1):(n1-1);    //         j1 = jj1   %i1   "   "    "   "   "      "
	  ay[j2][j1] += _jtj[i2][i1]*ax[k2][k1];   //         ay(j1,j2) = ay(j1,j2) + jtj(i1,i2) * Gm1(k1,k2)
	}                                          //                   = beta*G'DG*Gm1 + J'J*Gm1 = (J'J+beta*G'DG)*Gm1
      }                                            // FL: BUT WE WANT THE INVERSE!!!
    }

    private int n1,n2,m;
    private float _c,_beta;
    private Tensors2 _t;
    private float[][] _jtj, _s;

  }   //end class Multiply


  // Apply G'DG
  public static float[][] applyGDG(Tensors2 et2, float[][] x_in, float scale,float[][] meshs){
    int n1 = x_in[0].length;                          // [n1,n2] = size(Gm1)
    int n2 = x_in.length;                             //         = [mesh.m1,mesh.m2]
    float[][] y_out = new float[n2][n1];              // new array of size (m1,m2)

    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    //  apply(et2,c, s,ax,ay)
    ldk.apply(et2,scale,meshs,x_in,y_out);           // ay = ay+G'DG*ax = G'DG * Gm1

    return y_out;
  }


  // Apply I+G'DG
  public static float[][] applyIGDG(Tensors2 et2, float[][] x_in, float scale,float[][] meshs){
    int n1 = x_in[0].length;                  // [n1,n2] = size(Gm1)
    int n2 = x_in.length;                     //         = [mesh.m1,mesh.m2]
    float[][] y_out = new float[n2][n1];      // new array of size (m1,m2)
    copy(x_in,y_out);                         // y_out = x_in
 
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    //  apply(et2,c, s,ax,ay)
    ldk.apply(et2,scale,meshs,x_in,y_out);    // y = y + G'DG*x = (I+G'DG) * x

    return y_out;
  }


  public static float[][] getSmoothed(LocalSmoothingFilter lsf, Tensors2 et, float[][] a){
    lsf.apply(et,a,a);
    return a;
  }


  //FL: this does not really "get Cm" as a matrix. Instead, it applies it to a vector and output a vector...
  // It is now redundant with function applyGDG and applyIGDG, which enable more easily to choose to include
  // the identity matrix in the calculation or not.
  public static float[][] getCm(Tensors2 t, float c, float[][] s, float[][] q){
    LocalDiffusionKernel ldk = new LocalDiffusionKernel();
    ldk.apply(t,c,s,q,q);   // q = q+G'DG*q = (I+G'DG)*q with q=??
    return q;
  }


  //FL: another definition of Cm? (and its application to a vector q=?)
  //FL: actually, G'DG=Cm^-1, so this rather compute y = Cm^-1 q
  public static float[][] getCm(SmoothCovariance sc, Sampling s1, Sampling s2, Tensors2 t, float[][] q){
    sc.apply(s1,s2,t,q);
    return q;
  }


  // Apply the inverse of Cm to a vector q
  //FL: actually, G'DG=Cm^-1, so applying the inverse of G'DG is applying Cm and not its inverse...
  public static float[][] getCmInverse(SmoothCovariance sc, Sampling s1, Sampling s2, Tensors2 t, float[][] q){
    sc.applyInverse(s1,s2,t,q);
    return q;
  }


  public static float[][] gridForVariableTensorsI(SmoothCovariance scm, Sampling s1, Sampling s2, Tensors2 t, float[][] ap, float[][] jtj) {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float[][] aq = new float[n2][n1];
    VecArrayFloat2 vp = new VecArrayFloat2(ap);
    VecArrayFloat2 vq = new VecArrayFloat2(aq);
    Smooth s = new Smooth(s1,s2,t,scm,jtj);
    CgSolver cgs = new CgSolver(1e-4,10000);
    CgSolver.Info info = cgs.solve(s,vp,vq); // Solve(A a, Vec b, Vec x), solve Ax=b
    return aq;
  }


  // Apply (Inv(Cm)+J'J)
  private static class Smooth implements CgSolver.A {

    Smooth(Sampling s1, Sampling s2, Tensors2 t, SmoothCovariance cm, float[][] jtj)
    {
      n1 = s1.getCount();
      n2 = s2.getCount();
      m = jtj.length;
      _s1 = s1;
      _s2 = s2;
      _t = t;
      _cm = cm;
      _jtj = jtj;
    }

    public void apply(Vec x, Vec y) { 
      float[][] ax = ((VecArrayFloat2)x).getArray();   // ax = Gm1
      float[][] ay = ((VecArrayFloat2)y).getArray();
      copy(ax,ay);
      _cm.applyInverse(_s1,_s2,_t,ay);
      for (int i2=0; i2<m; ++i2) {                 // for i2 = 1 : mesh.num_param
	int k2 = i2/n1;                            //     k2 = ii2   %i2 index in model of size [m1,m2]
	int k1 = (i2%n1!=0)?(i2%n1-1):(n1-1);      //     k1 = ii1   %i1   "   "    "   "   "      "
	for (int i1=0; i1<m; ++i1) {               //     for i1 = 1 : mesh.num_param
	  int j2 = i1/n1;                          //         j2 = jj2   %i2 index in model of size [m1,m2]
	  int j1 = (i1%n1!=0)?(i1%n1-1):(n1-1);    //         j1 = jj1   %i1   "   "    "   "   "      "
	  ay[j2][j1] += _jtj[i2][i1]*ax[k2][k1];   //         ay(j1,j2) = ay(j1,j2) + J'J(i1,i2)*ax(k1,k2)
	}
      }
    }

    private int n1,n2,m;
    private Sampling _s1,_s2;
    private Tensors2 _t;
    private SmoothCovariance _cm;
    private float[][] _jtj;

  }   //end class Smooth

}   //end class Sm

