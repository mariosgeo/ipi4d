package interp;
import interp.SmoothCovariance.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;
import dnp.*;

public class Sm{

  public static float[][][] getEigenvalues(EigenTensors2 t, float[][] au, float[][] av){
    t.getEigenvalues(au,av);          // get eigenvalues of structure tensor mesh.et2 defined in initcm2.m
    return new float[][][]{au,av};
  }

  //            interp.Sm.getLaplacianI(  mesh.et2,         dm1,           J'J,  lagrn?, mesh.s=[]?,    lagrn?)
  public static float[][] getLaplacianI(Tensors2 t, float[][] x, float[][] jtj, float c,float[][] s,float beta){
    int n1 = x[0].length;                          // [n1,n2] = size(dm1)
    int n2 = x.length;                             //         = [mesh.m1,mesh.m2]
    float[][] y = new float[n2][n1];               // new array of size (m1,m2)
    VecArrayFloat2 vp = new VecArrayFloat2(x);     //  "    "   "   "      "
    VecArrayFloat2 vq = new VecArrayFloat2(y);     //  "    "   "   "      "
    Multiply m = new Multiply(x,t,jtj,c,s,beta);   //
    CgSolver cgs = new CgSolver(1e-4,10000);       // define new CG solver with params ?? 1e-4,1e4
    CgSolver.Info info = cgs.solve(m,vp,vq);       // Solve(A a, Vec b, Vec x), solve Ax=b, ie m*vq=dm1

    getCm Cm_out = new getCm(t,c,s,q)              // FL (TEST): compute Cm explicitly and output it

    return y;
  }


  // Apply (J'J+G'DG)
  private static class Multiply implements CgSolver.A {

    //               dm1,   mesh.et2,           J'J,  lagrn?,  mesh.s=[]?,     lagrn?
    Multiply(float[][] x, Tensors2 t, float[][] jtj, float c, float[][] s, float beta)
    {
      n1 = x[0].length;
      n2 = x.length;
      m = jtj.length;    // m = length(J'J) = mesh.num_param
      _t = t;
      _jtj = jtj;
      _c = c;
      _s = s;
      _beta = beta;
    }

    public void apply(Vec x, Vec y) { 
      float[][] ax = ((VecArrayFloat2)x).getArray();   // ax = dm1 = J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m
      float[][] ay = ((VecArrayFloat2)y).getArray();   // ay = output = G'DG*ax
      copy(ax,ay);

      //-- http://dhale.github.io/jtk/api/index.html doc for Class LocalDiffusionKernel --
      // Local diffusion kernel for use in anisotropic diffusion filtering.
      //   This kernel is a filter that computes y += G'DGx = y+G'DG*x where G is a gradient operator, G' is its adjoint,
      // and D is a local diffusion tensor field that determines for each image sample the filter coefficients.
      //   A local diffusion kernel is typically used in combinations with others. For example, the filter implied by
      // (I+G'DG)y = G'DGx acts as a notch filter. It attenuates features for which G'DG is zero while preserving other features.
      // The diffusion tensors in D control the width, orientation, and anisotropy of the spectral notch. Note that application
      // of this filter requires solution of a sparse symmetric positive-definite system of equations. An even simpler example is
      // the filter implied by (I+G'DG)y = x. This filter smooths features in the directions implied by the tensors D. Again, 
      // application of this filter requires solving a sparse symmetric positive-definite system of equations.
      //   The accumulation of the kernel output in y = y+G'DGx is useful when constructing such filter combinations. 
      // Given y = 0, this kernel computes y = G'DGx. Given y = x, it computes y = (I+G'DG)x.
      LocalDiffusionKernel ldk = new LocalDiffusionKernel();
      //  apply(et2,c, s,ax,ay)   with   Tensors2 t, float c, float[][] s, float[][] ax, float[][] ay,   ay = output
      ldk.apply(_t,_c,_s,ax,ay);                   // ay = ay+G'DG*ax = G'DG * dm1
      mul(_beta,ay,ay);                            // ay = beta * G'DG * dm1
      for (int i2=0; i2<m; ++i2) {                 // for i2 = 1 : mesh.num_param
	int k2 = i2/n1;                            //     k2 = ii2   %i2 index in model of size [m1,m2]
	int k1 = (i2%n1!=0)?(i2%n1-1):(n1-1);      //     k1 = ii1   %i1   "   "    "   "   "      "
	for (int i1=0; i1<m; ++i1) {               //     for i1 = 1 : mesh.num_param
	  int j2 = i1/n1;                          //         j2 = jj2   %i2 index in model of size [m1,m2]
	  int j1 = (i1%n1!=0)?(i1%n1-1):(n1-1);    //         j1 = jj1   %i1   "   "    "   "   "      "
	  ay[j2][j1] += _jtj[i2][i1]*ax[k2][k1];   //         ay(j1,j2) = ay(j1,j2) + jtj(i1,i2) * dm1(k1,k2)
	}                                          //                   = l*G'DG*dm1 + J'J*dm1 = (J'J+l*G'DG)*dm1
      }
    }

    private int n1,n2,m;
    private float _c,_beta;
    private Tensors2 _t;
    private float[][] _jtj, _s;

  }   //end class Multiply


  public static float[][] getSmoothed(LocalSmoothingFilter lsf, Tensors2 et, float[][] a){
    lsf.apply(et,a,a);
    return a;
  }


 // get Cm: to be output of getLaplacianI???
 public static float[][] getCm(Tensors2 t, float c, float[][] s, float[][] q){
   LocalDiffusionKernel ldk = new LocalDiffusionKernel();
      ldk.apply(t,c,s,q,q);   // q = q+G'DG*q = (I+G'DG)*q with q=??
    return q;
  }
 
  // or this one???
  public static float[][] getCm(SmoothCovariance sc, Sampling s1, Sampling s2, Tensors2 t, float[][] q){
    sc.apply(s1,s2,t,q);
    return q;
  }

  // get Cm^-1
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
      float[][] ax = ((VecArrayFloat2)x).getArray();   // ax = dm1
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


