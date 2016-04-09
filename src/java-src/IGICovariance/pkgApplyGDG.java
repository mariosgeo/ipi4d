/**
 * Package 'applyGDG': set of functions to apply the matrix-vector product G'DGx
 * for an input vector x, where G is the gradient operator and D the diffusion
 * tensor field contained in an input tensor array.
 * The matrix-vector product G'DG makes use of Dave Hale's routine LocalDiffusionKernel in
 * JTK library edu.mines.jtk.dsp. The input tensor array must be defined previously using JTK.
 * Comments refer to the use of the package in the frame of image-guided inversion.
 * <p>
 * Function 'getEigenvalues' gets the eigenvalues of tensors at all points of a 2D array
 * (while Dave Hale's routines only extract them at a given point).
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
 * <p>
 * Function 'applyIGDG(tensor T, float[][] x_in, float scale,float[][] meshs)'
 * computes the matrix-vector product y=(I+G'DG)x,
 *
 * @author Jieyi Zhou, Colorado School of Mines (getEigenvalues)
 *         Francois Lavoue, Colorado School of Mines (applyGDG, applyIGDG)
 * @version October 20, 2015
 */

package IGICovariance;

import dnp.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;


public class pkgApplyGDG{


  // Get eigenvalues of structure tensors in all points of the model
  // (while Dave Hale's routine only get them locally at a given point)
  public static float[][][] getEigenvalues(EigenTensors2 t, float[][] au, float[][] av){
    t.getEigenvalues(au,av);          // get eigenvalues of structure tensor mesh.et2 defined in initcm2.m
    return new float[][][]{au,av};
  }


  // Apply G'DG
  public static float[][] applyGDG(Tensors2 et2, float[][] x_in, float scale,float[][] meshs, String stencil){
    int n1 = x_in[0].length;                          // [n1,n2] = size(Gm1)
    int n2 = x_in.length;                             //         = [mesh.m1,mesh.m2]
    float[][] y_out = new float[n2][n1];              // new array of size (m1,m2)

    //LocalDiffusionKernel ldk = new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.stencil);
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D33);

    //  apply(et2,c, s,ax,ay)
    ldk.apply(et2,scale,meshs,x_in,y_out);           // ay = ay+G'DG*ax

    return y_out;
  }


  // Apply I+G'DG
  public static float[][] applyIGDG(Tensors2 et2, float[][] x_in, float scale,float[][] meshs, String stencil){
    int n1 = x_in[0].length;                  // [n1,n2] = size(Gm1)
    int n2 = x_in.length;                     //         = [mesh.m1,mesh.m2]
    float[][] y_out = new float[n2][n1];      // new array of size (m1,m2)
    copy(x_in,y_out);                         // y_out = x_in
 
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D33);

    //  apply(et2,c, s,ax,ay)
    ldk.apply(et2,scale,meshs,x_in,y_out);    // y = y + G'DG*x = (I+G'DG) * x

    return y_out;
  }


}   //end class pkgApplyGDG


