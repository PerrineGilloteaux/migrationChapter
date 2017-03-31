import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;

/**
 * v0.0 -> as submitted and published in JCS 
 * Starting point of the plugin: a stack on which one roi corresponding to the protrusion for each frame has need added
 * A lot of To be checked.
 * @author paul-gilloteaux-p
 * @TODO add the preprocessing step as part of the software (with back up of the segmentation, so that it can start from a movie alone.
 * @TODO remove the centroid parameter if the normal to the tengent method is used.
 * @TODO for now only ROI area would be processed (such as the one created when using Analyse particle on a binary mask. 
 * Could be more generic on any line roi
 */

public class Recruitment_Edge_Analysis
  implements PlugIn
{
  double[] profile1frame;
  double[] profile;
  double[] profilespeed;
  String method;
  boolean verbose;
  
  public Recruitment_Edge_Analysis() {}
  
  /**
   * ImageJ plugin run function (will be called when launching the plugin)
   * 
   * 
   */
  public void run(String arg)
  {
    ImagePlus imp = WindowManager.getCurrentImage();
    // Some checks:
    // an image should be open
    if (imp == null)
    {
      IJ.error("A stack and a roi manager containing the cell boundaries should be openned and active");
      return;
    }
    int numStacks = imp.getStackSize();
    ImageProcessor ip = imp.getProcessor();
    imp.setSlice(1);
    ip = imp.getProcessor();
    // this image should be a 2D movie
    if (numStacks == 1) {
      IJ.error("Sorry! This Plugin only works on Stacks, and it has to be the active window");
      return;
    }
    
    // Cell protrusion have to be presegmented
    RoiManager myroiManager = RoiManager.getInstance();
    if (myroiManager == null) {
      IJ.error("Roi manager should be openned and contains one roi per stack, obtained for exemple by analyse particle on the thresholded stack. Hint: You may want to use the size filtering in Analyze particle in order to get rid of small noisy area outside the cell. ");
      return;
    }
    // only one roi per stack in this version.
    Roi[] roilist = myroiManager.getRoisAsArray();
    if (roilist.length != numStacks) {
      IJ.error("Roi manager should be openned and contains one roi per stack, obtained for exemple by analyse particle on the thresholded stack. Hint: You may want to use the size filtering in Analyze particle in order to get rid of small noisy area outside the cell. ");
      return;
    }
    

    // Interface description
    int Ninit = 0;
    String[] items = new String[]{"Normal(perpendicular to tangent)", "Ray Casting"};
    ImageStatistics stats = ip.getStatistics();
    double CentroidX = stats.xCenterOfMass;
    double CentroidY = stats.yCenterOfMass;
    GenericDialog gd = new GenericDialog("Recruitment and Edge Dynamics Computation");
    gd.addMessage("Recruitment_Edge_dynamics takes as input a stack, for now already processed in order to get an ROI per frame. "
    		+ "It gives as an output both a kymograph of the speed and the intensity of the edge");
    gd.addMessage(" A reference point has to be chosen from which the ray will be casted (RAY CASTING OPTION) , and to sign the speed (all) ");
    gd.addNumericField("x", CentroidX, 0);
    gd.addNumericField("y", CentroidY, 0);
    gd.addNumericField("Sampling", 1.0D, 0);
    gd.addCheckbox(" Smooth boundary? running average", false);
    gd.addCheckbox(" Verbose (more log message)", false);
    gd.addChoice(" Method for band estimation", items, "Ray Casting");
    gd.addNumericField("Line width (the width of the band each side of the boundary of the recruitment", 3.0D, 0);
    gd.showDialog();
    if (gd.wasCanceled()) { return;
    }
    CentroidX = gd.getNextNumber();
    CentroidY = gd.getNextNumber();
    int sampling = (int)gd.getNextNumber();
    int linewidth = (int)gd.getNextNumber();
    boolean smooth = gd.getNextBoolean();
    this.verbose = gd.getNextBoolean();
    this.method = gd.getNextChoice();
    IJ.log("REFERENCE point " + IJ.d2s(CentroidX) + " " + IJ.d2s(CentroidY));
    IJ.log("Line width  " + IJ.d2s(linewidth));
    IJ.log("Method  " + this.method);
    
   
    
    
    //initialize LUT
    float[] reds = new float[32];
    float[] greens = new float[32];
    float[] blues = new float[32];
    int factor = ice(reds, greens, blues);
    //initialize draw
    Plot mydraw = new Plot("Matching points for indication", "x", "y");
    mydraw.setLimits(0.0D, ip.getWidth(), 0.0D, ip.getHeight());
    mydraw.setSize(ip.getWidth(), ip.getHeight());
    
    //initialize contour point list on first frame 
    //(under the strong assumption that roi list first item is the roi on first frame)
    FloatPolygon pinit = roilist[0].getInterpolatedPolygon(sampling, smooth);
    
    float[] intxinit = pinit.xpoints;
    float[] intyinit = pinit.ypoints;
    double[] xinittmp = new double[pinit.npoints];
    double[] yinittmp = new double[pinit.npoints];
    int j = 0;
    for (int i = 0; i < pinit.npoints; i++)
    {
      // remove points touching the image border (typically when the ROI comes from a zoom on a protrusion 	
      if ((intxinit[i] < 2.0F) || (intyinit[i] < 2.0F) || (intxinit[i] > ip.getWidth() - 2) || (intyinit[i] > ip.getHeight() - 2))
      {
        if (verbose) {
          IJ.log("Border point removed " + IJ.d2s(intxinit[i]) + IJ.d2s(intyinit[i]));
        }
      }
      else {
        xinittmp[j] = intxinit[i];
        yinittmp[j] = intyinit[i];
        j++;
      }
    }
    Ninit = j; // number of points sampled in first frame
    double[] xinit = new double[Ninit];
    double[] yinit = new double[Ninit];
    for (int i = 0; i < Ninit; i++)
    {
      xinit[i] = xinittmp[i];
      yinit[i] = yinittmp[i];
    }
    
    //plot the first sampled points (in first frame) in the draw.
    mydraw.setColor(Color.GREEN);
    mydraw.addPoints(xinit, yinit, 0);
    
    // 
    double[] pos1 = new double[2];
    double[] pos2 = new double[2];
    int fshift=3; // maximum shift? tbc
    for (int i = 0; i < Ninit; i++)
    {
      if (this.method.compareTo("Ray Casting") == 0)
      {
        pos1 = this.GetShiftedPositionRaycasting(-fshift, xinit[i], yinit[i], CentroidX, CentroidY);
        pos2 = this.GetShiftedPositionRaycasting(fshift, xinit[i], yinit[i], CentroidX, CentroidY);

      }
      else // compute by following the normal. No Centroid Used here.
      {
        pos1 = this.GetShiftedPosition(-fshift, xinit, yinit, i);
        pos2 = this.GetShiftedPosition(fshift, xinit, yinit, i);
      }
      // add lines along which the matching point will be found for all sampled point.
      mydraw.drawLine(pos1[0], pos1[1], pos2[0], pos2[1]);
    }
    
    // Then for each ROIs (aka each frame here)
    //       - interrupt the program if there is no ROI in this frame OR if it is not an area
    //       - get the intensity values along the sample point ?? Question why centroid?
    //		- for each frame n, every 10 sampled points: draw a line between the position 
    //		at frame n and te position at n+1
	//      for which color is proportional to the frame number
  //         - Draw the last position of sampled point in RED
    for (int r = 0; r < roilist.length; r++)
    {
      Roi roi = roilist[r];
      Roi roinext = roi;
      if (r < roilist.length - 1)
        roinext = roilist[(r + 1)];
      imp.setSlice(r + 1);
      ip = imp.getProcessor();
      // tbc: c'est pas l'inerse les conditions si nul ou si c'est pas une aire?
      //- interrupt the program if there is no ROI in this frame OR if it is not an area
      if (roi != null && !roi.isArea()) {
          IJ.error("Area selection required on each frame. Problem identified with frame"+r+1);
          return;
        
      }
      
      FloatPolygon p = roi.getInterpolatedPolygon(sampling, smooth);
      FloatPolygon pnext = roinext.getInterpolatedPolygon(sampling, smooth);
    //get the intensity values along the sample point in a band of linewidth, and keep the maximum value. Question why centroid?
      for (int shift = -linewidth; shift <= linewidth; shift++)
      { 
        double[] profile1frame1shift = this.getIrregularProfile(p, xinit, yinit, ip, shift, Ninit, CentroidX, CentroidY);
        if (shift == -linewidth)
        {
          this.profile1frame = profile1frame1shift;
        }
        else
        {
          for (int i = 0; i < Ninit; i++)
          {
            if (this.profile1frame[i] < profile1frame1shift[i]) {
              this.profile1frame[i] = profile1frame1shift[i];
            }
          }
        }
      }
     
      // concatene this into profile
      this.profile = this.concat(this.profile, this.profile1frame);
      double[] newpos = this.getnewpos(p, xinit, yinit, ip, Ninit, CentroidX, CentroidY);
      for (int i = 0; i < Ninit; i++)
      {
    	 //every 10 sampled points: draw a line between the position at frame n and te position at n+1
    	  //for which color is proportional to the frame number
        if (i % 10 == 0)
        {
          int indexedintensity2 = Math.round(r / numStacks * factor);
          Color rgb = new Color(reds[indexedintensity2], greens[indexedintensity2], blues[indexedintensity2]);
          mydraw.setColor(rgb);
          mydraw.drawLine(xinit[i], yinit[i], newpos[i], newpos[(i + Ninit)]);
        }
        //replace xinit by position at n+1
        xinit[i] = newpos[i];
        yinit[i] = newpos[(i + Ninit)];
      }
      //Draw the last position of sampled point in RED
      mydraw.setColor(Color.RED);
      if (r == roilist.length - 2) {
        mydraw.addPoints(xinit, yinit, 5);
      }
      //Compute the speed: should be between xinit and newpos, why pos next?? tbc (what is doing getnewpos excatly)
      double[] newposnext = this.getnewpos(pnext, xinit, yinit, ip, Ninit, CentroidX, CentroidY);
      double[] edgespeed = new double[Ninit];
      for (int n = 0; n < Ninit; n++) {
        double deltax = newposnext[n] - newpos[n];
        double deltay = newposnext[(n + Ninit)] - newpos[(n + Ninit)];
        
        edgespeed[n] = Math.sqrt(deltax * deltax + deltay * deltay);
        

        // sign the speed according to the fact it is going away or toward the centroid
        double distancefromcentroidx = newpos[n] - CentroidX;
        double distancefromcentroidy = newpos[(n + Ninit)] - CentroidY;
        double distanceold = Math.sqrt(distancefromcentroidx * distancefromcentroidx + distancefromcentroidy * distancefromcentroidy);
        double distancefromcentroidnextx = newposnext[n] - CentroidX;
        double distancefromcentroidnexty = newposnext[(n + Ninit)] - CentroidY;
        double distancenew = Math.sqrt(distancefromcentroidnextx * distancefromcentroidnextx + distancefromcentroidnexty * distancefromcentroidnexty);
        if (distancenew - distanceold < 0.0D) {
          edgespeed[n] = (-edgespeed[n]);
        }
      }
      
      this.profilespeed = this.concat(this.profilespeed, edgespeed);
    }
    

    // show speed image
    FloatProcessor nipspeed = new FloatProcessor(Ninit, numStacks, profilespeed);
    ImagePlus speed = new ImagePlus("Edge dynamics", nipspeed);
    // show membrane recruitment image
    FloatProcessor nip = new FloatProcessor(Ninit, numStacks, profile);
    ImagePlus kymo = new ImagePlus("Recruitment", nip);
    kymo.show();
    IJ.run("Fire");
    speed.show();
    IJ.run("Fire");
    // show our sheme @TODO add a legend
    mydraw.show();
  }
/**
 * getnewpos 
 * @param p the full list of point as sampled by ImageJ, with sampling and smoothing options at frame f
 * @TODO do not use sampling for frame>1 for a better curve tracking
 * @param xinit list of initial sampled x positions at frame f-1
 * @param yinit list of initial sampled y positions at frame f-1
 * @param ip Image?? why ??
 * @param nmax  number of initially sampled points (could be read as well from xinit size tbc
 * @param centroidX position of the refenece point X (used for ray casting only)
 * @param centroidY idem Y
 * @return the full list of new position in one array only [x1 x2.. xNmax y1...yNmax]
 */
  private double[] getnewpos(FloatPolygon p, double[] xinit, double[] yinit, ImageProcessor ip, int nmax, double centroidX, double centroidY)
  {
    double[] x = new double[nmax];
    double[] y = new double[nmax];
    for (int i = 0; i < nmax; i++)
    {
      double[] projectedpoint = this.GetProjection(p, xinit[i], yinit[i], centroidX, centroidY);
      x[i] = projectedpoint[0];
      y[i] = projectedpoint[1];
    }
    
    double[] values = new double[nmax * 2 + 1];
    for (int i = 0; i < nmax; i++) {
      values[i] = x[i];
      values[(i + nmax)] = y[i];
    }
    return values;
  }
  



/**
 * concatenate two arrays of doubles in one, A first.[a1 ...an1 b1...bn2]
 * @param A [a1 ...an1]
 * @param B [ b1...bn2]
 * @return [a1 ...an1 b1...bn2]
 */
  private double[] concat(double[] A, double[] B)
  {
    double[] C;
    


    if (A == null)
    {
       C = new double[B.length];
      System.arraycopy(B, 0, C, 0, B.length);
    }
    else
    {
      C = new double[A.length + B.length];
      System.arraycopy(A, 0, C, 0, A.length);
      System.arraycopy(B, 0, C, A.length, B.length);
    }
    return C;
  }
  
/**
 * 
 * @param p p the full list of point as sampled by ImageJ, with sampling and smoothing options at frame f
 * @param xinit list of  sampled x positions at frame f
 * @param yinit list of  sampled y positions at frame f
 * @param ip
 * @param shift
 * @param nmax
 * @param centroidx position of the reference point X (used for ray casting only)
 * @param centroidy position of the reference point Y (used for ray casting only)
 * @return an array containing the intensity value at the shifted position of shift in different sampled position for this frame
 */
  double[] getIrregularProfile(FloatPolygon p, double[] xinit, double[] yinit, ImageProcessor ip, int shift, int nmax, double centroidx, double centroidy)
  {
    double[] x = new double[nmax];
    double[] y = new double[nmax];
    double[] xshifted = new double[nmax];
    double[] yshifted = new double[nmax];
    for (int i = 0; i < nmax; i++)
    {
      double[] projectedpoint = this.GetProjection(p, xinit[i], yinit[i], centroidx, centroidy);
      x[i] = projectedpoint[0];
      y[i] = projectedpoint[1];
    }
    if (shift != 0)
    {
      double[] shiftedpoint = new double[2];
      for (int i = 0; i < nmax; i++)
      {
        if (method.compareTo("Ray Casting") == 0)
        {
          shiftedpoint = this.GetShiftedPositionRaycasting(shift, x[i], y[i], centroidx, centroidy);
        }
        else
        {
          shiftedpoint = this.GetShiftedPosition(shift, x, y, i);
        }
        xshifted[i] = shiftedpoint[0];
        yshifted[i] = shiftedpoint[1];
      }
    }
    else
    {
      for (int i = 0; i < nmax; i++)
      {
        xshifted[i] = x[i];
        yshifted[i] = y[i];
      }
    }
    
    double[] values = new double[nmax];
    for (int i = 0; i < nmax; i++)
    {
      values[i] = ip.getInterpolatedValue(xshifted[i], yshifted[i]);
    }
    
    return values;
  }
  
/**
 * getShiftedPositionRayCasting will look for a point by launching a ray from the centroid passing trough the currentx,y point,
 * and shifted by shift (used for profile band) 
 * @param shift
 * @param x
 * @param y
 * @param centroidx
 * @param centroidy
 * @return the position of a point on next frame curve
 */

  private double[] GetShiftedPositionRaycasting(int shift, double x, double y, double centroidx, double centroidy)
  {
    double cx = x - centroidx;
    double cy = y - centroidy;
    double[] projectedpoint = new double[2];
    
    double normc = Math.sqrt(cx * cx + cy * cy);
    cx /= normc;
    cy /= normc;
    projectedpoint[0] = x;
    projectedpoint[1] = y;
    
    projectedpoint[0] += shift * cx;
    projectedpoint[1] += shift * cy;
    return projectedpoint;
  }
  

  /**
   * getShiftedPosition will look for a point by launching a ray from the centroid passing trough the currentx,y point,
   * and shifted by shift (used for profile band) 
   * @param shift
   * @param x
   * @param y
   * @param i index of th epoint process in the table
   * @return the position of a point on next frame curve
   */
  private double[] GetShiftedPosition(int shift, double[] x, double[] y, int i)
  {
	  double racine2;
      double xmiddle;
      double a;
      double racine1;
      double ymiddle;
      double[] pos = new double[2];
      if (i > 0) {
          if (i < x.length - 1) {
              a = (x[i - 1] - x[i + 1]) / (y[i - 1] - y[i + 1]);
              xmiddle = (x[i + 1] - x[i - 1]) / 2.0 + x[i - 1];
              ymiddle = (y[i + 1] - y[i - 1]) / 2.0 + y[i - 1];
          } else {
              a = (x[i - 1] - x[i]) / (y[i - 1] - y[i]);
              xmiddle = (x[i] - x[i - 1]) / 2.0 + x[i - 1];
              ymiddle = (y[i] - y[i - 1]) / 2.0 + y[i - 1];
          }
      } else { //i==0
          a = (x[i] - x[i + 1]) / (y[i] - y[i + 1]);
          xmiddle = (x[i + 1] - x[i]) / 2.0 + x[i];
          ymiddle = (y[i + 1] - y[i]) / 2.0 + y[i];
      }
      a = -1.0 / a;
      double b = ymiddle - a * xmiddle;
      double newa = a * a + 1.0;
      double newb = 2.0 * a * b - 2.0 * y[i] * a - 2.0 * x[i];
      double newc = b * b - 2.0 * y[i] * b + y[i] * y[i] + x[i] * x[i] - (double)shift * (double)shift;
      double discriminant = newb * newb - 4.0 * newa * newc;
      if (discriminant > 0.0) {
          double bprime = newb / 2.0;
          double deltaprime = bprime * bprime - newa * newc;
          double q = - bprime + Math.signum(bprime) * Math.sqrt(deltaprime);
          racine1 = q / newa;
          racine2 = newc / q;
      } else if (discriminant == 0.0) {
          if (this.verbose) {
              IJ.log((String)("The discriminant was 0 for point " + IJ.d2s((double)i)));
          }
          racine2 = racine1 = (- newb) / (2.0 * newa);
      } else {
          if (this.verbose) {
              IJ.log((String)("The discriminant was negative (no solution) for point " + IJ.d2s((double)i)));
          }
          racine1 = x[i];
          racine2 = x[i];
      }
      pos[0] = shift <= 0 ? racine1 : racine2;
      pos[1] = a * pos[0] + b;
      return pos;
  }

/**
 * GetProjection only ray casting???
 * @param p the new area polygon at frame f
 * @param x position on f-1
 * @param y position on f-1 (on previous curve)
 * @param centroidx
 * @param centroidy
 * @return the projection of the point x y 
 */

  private double[] GetProjection(FloatPolygon p, double x, double y, double centroidx, double centroidy)
  {
    double cx = x - centroidx;
    double cy = y - centroidy;
    double[] projectedpoint = new double[2];
    
    double normc = Math.sqrt(cx * cx + cy * cy);
    cx /= normc;
    cy /= normc;
    projectedpoint[0] = x;
    projectedpoint[1] = y;
    int count = 0;
    if (p.contains((float)x, (float)y)) {
      while (p.contains((float)projectedpoint[0], (float)projectedpoint[1]))
      {
        projectedpoint[0] += 0.1D * cx;
        projectedpoint[1] += 0.1D * cy;
      }
    } else {
      while (!p.contains((float)projectedpoint[0], (float)projectedpoint[1]))
      {
        projectedpoint[0] -= 0.1D * cx;
        projectedpoint[1] -= 0.1D * cy;
        count++;
        if (count > 100){
        	if (this.verbose){
            	IJ.log("did not find the contour!!");
            }
          return projectedpoint;
          }
      }
    }
    return projectedpoint;
  }
  

/**
 * just an help message
 */
  void showAbout()
  {
    IJ.showMessage("About Recruitment Edge Dynamic...", 
      "This plugin aims at studying the boundary of a structure evolving over time, for exemple a protrusion by measuring both its speed and intensity. Authors: Perrine Paul-Gilloteaux, Maria Carla Parrini");
  }
  





/**
 * 
 * @param reds
 * @param greens
 * @param blues
 * @return the LUT Fire in 3 arrays reds, greens and blues, and return the length of it for indication for any rescale.
 */

  int ice(float[] reds, float[] greens, float[] blues)
  {
    int[] r = { 0, 0, 0, 0, 0, 0, 19, 29, 50, 48, 79, 112, 134, 158, 186, 201, 217, 229, 242, 250, 250, 250, 250, 251, 250, 250, 250, 250, 251, 251, 243, 230 };
    int[] g = { 156, 165, 176, 184, 190, 196, 193, 184, 171, 162, 146, 125, 107, 93, 81, 87, 92, 97, 95, 93, 93, 90, 85, 69, 64, 54, 47, 35, 19, 0, 4, 0 };
    int[] b = { 140, 147, 158, 166, 170, 176, 209, 220, 234, 225, 236, 246, 250, 251, 250, 250, 245, 230, 230, 222, 202, 180, 163, 142, 123, 114, 106, 94, 84, 64, 26, 27 };
    for (int i = 0; i < r.length; i++)
    {
      reds[i] = (r[i] / 255.0F);
      greens[i] = (g[i] / 255.0F);
      blues[i] = (b[i] / 255.0F);
    }
    return r.length;
  }
}
