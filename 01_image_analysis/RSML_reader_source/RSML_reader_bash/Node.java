

import java.awt.geom.*;

/**
 * @author Xavier Draye - Universit� catholique de Louvain
 * @author Guillaume Lobet - Universit� de Li�ge 
 * 
 * Node class.
 * 
 */


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/**
 * Constructor
 * @author guillaumelobet
 *
 */
class Node {
   float x, y, z, theta, length, cLength, diameter, birthTime;
   
   // length and cLength are in pixels
   Node child;
   Node parent;
   boolean needsRefresh;
   boolean bCross01 = false, bCross23 = false;
   boolean pCross01 = false, pCross23 = false;   


  /**
   * Constructor
   * @param x
   * @param y
   * @param d
   * @param n
   * @param after
   */
   public Node(float x, float y, float z, float d, Node n, boolean after) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.diameter = d;
      if (after) {
         parent = n;
         if (parent != null) parent.child = this;
         child = null;
         }
      else {
         child = n;
         if (child != null) child.parent = this;
         parent = null;
         }
      needsRefresh = true;      
      }

   /**
    * 
    * @param x
    * @param y
    * @param n
    * @param after
    */
   public Node(float x, float y, float z, Node n, boolean after) {
      this(x, y, z, 0f, n, after);
      }

   /**
    * Build the node
    */
   public void buildNode() {
      if (diameter < 1.0f) diameter = 1.0f;

      // calculate length and theta where required
      if (parent != null) {
         float dx = x - parent.x;
         float dy = y - parent.y;
         float dz = z - parent.z;
         parent.theta = vectToTheta(dx, dy, dz);
         parent.length = norm(dx, dy, dz);
         }
      if (child != null) {
         float dx = child.x - x;
         float dy = child.y - y;
         float dz = child.z - z;
         theta = vectToTheta(dx, dy, dz);
         length = norm(dx, dy, dz);
         }
      }

   /**
    * 
    * @param dx
    * @param dy
    * @return
    */
   public static float norm(float dx, float dy, float dz) {
	      return (float) Math.sqrt(dx * dx + dy * dy + dz * dz); 
	      }
  
   /**
    * 
    * @param n
    */
   public void copy(Node n) {
      x = n.x;
      y = n.y;
      theta = n.theta;
      length = n.length;
      cLength = n.cLength;
      diameter = n.diameter;
      if (parent != null) parent.needsRefresh = true;
      if (child != null) child.needsRefresh = true;
      needsRefresh = true;
      }

   /**
    * Move the node of a given x and y deviation
    * @param dx
    * @param dy
    */
   public void translate(float dx, float dy) {
      x += dx;
      y += dy;    
      needsRefresh = true;
      }

 

   /**
    * Get the distance between this node and an other one
    * @param n
    * @return
    */
   public float getDistanceTo (Node n) {
      float d = 0.0f;
      for (Node n1 = this; n1 != n; n1 = n1.child) { 
         if (n1 == null) return 0.0f;
         d += (float) Point2D.distance((double) n1.x, (double) n1.y,
                                       (double) n1.child.x, (double) n1.child.y);
         }
      return d;
      }

   /**
    * Compute the length between the base of the root and the node
    */
   public void calcCLength() {
      calcCLength(this.cLength);
      }

	/**
	 * Compute the length between the base of the root and the node
	 * @param startValue
	 */
   public void calcCLength(float startValue) {
      this.cLength = startValue;
      Node n = this;
      while (n.child != null) {
         n.child.cLength = n.cLength + n.length;
         n = n.child;
         }
      }
   
   /**
    * Read the node information from and RSML file
    * @param parentDOM the xml elemnt containg the x/y coordinates
    * @param diamDOM the xml element contining the diameter elements
    * @param dpi
    */
   public void readRSML(org.w3c.dom.Node parentDOM, org.w3c.dom.Node diamDOM, float dpi) {
	   
		  org.w3c.dom.Node nn = parentDOM.getAttributes().getNamedItem("x");
		  if (nn != null) x = Float.valueOf(nn.getNodeValue()).floatValue() * dpi;
		  nn = parentDOM.getAttributes().getNamedItem("y");
		  if (nn != null) y = Float.valueOf(nn.getNodeValue()).floatValue() * dpi;
		  nn = parentDOM.getAttributes().getNamedItem("z");
		  if (nn != null) z = Float.valueOf(nn.getNodeValue()).floatValue() * dpi;
		  if(diamDOM != null){
			  diameter = Float.valueOf(diamDOM.getFirstChild().getNodeValue()).floatValue() * dpi;
		  }
		  else{
			  diameter = 2;
		  }
		  
		  
	      if (parent != null) {
	          float dx = x - parent.x;
	          float dy = y - parent.y;
	          float dz = z - parent.z;
	          parent.theta = vectToTheta(dx, dy, dz);
	          parent.length = (float) Math.sqrt(dx * dx + dy * dy);
//	          parent.calcBorders();
	          //parent.calcPoles(0);
	         }
	      needsRefresh = true;
	      }
   
   
   /**
    * Convert a vector to an ange
    * @param dirX
    * @param dirY
    * @return
    */
   public static float vectToTheta (float dirX, float dirY, float dirZ) {
	      float norm = (float) Math.sqrt(dirX * dirX + dirY * dirY + dirZ * dirZ);
	      return (float) (dirY <= 0 ? Math.acos(dirX / norm) : 2.0 * Math.PI - Math.acos(dirX / norm));      

	 }
}


