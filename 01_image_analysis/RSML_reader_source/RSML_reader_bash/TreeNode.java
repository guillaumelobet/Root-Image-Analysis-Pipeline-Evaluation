
import ij.IJ;

import java.util.ArrayList;



class TreeNode {
	
	TreeNode parent;
	public ArrayList<TreeNode> childList = new ArrayList<TreeNode>();
	float length;
	int order;
	int id;
	int magnitude;
	static int altitude;
	static int childNumber;
	static int treeLength;
	int pathLength;
	
	
	public TreeNode(TreeNode p, float l, int o, int i){
		parent = p;
		length = l;
		order = o;
		id = i;
		magnitude = 1;
		pathLength = 0;
		if(parent != null) pathLength = parent.pathLength++;
	}
	
	public void addChild(TreeNode c){
		childList.add(c);
	}
	
	/**
	 * Get the total length of the network 
	 * @param t
	 */
	public float getLength(boolean start) {
		if(start) treeLength = 0;
		if (childList != null) {
			for (TreeNode resFils : childList) {
				resFils.getLength(false);
			}
		}
		if(childList != null) treeLength += length;
		return treeLength;
	}
	
	
	
	/**
	 * Get the number of children of the network
	 * @return the number of children
	 */
	public int getNumberOfChildren(boolean start){
		if(start) childNumber = 0;
		if (childList != null) {
			for (TreeNode resFils : childList) {
				resFils.getNumberOfChildren(false);
			}
		}
		if(childList != null) childNumber += childList.size();
		return childNumber;
	}
	
	
	/**
	 * Get the topological index of the network (Fitter 87)
	 * TODO
	 */
	public int getTopologicalIndex(boolean start){
		return getMagnitude() / getAltitude(true);
	}	
	/**
	 * Get the magnitude of the network (Fitter 87)
	 * @return the number of children
	 * TODO
	 */
	public void setMagnitude(){
		
		// Navigate to the leaves
		if (childList != null) {
			for (TreeNode resFils : childList) {
				resFils.setMagnitude();
			}
		}

		// Compute the magnitude
		if (childList.size() > 0) {
			magnitude = 0;
			for (TreeNode resFils : childList) {
				magnitude += resFils.magnitude;
			}
		}
	}
	
	public int getMagnitude(){
		return magnitude;
	}
	
	
	/**
	 * Get the altitude of the network (Fitter 87)
	 * @return the number of children
	 * TODO
	 */
	public int getAltitude(boolean start){
		if(start) altitude = 0;
		if (childList != null) {
			for (TreeNode resFils : childList) {
				resFils.getAltitude(false);
			}
		}
		altitude = Math.max(pathLength, altitude);
		return altitude;
	}	
	
	/**
	 * Get the external path length of the network (Fitter 87)
	 * @return the number of children
	 * TODO
	 */
	public int getExternalPathLength(boolean start){
		if(start) childNumber = 0;
		if (childList != null) {
			for (TreeNode resFils : childList) {
				resFils.getNumberOfChildren(false);
			}
		}
		if(childList != null) childNumber += childList.size();
		return childNumber;
	}
	

}
