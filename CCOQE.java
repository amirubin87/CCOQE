package CCOQE;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class CCOQE {

	public final static Integer C = 65650;
	
	public static void main(String[] args) throws IOException {
		//      Read input
        /*String GraphFile = args[0];//"C:/cygwin64/home/t-amirub/binary_networks/network.dat"; // 
        String OutputDir =  args[1];//"C:/cygwin64/home/t-amirub/binary_networks/res/";// 
        double[] bettas = parseBettas(args[2]);//{1.1};// 
        double alpha = Double.parseDouble(args[3]);//0.5;//*/
        
		String GraphFile = "C:/Users/t-amirub/Desktop/net/Les Miserables/net.txt";
		String OutputDir = "C:/Users/t-amirub/Desktop/net/Les Miserables/CCOQE/";
		double[] bettas = parseBettas("1.1");//,1.04,1.06,1.08,1.1,1.2,1.3,1.4        
        double alpha = Double.parseDouble("0.8");
        
        int iteratioNumToStartMerge = Integer.parseInt("0");
        int maxIterationsToRun = Integer.parseInt("30");
                        
        // read graph file
        int N =0;
        int E = 0;
        
        Map<Integer, List<Integer>> A= new HashMap<Integer, List<Integer>> ();        
        Map<Integer,Integer> ranks = new HashMap<Integer,Integer>();
        Map<Integer, List<Integer>> node2comms= new HashMap<Integer, List<Integer>> ();
        Map<Integer, List<Integer>> comm2nodes= new HashMap<Integer, List<Integer>> ();
        //TODO- use this for perf issues
        Map<Integer, Map<Integer, Double>> node2Comm2Qe= new HashMap<Integer, Map<Integer, Double>>();
        //Map<Integer, Map<Integer,Integer>> Intersection_c1_c2 = new HashMap<Integer, Map<Integer,Integer>> ();
        
        	System.out.println("Reading edges");
        	
		    List<String> edges = Files.readAllLines(Paths.get(GraphFile), Charset.defaultCharset());
		    System.out.println("Done Read edges");
		    E = edges.size();
		    System.out.println("Iterating edges");
		    int commId = 0;
		    for(String edge: edges){
		        String[] nodes = edge.split("\\s");
		        Integer from = Integer.parseInt(nodes[0]);
		        Integer to = Integer.parseInt(nodes[1]);		
		        if(from == to){
		        	E--;    
		        	continue;		                
		        }
		        
		     // update A		        
		        if(! UpdateA(A, from, to)){
		        	E--;
		        	continue;
		        }
		        UpdateA(A,to, from);
		        
		        //update data struture
		        commId = AddNodeToDataStructure(node2comms, comm2nodes, from, commId);
		        commId = AddNodeToDataStructure(node2comms, comm2nodes, to, commId); 		        
		        
		        // update ranks
		        ranks.put(from, ranks.getOrDefault(from, 0) + 1);
		        ranks.put(to, ranks.getOrDefault(to, 0) + 1);
		    }

        N = ranks.keySet().size();
        System.out.println("Done processing graph");
        System.out.println("N = " + N +". E = " + E);        
        

        Map<Integer, List<Integer>> Tnode2comms;
        Map<Integer, List<Integer>> Tcomm2nodes;
        Map<Integer, Map<Integer, Double>> Tnode2Comm2Qe;
        Collection<List<Integer>> comms;
        for (double betta:bettas){
        	System.out.println("pathToGraph:             "+GraphFile);      
			System.out.println("outputPath:              "+OutputDir);			
			System.out.println("alpha:                   "+alpha);
			System.out.println("iteratioNumToStartMerge: "+iteratioNumToStartMerge);
			System.out.println("maxIterationsToRun:      "+maxIterationsToRun);
        	 System.out.println("                Betta  " + betta);
        	 
        	 Tnode2comms = DeepCopy(node2comms);
        	 Tcomm2nodes = DeepCopy(comm2nodes);
        	 comms = FindCommunities(Tnode2comms, Tcomm2nodes, betta, alpha, A, ranks, N, E, iteratioNumToStartMerge, maxIterationsToRun);
        	 WriteCommsToFile(comms, OutputDir, betta);
        	 //System.out.println(comms);
        }
    }
 
	private static Collection<List<Integer>> FindCommunities(
			Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes, double betta,
			double alpha, Map<Integer, List<Integer>> a,
			Map<Integer, Integer> ranks, int N, int E, int iteratioNumToStartMerge, int maxIterationsToRun) {
    	
    	int numOfIterations = 0;
    	int amountOfDone = 0;
    	boolean haveMergedComms = false; 	
    	    	
		int nonEmptyComms = comm2nodes.keySet().size();
    	while(numOfIterations<2 ||( amountOfDone < nonEmptyComms && numOfIterations < maxIterationsToRun)){
    		
    		System.out.println(""); 
    		System.out.println("numOfIterations: " + numOfIterations);
    		System.out.println("amount of communities: "+ nonEmptyComms);
    		System.out.print("progress in current iteration: ");
    		nonEmptyComms = 0;
    		amountOfDone = 0;
    		numOfIterations++;
    		// Go over all comms
    		List<Integer> comms = new ArrayList(comm2nodes.keySet());
    		int n = comms.size();
    	    int tenPercent = n/10+1;
    	    int commCounter=0;
    		for(Integer comm : comms){
    			commCounter++;
    		
    			if(comm2nodes.get(comm)== null || comm2nodes.get(comm).size()<1){
    				continue;
    			}
    			nonEmptyComms++;
    			if ((commCounter%tenPercent) == 0){
	        		System.out.print(commCounter/tenPercent*10 + "%  ");
	        	}
    			
        		Map<Integer, Double> nodes_inc = FindNewNodes(node2comms, comm2nodes, a, ranks, E, comm, betta);
        		if(nodes_inc.size()==0)
        			continue;
        		double maxInc = Collections.max(nodes_inc.values());
        		List<Integer> newNodes =FilterBestNodes(nodes_inc, betta);
        		boolean nodesHaveChanged = CheckIfOtherNodesAreToBeRemovedFromComm(newNodes, comm, node2comms, comm2nodes, a, ranks, E, betta, maxInc);
        		
    			Map<Integer,Map<Integer, Double>> commsCouplesIntersectionRatio = AddNodesToComm(node2comms, comm2nodes, comm, newNodes, alpha);    			
    			if(numOfIterations>iteratioNumToStartMerge){
    				haveMergedComms = FindAndMergeComms(node2comms, comm2nodes, commsCouplesIntersectionRatio);   			
    			}    			
	            
    			if (!haveMergedComms && !nodesHaveChanged ){
    				amountOfDone++;
    			}
    		}
    	}
    	System.out.println("numOfIterations: " + numOfIterations);
		return comm2nodes.values();
	}
	
	private static Map<Integer,Map<Integer, Double>> AddNodesToComm(
			Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes, 
			Integer comm,
			List<Integer> newNodes, double alpha) {
		Map<Integer,Map<Integer, Double>> ans = new HashMap<>();
		List<Integer> commNodes = comm2nodes.get(comm);
		for (Integer node : newNodes){
			//add to community
			commNodes.add(node);
			// add to nodes list
			node2comms.get(node).add(comm);
						
			for(Integer otherComm : node2comms.get(node)){
				if(otherComm!=comm){
					Integer highComm = comm;
					Integer lowComm = otherComm;
					if(highComm<lowComm){
						highComm = otherComm;
						lowComm = comm;
					}
				    List<Integer> comm1 = commNodes;
				    List<Integer> comm2 = comm2nodes.get(otherComm);				
		            int intersectionSize = intersectionSize(comm1,comm2);
		            if (intersectionSize < 3){
		            	intersectionSize= 0;
		            }	            
		            
		            double intersectionRatio = (double)intersectionSize/(double)Math.min(comm1.size(), comm2.size());
		            if (intersectionRatio > alpha){
			            Map<Integer,Double> lowMap = ans.getOrDefault(lowComm, new HashMap<Integer,Double>());
			            lowMap.put(highComm,intersectionRatio);
			            ans.put(lowComm, lowMap);
		            }
				}
			}
		}	
		
		return ans;
	}
	
	private static boolean CheckIfOtherNodesAreToBeRemovedFromComm(List<Integer> newNodes, Integer comm,
			Map<Integer, List<Integer>> node2comms, Map<Integer, List<Integer>> comm2nodes,
			Map<Integer, List<Integer>> a, Map<Integer, Integer> ranks, int e, double betta, double maxInc) {
		boolean ans = false;
		List<Integer> nodesToRemove = new ArrayList<Integer>();
		if(comm2nodes.get(comm).size()<3){
			return false;
		}
		//go over all nodes in comm
		for(Integer node :comm2nodes.get(comm)){
			//only check for other nodes
			if(!newNodes.contains(node)){
				double qe= CalcQeImprovement(node2comms, comm2nodes, a, ranks, e, comm, node);
				if (qe*betta<maxInc){
					ans = true;
					nodesToRemove.add(node);					
				}
			}
		}
		for(Integer node : nodesToRemove){
			RemoveNodeFromComm(node2comms, comm2nodes, node, comm);
		}
		return ans;
	}

	private static void RemoveNodeFromComm(
			Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes, Integer node, Integer comm) {
		//Remove from comm	
		comm2nodes.get(comm).remove((Integer)node);			
		//Remove from nodes list of comms
		node2comms.get(node).remove((Integer)comm);	
						
	}
    
	 private static Map<Integer, List<Integer>> DeepCopy(
				Map<Integer, List<Integer>> map) {
	 	Map<Integer, List<Integer>> ans = new HashMap<Integer, List<Integer>>();
	 	for (Integer k : map.keySet()){
	 		List<Integer> clone = cloneList(map.get(k));
	 		ans.put(k, clone);
	 	}
			return ans;
		}
	 
	 public static List<Integer> cloneList(List<Integer> list) {
	     List<Integer> clone = new ArrayList<Integer>(list.size());
	     for(Integer item: list) clone.add(item);
	     return clone;
	 }
 
	private static Map<Integer,Double> FindNewNodes(Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes, Map<Integer, List<Integer>> a, Map<Integer, Integer> ranks, int e,
			int comm, double betta) {
		Map<Integer,Double> nodes_inc = new HashMap<Integer,Double>(); 
		// Go over all border nodes
		List<Integer> borderNodes = FindBorderNodes(comm,node2comms, comm2nodes, a);
        if (borderNodes.size() == 0){
        	return new HashMap<Integer,Double>();
        }
        for (int borderNode : borderNodes){
        	double inc= CalcQeImprovement(node2comms, comm2nodes, a, ranks , e, comm, borderNode);            
        	nodes_inc.put(borderNode, inc);
        }		
        return nodes_inc;
	}
	
	private static boolean FindAndMergeComms(
			Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes,
			Map<Integer, Map<Integer, Double>> commsCouplesIntersectionRatio) {
		boolean haveMergedComms = false;
		Map<Integer, Double> lowMap;
		for (Integer lowK : commsCouplesIntersectionRatio.keySet()){
			lowMap = commsCouplesIntersectionRatio.get(lowK);
			for (Integer highK : lowMap.keySet()){				
				haveMergedComms = MergeComms(node2comms, comm2nodes, lowK, highK) || haveMergedComms;				
			}
		}
		return haveMergedComms;
	}
	
	private static Boolean MergeComms(Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes, Integer lowK, Integer highK) {
		List<Integer> commHigh = comm2nodes.getOrDefault(highK, new ArrayList<Integer>());
		List<Integer> commsLow = comm2nodes.getOrDefault(lowK, new ArrayList<Integer>());
		if(commHigh.size()==0 || commsLow.size()==0){
			return false;
		}
		for (Integer node : commHigh){
			List<Integer> nodeComm = node2comms.get(node);
			
			// remove comm from nodes comm
			nodeComm.remove((Integer)highK);
			
			//adds low comm to node
			if (!nodeComm.contains((Integer)lowK)){
				nodeComm.add((Integer)lowK);
			}
			
			// adds node to low comm
			if (!commsLow.contains((Integer)node)){
				commsLow.add((Integer)node);
			}
			
		}
		// removes all nodes from comm
		comm2nodes.remove((Integer)highK);
		return true;
		
	}

	public static <T> int intersectionSize(List<T> list1, List<T> list2) {
        int ans = 0;

        for (T t : list1) {
            if(list2.contains(t)) {
                ans++;
            }
        }

        return ans;
    }

	private static List<Integer> FilterBestNodes(
			Map<Integer, Double> nodes_inc, double betta) {
		List<Integer> ans = new ArrayList<Integer>(); 
		double max = Collections.max(nodes_inc.values());
		if (max <=0){
			return ans;
		}
		
		for (Integer k : nodes_inc.keySet()){
			if(nodes_inc.get(k)*betta > max){
				ans.add(k);
			}
		}
		return ans;
	}

	private static double CalcQeImprovement(
			Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes,
			Map<Integer, List<Integer>> a,
			Map<Integer, Integer> ranks, int e, Integer neighborComm, Integer node) {
		double Qe = 0;
		//int Onode = node2comms.get(node).size() + 1;
		int Knode = ranks.getOrDefault(node, 0);
		List<Integer> neis = a.get(node);
		int Anodei = 0;
		int Ki = 0;
		int Oi = 0;
		for (Integer i : comm2nodes.getOrDefault(neighborComm, new ArrayList<Integer>())){
			Anodei = 0;
			if (neis.contains(i)){
				Anodei = 1;
			}
			Ki = ranks.getOrDefault(i, 0);
			Oi = node2comms.get(i).size();
			Qe = Qe + (
						(
							Anodei - 
								(
									(Knode*Ki)
									/((double)2*e)
								)
						)
						/(double)Oi
					   );
			
		}
		return Qe ;
	}

	private static List<Integer> FindBorderNodes(
			Integer comm, Map<Integer, List<Integer>> node2comms, 
			Map<Integer, List<Integer>> comm2nodes,
			Map<Integer, List<Integer>> a
			) {
		Set<Integer> ans = new HashSet<>();
		List<Integer> nodes = comm2nodes.get(comm);
		for(Integer node:nodes){
			List<Integer> neis = a.get(node);
			for (int nei : neis){
				if(!nodes.contains(nei)){
					ans.add(nei);
				}
			}
		}
		List<Integer> ret = new ArrayList<>();
		ret.addAll(ans);
		return ret;
	}

	private static void WriteCommsToFile(Collection<List<Integer>> comms,
			String outputDir, double betta) throws FileNotFoundException, UnsupportedEncodingException {
		
		PrintWriter writer = new PrintWriter(outputDir + "Betta-" + betta + ".txt", "UTF-8");
		for (List<Integer> comm : comms){
			if(comm.size()<3)
				continue;
			String line ="";
			for (Integer node : comm)
			{
				line += node + " ";
			}
			if (line.length() > 0)
				writer.println(line);
		}
		writer.close();		
	}

	private static boolean UpdateA(Map<Integer, List<Integer>> A, Integer from, Integer to) {
		List<Integer> oldNei = A.get(from);
        if(oldNei == null){
        	oldNei = new ArrayList<Integer>();
        	A.put(from, oldNei);
        }
        
        else if(oldNei.contains(to)){        	
        	return false;        	
        }
        
        oldNei.add(to);
        
        return true;
	}
	
	private static int AddNodeToDataStructure(
			Map<Integer, List<Integer>> node2comms,
			Map<Integer, List<Integer>> comm2nodes,
			Integer node, Integer commId) {
			List<Integer> oldComms = node2comms.get(node);
			//We only do this if no comm has been set to this node
	        if(oldComms == null){
	        	oldComms = new ArrayList<Integer>();
	        	oldComms.add(commId);
	        	node2comms.put(node, oldComms);
	        	List<Integer> oldnodes = new ArrayList<Integer>();
	        	oldnodes.add(node);
	        	comm2nodes.put(commId, oldnodes);	        	
		        commId++;
	        }
	        return commId;		
	}

	private static double[] parseBettas(String iBettas) {
		String[] bettas = iBettas.split(",");
		double[] ans = new double[bettas.length];
		int i = 0;
		for(String betta : bettas ){
			ans[i] = Double.parseDouble(betta);
			i++;
		}
		return ans;
	}
}

