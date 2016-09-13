package hotnet2;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

import org.junit.Test;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

public class RandomGraphGenerator{
	
	public RandomGraphGenerator(){
		
	}
	
	@Test
	public void Test() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";

/*		//Create the largest component using the whole ReactomeFI network graph
		GraphUtils gu = new GraphUtils(); 
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraphWrapper(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraphWrapper(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		Set<String> geneSet = gu.getGeneGraphSetWrapper(largestComponent);
*/
		testRandomNetworkPermutationSwitch(directory, 6000, 5);
	}
	
	/**
	 * Creates a random connected network using the number of provided nodes.
	 * @param numNodes - Number of nodes that should be in graph.
	 * @return a connected graph with the provided number of nodes. 
	 */
	private Graph<String, String> createRandomConnectedGraphTest (int numNodes){
		//Create and add nodes to graph and obtain a list of the degree each node should have
		Graph<String, String> graph = new UndirectedSparseGraph<String, String>();
		List<Integer> randomDegreeList = new ArrayList<Integer>(); 
		for (int i=0; i<numNodes; i++){
			graph.addVertex(""+i);
			Random randomGen = new Random();
			int random = randomGen.nextInt(numNodes-1)+1;
			randomDegreeList.add(random);
		}
		
		//Create and add random edges to the graph.  Self-edges are not permitted.  
		List<String> nodeList = new ArrayList<String>(graph.getVertices());
		List<String> nodeCheckList = new ArrayList<String>(nodeList);
		for (int i=0; i<nodeList.size(); i++){
			for (int j=0; j<randomDegreeList.get(j); j++){
				String source = nodeList.get(j);
				Random randomGen = new Random();
				int random = randomGen.nextInt(numNodes-1);
				String target = nodeList.get(random);
				if (source != target){
					graph.addEdge(source + "\t" + target, source, target);
					nodeCheckList.remove(source);
					nodeCheckList.remove(target);
				} else {
					if (random+1 < numNodes){
						String newTarget = nodeList.get(random+1);
						graph.addEdge(source + "\t" + newTarget, source, newTarget);
						nodeCheckList.remove(source);
						nodeCheckList.remove(target);
					} else {
						String newTarget = nodeList.get(random-1);
						graph.addEdge(source + "\t" + newTarget, source, newTarget);
						nodeCheckList.remove(source);
						nodeCheckList.remove(target);
					}
				}
			}
		}
		//Nodes that weren't included in network will be added to the 1st node to ensure only 1 component exists
		if (!nodeCheckList.isEmpty()){
			String source = nodeList.get(0);
			for (int i=0; i<nodeCheckList.size(); i++){
				String target = nodeCheckList.get(i);
				graph.addEdge(source + "\t" + target, source, target);
			}
			
		}

		//Checks network is connected 
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		if (components.size() != 1){
			System.err.println("Provided graph is not connected");
			System.err.println("components: " + components);
			System.exit(0);
		}
		return graph; 
	}

	/**
	 * Checks that the connected permuted graph has the same genes and same degrees as the provided graph.
	 * @param directory - Directory containing permuted graph's edge list.
	 * @param permuteEdgeListFile - File with permuted graph's edge list.
	 * @param graph - Original graph used to create permuted graph.
	 * @return true if graph has the same gene and degrees, else false.
	 * @throws IOException
	 */
	private Boolean checkPermutation(String directory, String permuteEdgeListFile, Graph<String, String> graph) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		//Create permuted graph from file and check there is only 1 component
		Path permuteFilePath = Paths.get(directory, permuteEdgeListFile);
		Set<String> allGenesPermuted = fu.getAllGenesWrapper(permuteFilePath);
		Set<String> allPairsPermuted = fu.getAllInteractionPairsWrapper(permuteFilePath);
		Graph<String, String> permutedGraph = gu.createGraphWrapper(allGenesPermuted, allPairsPermuted);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(permutedGraph);
		if (components.size() != 1){
			throw new IllegalArgumentException("Provided permuted graph " + permuteEdgeListFile + " is not connected, it has " + components + " components.");
		}
		
		//Compare if the 2 graphs have the same genes and the gene's have the same degree
		Set<String> permutedGeneSet = gu.getGeneGraphSetWrapper(permutedGraph);
		List<Integer> permutedGraphDegree = gu.getGeneDegreeListWrapper(permutedGraph, permutedGeneSet);
		Set<String> originalGeneSet = gu.getGeneGraphSetWrapper(graph);
		List<Integer> originalGraphDegree = gu.getGeneDegreeListWrapper(graph, originalGeneSet);
		if (originalGeneSet.containsAll(permutedGeneSet) && permutedGeneSet.containsAll(originalGeneSet)){
			if (originalGraphDegree.containsAll(permutedGraphDegree) && permutedGraphDegree.containsAll(originalGraphDegree))
				return true;
			else{
				throw new IllegalArgumentException("Provided permuted graph " + permuteEdgeListFile + " does not have the same degrees as the provided graph");
/*				System.err.println("Permuted graph '" + permuteEdgeListFile + " does not have the same degrees as the provided graph");
				System.err.println("original Degrees: " + originalGraphDegree);
				System.err.println("permuted Degrees: " + permutedGraphDegree);
				return false;
*/				
			}
		} else {
			throw new IllegalArgumentException("Provided permuted graph " + permuteEdgeListFile + " does not have the same degrees as the provided graph");
/*			System.err.println("Permuted graph '" + permuteEdgeListFile + "' does not have the same genes as the provided graph");
			System.err.println("original GeneSet: " + originalGeneSet);
			System.err.println("permuted GeneSet: " + permutedGeneSet);
			return false;
*/			
		}
	}
	
	/**
	 * Creates a random network and use it to generate random permutations using the switch algorithm.  All graphs are saved.
	 * <p>
	 * <b>Note:</b>  this is for testing purposes to try the switch algorithm on different networks.
	 * @param directory - Directory to save files in.
	 * @param numNodes - Number of nodes a graph should have.
	 * @throws IOException
	 */
	private void testRandomNetworkPermutationSwitch(String directory, int numNodes, int numPermutations) throws IOException{
		GraphUtils gu = new GraphUtils();
		FileUtils fu = new FileUtils();
		Graph<String, String> graph = createRandomConnectedGraphTest(numNodes);
		Set<String> pairs = gu.getPairsWrapper(graph);
		fu.saveSetToFile(directory+"/output/permutations/", "randomGraphForPermutation", pairs);
		
		System.out.println("Graph Size: "+ "\t nodes: " + graph.getVertexCount() + " edges: " + graph.getEdgeCount());
		
/*		System.out.println("\nnetwork config");
		String permuteEdgeListFile2 = "edgeListPermuted";
		permuteNetworkConfig(directory, permuteEdgeListFile2, numPermutations, graph);
*/		
		String permuteEdgeListFile = "edgeListPermuted";
		System.out.println("\n-------------------------------------Switch original-----------------------------------------");
//		List<Integer> numSwaps = permuteNetworkSwitch(directory, permuteEdgeListFile, numPermutations, graph);
		List<Integer> numSwapsRedo = permuteNetworkSwitchRedo(directory, permuteEdgeListFile, numPermutations, graph);
		System.out.println("\n-------------------------------------Switch Redo-----------------------------------------");
//		List<Integer> numSwapsRedo = permuteNetworkSwitchRedo(directory, permuteEdgeListFile, numPermutations, graph);
		List<Integer> numSwaps = permuteNetworkSwitch(directory, permuteEdgeListFile, numPermutations, graph);
	}

	/**
	 * Permutes network n times and save each permutation in a file.\
	 * <p>
	 * Based on Python HotNet2's permuteNetwork.py permute_network().
	 * @param directory - Directory to store files. 
	 * @param fileName - File name to store permutated graphs' edge lists in.
	 * @param numPermutations - Number of permutated graphs to create. 
	 * @param graph - Graph used for permutations.
	 * @return a list containing the number of swaps performed for each permutation. 
	 * @throws IOException
	 */
	private List<Integer> permuteNetworkSwitch(String directory, String fileName, int numPermutations, Graph<String, String> graph) throws IOException{
		System.out.println("\t---------- Creating permuted Networks ----------\t original");
		
		GraphUtils gu = new GraphUtils();
		
//		int edgeSwapConstant = 10;//115;
//		int numEdges = graph.getEdgeCount();
//		int numSwap = edgeSwapConstant * numEdges;
		int numSwap = 1; 
		int windowThreshold = 3;
		
		HashMap<String, Integer> degrees = gu.getDegreeMapWrapper(graph);
		List<String> degreeKeys = new ArrayList<String>(degrees.keySet());	
		Set<String> geneSet = gu.getGeneGraphSetWrapper(graph);
		List<Integer> distribution = gu.getGeneDegreeListWrapper(graph, geneSet);
		List<Double> cdf = calculateCumulativeDistribution(distribution);
		
		List<Integer> swaps = new ArrayList<Integer>();
		List<String> timePermutations = new ArrayList<String>();
		String newDir = directory+"/output/permutations/";			
		for (int i=1; i<numPermutations+1; i++){
			Graph<String, String> tempGraph = new UndirectedSparseGraph<String, String>();
			for (String v: graph.getVertices())
				tempGraph.addVertex(v);
			for (String e: graph.getEdges())
				tempGraph.addEdge(e, graph.getIncidentVertices(e));
			String newFile = fileName+"_"+i;
		
			long start1 = System.currentTimeMillis();
//			swaps.add(switchAlgorithm(newDir, newFile, tempGraph, numSwap, windowThreshold, degreeKeys, cdf));
			swaps.add(switchAlgorithmRedo2(newDir, newFile, tempGraph, numSwap, windowThreshold, degreeKeys));
			long end1 = System.currentTimeMillis();
			
			//Check that permutated graph has the same nodes as original graph and the same number of degrees <---- not sure if needed, switchAlgorithm should ensure this
/*			Graph<String, String> tempGraph2 = graph;
			checkPermutation(newDir, newFile, tempGraph2);	
*/
			int numVertex = tempGraph.getVertexCount();
			int numEdge = tempGraph.getEdgeCount();
//			System.out.println("original\t network " + i + "\\" +  numPermutations + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwaps: " + swaps.get(i-1) + "\tWindowThreshold: " + windowThreshold + "\tSwitch Algorithm Time: " + ((end1 - start1) / 1000) + " seconds");
			System.out.println("original\t network " + i + "\\" +  numPermutations + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwaps: " + swaps.get(i-1) + "\tWindowThreshold: " + windowThreshold + "\tSwitch Algorithm Time: " + ((end1 - start1)) + " ms");
			timePermutations.add(i + "\t" + numPermutations + "\t" + numVertex + "\t" + numEdge + "\t" + swaps.get(i-1) + "\t" + windowThreshold + "\t" + ((end1 - start1)/1000));
		}
		String permuteTimeFile = "permuteTime_" + graph.getVertexCount(); 
		savePermuteTimeListToFile(newDir, permuteTimeFile, timePermutations);
		
		return swaps; 
	}
	 
	/**
	 * Performs switch algorithm for generating and saving a random network.
	 * <p>
	 * Based on Python NetworkX's <a href="https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/algorithms/swap.html#connected_double_edge_swap">connected_double_edge_swap()</a>.
	 * @param directory - Directory to save random network file in.
	 * @param fileName - File to save random network in.
	 * @param graph - Undirected graph.
	 * @param numSwap - Number of swaps to perform.
	 * @param windowThreshold - Number of times to swap or attempt swaps before checking if graph is still connected.
	 * @return number of swaps used to generate random network.
	 * @throws IOException 
	 */
	private int switchAlgorithm(String directory, String fileName, Graph<String, String> graph, int numSwap, int windowThreshold, List<String> degreeKeys, List<Double>cdf) throws IOException{	
		long time = 0; 
		long start1 = System.currentTimeMillis();
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		long end1 = System.currentTimeMillis();
		time += (end1 - start1); 
		
		if (components.size() != 1){
			System.err.println("Provided graph is not connected");
			System.exit(0);
		}
		if (graph.getVertexCount() < 4){
			System.err.println("Provided graph has less than 4 nodes");
			System.exit(0);
		}
		
		int n = 0;
		int swapcount = 0;
		int window = 1;

		while(n < numSwap){
			int wcount = 0;
			Queue<String> swapped = new LinkedList<String>(); 
			if(window < windowThreshold){
				Boolean fail = false; 
				while(wcount<window && n<numSwap){
					//pick 2 random edges
					List<Integer> random2Edges = calculateDiscreteSequence(2, cdf);
					int ui = random2Edges.get(0);
					int xi = random2Edges.get(1);
					//skip if same source node indices 
					if (ui == xi)
						continue;
					//convert indices to source nodes
					String u = degreeKeys.get(ui);
					String x = degreeKeys.get(xi);
					//Choose targets nodes from source's neighboring nodes
					List<String> uList = new ArrayList<String>();
					uList.addAll(graph.getNeighbors(u));				
					int random = getRandIndexTargetNode(uList);
					String v = uList.get(random);
					List<String> xList = new ArrayList<String>();
					xList.addAll(graph.getNeighbors(x));
					random = getRandIndexTargetNode(xList);
					String y = xList.get(random);
					//skip if same target nodes
					if (v == y)
						continue;
					//swap if x is not a neighbor of u and y is not a neighbor of v
					List<String> vList = new ArrayList<String>();
					vList.addAll(graph.getNeighbors(v));
					if (!uList.contains(x) && !vList.contains(y)){
						if (graph.containsEdge(u + "\t" + v))
							graph.removeEdge(u + "\t" +v);
						else
							graph.removeEdge(v + "\t" + u);
						if (graph.containsEdge(x + "\t" + y))
							graph.removeEdge(x + "\t" + y);
						else
							graph.removeEdge(y + "\t" + x);
						graph.addEdge(u + "\t" + x, u, x);
						graph.addEdge(v + "\t" + y, v, y);
						swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
						swapcount += 1;
					}
					n += 1; 
					//if graph is connected, increase window size
					start1 = System.currentTimeMillis();
					components = wcc.transform(graph);
					end1 = System.currentTimeMillis();
					time += ((end1 - start1));
					if (components.size() == 1)
						wcount += 1; 
					else {
						//undo swap, graph is disconnected
						graph.addEdge(u + "\t" + v, u, v);
						graph.addEdge(x + "\t" + y, x, y);
						graph.removeEdge(u + "\t" + x);
						graph.removeEdge(v + "\t" + y);
						swapcount -= 1; 
						fail = true; 
					}
				}
				//if a swap fails, reduce window size
				if (fail)
					window = (int) Math.ceil(window/2);
				else
					window += 1; 
			}
			//a large window indicates many swaps should work, so do all swaps then check connection
			else{
				while(wcount<window && n<numSwap){
					//pick 2 random edges
					List<Integer> random2Edges = calculateDiscreteSequence(2, cdf);
					int ui = random2Edges.get(0);
					int xi = random2Edges.get(1);
					//skip same source nodes
					if (ui == xi)
						continue;
					//convert index to node label
					String u = degreeKeys.get(ui);
					String x = degreeKeys.get(xi);
					//Choose targets from neighboring nodes
					List<String> uList = new ArrayList<String>();
					uList.addAll(graph.getNeighbors(u));
					int random = getRandIndexTargetNode(uList);
					String v = uList.get(random);
					List<String> xList = new ArrayList<String>();
					xList.addAll(graph.getNeighbors(x));
					random = getRandIndexTargetNode(xList);
					String y = xList.get(random);
					//skip if same target nodes
					if (v == y)
						continue;
					//swap if x is not a neighbor of u and y is not a neighbor of v
					List<String> vList = new ArrayList<String>();
					vList.addAll(graph.getNeighbors(v));
					if (!uList.contains(x) && !vList.contains(y)){
						if (graph.containsEdge(u + "\t" + v))
							graph.removeEdge(u + "\t" +v);
						else
							graph.removeEdge(v + "\t" + u);
						if (graph.containsEdge(x + "\t" + y))
							graph.removeEdge(x + "\t" + y);
						else
							graph.removeEdge(y + "\t" + x);
						graph.addEdge(u + "\t" + x, u, x);
						graph.addEdge(v + "\t" + y, v, y);						
						swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
						swapcount += 1;
					}
					n += 1; 
					wcount += 1; 
					start1 = System.currentTimeMillis();
					components = wcc.transform(graph);
					end1 = System.currentTimeMillis();
					time += ((end1 - start1));
					if (components.size() == 1)
						window += 1; 
					else {
						//undo swaps from previous window and reduce window size, because the graph is disconnected
						while (!swapped.isEmpty()){
							String fourNodes = swapped.poll();
							String[] nodes = fourNodes.split("\t");
							u = nodes[0];
							v = nodes[1];
							x = nodes[2];
							y = nodes[3];
							graph.addEdge(u + "\t" + v, u, v);
							graph.addEdge(x + "\t" + y, x, y);
							graph.removeEdge(u + "\t" + x);
							graph.removeEdge(v + "\t" + y);
							swapcount -= 1; 
						}
						window = (int) Math.ceil(window/2);
					}
				}
			}
		}
		//Save the interaction pairs in the permuted network and return the number of swaps
		GraphUtils gu = new GraphUtils();
		FileUtils fu = new FileUtils();
		Set<String> pairs = gu.getPairsWrapper(graph);
		fu.saveSetToFile(directory, fileName, pairs);
		System.out.println("\t\t\t\t\t\t\t\t\t\t\t\tTotal component check time: " +  (time / 1000) + " seconds");
		return swapcount;
	}
	
	//**************************************************************************************************************asdfasdf 

	private void testRandomNetworkPermutationSwitchRedo(String directory, int numNodes, int numPermutations) throws IOException{
		GraphUtils gu = new GraphUtils();
		FileUtils fu = new FileUtils();
		Graph<String, String> graph = createRandomConnectedGraphTest(numNodes);
		Set<String> pairs = gu.getPairsWrapper(graph);
		fu.saveSetToFile(directory+"/output/permutations/", "randomGraphForPermutation", pairs);
		
		String permuteEdgeListFile = "edgeListPermuted";
		List<Integer> numSwaps = permuteNetworkSwitchRedo(directory, permuteEdgeListFile, numPermutations, graph);
	}
		
	private List<Integer> permuteNetworkSwitchRedo(String directory, String fileName, int n, Graph<String, String> graph) throws IOException{
		System.out.println("---------- Creating permuted Networks ----------\t redo");
		
		GraphUtils gu = new GraphUtils();
		
//		int edgeSwapConstant = 10;//115;
//		int numEdges = graph.getEdgeCount();
//		int numSwap = edgeSwapConstant * numEdges;
		int numSwap = 1; 
		int windowThreshold = 3;
		
		HashMap<String, Integer> degrees = gu.getDegreeMapWrapper(graph);
		List<String> degreeKeys = new ArrayList<String>(degrees.keySet());
		Set<String> geneSet = gu.getGeneGraphSetWrapper(graph);
//		List<Integer> distribution = getGeneDegreeList(graph, geneSet);
//		List<Double> cdf = calculateCumulativeDistribution(distribution);
		
		List<Integer> swaps = new ArrayList<Integer>();
		List<String> timePermutations = new ArrayList<String>();
		String newDir = directory+"/output/permutations2/";
		for (int i=1; i<n+1; i++){
			Graph<String, String> tempGraph = new UndirectedSparseGraph<String, String>();
			for (String v: graph.getVertices())
				tempGraph.addVertex(v);
			for (String e: graph.getEdges())
				tempGraph.addEdge(e, graph.getIncidentVertices(e));
			String newFile = fileName+"_"+i;
			
			long start1 = System.currentTimeMillis();
			swaps.add(switchAlgorithmRedo(newDir, newFile, tempGraph, numSwap, windowThreshold, degreeKeys));
			long end1 = System.currentTimeMillis();
			//Check that permutated graph has the same nodes as original graph and the same number of degrees <---- not sure if needed, switchAlgorithm should ensure this
/*			Graph<String, String> tempGraph2 = graph;
			checkPermutation(newDir, newFile, tempGraph2);	
*/			
			int numVertex = tempGraph.getVertexCount();
			int numEdge = tempGraph.getEdgeCount();
//			System.out.println("redo\t network " + i + "\\" +  n + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwaps: " + swaps.get(i-1) + "\tWindowThreshold: " + windowThreshold + "\tSwitch Algorithm Time: " + ((end1 - start1) / 1000) + " seconds");
			System.out.println("redo\t network " + i + "\\" +  n + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwaps: " + swaps.get(i-1) + "\tWindowThreshold: " + windowThreshold + "\tSwitch Algorithm Time: " + ((end1 - start1)) + " ms");
			timePermutations.add(i + "\t" + n + "\t" + numVertex + "\t" + numEdge + "\t" + swaps.get(i-1) + "\t" + windowThreshold + "\t" + ((end1 - start1)/1000));
		}
		String permuteTimeFile = "permuteTime_" + graph.getVertexCount(); 
		savePermuteTimeListToFile(newDir, permuteTimeFile, timePermutations);
		
		return swaps; 
	}
		
	private int switchAlgorithmRedo(String directory, String fileName, Graph<String, String> graph, int numSwap, int windowThreshold, List<String> degreeKeys) throws IOException{	
		long time = 0; 
		long start1 = System.currentTimeMillis();
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		long end1 = System.currentTimeMillis();
		time += (end1 - start1);
		
		if (components.size() != 1){
			System.err.println("Provided graph is not connected");
			System.exit(0);
		}
		if (graph.getVertexCount() < 4){
			System.err.println("Provided graph has less than 4 nodes");
			System.exit(0);
		}
		int n = 0;
		int swapcount = 0;
		int window = 1; 
		while(n < numSwap){
			int wcount = 0;
			Queue<String> swapped = new LinkedList<String>(); 
			if(window < windowThreshold){
				Boolean fail = false; 
				while(wcount<window && n<numSwap){
					//pick 2 random edges
					int ui = getRandomIndex(degreeKeys.size());
					int xi = getRandomIndex(degreeKeys.size());
					
					//skip if same source node indices 
					if (ui == xi)
						continue;
					//convert indices to source nodes
					String u = degreeKeys.get(ui);
					String x = degreeKeys.get(xi);
					//Choose targets nodes from source's neighboring nodes
					List<String> uList = new ArrayList<String>();
					uList.addAll(graph.getNeighbors(u));
					int random = getRandIndexTargetNode(uList);
					String v = uList.get(random);
					List<String> xList = new ArrayList<String>();
					xList.addAll(graph.getNeighbors(x));
					random = getRandIndexTargetNode(xList);
					String y = xList.get(random);
					//skip if same target nodes
					if (v == y)
						continue;
					//swap if x is not a neighbor of u and y is not a neighbor of v
					List<String> vList = new ArrayList<String>();
					vList.addAll(graph.getNeighbors(v));
					if (!uList.contains(x) && !vList.contains(y)){
						if (graph.containsEdge(u + "\t" + v))
							graph.removeEdge(u + "\t" +v);
						else
							graph.removeEdge(v + "\t" + u);
						if (graph.containsEdge(x + "\t" + y))
							graph.removeEdge(x + "\t" + y);
						else
							graph.removeEdge(y + "\t" + x);
						graph.addEdge(u + "\t" + x, u, x);
						graph.addEdge(v + "\t" + y, v, y);
						swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
						swapcount += 1;
					}
					n += 1; 
					//if graph is connected, increase window size
					start1 = System.currentTimeMillis();
					components = wcc.transform(graph);
					end1 = System.currentTimeMillis();
					time += ((end1 - start1));
					if (components.size() == 1)
						wcount += 1; 
					else {
						//undo swap, graph is disconnected
						graph.addEdge(u + "\t" + v, u, v);
						graph.addEdge(x + "\t" + y, x, y);
						graph.removeEdge(u + "\t" + x);
						graph.removeEdge(v + "\t" + y);
						swapcount -= 1; 
						fail = true; 
					}
				}
				//if a swap fails, reduce window size
				if (fail)
					window = (int) Math.ceil(window/2);
				else
					window += 1; 
			}
			//a large window indicates many swaps should work, so do all swaps then check connection
			else{
				while(wcount<window && n<numSwap){
					//pick 2 random edges
					int ui = getRandomIndex(degreeKeys.size());
					int xi = getRandomIndex(degreeKeys.size());
					
					//skip same source nodes
					if (ui == xi)
						continue;
					//convert index to node label
					String u = degreeKeys.get(ui);
					String x = degreeKeys.get(xi);
					//Choose targets from neighboring nodes
					List<String> uList = new ArrayList<String>();
					uList.addAll(graph.getNeighbors(u));
					int random = getRandIndexTargetNode(uList);
					String v = uList.get(random);
					List<String> xList = new ArrayList<String>();
					xList.addAll(graph.getNeighbors(x));
					random = getRandIndexTargetNode(xList);
					String y = xList.get(random);
					//skip if same target nodes
					if (v == y)
						continue;
					//swap if x is not a neighbor of u and y is not a neighbor of v
					List<String> vList = new ArrayList<String>();
					vList.addAll(graph.getNeighbors(v));
					if (!uList.contains(x) && !vList.contains(y)){
						if (graph.containsEdge(u + "\t" + v))
							graph.removeEdge(u + "\t" +v);
						else
							graph.removeEdge(v + "\t" + u);
						if (graph.containsEdge(x + "\t" + y))
							graph.removeEdge(x + "\t" + y);
						else
							graph.removeEdge(y + "\t" + x);
						graph.addEdge(u + "\t" + x, u, x);
						graph.addEdge(v + "\t" + y, v, y);						
						swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
						swapcount += 1;
					}
					n += 1; 
					wcount += 1;
					start1 = System.currentTimeMillis();
					components = wcc.transform(graph);
					end1 = System.currentTimeMillis();
					time += ((end1 - start1));
					if (components.size() == 1)
						window += 1; 
					else {
						//undo swaps from previous window and reduce window size, because the graph is disconnected
						while (!swapped.isEmpty()){
							String fourNodes = swapped.poll();
							String[] nodes = fourNodes.split("\t");
							u = nodes[0];
							v = nodes[1];
							x = nodes[2];
							y = nodes[3];
							graph.addEdge(u + "\t" + v, u, v);
							graph.addEdge(x + "\t" + y, x, y);
							graph.removeEdge(u + "\t" + x);
							graph.removeEdge(v + "\t" + y);
							swapcount -= 1; 
						}
						window = (int) Math.ceil(window/2);
					}
				}
			}
		}
		//Save the interaction pairs in the permuted network and return the number of swaps
		GraphUtils gu = new GraphUtils();
		FileUtils fu = new FileUtils();
		Set<String> pairs = gu.getPairsWrapper(graph);
		fu.saveSetToFile(directory, fileName, pairs);
		System.out.println("\t\t\t\t\t\t\t\t\t\t\t\tTotal component check time: " + (time/1000) + " seconds");
		return swapcount;
	}

	
	
	
	
	
	private int switchAlgorithmRedo2(String directory, String fileName, Graph<String, String> graph, int numSwap, int windowThreshold, List<String> degreeKeys) throws IOException{	
		long time = 0; 
		long start1 = System.currentTimeMillis();
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		long end1 = System.currentTimeMillis();
		time += (end1 - start1);
		
		if (components.size() != 1){
			System.err.println("Provided graph is not connected");
			System.exit(0);
		}
		if (graph.getVertexCount() < 4){
			System.err.println("Provided graph has less than 4 nodes");
			System.exit(0);
		}
		int n = 0;
		int swapcount = 0;
		int window = 1; 
		while(n < numSwap){
			int wcount = 0;
			Queue<String> swapped = new LinkedList<String>(); 
			if(window < windowThreshold){
				Boolean fail = false; 
				while(wcount<window && n<numSwap){
					//pick 2 random edges
					int ui = getRandomIndex(degreeKeys.size());
					int xi = getRandomIndex(degreeKeys.size());
					
					//skip if same source node indices 
					if (ui == xi)
						continue;
					//convert indices to source nodes
					String u = degreeKeys.get(ui);
					String x = degreeKeys.get(xi);
					//Choose targets nodes from source's neighboring nodes
					List<String> uList = new ArrayList<String>();
					uList.addAll(graph.getNeighbors(u));
					int random = getRandIndexTargetNode(uList);
					String v = uList.get(random);
					List<String> xList = new ArrayList<String>();
					xList.addAll(graph.getNeighbors(x));
					random = getRandIndexTargetNode(xList);
					String y = xList.get(random);
					//skip if same target nodes
					if (v == y)
						continue;
					//swap if x is not a neighbor of u and y is not a neighbor of v
					List<String> vList = new ArrayList<String>();
					vList.addAll(graph.getNeighbors(v));
					if (!uList.contains(x) && !vList.contains(y)){
						if (graph.containsEdge(u + "\t" + v))
							graph.removeEdge(u + "\t" +v);
						else
							graph.removeEdge(v + "\t" + u);
						if (graph.containsEdge(x + "\t" + y))
							graph.removeEdge(x + "\t" + y);
						else
							graph.removeEdge(y + "\t" + x);
						graph.addEdge(u + "\t" + x, u, x);
						graph.addEdge(v + "\t" + y, v, y);
						swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
						swapcount += 1;
					}
					n += 1; 
					//if graph is connected, increase window size
					start1 = System.currentTimeMillis();
					boolean connected = isConnectedWrapper(u, x, v, y, graph);
					end1 = System.currentTimeMillis();
					time += ((end1 - start1));
					if (connected == true)
						wcount += 1; 
					else {
						//undo swap, graph is disconnected
						graph.addEdge(u + "\t" + v, u, v);
						graph.addEdge(x + "\t" + y, x, y);
						graph.removeEdge(u + "\t" + x);
						graph.removeEdge(v + "\t" + y);
						swapcount -= 1; 
						fail = true; 
					}
				}
				//if a swap fails, reduce window size
				if (fail)
					window = (int) Math.ceil(window/2);
				else
					window += 1; 
			}
			//a large window indicates many swaps should work, so do all swaps then check connection
			else{
				while(wcount<window && n<numSwap){
					//pick 2 random edges
					int ui = getRandomIndex(degreeKeys.size());
					int xi = getRandomIndex(degreeKeys.size());
					
					//skip same source nodes
					if (ui == xi)
						continue;
					//convert index to node label
					String u = degreeKeys.get(ui);
					String x = degreeKeys.get(xi);
					//Choose targets from neighboring nodes
					List<String> uList = new ArrayList<String>();
					uList.addAll(graph.getNeighbors(u));
					int random = getRandIndexTargetNode(uList);
					String v = uList.get(random);
					List<String> xList = new ArrayList<String>();
					xList.addAll(graph.getNeighbors(x));
					random = getRandIndexTargetNode(xList);
					String y = xList.get(random);
					//skip if same target nodes
					if (v == y)
						continue;
					//swap if x is not a neighbor of u and y is not a neighbor of v
					List<String> vList = new ArrayList<String>();
					vList.addAll(graph.getNeighbors(v));
					if (!uList.contains(x) && !vList.contains(y)){
						if (graph.containsEdge(u + "\t" + v))
							graph.removeEdge(u + "\t" +v);
						else
							graph.removeEdge(v + "\t" + u);
						if (graph.containsEdge(x + "\t" + y))
							graph.removeEdge(x + "\t" + y);
						else
							graph.removeEdge(y + "\t" + x);
						graph.addEdge(u + "\t" + x, u, x);
						graph.addEdge(v + "\t" + y, v, y);						
						swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
						swapcount += 1;
					}
					n += 1; 
					wcount += 1;
					start1 = System.currentTimeMillis();
					boolean connected = isConnectedWrapper(u, x, v, y, graph);
					end1 = System.currentTimeMillis();
					time += ((end1 - start1));
					if (connected == true)
						window += 1;
					else {
						//undo swaps from previous window and reduce window size, because the graph is disconnected
						while (!swapped.isEmpty()){
							String fourNodes = swapped.poll();
							String[] nodes = fourNodes.split("\t");
							u = nodes[0];
							v = nodes[1];
							x = nodes[2];
							y = nodes[3];
							graph.addEdge(u + "\t" + v, u, v);
							graph.addEdge(x + "\t" + y, x, y);
							graph.removeEdge(u + "\t" + x);
							graph.removeEdge(v + "\t" + y);
							swapcount -= 1; 
						}
						window = (int) Math.ceil(window/2);
					}
				}
			}
		}
		//Save the interaction pairs in the permuted network and return the number of swaps
		GraphUtils gu = new GraphUtils();
		FileUtils fu = new FileUtils();
		Set<String> pairs = gu.getPairsWrapper(graph);
		fu.saveSetToFile(directory, fileName, pairs);
		System.out.println("\t\t\t\t\t\t\t\t\t\t\t\tTotal component check time: " + (time/1000) + " seconds");
		return swapcount;
	}

	
	
	
	
	
	
	
	private List<Integer> get2RandomEdgeIndices(int max){
		List<Integer> randomEdges = new ArrayList<Integer>();
		Random randomGen = new Random();
		int  random = randomGen.nextInt(max);
		randomEdges.add(random);
		int random2 = randomGen.nextInt(max);
		randomEdges.add(random2);
		return randomEdges;
	}
		
	
	/** 
	 * Gets random index of a list. 
	 * @param max - Size of list.
	 * @return a random index of a list.
	 */
	private int getRandomIndex(int max){ 
		Random randomGen = new Random();
		int  random = randomGen.nextInt(max);
		return random; 
	}

	/**
	 * Gets the random index from a list containing neighboring nodes.
	 * <p>
	 * <b>Note:</b> used in the switchAlgorithm() to obtain random target nodes for swaps. 
	 * @param neighborList- list for obtaining random index from.
	 * @return a random index of a list.
	 */
	private int getRandIndexTargetNode (List<String> neighborList){		
		Random randomGen = new Random();
		int maxRange = neighborList.size()-1; 
		if (maxRange == 0)
			return 0;
		int random = randomGen.nextInt(neighborList.size()-1);
		return random; 
	}
	
	/**
	 * Calculates the normalized cumulative distribution.
	 * <p>
	 * Based on Python NetworkX's <a href="https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/utils/random_sequence.html#cumulative_distribution">cumulative_distribution()</a>.
	 * @param distribution - used to create normalized cumulative distribution.
	 * @return the normalized cumulative distribution.
	 */
	private List<Double> calculateCumulativeDistribution(List<Integer> distribution){
		List<Double> cdf = new ArrayList<Double>();
		cdf.add(0.0);
		double psum = 0;
		for(int i=0; i<distribution.size(); i++)
			psum += distribution.get(i);
		for (int i=0; i<distribution.size(); i++)
			cdf.add(cdf.get(i)+distribution.get(i)/psum);		
		return cdf; 
	}
	
	/**
	 * Calculates the discrete sequence.
	 * <p>
	 * Based on Python NetworkX's <a href="https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/utils/random_sequence.html#discrete_sequence">discrete_sequence()</a>.			
	 * @param n - Length of the sequence.
	 * @param cumulativeDistribution - Normalized cumulative distribution from {@link #calculateCumulativeDistribution(List)}.
	 * @return sequence of length n from the normalized cumulative distribution.
	 */
	private List<Integer> calculateDiscreteSequence(int n, List<Double> cumulativeDistribution){
		List<Double> inputseq = new ArrayList<Double>();
		for (int i=0; i<n; i++){
			Random randomGen = new Random();
			Double random = randomGen.nextDouble();
			inputseq.add(random);
		}
		List<Integer> seq = new ArrayList<Integer>();
		for (Double s: inputseq)
			seq.add(calculateBisectLeft(cumulativeDistribution, s)-1);
		return seq; 
	}
		
	/**
	 * Calculates where to insert a sequence so that the cumulative distribution list will still be sorted.
	 * <p>
	 * Based on Python's <a href="https://hg.python.org/cpython/file/2.7/Lib/bisect.py">bisect_left()</a>.
	 * @param cumulativeDistribution - Normalized cumulative distribution from {@link #calculateCumulativeDistribution(List)}.
	 * @param seq - Value to be inserted in cumulative distribution list.
	 * @return the index where seq should be inserted in the sorted cumulativeDistribution list.
	 */
	private int calculateBisectLeft(List<Double> cumulativeDistribution, Double seq){
		int lo = 0;
		int hi = cumulativeDistribution.size();
		while (lo < hi){
			int mid = Math.floorDiv(lo+hi, 2);
			if (cumulativeDistribution.get(mid) < seq)
				lo = mid + 1;
			else
				hi = mid;
		}
		return lo; 
	}

	/**
	 * Saves a list with permutation times into a file.
	 * <p>
	 * <b>Note:</b> details of list will be saved in the following order: permutation, numPermutes, numVertex, numEdge, networkSize, time(seconds)
	 * @param directory - Directory to save the file in.	
	 * @param fileName - File name to save list in.
	 * @param listSave - List with permutation times to be saved.
	 * @throws IOException
	 */
	private void savePermuteTimeListToFile(String directory, String fileName, List<String> listSave) throws IOException{
		File file = new File(directory+ fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + fileName;
		PrintWriter pw = new PrintWriter(filePath);		
		pw.println("permutation\t" + "numPermutes\t" + "numVertex\t" + "numEdge\t" + "numSwaps\t" + "time(seconds)");
		for (int i=0; i<listSave.size(); i++){
			List<String> sub = listSave.subList(i, Math.min(listSave.size(), i+1));
			for (String s: sub)
				pw.print(s + "\t");
			pw.println();
		}
		pw.close();
	}
	
	@Test
	public void testReactomeConnectionTime() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
//		String fiFile = "FIsInGene_temp.txt";
		String fiFile = "FIsInGene_031516_with_annotations.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		GraphUtils gu = new GraphUtils(); 
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraphWrapper(directory, fiFile);
		Graph<String, String> graph = gu.createLargestComponentGraphWrapper(allGenesGraph);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		System.out.println("components: " + components.size() +" " + components);
		System.out.println("size: " + graph.getVertices().size());

		Set<String> connectTimeList = new HashSet<String>();
		for (int i=0; i<1000; i++){
			Graph<String, String> tempGraph = new UndirectedSparseGraph<String, String>();
			for (String v: graph.getVertices())
				tempGraph.addVertex(v);
			for (String e: graph.getEdges())
				tempGraph.addEdge(e, graph.getIncidentVertices(e));

			String connectTime = testReactomeRepeat(tempGraph);
			connectTimeList.add(connectTime);
		}
		
		FileUtils fu = new FileUtils();
		String fileName = "conectionTimes.txt";
		fu.saveSetToFile(directory+"/output/", fileName, connectTimeList);
		System.out.println(directory+"/output/"+fileName);
	}

	/**
	 * Performs different connection tests to determine if the Reactome FI network is still connected after swapping 2 edges.
	 * @param graph - Graph to perform swapping on
	 */
	private String testReactomeRepeat(Graph<String, String> graph){		
		List<String> edges = new ArrayList<String>(graph.getEdges());
		List<Integer> random2Edges = new ArrayList<Integer>();
		for (int i=0; i<2; i++){
			Random randomGen = new Random();
			int random = randomGen.nextInt(graph.getVertexCount()-1)+1;
			random2Edges.add(random);
		}
		String[] removeEdge1 = edges.get(random2Edges.get(0)).split("\t");
		String u = removeEdge1[0];
		String v = removeEdge1[1];
		String[] removeEdge2 = edges.get(random2Edges.get(1)).split("\t");
		String x = removeEdge2[0];
		String y = removeEdge2[1];
		
		graph.removeEdge(u + "\t" + v);
		graph.removeEdge(x + "\t" + y);
		graph.addEdge(u + "\t" + x, u, x);
		graph.addEdge(v + "\t" + y, v, y);
		
//		System.out.println("u: " + u + " x: " + x + " v: " + v + " y: " + y);
		
		//Use u and x to find if there is a connection to v and y
		long start1 = System.currentTimeMillis();
		boolean value = isConnectedWrapper(u, x, v, y, graph);
		long end1 = System.currentTimeMillis();
		long bfsTime = end1 - start1;
//		System.out.println(value);

		//Use u and x to see if there is connection to v and y and use v and y to see if there is a connection to u and x
		long start2 = System.currentTimeMillis();
		boolean value2 = isConnectedWrapper2(u, x, v, y, graph);
		long end2 = System.currentTimeMillis();
		long bfsTime2 = end2 - start2;
//		System.out.println(value);		
		
		//Use u and see if there is a connection to v and y
		long start5 = System.currentTimeMillis();
		boolean value5 = isConnectedWrapper3(u, x, v, y, graph);
		long end5 = System.currentTimeMillis();
		long bfsTime5 = end5 - start5;
//		System.out.println(value);
		
		//Use JUNG to see if graph is connected
		long start3 = System.currentTimeMillis(); 
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		if (components.size() != 1){
			throw new IllegalArgumentException("Provided graph is not connected. Components: " + components);
		}
		long end3 = System.currentTimeMillis();
		long wccTime = end3 - start3;
				
		//Use shortest path to see if graph is connected. This bad idea, slowest.
		long start4 = System.currentTimeMillis();
		DijkstraShortestPath<String, String> dsp = new DijkstraShortestPath<String, String> (graph);
		List<String> path = dsp.getPath(u, x);
		long end4 = System.currentTimeMillis();
		long pathTime = end4 - start4;
		
//		System.out.println("Time BFS1:\t" + bfsTime + "\tms\t BFS2:\t"+ bfsTime2 + "\tms\t WCC:\t" + wccTime + "\tms\t PATH:\t" + pathTime +"\tms\t" + "u: " + u + " x: " + x + " v: " + v + " y: " + y);
//		System.out.println("Time BFS1:\t" + bfsTime + "\tms\t BFS2:\t"+ bfsTime2 + "\tms\t WCC:\t" + wccTime + "\tms\t" + "u: " + u + " x: " + x + " v: " + v + " y: " + y);
		return("Time BFS1:\t" + bfsTime + "\tms\t BFS2:\t"+ bfsTime2 + "\tms\t BFS3:\t"+ bfsTime5 + "\tms\t WCC:\t" + wccTime + "\tms\t PATH:\t" + pathTime + "\tms\t" + "u: " + u + " x: " + x + " v: " + v + " y: " + y);
	}
	
	
	@Test
	public void testConnection(){
		GraphUtils gu = new GraphUtils();
		
		Set<String> genes = new HashSet();
		for (int i=1; i<21; i++){
			genes.add("A"+i);
		}
		Set<String> pairs = new HashSet();
		pairs.add("A1	A2");
		pairs.add("A2	A3");
		pairs.add("A3	A4");
		pairs.add("A5	A4");
		pairs.add("A6	A4");
		pairs.add("A6	A7"); //comment out for disconnected graph A1-A2 to A7-A9
		pairs.add("A8	A7");
		pairs.add("A9	A7");	
		pairs.add("A3	A10");
		pairs.add("A10	A11");
		pairs.add("A11	A12");
		pairs.add("A12	A14");
		pairs.add("A12	A13");
		pairs.add("A13	A15");
		pairs.add("A16	A15");
		pairs.add("A16	A17");
		pairs.add("A18	A17");
		pairs.add("A18	A19");
		pairs.add("A20	A19");
		pairs.add("A1	A15");
		Graph<String, String> graph = gu.createGraphWrapper(genes, pairs);
		boolean value = isConnectedWrapper("A1", "A2", "A9", "A7", graph);
		System.out.print(value);
		boolean value2 = isConnectedWrapper("A1", "A2", "A9", "A7", graph);
		System.out.print(value2);
	}
	
	@Test
	public void testRandomNetworkConnection(){
		GraphUtils gu = new GraphUtils();
		int numNodes = 12000; 
		Graph<String, String> graph = createRandomConnectedGraphTest(numNodes);
		System.out.println("nodes: " + graph.getVertexCount() + " edges:" + graph.getEdgeCount());
		List<String> edges = new ArrayList<String>(graph.getEdges());
		List<Integer> random2Edges = new ArrayList<Integer>();
		for (int i=0; i<2; i++){
			graph.addVertex(""+i);
			Random randomGen = new Random();
			int random = randomGen.nextInt(numNodes-1)+1;
			random2Edges.add(random);
		}
		String[] removeEdge1 = edges.get(random2Edges.get(0)).split("\t");
		String u = removeEdge1[0];
		String v = removeEdge1[1];
		String[] removeEdge2 = edges.get(random2Edges.get(1)).split("\t");
		String x = removeEdge2[0];
		String y = removeEdge2[1];
		
		graph.removeEdge(u + "\t" + v);
		graph.removeEdge(x + "\t" + y);
		graph.addEdge(u + "\t" + x, u, x);
		graph.addEdge(v + "\t" + y, v, y);
		
		System.out.println("u: " + u + " x: " + x + " v: " + v + " y: " + y);
		
		long start1 = System.currentTimeMillis();
		boolean value = isConnectedWrapper(u, x, v, y, graph);		
		long end1 = System.currentTimeMillis();
		long bfsTime1 = (end1 - start1);
		System.out.println("Time BFS1_: " + ((end1 - start1)) + " ms");
		
		start1 = System.currentTimeMillis();
		value = isConnectedWrapper(u, x, v, y, graph);		
		end1 = System.currentTimeMillis();
		long bfsTime2 = (end1 - start1);
		System.out.println("Time BFS1_: " + ((end1 - start1)) + " ms");

		start1 = System.currentTimeMillis();
		value = isConnectedWrapper2(u, x, v, y, graph);		
		end1 = System.currentTimeMillis();
		long bfsTimeW2 = (end1 - start1);
		System.out.println("Time BFS_2: " + ((end1 - start1)) + " ms");
		
		start1 = System.currentTimeMillis();
		value = isConnectedWrapper3(u, x, v, y, graph);		
		end1 = System.currentTimeMillis();
		long bfsTimeW3 = (end1 - start1);
		System.out.println("Time BFS_3: " + ((end1 - start1)) + " ms");
		
		start1 = System.currentTimeMillis(); 
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		if (components.size() != 1){
			throw new IllegalArgumentException("Provided graph is not connected. Components: " + components);
		}
		end1 = System.currentTimeMillis();
		long wccTime = (end1 - start1);
		System.out.println("Time WCC: " + ((end1 - start1)) + " ms");
		
		System.out.println("Time BFS1:\t" + bfsTime1 + "\tms\t BFS1:\t" + bfsTime2 + "\tms\t BFS2:\t" + bfsTimeW2+ "\tms\t BFS3:\t" + bfsTimeW3 + "\tms\t WCC:\t"+ wccTime +"\tms\t" + "u: " + u + " x: " + x + " v: " + v + " y: " + y);
	}
	
	/**
	 * Checks if swapped edges are connected by traversing from edge1's nodes towards edge2's nodes. 
	 * @param geneU	- Newly formed edge1's node1.
	 * @param geneX - Newly formed edge1's node2.
	 * @param geneV - Newly formed edge2's node1. 
	 * @param geneY - Newly formed edge2's node2. 
	 * @param graph - Graph used to check connection. 
	 * @return true if the graph is connected and false if the graph is not connected.
	 */
	private boolean isConnectedWrapper(String geneU, String geneX, String geneV, String geneY, Graph<String, String> graph){
		Set<String> neighbors = new HashSet<String>();
		neighbors.addAll(graph.getNeighbors(geneU));
		neighbors.addAll(graph.getNeighbors(geneX));
		Set<String> checked = new HashSet<String>();
		Set<String> result = isConnected(geneU, geneX, geneV, geneY, graph, neighbors, checked); 
		
		if (result.contains("true"))
			return true;
		return false; 
	}

	/**
	 * Checks if swapped edges are connected by traversing from one of edge1's nodes towards edge2's nodes. 
	 * @param geneU	- Newly formed edge1's node1.
	 * @param geneX - Newly formed edge1's node2.
	 * @param geneV - Newly formed edge2's node1. 
	 * @param geneY - Newly formed edge2's node2. 
	 * @param graph - Graph used to check connection. 
	 * @return true if the graph is connected and false if the graph is not connected.
	 */
	private boolean isConnectedWrapper3(String geneU, String geneX, String geneV, String geneY, Graph<String, String> graph){
		Set<String> neighbors = new HashSet<String>();
		neighbors.addAll(graph.getNeighbors(geneU));
		Set<String> checked = new HashSet<String>();
		Set<String> result = isConnected(geneU, geneX, geneV, geneY, graph, neighbors, checked); 
		
		if (result.contains("true"))
			return true;
		return false; 
	}

	/**
	 * Checks if the graph is connected.
	 * <p>
	 * Edges u-x and v-y are newly formed from edges u-v and x-y. This checks that u-x is able to make a path to v-y using recursion. 
	 * Example of swapped edges:<br>
	 * u - v &emsp;&emsp;&emsp;&emsp; u &emsp;v <br>
	 * &emsp;&emsp;&emsp;---> &emsp;&nbsp; | &emsp; | <br>
	 * x - y &emsp;&emsp;&emsp;&emsp; x	 &emsp; y
	 * @param geneU	- Newly formed edge1's node1.
	 * @param geneX - Newly formed edge1's node2.
	 * @param geneV - Newly formed edge2's node1. 
	 * @param geneY - Newly formed edge2's node2. 
	 * @param graph - Graph used to check connection. 
	 * @param neighbors - Set of direct neighbors from edge1.
	 * @param checked - Set of genes checked for connection to edge2.
	 * @return true if the graph is connected and false if the graph is not connected.
	 */
	private Set<String> isConnected(String geneU, String geneX, String geneV, String geneY, Graph<String, String> graph, Set<String> neighbors, Set<String> checked){
		if (checked.contains(geneV) || checked.contains(geneY)){
//			System.out.println("-- connected --");
			Set<String> returnSet = new HashSet<String>();
			returnSet.add("true");
			return returnSet;
		}
		if (checked.containsAll(neighbors)){
//			System.out.println("-- unconnected - checked has all neighbors --");
			Set<String> returnSet = new HashSet<String>();
			returnSet.add("false");
			return returnSet;
		} 
		else {
			HashSet<String> tempNeighbors = new HashSet<String>();
			for (String g: neighbors){
				if (! checked.contains(g)){
					tempNeighbors.addAll(graph.getNeighbors(g));
					checked.add(g);
				}
			}
			neighbors.addAll(tempNeighbors);
			return isConnected(geneU, geneX, geneV, geneY, graph, neighbors, checked);
		}
	}

	/**
	 * Checks if swapped edges are connected by traversing from edge1's nodes towards edge2's nodes and edge2's nodes towards edge1's nodes. 
	 * @param geneU	- Newly formed edge1's node1.
	 * @param geneX - Newly formed edge1's node2.
	 * @param geneV - Newly formed edge2's node1. 
	 * @param geneY - Newly formed edge2's node2. 
	 * @param graph - Graph used to check connection. 
	 * @return true if the graph is connected and false if the graph is not connected.
	 */
	private boolean isConnectedWrapper2(String geneU, String geneX, String geneV, String geneY, Graph<String, String> graph){
		Set<String> neighborsUX = new HashSet<String>();
		neighborsUX.addAll(graph.getNeighbors(geneU));
		neighborsUX.addAll(graph.getNeighbors(geneX));
		Set<String> checkedUX = new HashSet<String>();
		Set<String> neighborsVY = new HashSet<String>();
		neighborsVY.addAll(graph.getNeighbors(geneV));
		neighborsVY.addAll(graph.getNeighbors(geneY));
		Set<String> checkedVY = new HashSet<String>();
		Set<String> result = isConnected2(geneU, geneX, geneV, geneY, graph, neighborsUX, checkedUX, neighborsVY, checkedVY); 
		
		if (result.contains("true"))
			return true;
		return false; 
	}

	/**
	 * Checks if the graph is connected by traversing from edge1's nodes towards edge2's nodes and edge2's nodes towards edge1's nodes.
	 * @param geneU	- Newly formed edge1's node1.
	 * @param geneX - Newly formed edge1's node2.
	 * @param geneV - Newly formed edge2's node1. 
	 * @param geneY - Newly formed edge2's node2. 
	 * @param neighborsUX - Set of direct neighbors from edge1.
	 * @param checkedUX - Set of genes checked for connection to edge2.
	 * @param neighborsVY - Set of direct neighbors from edge2.
	 * @param checkedVY - Set of genes checked for connection to edge1.
	 * @return true if the graph is connected and false if the graph is not connected.
	 */
	private Set<String> isConnected2(String geneU, String geneX, String geneV, String geneY, Graph<String, String> graph, Set<String> neighborsUX, Set<String> checkedUX, Set<String> neighborsVY, Set<String> checkedVY){
		if (checkedUX.contains(checkedVY)){
//			System.out.println("-- connected --");
			Set<String> returnSet = new HashSet<String>();
			returnSet.add("true");
			return returnSet;
		}
		if (checkedUX.containsAll(neighborsUX)){
//			System.out.println("-- unconnected - checked has all neighbors --");
			Set<String> returnSet = new HashSet<String>();
			returnSet.add("false");
			return returnSet;
		} 
		else {
			HashSet<String> tempNeighborsUX = new HashSet<String>();
			for (String g: neighborsUX){
				if (! checkedUX.contains(g)){
					tempNeighborsUX.addAll(graph.getNeighbors(g));
					checkedUX.add(g);
				}
			}
			neighborsUX.addAll(tempNeighborsUX);
			HashSet<String> tempNeighborsVY = new HashSet<String>();
			for (String g: neighborsVY){
				if (! checkedVY.contains(g)){
					tempNeighborsVY.addAll(graph.getNeighbors(g));
					checkedVY.add(g);
				}
			}
			neighborsVY.addAll(tempNeighborsVY);
			return isConnected2(geneU, geneX, geneV, geneY, graph, neighborsUX, checkedUX, neighborsVY, checkedVY);
		}
	}
}