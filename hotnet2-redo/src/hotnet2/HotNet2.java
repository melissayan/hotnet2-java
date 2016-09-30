package hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.DoubleStream;

import org.junit.*;

import edu.uci.ics.jung.graph.*;//Graph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.algorithms.cluster.*;//WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.filters.*;//FilterUtils;
import edu.uci.ics.jung.algorithms.importance.*;//BetweennessCentrality;
import edu.uci.ics.jung.algorithms.matrix.*;
import edu.uci.ics.jung.graph.*;

import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SparseFieldMatrix;
import org.apache.commons.math3.stat.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;
import org.ojalgo.*;
import org.ojalgo.OjAlgoUtils;
import org.ojalgo.access.Access2D.Builder;
import org.ojalgo.matrix.*;
import org.ojalgo.matrix.jama.*;
import org.ojalgo.matrix.store.PrimitiveDenseStore;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;

public class HotNet2 {
	
	public static void main(String args[])throws IOException {
		System.out.println("--------------------------- START ---------------------------"); 
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String allGenesFile = "allGenes.txt";
		String allPairsFile = "allPairs.txt";
		String largestCompFile = "largestComponentPairs.txt";
//		String refLargestCompFile = "FIsInGene_031516_BigComp_temp.txt";
//		String refLargestCompFile = "FIsInGene_031516_BigComp.txt";
		String betweennessScoreFile = "betweennessScore.txt";
		String influenceFile = "influence.txt";
		String matrixOrderFile = "matrixOrder.txt";
		String heatScoreFile = "mutation_frequency_expr_filtered.txt";
		String geneMissingHeatFile = "geneMissingHeat.txt";
		String geneIndexFile = "geneIndexReactome.txt";
		String edgeListFile = "edgeListReactome.txt";
		String permuteEdgeListFile = "edgeListPermuted";
		HotNet2 obj = new HotNet2();
		DeltaSelection ds = new DeltaSelection();

		if (args.length==0){
			System.err.println("Arg[0] must indicate one of the ReactomeFI files below using a number:");
			System.err.println("\t 1: FIsInGene_temp.txt\n\t 2: FIsInGene_031516_with_annotations.txt");
			System.exit(0);
		}
		if (Integer.parseInt(args[0])==1){
			String fiFile = "FIsInGene_temp.txt";
			Path fiFilePath = Paths.get(directory, fiFile);
			obj.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile);
		}
		if (Integer.parseInt(args[0])==2){
			String fiFile = "FIsInGene_031516_with_annotations.txt";
			Path fiFilePath = Paths.get(directory, fiFile);
//			obj.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile);
			
			
			Graph<String, String> allGenesGraph = obj.createReactomeFIGraph(directory, fiFile);
			Graph<String, String> largestComponent = obj.createLargestComponentGraph(allGenesGraph);
			System.out.println("size: " + largestComponent.getVertices().size());
			List<Integer> numSwaps = obj.permuteNetworkSwitch(directory, permuteEdgeListFile, 5, largestComponent);
		}
		if (Integer.parseInt(args[0])==3){
			String fiFile = "iref_edge_list";
//			String fiFile = "iref_edge_list_temp";
			fiFile = obj.modifyEdgeFileToTabDelimited(directory, fiFile, false, " ");
			System.out.println(fiFile);
			Path fiFilePath = Paths.get(directory, fiFile);
			obj.selectBetaForIrefindex(directory, fiFile, betweennessScoreFile, influenceFile);
		}
		if (Integer.parseInt(args[0])==4){
			if (Integer.parseInt(args[1])==1){
				System.out.println("Network with 1000 edges permuted 10 times");
				obj.testRandomNetworkPermutationSwitch(directory, 1000, 10);			
			}
			if (Integer.parseInt(args[1])==2){
				System.out.println("Network with 10000 edges permuted 10 times");
				obj.testRandomNetworkPermutationSwitch(directory, 10000, 10);			
			}
		}
		
		//DeltaSelection
		if (Integer.parseInt(args[0])==5){
			System.out.println("Reactome FI network");
			BigDecimal tempBeta = new BigDecimal (args[1]);
			ds.selectDeltaWrapper(directory, tempBeta);			
		}
		if (Integer.parseInt(args[0])==6){
			System.out.println("iRefIndex network");
			BigDecimal tempBeta = new BigDecimal ("0.45");
			ds.selectDeltaForIrefindexWrapper(directory, tempBeta);			
		}
		
		
		
		System.out.println("---------------------------- END ----------------------------");
	}
	
	/**
	 * Modifies a file containing gene pairs by saving it as a new tab delimited file
	 * 		- used to test out interactions from the python implementation
	 * @param directory		directory containing file to be modified
	 * @param fileToModify	the file that will be modified into a tab delimited file
	 * @param header		Boolean to determine if the file has a header (true) or no header (false)
	 * @param delimiter		the delimiter in the file that will be modified
	 * @return				the name of the newly modified file
	 * @throws IOException
	 */
	private String modifyEdgeFileToTabDelimited(String directory, String fileToModify, Boolean header, String delimiter) throws IOException{
		//Set up for file that will get modified
		Path path = Paths.get(directory, fileToModify);
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		if (header == true)
			br.readLine();
		//Set up for modified file
		String modifiedFileName = "modified_" + fileToModify; 
		File file = new File(modifiedFileName);		
		if (!file.exists())
			file.createNewFile();
		String filePath = directory +"/" + modifiedFileName;
		PrintWriter pw = new PrintWriter(filePath);
		pw.println("gene1 \tgene2");
		//Read file that 
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene1 = row[0];
			String gene2 = row[1];
			pw.println(gene1 + "\t" + gene2);
		}
		br.close();
		pw.close();
		return modifiedFileName;
	}
	
	
	
	@Test
	public void Test() throws IOException{
	
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
//		String fiFile = "FIsInGene_temp.txt";
		String fiFile = "FIsInGene_031516_with_annotations.txt";
		Path fiFilePath = Paths.get(directory, fiFile);
		String allGenesFile = "allGenes.txt";
		String allPairsFile = "allPairs.txt";
		String largestCompFile = "largestComponentPairs.txt";
//		String refLargestCompFile = "FIsInGene_031516_BigComp_temp.txt";
//		String refLargestCompFile = "FIsInGene_031516_BigComp.txt";
		String betweennessScoreFile = "betweennessScore.txt";
		String influenceFile = "influence.txt";
		String matrixOrderFile = "matrixOrder.txt";
//		String heatScoreFile = "mutation_frequency_expr_filtered.txt";
		String heatScoreFile = "heatScore.txt";
		String geneIndexFile = "geneIndexReactome.txt";
		String edgeListFile = "edgeListReactome.txt";
		String permuteEdgeListFile = "edgeListPermuted";
		
		//Save ReactomeFI network files related to genes and gene interaction pairs
//		saveReactomeFIGraphFiles(directory, fiFile, allGenesFile, allPairsFile, largestCompFile);

		//Convert the largest component file into files to run the Python version of HotNet2 - requires saveReactomeFIGraphFiles() to work
//		convertFileForPythonHotNet2(directory, largestCompFile, geneIndexFile, edgeListFile);
//		convertFileForPythonHotNet2(directory, allPairsFile, geneIndexFile, edgeListFile);
		
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = createLargestComponentGraph(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());
		
		//Calculate and save betweenness centrality for all genes in largest component
//		saveBetweennessCentrality(directory, betweennessScoreFile, largestComponent);
		
		//Get 5 source proteins from betweennness centrality scores
//		getSourceProteins(directory, betweennessScoreFile);		
		
		
//		selectBeta(directory, fiFile, betweennessScoreFile, influenceFile);
		
		//Set of genes in the largest component, ordering determines matrix content ordering
		Set<String> geneSet = getGeneGraphSet(largestComponent);
		saveSetToFile(directory+"/output/", matrixOrderFile, geneSet);
		
		//Saves of set of genes with missing heat scores
		getMissingHeatGenes(directory, heatScoreFile, largestComponent);


/*		//HotNet2 Algorithm with Apache Commons Math for all steps
		BigDecimal tempBeta = new BigDecimal("0.5");	//get beta from user
		double beta = tempBeta.doubleValue();
		BigDecimal tempDelta = new BigDecimal("0.05");	//get delta from user
		double delta = tempDelta.doubleValue();
		RealMatrix F = createDiffusionMatrix(largestComponent, geneSet, beta);
//		createTempHeatScoreFile(directory, fiFile, heatScoreFile); //REMOVE THIS and use real heat scores, this is just for testing
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		RealMatrix E = createExchangedHeatMatrix(F, geneSet, heatScoreMap);
		RealMatrix H = identifyHotSubnetworks(E, delta);
//		String diffusionTempFile = "diffusion.txt";
//		saveMatrix(directory+"/output/", diffusionTempFile, F);
*/
		
/*		//HotNet2 Algorithm with ojAlgo diffusion
		BigDecimal tempBeta = new BigDecimal("0.5");	//get beta from user
		double beta = tempBeta.doubleValue();
		BigDecimal tempDelta = new BigDecimal("0.05");	//get delta from user
		double delta = tempDelta.doubleValue();
		PrimitiveMatrix tempF = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix F = convertOJAToACM(tempF); 
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		RealMatrix E = createExchangedHeatMatrix(F, geneSet, heatScoreMap);
		RealMatrix H = identifyHotSubnetworks(E, delta);		
		obtainSubnetworkSet(H, geneSet);
*/		
		
		//Converts OJA to ACM
/*		PrimitiveMatrix ojaM = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix acmM = convertOJAToACM(ojaM);
		System.out.println("ojaM: ------------\n" + ojaM);
		System.out.println("acmM: ------------\n" + acmM);
*/		
		
//		List<Integer> numSwaps = permuteNetwork(directory, permuteEdgeListFile, 5, largestComponent);
		
		
/*		String pythonEdgeListFile = "randomGraph_5.txt";
		testRandomNetwork(directory, pythonEdgeListFile);
*/		

		//Testing Random Networks
//		testRandomNetworkPermutationSwitch(directory, 1000, 5);
//		testRandomNetworkPermutationConfig(directory, 10000, 5);
		
		testRandomNetworkPermutationSwitch(directory, 1000, 5);
//		testRandomNetworkPermutationSwitchRedo(directory, 1000, 5);

		
		System.out.println("-----");
	}
	
	
	
	
	
	/**
	 * Create a random connected network using the number of provided nodes
	 * @param numNodes	number of nodes that should be in graph
	 * @return			a connected graph with the provided number of nodes 
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
	 * Uses a random network generated by Python
	 * 		--Python code to generate random graph with a given numNodes and save is:
	 * 			randomGraph = nx.complete_graph(numNodes);
	 * 			nx.write_edgelist(randomGraph, "randomGraph.txt")
	 * @param directory		directory where random network file is
	 * @param edgeListFile	name of random network file used for permutations
	 * @throws IOException
	 */
	private void testRandomPyNetwork(String directory, String edgeListFile) throws IOException{
		//Create graph from file and check there is only 1 component
		Path filePath = Paths.get(directory, edgeListFile);
		Set<String> allGenes = getAllGenes(filePath, " ");
		Set<String> allPairs = getAllInteractionPairs(filePath, " ");
		Graph<String, String> graph = createGraph(allGenes, allPairs);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		if (components.size() != 1){
			System.err.println("Provided edgeListFile " + edgeListFile + "'s graph is not connected");
			System.err.println("components: " + components);
			System.exit(0);
		}

		String permuteEdgeListFile = "edgeListPermuted";
		List<Integer> numSwaps = permuteNetworkSwitch(directory, permuteEdgeListFile, 10, graph);
	}
	
	/**
	 * Checks that the connected permuted graph has the same genes and same degrees as the provided graph
	 * @param directory				directory containing permuted graph's edge list
	 * @param permuteEdgeListFile	file with permuted graph's edge list
	 * @param graph					original graph used to create permuted graph
	 * @return						true if graph has the same gene and degrees, else false
	 * @throws IOException
	 */
	private Boolean checkPermutation(String directory, String permuteEdgeListFile, Graph<String, String> graph) throws IOException{
		//Create permuted graph from file and check there is only 1 component
		Path permuteFilePath = Paths.get(directory, permuteEdgeListFile);
		Set<String> allGenesPermuted = getAllGenes(permuteFilePath);
		Set<String> allPairsPermuted = getAllInteractionPairs(permuteFilePath);
		Graph<String, String> permutedGraph = createGraph(allGenesPermuted, allPairsPermuted);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(permutedGraph);
		if (components.size() != 1){
			System.err.println("Provided permuted graph " + permuteEdgeListFile + " is not connected");
			System.err.println("permutedGraph: " + permutedGraph);
			System.err.println("components: " + components);
			System.exit(0);
		}
		
		//Compare if the 2 graphs have the same genes and the gene's have the same degree
		Set<String> permutedGeneSet = getGeneGraphSet(permutedGraph);
		List<Integer> permutedGraphDegree = getGeneDegreeList(permutedGraph, permutedGeneSet);
		Set<String> originalGeneSet = getGeneGraphSet(graph);
		List<Integer> originalGraphDegree = getGeneDegreeList(graph, originalGeneSet);
		
		if (originalGeneSet.containsAll(permutedGeneSet) && permutedGeneSet.containsAll(originalGeneSet)){
			if (originalGraphDegree.containsAll(permutedGraphDegree) && permutedGraphDegree.containsAll(originalGraphDegree))
				return true;
			else{
				System.err.println("Permuted graph '" + permuteEdgeListFile + "' does not have the same degrees as the provided graph");
				System.err.println("original Degrees: " + originalGraphDegree);
				System.err.println("permuted Degrees: " + permutedGraphDegree);
				return false;
			}
		} else {
			System.err.println("Permuted graph '" + permuteEdgeListFile + "' does not have the same genes as the provided graph");
			System.err.println("original GeneSet: " + originalGeneSet);
			System.err.println("permuted GeneSet: " + permutedGeneSet);
			return false;
		}
	}
	
	/**
	 * Creates a random network and use it to generate random permutations using the switch algorithm.  All graphs are saved.
	 * 		--note: this is for testing purposes to try the switch algorithm on different networks
	 * @param directory		directory to save files in
	 * @param numNodes		number of nodes a graph should have
	 * @throws IOException
	 */
	private void testRandomNetworkPermutationSwitch(String directory, int numNodes, int numPermutations) throws IOException{
		Graph<String, String> graph = createRandomConnectedGraphTest(numNodes);
		Set<String> pairs = getPairs(graph);
		saveSetToFile(directory+"/output/permutations/", "randomGraphForPermutation", pairs);
		
/*		System.out.println("\nnetwork config");
		String permuteEdgeListFile2 = "edgeListPermuted";
		permuteNetworkConfig(directory, permuteEdgeListFile2, numPermutations, graph);
*/		
		String permuteEdgeListFile = "edgeListPermuted";
		System.out.println("\n-------------------------------------Switch original-----------------------------------------");
		List<Integer> numSwaps = permuteNetworkSwitch(directory, permuteEdgeListFile, numPermutations, graph);
		System.out.println("\n-------------------------------------Switch Redo-----------------------------------------");
		List<Integer> numSwapsRedo = permuteNetworkSwitchRedo(directory, permuteEdgeListFile, numPermutations, graph);		
	}

	/**
	 * Permute network n times and save each permutation in a file 
	 * 		-- based on Python HotNet2's permuteNetwork.py permute_network()
	 * @param graph
	 * @param edgeSwapConstant
	 * @param directory
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	private List<Integer> permuteNetworkSwitch(String directory, String fileName, int numPermutations, Graph<String, String> graph) throws IOException{
		System.out.println("\t---------- Creating permuted Networks ----------\t original");
		
		int edgeSwapConstant = 10;//115;
		int numEdges = graph.getEdgeCount();
		int numSwap = edgeSwapConstant * numEdges; 
		int windowThreshold = 3;
		
		HashMap<String, Integer> degrees = getDegreeMap(graph);
		List<String> degreeKeys = new ArrayList<String>(degrees.keySet());	
		Set<String> geneSet = getGeneGraphSet(graph);
		List<Integer> distribution = getGeneDegreeList(graph, geneSet);
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
			swaps.add(switchAlgorithm(newDir, newFile, tempGraph, numSwap, windowThreshold, degreeKeys, cdf));
			long end1 = System.currentTimeMillis();
			
			//Check that permutated graph has the same nodes as original graph and the same number of degrees <---- not sure if needed, switchAlgorithm should ensure this
/*			Graph<String, String> tempGraph2 = graph;
			checkPermutation(newDir, newFile, tempGraph2);	
*/
			int numVertex = tempGraph.getVertexCount();
			int numEdge = tempGraph.getEdgeCount();
			System.out.println("original\t network " + i + "\\" +  numPermutations + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwaps: " + swaps.get(i-1) + "\tWindowThreshold: " + windowThreshold + "\tSwitch Algorithm Time: " + ((end1 - start1) / 1000) + " seconds");
			timePermutations.add(i + "\t" + numPermutations + "\t" + numVertex + "\t" + numEdge + "\t" + swaps.get(i-1) + "\t" + windowThreshold + "\t" + ((end1 - start1)/1000));
		}
		String permuteTimeFile = "permuteTime_" + graph.getVertexCount(); 
		savePermuteTimeListToFile(newDir, permuteTimeFile, timePermutations);
		
		return swaps; 
	}
	 
	/**
	 * Switch algorithm for generating and saving a random network
	 * 		-- based on Python networkx's connected_double_edge_swap()
	 * 			https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/algorithms/swap.html#connected_double_edge_swap
	 * @param directory			directory to save random network file in
	 * @param fileName			file to save random network in 
	 * @param graph				undirected graph
	 * @param numSwap			number of swaps to perform
	 * @param windowThreshold	number of times to swap or attempt swaps before checking if graph is still connected
	 * @return					number of swaps used to generate random network			
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
		Set<String> pairs = getPairs(graph);
		saveSetToFile(directory, fileName, pairs);
		System.out.println("\t\t\t\t\t\t\t\t\t\t\t\tTotal component check time: " +  (time / 1000) + " seconds");
		return swapcount;
	}
	

	
	
	
	
	
	
//**************************************************************************************************************asdfasdf 

	private void testRandomNetworkPermutationSwitchRedo(String directory, int numNodes, int numPermutations) throws IOException{
		Graph<String, String> graph = createRandomConnectedGraphTest(numNodes);
		Set<String> pairs = getPairs(graph);
		saveSetToFile(directory+"/output/permutations/", "randomGraphForPermutation", pairs);
		
		String permuteEdgeListFile = "edgeListPermuted";
		List<Integer> numSwaps = permuteNetworkSwitchRedo(directory, permuteEdgeListFile, numPermutations, graph);
	}
	
	private List<Integer> permuteNetworkSwitchRedo(String directory, String fileName, int n, Graph<String, String> graph) throws IOException{
		System.out.println("---------- Creating permuted Networks ----------\t redo");
		
		int edgeSwapConstant = 10;//115;
		int numEdges = graph.getEdgeCount();
		int numSwap = edgeSwapConstant * numEdges;
		int windowThreshold = 3;
		
		HashMap<String, Integer> degrees = getDegreeMap(graph);
		List<String> degreeKeys = new ArrayList<String>(degrees.keySet());
		Set<String> geneSet = getGeneGraphSet(graph);
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
			System.out.println("redo\t network " + i + "\\" +  n + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwaps: " + swaps.get(i-1) + "\tWindowThreshold: " + windowThreshold + "\tSwitch Algorithm Time: " + ((end1 - start1) / 1000) + " seconds");
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
		Set<String> pairs = getPairs(graph);
		saveSetToFile(directory, fileName, pairs);
		System.out.println("\t\t\t\t\t\t\t\t\t\t\t\tTotal component check time: " + (time/1000) + " seconds");
		return swapcount;
	}
/*	
	private int switchAlgorithmRedo2(String directory, String fileName, Graph<String, String> graph, int numSwap, int windowThreshold, HashMap<String,Integer> degrees, List<String> degreeKeys, List<Double> cdf) throws IOException{	
//		List<String> edgeList = new ArrayList<String>(graph.getEdges());
		
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		if (components.size() != 1){
			System.err.println("Provided graph is not connected");
			System.out.println(graph);
			System.exit(0);
		}
		if (graph.getVertexCount() < 4){
			System.err.println("Provided graph has less than 4 nodes");
			System.exit(0);
		}
		
		HashMap<String, Integer> geneDegreeMap = getDegreeMap(graph);
		List<String> genes = new ArrayList<String>(graph.getVertices());
		HashMap<String, Set<String>> geneNeighbors = new HashMap<String, Set<String>> ();
		for (String g: graph.getVertices()){
			Set<String> neighbors = new HashSet<String> (graph.getNeighborCount(g));
			geneNeighbors.put(g, neighbors);
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
					int ui = getRandomIndex(genes.size());
					int xi = getRandomIndex(genes.size());
					
					
					int vi = getRandomIndex(genes.size());
					
					
					int yi = getRandomIndex(genes.size());
					if (ui==xi || ui==yi ||  vi==xi || vi==yi ){
						System.out.println("---samei:\t" + ui +"\t"+ vi +"\t"+ xi +"\t"+ yi);
						continue;
					}
					String u = genes.get(ui);
					List<String> uList = new ArrayList<String>(graph.getNeighbors(u));
					String v = uList.get(getRandomIndex(uList.size()));
					
					
					//pick 2 random edges
					List<Integer> randomEdgeIndex = get2RandomEdgeIndices(edgeList.size());
					String edge1 = edgeList.get(randomEdgeIndex.get(0));
					String edge2 = edgeList.get(randomEdgeIndex.get(1));
					//if same edges, skip
					if (edge1 == edge2)
						continue;
//					edgeList.remove(randomEdgeIndex.get(0));
//					edgeList.remove(randomEdgeIndex.get(1));
					//get 4 different nodes from the 2 edges
					String[] nodes1 = edge1.split("\t");
					String[] nodes2 = edge2.split("\t");
					String u = nodes1[0];
					String v = nodes1[1];
					String x = nodes2[0];
					String y = nodes2[1];
					//if any nodes are the same, skip
					if (u==x || u==y ||  v==x || v==y ){
						System.out.println("---same:\t" + u +"\t"+ v +"\t"+ x +"\t"+ y);
						continue;
					}
					//swap edges 
					graph.removeEdge(edge1);
					graph.removeEdge(edge2);
					graph.addEdge(u + "\t" + x, u, x);
					graph.addEdge(v + "\t" + y, v, y);		
					swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
					swapcount += 1;
					n +=1;
					//if graph is connected, increase window size
					components = wcc.transform(graph);
					if (components.size() == 1)
						wcount += 1; 
					else {
						graph.removeEdge(u + "\t" + x);
						graph.removeEdge(v + "\t" + y);
						graph.addEdge(u + "\t" + v, u, v);
						graph.addEdge(x + "\t" + y, x, y);
//						edgeList.add(edge1);
//						edgeList.add(edge2);
						swapcount -=1; 
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
					List<Integer> randomEdgeIndex = get2RandomEdgeIndices(edgeList.size());
					String edge1 = edgeList.get(randomEdgeIndex.get(0));
					String edge2 = edgeList.get(randomEdgeIndex.get(1));
					//if same edges, skip
					if (edge1 == edge2)
						continue;
//					edgeList.remove(randomEdgeIndex);
					//get 4 different nodes from the 2 edges
					String[] nodes1 = edge1.split("\t");
					String[] nodes2 = edge2.split("\t");
					String u = nodes1[0];
					String v = nodes1[1];
					String x = nodes2[0];
					String y = nodes2[1];
					//if any nodes are the same, skip
					if (u==x || u==y ||  v==x || v==y ){
						System.out.println("---same:\t" + u +"\t"+ v +"\t"+ x +"\t"+ y);
						continue;
					}
					//swap edges 
					graph.removeEdge(edge1);
					graph.removeEdge(edge2);
					graph.addEdge(u + "\t" + x, u, x);
					graph.addEdge(v + "\t" + y, v, y);		
					swapped.add(u + "\t" + v + "\t" + x + "\t" + y);
					swapcount += 1;
					n +=1;
					wcount += 1;					
					//if graph is connected, increase window size
					components = wcc.transform(graph);
					if (components.size() == 1)
						window += 1; 
					else {
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
//							edgeList.add(edge1);
//							edgeList.add(edge2);
							swapcount -= 1;
						}
						window = (int) Math.ceil(window/2); 
					}
				}
			}	
		}
		//Save the interaction pairs in the permuted network and return the number of swaps
		Set<String> pairs = getPairs(graph);
		saveSetToFile(directory, fileName, pairs);
		return swapcount; 
	}
*/	
	
	
	
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
	 * Get random index of list 
	 * @param max	size of list
	 * @return		random index
	 */
	private int getRandomIndex(int max){ 
		Random randomGen = new Random();
		int  random = randomGen.nextInt(max);
		return random; 
	}

	/**
	 * Get the random index from a list containing neighboring nodes
	 *		*note: used in the switchAlgorithm() to obtain random target nodes for swaps
	 * @param neighborList	list for obtaining random index from
	 * @return				random index of list
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
	 * Get a degree hash map of a graph using gene as key and degree as value
	 * 		-- based on python networkx's Graph.degree()
	 * 			https://networkx.readthedocs.io/en/stable/_modules/networkx/classes/graph.html#Graph.degree
	 * @param graph		used to create hashmap
	 * @return			a map with gene as key and gene's degree as value
	 */
	private HashMap<String, Integer> getDegreeMap(Graph<String,String> graph){
		HashMap<String, Integer> degreeMap = new HashMap<String, Integer>();
		for (String v: graph.getVertices())
			degreeMap.put(v, graph.getOutEdges(v).size());
		return degreeMap;
	}
	
	/**
	 * Calculate normalized cumulative distribution 
	 * 		-- based on python networkx's cumulative_distribution()
	 * 			https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/utils/random_sequence.html#cumulative_distribution
	 * @param distribution		used to create normalized cumulative distribution
	 * @return					normalized cumulative distribution
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
	 * 
	 * 			https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/utils/random_sequence.html#discrete_sequence
	 * @param n
	 * @param cumulativeDistribution
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
	 * Calculate where to insert a sequence so that the cumulative distribution list will still be sorted.
	 * 		-- based on python's bisect_left() 
	 * 			https://hg.python.org/cpython/file/2.7/Lib/bisect.py
	 * @param cumulativeDistribution	
	 * @param seq
	 * @return 
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
	 * Creates a random network and use it to generate random permutations using the configuration model algorithm.  All graphs are saved.
	 * 		--note: this is for testing purposes to try the switch algorithm on different networks
	 * @param directory		directory to save files in
	 * @param numNodes		number of nodes a graph should have
	 * @throws IOException
	 */
	private void testRandomNetworkPermutationConfig(String directory, int numNodes, int numPermutations) throws IOException{
		Graph<String, String> graph = createRandomConnectedGraphTest(numNodes);
		Set<String> pairs = getPairs(graph);
		saveSetToFile(directory+"/output/permutations/", "randomGraphForPermutation", pairs);
		
		String permuteEdgeListFile = "edgeListPermuted";
		permuteNetworkConfig(directory, permuteEdgeListFile, numPermutations, graph);
	}
	/**
	 * Permute network n times and save each permutation in a file 
	 * @param directory
	 * @param fileName
	 * @param graph
	 * @return
	 * @throws IOException
	 */
	private void permuteNetworkConfig(String directory, String fileName, int n, Graph<String, String> graph) throws IOException{
		System.out.println(" ---------- Creating permuted Networks ---------- ");
		
		String newDir = directory+"/output/permutations/";
		List<String> timePermutations = new ArrayList<String>();
		for (int i=1; i<n+1; i++){
			Graph<String, String> tempGraph1 = graph;
			String newFile = fileName+"_"+i;
			
			long start1 = System.currentTimeMillis();
			configurationModelAlgorithm(newDir, newFile, tempGraph1);
			long end1 = System.currentTimeMillis();
			
			int numVertex = tempGraph1.getVertexCount();
			int numEdge = tempGraph1.getEdgeCount();
			System.out.println("\t network " + i + "\\" +  n + "\tVertices: " + numVertex + "\tEdges: " + numEdge + "\tSwitch Algorithm Time: " + ((end1 - start1) / 1000) + " seconds");
			timePermutations.add(i + "\t" + n + "\t" + numVertex + "\t" + numEdge  + "\t" + ((end1 - start1)/1000));
		}
		String permuteTimeFile = "permuteTime_" + graph.getVertexCount(); 
		savePermuteTimeListToFile(newDir, permuteTimeFile, timePermutations); 
	}
	
	/**
	 * Configuration Model algorithm for generating and saving a random network
	 * 		-- based on Python networkx's configuration_model()
	 * 			https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/generators/degree_seq.html#configuration_model
	 * @param directory	directory to save random graph file 
	 * @param fileName	name of file to store edgelists of random network
	 * @param graph		graph used to generate a random permutation
	 * @throws IOException 
	 */
	private void configurationModelAlgorithm(String directory, String fileName, Graph<String, String> graph) throws IOException{
		Set<String> geneSet = getGeneGraphSet(graph);
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> degreeList = getGeneDegreeList(graph, geneSet);
		
		//Check that there's even number of degrees so no node will be leftover
		double degSum = 0;
		for(int i=0; i<degreeList.size(); i++)
			degSum += degreeList.get(i);
		if (degSum%2 != 0){
			System.err.println("Sum of degrees must be even");
			System.exit(0);
		}

		//number of nodes in graph
		int N = degreeList.size();

		//create empty list with nodes
		Set<String> nopairs = new HashSet<String>();
		Graph<String, String> randomGraph = createGraph(geneSet, nopairs);	

		//Done, there are no edges to add
		if (N==0 || degSum==0)
			return ;
	
		//Build stublist. For each node in geneSet, 
		List<String> tempStublist = new ArrayList<String>(); 
		for (int n=0; n<geneList.size(); n++){
			for (int i=0; i<degreeList.get(n); i++)
				tempStublist.add(geneList.get(n));
		}
		
		//Shuffle stublist and assign pairs by remove 2 items at a time
		Collections.shuffle(tempStublist);
		Queue<String> stublist = new LinkedList<String>(tempStublist);
		while (!stublist.isEmpty()){
			String n1 = stublist.poll();
			String n2 = stublist.poll();		
			randomGraph.addEdge(n1 + "\t" + n2 , n1 , n2);
		}
		
		//Save the interaction pairs in the permuted network and return the number of swaps
		Set<String> pairs = getPairs(randomGraph);
		saveSetToFile(directory, fileName, pairs);
		return; 
	}
	
	
	
	/**
	 * Convert the largest component file into files suitable for the Python version of HotNet2
	 * 		*note: HotNet2 Python version requires a edge list file and a gene index file
	 * @param directory			directory where the largest component of the graph is stored as a file of gene pairs
	 * @param largestCompFile	name of the file containing interaction pairs from the largest component
	 * @param geneIndexFile		name of the file where the index and gene will be saved
	 * @param edgeListFile		name of the file where the gene pairs will be saved
	 * @throws IOException
	 */
	private void convertFileForPythonHotNet2(String directory, String largestCompFile, String geneIndexFile, String edgeListFile) throws IOException{
		Path largeCompFilePath = Paths.get(directory+"/output/", largestCompFile);

		//Get all genes and all interaction pairs from largest Component file
		Set<String> allGenes = getAllGenes(largeCompFilePath);
		Set<String> allPairs = getAllInteractionPairs(largeCompFilePath);	
		
		//Create graph and get the largest component (this checks to ensure file only has 1 component)
		Graph<String, String> allGenesGraph = createGraph(allGenes, allPairs);
		Graph<String, String> largestComponent = createLargestComponentGraph(allGenesGraph);
		Set<String> geneSet = getGeneGraphSet(largestComponent);
		Set<String> genePairs = getPairs(largestComponent);
		
		//Create gene index file and a hash map with the gene and index for use in creating the edge list file
		File geneFile = new File(directory + "/output/" + geneIndexFile);
		geneFile.getParentFile().mkdirs();
		if (!geneFile.exists())
			geneFile.createNewFile();
		String geneFilePath = directory + "/output/" + geneIndexFile;
		PrintWriter gPw = new PrintWriter(geneFilePath);
		int i=1;
		Map<String, Integer> geneMap = new HashMap<String, Integer>();
		for(String gs: geneSet){
			gPw.println(i + "\t" + gs);
			geneMap.put(gs, i);
			i++;
		}
		gPw.close();
		
		//Create edge list file
		File edgeFile = new File(directory + "/output/" + edgeListFile);
		edgeFile.getParentFile().mkdirs();
		if (!edgeFile.exists())
			edgeFile.createNewFile();
		String edgeFilePath = directory + "/output/" + edgeListFile;
		PrintWriter ePw = new PrintWriter(edgeFilePath);
		for(String gp: genePairs){
			String[] pair = gp.split("\t");
			String gene1 = pair[0];
			String gene2 = pair[1];
			ePw.println(geneMap.get(gene1) + "\t" + geneMap.get(gene2));
		}
		ePw.close();
	}
	
	/**
	 * Save 3 different files related to the ReactomeFI network:
	 * 		- all genes from the ReactomeFI network
	 * 		- all interaction pairs from the ReactomeFI network
	 * 		- all interaction pairs from the largest component in the ReactomeFI network
	 * 		*note: this can also check that the largest component created matches the largest component in a reference file 
	 * @param directory			directory ReactomeFI network file is located and place to save new files	
	 * @param fiFile			name of the ReactomeFI network file
	 * @param allGenesFile		name of the file for genes in the ReactomeFI network 
	 * @param allPairsFile		name of the file for all interaction pairs in the ReactomeFI network
	 * @param largestCompFile	name of the file to save interaction pairs from the largest component
	 * @throws IOException
	 */
	private void saveReactomeFIGraphFiles(String directory, String fiFile, String allGenesFile, String allPairsFile, String largestCompFile) throws IOException{
		Path fiFilePath = Paths.get(directory, fiFile);

		//Get all genes from interaction file and save
		Set<String> allGenes = getAllGenesReactome(fiFilePath);
		saveSetToFile(directory+"/output/", allGenesFile, allGenes);
		
		//Get all interaction pairs from interaction file and save
		Set<String> allPairs = getAllInteractionPairsReactome(fiFilePath);
		saveSetToFile(directory+"/output/", allPairsFile, allPairs);
		
		//Create graph, get largest component, and save the interaction pairs in the largest component
		Graph<String, String> allGenesGraph = createGraph(allGenes, allPairs);
		Graph<String, String> largestComponent = createLargestComponentGraph(allGenesGraph);
		Set<String> genePairs = getPairs(largestComponent);
		saveSetToFile(directory+"/output/", largestCompFile, genePairs);
		
		//Check that interaction pairs in largest component matches reference file
		String refLargestCompFile = "FIsInGene_031516_BigComp.txt";
		compareFiles(directory, refLargestCompFile, "/output/"+largestCompFile);
	}
	
	/**
	 * Create a graph of the whole ReactomeFI network
	 * @param directory		directory ReactomeFI network file is located and place to save new files
	 * @param fiFile		name of the ReactomeFI network file
	 * @return				a graph of the whole ReactomeFI network
	 * @throws IOException
	 */
	private Graph<String, String> createReactomeFIGraph(String directory, String fiFile) throws IOException{
		Path fiFilePath = Paths.get(directory, fiFile);
		Set<String> allGenes = getAllGenesReactome(fiFilePath);
		Set<String> allPairs = getAllInteractionPairsReactome(fiFilePath);
		Graph<String, String> allGenesGraph = createGraph(allGenes, allPairs);
		return allGenesGraph; 
	}
	
	/**
	 * Save the results from calculating the betweenness centrality for all nodes in a graph
	 * @param directory		directory to save the file in
	 * @param fileName		name of file to save the betweenness centrality results in
	 * @param graph			graph used for calculating betweenness centrality 
	 * @throws IOException 
	 */
	private void saveBetweennessCentrality(String directory, String fileName, Graph<String, String> graph) throws IOException{
		File file = new File(directory + "/output/" + fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + "/output/" + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		BetweennessCentrality<String, String> bc = new BetweennessCentrality<String, String>(graph);
		bc.setRemoveRankScoresOnFinalize(false);
		bc.evaluate();
		for (String vertex: graph.getVertices())
			pw.println(bc.getVertexRankScore(vertex) + "\t" + vertex);
		pw.close();
	}

	/**
	 * Get 5 quartile source proteins based results in the betweenness centrality file  
	 * @param directory		directory where the betweenness centrality score file is located
	 * @param fileName		name of the betweenness centrality score file
	 * @return				a set containing the genes from the minimum, 25% quartile, median, 75% quartile, and maximum betweenness centrality 
	 * @throws IOException
	 */
	private List<String> getSourceProteins(String directory, String fileName) throws IOException{
		File file = new File(directory + "/" + fileName);
		if (!file.exists()){
			System.err.println("betweennessScore.txt file with betweenness scores for all genes in largest network does not exist");
			System.exit(0);
		}

		//read in scores for finding genes in quartile
		Path path = Paths.get(directory, fileName);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);

		List<Double> bcScores = new ArrayList<Double>(); 
		Map<Double, List<String>> bcMap = new HashMap<Double, List<String>>();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			Double score = Double.parseDouble(row[0]);
			String gene = row[1];
			bcScores.add(score);
			List<String> geneList = bcMap.get(score);
			if (geneList == null){
				geneList = new ArrayList<String>();
				geneList.add(gene);
				bcMap.put(score, geneList);
			}
			else {
				geneList.add(gene);
				bcMap.put(score, geneList);
			}
		}
		br.close();
		Collections.sort(bcScores);
		
		List<Double> quartile = new ArrayList<Double>();
		List<String> quartileGenes = new ArrayList<String>();
		int min = 0;  
		int q25 = (int) Math.ceil(bcScores.size()*0.25);
		int med = (int) Math.ceil(bcScores.size()*0.50);
		int q75 = (int) Math.ceil(bcScores.size()*0.75);
		int max = (int) Math.ceil(bcScores.size()-1);
		quartile.add(bcScores.get(min));
		quartile.add(bcScores.get(q25));
		quartile.add(bcScores.get(med));
		quartile.add(bcScores.get(q75));
		quartile.add(bcScores.get(max));
		for (double q: quartile){
			List<String> geneList = bcMap.get(q);
			Integer size = geneList.size();
			Random randomGen = new Random();
			Integer random = randomGen.nextInt(size);
			quartileGenes.add(geneList.get(random));
		}
		System.out.println("Source Proteins: " + quartileGenes);
		return quartileGenes;
	}
	
	/**
	 * Get a ordered set of genes in the graph
	 * @param graph			graph used to obtain genes from 
	 * @return				a ordered set of genes from the graph
	 */
	private Set<String> getGeneGraphSet (Graph<String, String> graph){
		Set<String> geneSet = new TreeSet<String>();
		for (String vertex: graph.getVertices())
			geneSet.add(vertex);
		return geneSet;
	}
	
	/**
	 * Create Diffusion Matrix F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1)
	 * @param graph			graph used to create matrix
	 * @param geneSet		set of genes in the graph used to determine matrix order
	 * @param beta			fraction of own heat each gene retains
	 * @return				a diffusion matrix for HotNet2
	 */
	private RealMatrix createDiffusionMatrix (Graph<String, String> graph, Set<String> geneSet, double beta){
		RealMatrix m = createNormAdjMatrix(graph, geneSet);	
		m = m.scalarMultiply(1-beta);
		int dim = geneSet.size();
		RealMatrix identityMatrix = MatrixUtils.createRealIdentityMatrix(dim);
		m = identityMatrix.subtract(m);
		m = MatrixUtils.inverse(m);
		m = m.scalarMultiply(beta);
		return m;
	}
	
	/**
	 * Create a normalized adjacency matrix
	 * @param graph			graph used to generate adjacency matrix
	 * @param geneSet		set of genes used to determine the order and degree of adjacency matrix elements
	 * @return				a normalized adjacency matrix
	 */
	private RealMatrix createNormAdjMatrix(Graph<String, String> graph, Set<String> geneSet) {
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> degreeList = getGeneDegreeList(graph, geneSet);		
		int dim = geneList.size();
		RealMatrix m = MatrixUtils.createRealMatrix(dim, dim);
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				String elm1 = geneList.get(i);
				String elm2 = geneList.get(j);				
				if (graph.containsEdge(elm1 + "\t" + elm2) || graph.containsEdge(elm2 + "\t" + elm1)){ 	//CHECK this
					double value = 1.0/degreeList.get(i);
					m.setEntry(i,j,value);
				}
			}
		}
		return m; 
	}

	/**
	 * Get a list of degree for each gene in the provided set  
	 * @param graph			graph used to obtain degree from for each gene
	 * @param geneSet		set of genes from the graph used to determine the list order of gene degree  
	 * @return				a list containing the degree for each gene in the provided set
	 */
	private List<Integer> getGeneDegreeList(Graph<String,String> graph, Set<String> geneSet){
		List<Integer> degreeList = new ArrayList<Integer>();
		List<String> geneList = new ArrayList<String>(geneSet);
		for (String g: geneList){
			graph.getOutEdges(g);
			degreeList.add(graph.getOutEdges(g).size());
		}
		return degreeList;
	}
	
	/**
	 * Create Exchanged Heat Matrix E = F * diagonal heat Matrix
	 * @param diffusionMatrix	diffusion matrix created by createDiffusionMatrix()
	 * @param geneSet			set of genes in the graph used to determine matrix order
	 * @param heatScoreMap		HashMap with gene as key and gene's heat score as value
	 * @return					a exchanged heat matrix for HotNet2
	 */	
	private RealMatrix createExchangedHeatMatrix(RealMatrix diffusionMatrix, Set<String> geneSet, HashMap<String, Double> heatScoreMap) {
		int dim = diffusionMatrix.getColumnDimension();
		RealMatrix heatMatrix = MatrixUtils.createRealMatrix(dim, dim);
		List<String> geneList = new ArrayList<String>(geneSet);
		for (int i=0; i<dim; i++){
			String gene = geneList.get(i);
			if (heatScoreMap.get(gene) != null)
				heatMatrix.setEntry(i, i, heatScoreMap.get(gene));
		}
		heatMatrix = diffusionMatrix.multiply(heatMatrix);
		return heatMatrix; 
	}

	/**
	 * Store heat scores of genes from a file into a HashMap, only genes in the provided set are included. 
	 * @param directory		directory where the file is located
	 * @param heatScoreFile name of the file with no header containing genes and heat scores separated a space 
	 * @param geneSet		set of genes that should be included in result
	 * @return				a HashMap with the gene as key and it's heat score as value
	 * @throws IOException
	 */
	private HashMap<String, Double> getHeatScoreMap(String directory, String heatScoreFile, Set<String> geneSet) throws IOException {
		HashMap<String, Double> geneScores = new HashMap<String, Double>();
		Path filePath = Paths.get(directory, heatScoreFile);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);		
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene = row[0];
			if (geneSet.contains(gene)){
				Double score = Double.parseDouble(row[1]);
				geneScores.put(gene, score);
			}
		}
		br.close();
		return geneScores;
	}
	
	/**
	 * Creates a file containing random heat scores for each gene in the ReactomeFI network
	 * 			----- this is for testing purposes only.  
	 * @param directory		directory where files are located and will be stored
	 * @param fiFile		name of the ReactomeFI network file
	 * @param heatScoreFile	name of the individual gene heat score file
	 * @throws IOException
	 */
	private void createTempHeatScoreFile(String directory, String fiFile, String heatScoreFile) throws IOException{
		Path fiFilePath = Paths.get(directory, fiFile);
		Set<String> allGenes = getAllGenesReactome(fiFilePath);
		String filePath = directory + "/" + heatScoreFile;
		PrintWriter pw = new PrintWriter(filePath);
		Random randomGen = new Random();
		for(String g: allGenes){
			Double random = randomGen.nextDouble();
			pw.println(g + " " + random);
		}
		pw.close();
	}
 	
	/** Identify Hot Subnetworks by removing values in exchanged heat matrix below the provided threshold 
	 * 		-- modified walkInOptimizedOrder to replace values below the threshold as 0
	 * @param exchangedHeatMatrix	exchanged heat matrix	
	 */
	private RealMatrix identifyHotSubnetworks(RealMatrix exchangedHeatMatrix, double edgeWeight){
		final double weight= edgeWeight;
		exchangedHeatMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
			public double visit(int row, int column, double value){	
				if (value <= weight)
					return 0;	
				return value; 
			}
		});	
		return exchangedHeatMatrix;
	}
	
	
	
//**************************************************** CHANGE THIS	
	/**
	 *  Obtain a set of subnetwork genes
	 * @param matrix	matrix to extract subnetworks from
	 * @param geneSet	set of genes used to determine matrix order
	 * @return			a set containing sets of genes that form subnetworks
	 */
	private Set<Set<String>> obtainSubnetworkSet (RealMatrix matrix, Set<String> geneSet){
		Graph<String, String> graph = covertMatrixToGraph(matrix, geneSet);

		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		System.out.println("There was " + components.size() + " subnetwork(s) identified");
	
		List<Set<String>> componentList = new ArrayList<Set<String>>(components);	
		for (int i=0; i<components.size(); i++){
			List<String> component = new ArrayList(componentList.get(i));
			String gene = component.get(0);
			
			graph.getOutEdges(gene);//Edges for a specific component
			//********** save each subnetwork into a edge list file?
		}
		
		
		
		return components;
	}
	
	/**
	 * Converts the provided RealMatrix into a graph 
	 * @param matrix	matrix to convert into graph
	 * @param geneSet	set of genes used to determine matrix order
	 * @return			a graph generated from the matrix
	 */
	private Graph<String, String> covertMatrixToGraph (RealMatrix matrix, Set<String> geneSet){
		List<String> geneList = new ArrayList<String>(geneSet);
		Set<String> genes = new TreeSet<String>();
		Set<String> pairs = new TreeSet<String>();
		int dim = matrix.getRowDimension();
		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				double value = matrix.getEntry(i, j);
				if (value != 0){
					String gene1 = geneList.get(i);
					String gene2 = geneList.get(j);
					genes.add(gene1);
					genes.add(gene2);
					pairs.add(gene1 + "\t" + gene2);
				}
			}
		}
		Graph<String, String> graph = createGraph(genes, pairs); 
		return graph; 
	}
	
	
	/**
	 * Save specified matrix into a tab delimited file
	 * @param directory		directory to store file
	 * @param fileName		name of file to store matrix contents
	 * @param matrix		matrix to be stored
	 * @throws IOException
	 */
	private void saveMatrix(String directory, String fileName, RealMatrix matrix) throws IOException{
		File file = new File(directory+"/output/"+ fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + "/output/" + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		int matrixDim = matrix.getRowDimension();
		for (int i=0; i<matrixDim; i++){
			double[] rowValues = matrix.getRow(i);
			for (double r: rowValues)
				pw.print(r + "\t");
			pw.println();
		}
		pw.close();
	}
	

	/**
	 * Select the beta parameter to assign an amount of heat retained by each gene for creating the diffusion matrix
	 * @param directory				directory of file
	 * @param fiFile
	 * @param betweennessScoreFile
	 * @param influenceFile
	 * @throws IOException
	 */
	private void selectBeta(String directory, String fiFile, String betweennessScoreFile, String influenceFile) throws IOException{
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = createLargestComponentGraph(allGenesGraph);

		//Calculate and save betweenness centrality for all genes in largest component
//		saveBetweennessCentrality(directory, betweennessScoreFile, largestComponent);
		
		//Get 5 source proteins from betweennness centrality scores
		String newDirectory = directory + "/output/"; 
		List<String> sourceProteins = getSourceProteins(newDirectory, betweennessScoreFile);

/*		//Get TP53 as a source protein
		List<String> sourceProteins =  new ArrayList<String>();
		sourceProteins.add("TP53");
*/		
		//Generate 20 diffusion matrices with different beta ranging from: 0.05, 0.10, 0.15,..., 1.00
		Set<String> geneSet = getGeneGraphSet(largestComponent);
		for(int i=1; i<21; i++){
			BigDecimal tempBeta = new BigDecimal("0.05");
			tempBeta = tempBeta.multiply(new BigDecimal(i));
			double beta = tempBeta.doubleValue();
			System.out.println("----beta: "+ beta);
			
			long start1 = System.currentTimeMillis();
			RealMatrix diffusionMatrix = createDiffusionMatrix(largestComponent, geneSet, beta);
			long end1 = System.currentTimeMillis();
			System.out.println("\tDiffusion Matrix Time Taken: " + ((end1 - start1) / 1000) + " seconds");
//			String newDirectory = directory +"/output/inflectionPoint"; 
//			String newFileName = tempBeta + "_diffusion.txt";
//			saveMatrix(newDirectory, newFileName, diffusionMatrix);

			saveInfluenceGeneCountRange(newDirectory, influenceFile, largestComponent, sourceProteins, diffusionMatrix, geneSet, tempBeta);			
		}
	
	}

	
	private void selectBetaForIrefindex(String directory, String fiFile, String betweennessScoreFile, String influenceFile) throws IOException{
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = createLargestComponentGraph(allGenesGraph);
		
		//Get TP53 source protein from iref_edge_list
		List<String> sourceProteins =  new ArrayList<String>();
		sourceProteins.add("10922");
//		sourceProteins.add("6911"); //source protein for iref_edge_list_temp
		
		//Generate diffusion matrix with beta: 0.45
		Set<String> geneSet = getGeneGraphSet(largestComponent);
		for(int i=1; i<2; i++){
			BigDecimal tempBeta = new BigDecimal("0.45");
			tempBeta = tempBeta.multiply(new BigDecimal(i));
			double beta = tempBeta.doubleValue();
			System.out.println("----beta: "+ beta);
			
			long start1 = System.currentTimeMillis();
			RealMatrix diffusionMatrix = createDiffusionMatrix(largestComponent, geneSet, beta);
			long end1 = System.currentTimeMillis();
			System.out.println("\tDiffusion Matrix Time Taken: " + ((end1 - start1) / 1000) + "seconds");

			saveInfluenceGeneCountRange(directory, influenceFile, largestComponent, sourceProteins, diffusionMatrix, geneSet, tempBeta);			
		}
	}

	
	/**
	 * Save the number of proteins with at least a diffusion value of theta influence for the 3 categories below related to the source proteins:
	 * 		- direct neighbors (source protein's direct neighbor)
	 * 		- secondary neighbors (source protein's secondary neighbor.  nodes distance 2 away from source protein)
	 * 		- all genes (all genes in network)
	 * @param directory			directory containing file
	 * @param fileName			name of file to save number of genes with at least theta influence for the 3 categories in
	 * @param graph				a graph used to find neighbors 
	 * @param sourceProteins	a list of 5 source proteins obtained from getSourceProteins()
	 * @param diffusionMatrix	a diffusion matrix who's contents are compared with the theta influence value
	 * @param geneSet			a set of genes that determines matrix ordering 
	 * @param beta				a parameter that determines the fraction of own heat each gene retains
	 * @throws IOException 
	 */
	private void saveInfluenceGeneCountRange(String directory, String fileName, Graph<String,String> graph, List<String> sourceProteins, RealMatrix diffusionMatrix, Set<String> geneSet, BigDecimal tempBeta) throws IOException{
		Set<String> influenceSet = new TreeSet<String>();
		for (int i=0; i<1001; i++){
			BigDecimal tempPt = new BigDecimal("0.0001");
			tempPt = tempPt.multiply(new BigDecimal(i));
			double influence = tempPt.doubleValue();
			int directNeighborsNum = 0;
			int secondaryNeighborsNum = 0;
			int allGenesNum = 0;
			
			System.out.println(sourceProteins + "\n" + graph + "\n" + diffusionMatrix + "\n" + directory + "\n" + fileName + "\n---------");
			for(String sp: sourceProteins){
				//Obtain the source protein's direct neighbors, secondary neighbors, and all proteins in the network 
				Set<String> directNeighbors = new TreeSet<String>();
				Set<String> secondaryNeighbors = new TreeSet<String>();
				Set<String> allGenes = new TreeSet<String>();
				directNeighbors.addAll(graph.getNeighbors(sp));
				for (String d: directNeighbors)
					secondaryNeighbors.addAll(graph.getNeighbors(d));
				secondaryNeighbors.remove(sp);
				secondaryNeighbors.removeAll(directNeighbors);
				allGenes.addAll(graph.getVertices());
				
				//Get distribution of proteins for the 3 influence point value results
				directNeighborsNum += calculateInfluenceQuantity(diffusionMatrix, geneSet, directNeighbors, sp, influence);
				secondaryNeighborsNum += calculateInfluenceQuantity(diffusionMatrix, geneSet, secondaryNeighbors, sp, influence);
				allGenesNum += calculateInfluenceQuantity(diffusionMatrix, geneSet, allGenes, sp, influence);
				
				//Store the 3 different distribution of proteins
				influenceSet.add(tempBeta +  "\t" + tempPt + "\t" + sp + "\t" + directNeighborsNum + "\t" + secondaryNeighborsNum + "\t" + allGenesNum);
			}
			String newDirectory = directory +"/inflectionPoint/"; 
			String newFileName = tempBeta + "_" + fileName;
			saveSetToFile(newDirectory, newFileName, influenceSet);
		}
	}
	
	/**
	 * Calculate the number of proteins greater than or equal to the theta influence point
	 * @param diffusionMatrix		diffusion matrix obtained from createDiffusionMatrix() 
	 * @param geneSet				a set of genes that determines matrix ordering
	 * @param specificGeneSet		a set of genes who must be found in a matrix
	 * @param sourceProtein			the source protein used to find matrix value for comparison against the theta influence
	 * @param thetaInfluence		used to determine number of genes with at least this value
	 * @return						the number of genes greater than or equal to the theta influence
	 */
	private int calculateInfluenceQuantity(RealMatrix diffusionMatrix, Set<String> geneSet, Set<String> specificGeneSet, String sourceProtein, double thetaInfluence){
		int quantity = 0; 
		List<Integer> proteinList = getSpecificGeneIndexInSet(geneSet, specificGeneSet, sourceProtein);
		for (Integer p: proteinList){
			if (proteinList.get(0) != p){
				if (thetaInfluence <= diffusionMatrix.getEntry(proteinList.get(0),p))
					quantity+=1;
			}
		}	
		return quantity; 
	}
	
	/**
	 * Get the indices of specific genes within a matrix
	 * 		*note: only used in selectBeta() and matrix is ordered by geneSet
	 * @param geneSet			a set of genes that determines matrix ordering 
	 * @param specificGeneSet	a set of genes whose indices must be found in a matrix 
	 * @param sourceProtein		the source protein used to find matrix value for comparison against the theta influence
	 * @return					a list of matrix indicies needed for calculateInfleuenceQuantity()
	 */
	private List<Integer> getSpecificGeneIndexInSet(Set<String> geneSet, Set<String> specificGeneSet, String sourceProtein) {
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> indexList = new ArrayList<Integer>();
		indexList.add(geneList.indexOf(sourceProtein));
		for (String s: specificGeneSet){
			indexList.add(geneList.indexOf(s));}
		return indexList;
	}
	
	

	private void selectDelta(){
		Path currentPath = Paths.get("");
		
	}
	

	
	
	
	/**
	 * Get all genes from network text file without a header
	 * @param path			path to file location
	 * @return				a set containing all genes
	 * @throws IOException
	 */
	private Set<String> getAllGenes(Path path) throws IOException{
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			genes.add(gene1);
			String gene2 = row[1];
			genes.add(gene2);
		}
		br.close();
		return genes;
	}
	
	/**
	 * Get all interaction pairs from network text file without a header
	 * @param path			path to the file location
	 * @return				a set containing all interaction pairs
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairs(Path path) throws IOException{
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			String gene2 = row[1];
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}


	/**
	 * Get all genes from the Reactome FI network text file
	 * @param path			path to file location
	 * @return				a set containing all genes
	 * @throws IOException
	 */
	private Set<String> getAllGenesReactome(Path path) throws IOException{
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			genes.add(gene1);
			String gene2 = row[1];
			genes.add(gene2);
		}
		br.close();
		return genes;
	}
	
	/**
	 * Get all interaction pairs from the Reactome FI network text file
	 * @param path			path to the file location
	 * @return				a set containing all interaction pairs
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairsReactome(Path path) throws IOException{
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			String gene2 = row[1];
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}

	/**
	 * Get all genes from any network text file with specified delimiter and header
	 * @param path			path to file location
	 * @param delimiter		delimiter that separates columns in file
	 * @return				a set containing all genes
	 * @throws IOException
	 */
	private Set<String> getAllGenes(Path path, String delimiter) throws IOException{
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene1 = row[0];
			genes.add(gene1);
			String gene2 = row[1];
			genes.add(gene2);
		}
		br.close();
		return genes;
	}
	
	/**
	 * Get all interaction pairs from any network text file with specified delimiter and header 
	 * @param path			path to the file location
	 * @param delimiter		delimiter that separates columns in file
	 * @return				a set containing all interaction pairs
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairs(Path path, String delimiter) throws IOException{
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene1 = row[0];
			String gene2 = row[1];
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}
	
	
	/**
	 * Save a set into a file
	 * @param directory		directory to save the file in	
	 * @param fileName		name of file to save set in
	 * @param setSave		a set to be saved 
	 * @throws IOException 
	 */
	private void saveSetToFile(String directory, String fileName, Set<String> setSave) throws IOException{
		File file = new File(directory + fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		for(String s: setSave){
			pw.println(s);
		}
		pw.close();
	}
	
	/**
	 * Create a graph 
	 * @param genes			set of genes in the graph
	 * @param pairs			set of interaction pairs in the graph
	 * @return				a undirected graph 
	 */
	private Graph<String,String> createGraph (Set<String> genes, Set<String> pairs){
		Graph<String, String> graph = new UndirectedSparseGraph<String, String>();	
		for (String g: genes)
			graph.addVertex(g);
		if(!pairs.isEmpty()){
			for (String p: pairs){ 
				String[] row = p.split("\t");
				String gene1 = row[0];
				String gene2 = row[1];
				graph.addEdge(p ,gene1, gene2); 
			}			
		}
		return graph;
	}
	
	/**
	 * Get the largest component in a graph
	 * @param graph			graph used in finding the largest component
	 * @return				a ordered set of genes in the largest component
	 */
	private Set<String> getLargestComponent(Graph<String, String> graph){
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		Set<String> largestComponent = new HashSet<String>();
		for (Set<String> c: components){	
			if(c.size() > largestComponent.size()){
				largestComponent = c; 
			}
		}
		return largestComponent; 
	}
	
	/**
	 * Create a graph of only the largest component
	 * @param graph			graph of whole network
	 * @return				a graph of only the largest component
	 */
	private Graph<String, String> createLargestComponentGraph (Graph<String, String> graph){
		Set<String> genes = getLargestComponent(graph);
		Graph<String, String> largestComponentGraph = FilterUtils.createInducedSubgraph(genes,graph);	
		return largestComponentGraph; 
	}
	
	/**
	 * Get interaction pairs from a graph
	 * @param graph		graph for extracting edges
	 * @return			a set containing the graph's interaction pairs   
	 */
	private Set<String> getPairs(Graph<String, String> graph){
		Set<String> pairs = new HashSet<String>();
		for (String v: graph.getVertices()){
			for (String e: graph.getEdges()){
				//note: graph.Endpoints(e) stores edges as Pair<L,R> (ex. <gene1, gene2>), so the first "<", ", ", and the last ">" must be removed
				String tempString = graph.getEndpoints(e).toString().replaceFirst("<", "").replace(", ", "\t").replaceFirst("[>]$","");	
				pairs.add(tempString);
			}
			break;
		}
		return pairs; 
	}
	
	/**
	 * Get file contents and store into a set
	 * @param directory		directory containing file	
	 * @param fileName		name of file to read in
	 * @return				a ordered set of file line contents
	 * @throws IOException
	 */
	private Set<String> getFileContents (String directory, String fileName) throws IOException{
		Set<String> contents = new TreeSet<String>();
		Path path = Paths.get(directory, fileName);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		for (String line = null; (line = br.readLine()) != null;){
			contents.add(line);
		}
		br.close();
		return contents;
	}
	
	/**
	 * Compare 2 different files for similarities and differences 
	 * @param directory		directory containing files
	 * @param file1			name of the first file to read in
	 * @param file2			name of the second file to read in
	 * @throws IOException
	 */
	private void compareFiles (String directory, String file1, String file2) throws IOException{
		Set<String> content1 = getFileContents(directory, file1);
		Set<String> content2 = getFileContents(directory, file2);
		if (content1.equals(content2) == true)
			System.out.println(file1 + " and " + file2 + "have the same contents");
		else{
			System.out.println("The following are in " + file1 + " but not in " + file2 +":");
			Set<String> temp = content1; 
			temp.removeAll(content2);
			Integer i = 1;
			for (String t: temp){
				System.out.println(i + ".\t" + t);
				i++; 
			}	
		}
	}

	/**
	 * Get genes with missing heat scores
	 * @param directory
	 * @param heatScoreFile
	 * @param graph
	 * @throws IOException 
	 */
	private void getMissingHeatGenes(String directory, String heatScoreFile, Graph<String,String> graph) throws IOException{
		String geneMissingHeatFile = "geneMissingHeat.txt";
		Set<String> geneSet = getGeneGraphSet(graph);
		Set<String> geneMissingHeat = getMissingHeatScoreGenes(directory, heatScoreFile, geneSet);
		saveSetToFile(directory+"/output/", geneMissingHeatFile, geneMissingHeat);
	}

	/**
	 * Compare heat score file genes with genes from a set
	 * @param directory		directory where the file is located
	 * @param heatScoreFile name of the file containing genes and heat scores separated a tab
	 * @param geneSet		set of genes that should be included in result
	 * @return				a set of genes and scores ordered by gene name
	 * @throws IOException
	 */
	private Set<String> getMissingHeatScoreGenes(String directory, String heatScoreFile, Set<String> geneSet) throws IOException {
		Set<String> heatScoreGenes = new HashSet<String>();
		Path filePath = Paths.get(directory, heatScoreFile);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene = row[0];
			heatScoreGenes.add(gene);
		}
		br.close();
		geneSet.removeAll(heatScoreGenes);
		return geneSet;
	}
	
	/**
	 * Save a list with diffusion times into a file
	 * 		---*note: 7 times from list will be saved into a file line in the following order:
	 * 				 normalized, scalar, identity, subtract, inverse, scalar, diffusion
	 * @param directory		directory to save the file in	
	 * @param fileName		name of file to save list in
	 * @param listSave		a list to be saved 
	 * @throws IOException 
	 */
	
	/**
	 * Save a list with permutation times into a file
	 * ---*note: 6 times from list will be saved into a file line in the following order:
	 * 				 permutation, numPermutes, numVertex, numEdge, networkSize, time(seconds)
	 * @param directory		directory to save the file in	
	 * @param fileName		name of file to save list in
	 * @param listSave		a list with permutation times to be saved
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
	
	
	
	
	/**
	 * Create Diffusion Matrix F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1)
	 * @param graph			graph used to create matrix
	 * @param geneSet		set of genes in the graph used to determine matrix order
	 * @param beta			fraction of own heat each gene retains
	 * @return				a diffusion matrix for HotNet2
	 */
	private PrimitiveMatrix createDiffusionMatrixOJA(Graph<String, String> graph, Set<String> geneSet, double beta){
		PrimitiveMatrix identityMatrix = PrimitiveMatrix.FACTORY.makeEye(geneSet.size(),geneSet.size());
		double[][] createdAdjMatrix = createNormAdjMatrixOJA(graph, geneSet);
        PrimitiveMatrix adjMatrix = PrimitiveMatrix.FACTORY.rows(createdAdjMatrix);
		PrimitiveMatrix m = adjMatrix.multiply(1-beta);
		m = identityMatrix.subtract(m);
		m = m.invert();
		m = m.multiply(beta);
		return m; 
	}
	
	/**
	 * Create a normalized adjacency matrix
	 * @param graph			graph used to generate adjacency matrix
	 * @param geneSet		set of genes used to determine the order and degree of adjacency matrix elements
	 * @return				a normalized adjacency matrix
	 */
	private double[][] createNormAdjMatrixOJA(Graph<String, String> graph, Set<String> geneSet) {
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> degreeList = getGeneDegreeList(graph, geneSet);		
		int dim = geneList.size();
		double[][] m = new double [dim][dim];
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				String elm1 = geneList.get(i);
				String elm2 = geneList.get(j);				
				if (graph.containsEdge(elm1 + "\t" + elm2) || graph.containsEdge(elm2 + "\t" + elm1)){ 	//CHECK this
					double value = 1.0/degreeList.get(i);
					m[i][j] = value;
				}
			}
		}
		return m; 
	}
	
	private RealMatrix convertOJAToACM(PrimitiveMatrix ojaMatrix){
		int dim = (int) ojaMatrix.countRows();
		RealMatrix m = MatrixUtils.createRealMatrix(dim, dim);
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				double value = ojaMatrix.get(i, j);
				m.setEntry(i, j, value);
			}
		}
		return m; 
	}
}


