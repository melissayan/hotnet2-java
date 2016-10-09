package edu.ohsu.hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class DeltaSelection{
	private Path currentPath;
	private String directory;

	public DeltaSelection(){
		this.currentPath = Paths.get("");
		this.directory = currentPath.toAbsolutePath().toString();
	}
	
	@Test
	public void testReactomeDelta() throws IOException{
		double beta = 0.25; //get beta from user --- should be what was obtained from BetaSelection for Reactome FI (0.25 or 0.30)
		int numPermutation = 100;
		selectDeltaByCompSize(directory, beta, numPermutation);
	//	selectDeltaByHeatScore (directory, beta, numPermutation);
	}
	
	@Test
	public void testIrefIndexDelta() throws IOException{
		double beta = 0.45;
		int numPermutation = 2; 
		selectDeltaForIrefindexByCompSize(directory, beta, numPermutation);
	}
	
	/**
	 * Selects the delta parameter for the Reactome FI network based on random permutation networks.
	 * <p>
	 * <b>Note:</b> results are saved in a textfile for processing in excel.
	 * @param directory - Directory of Reactome FI network file and place to save files. 
	 * @param beta - Fraction of own heat each gene retains. Should be beta=0.25 or 0.30.
	 * @throws IOException
	 */
	public void selectDeltaByCompSize(String directory, double beta, int numPermutation) throws IOException{
		System.out.println("selectDelta Function");
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		String geneIndexDirectory = directory + "/PythonPermutation/";
		String geneIndexFile = "geneIndexReactome.txt";
		Set<String> genes = fu.getAllGenesPY(geneIndexDirectory, geneIndexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, "\t");
		//List for storing delta values for the maximum component size 10
		List<Double> maxCompSize5 = new ArrayList<Double>();
		List<Double> maxCompSize10 = new ArrayList<Double>();
		List<Double> maxCompSize15 = new ArrayList<Double>();
		List<Double> maxCompSize20 = new ArrayList<Double>();
		//Use random permutated networks to obtain delta
		String edgeListDirectory = directory + "/PythonPermutation/permutations/";
		File folder = new File(edgeListDirectory);
		File[] listOfFiles = folder.listFiles();
		int iterations = 0;
		int i = 0;
		System.out.println("Start reading files");
		for (File file: listOfFiles){
			long start1 = System.currentTimeMillis();
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
			String heatScoreFile = "heatScore.txt";
			//Create graph from file and check there is only 1 component
			Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, " ");
			Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				//throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
				System.out.println("Provided permuted graph " + edgeListFile + " is not connected, it has " + components.size() + " components.");
			}
			//Store smallest delta value corresponding with maximum component sizes of 5, 10, 15, and 20
			HashMap<Integer, Double> compSizeToDeltaMap = getByCompSizeMap(heatScoreDirectory, heatScoreFile, beta, largestComponent);
			maxCompSize5.add(compSizeToDeltaMap.get(5));
			maxCompSize10.add(compSizeToDeltaMap.get(10));
			maxCompSize15.add(compSizeToDeltaMap.get(15));
			maxCompSize20.add(compSizeToDeltaMap.get(20));
			long end1 = System.currentTimeMillis();
			System.out.println(iterations + "\t" + edgeListFile + "\t" + (end1-start1) + "\tms");
			iterations += 1;
			 if (++i > numPermutation)
				 break;
		}
		//Save deltas for each maximum component size
		String deltaCompDirectory = directory + "/outcome/deltaSelection/";
		String size5File = "deltaMaxComp5.txt";
		String size10File = "deltaMaxComp10.txt";
		String size15File = "deltaMaxComp15.txt";
		String size20File = "deltaMaxComp20.txt";
		Collections.sort(maxCompSize5);
		fu.saveListToFile(deltaCompDirectory, size5File, maxCompSize5);
		Collections.sort(maxCompSize10);
		fu.saveListToFile(deltaCompDirectory, size10File, maxCompSize10);
		Collections.sort(maxCompSize15);
		fu.saveListToFile(deltaCompDirectory, size15File, maxCompSize15);
		Collections.sort(maxCompSize20);
		fu.saveListToFile(deltaCompDirectory, size20File, maxCompSize20);
		return;
	}
	
	HashMap<Integer, Double> sample(String heatScoreDirectory, String heatScoreFile, double beta, Graph<String,String> graph) throws IOException{
		HashMap<Integer, Double> hm = new HashMap<Integer, Double>();
		hm.put(5, 0.2739485);
		hm.put(10, 0.2739485);
		hm.put(15, 0.2739485);
		hm.put(20, 0.2739485);
		return hm;
	}
	
	
	/**
	 * Selects the delta parameter for the iRefIndex network based on random permutation networks.
	 * <p>
	 * <b>Note:</b> to compare implementation matches HotNet2 Supplementary Figure 25 for iRefIndex network with beta=0.45 and maxCompSize L_MAX=10.
	 * The 100 random permutations are from <a href="http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/iref9.tar">Raphael-Group Hotnet2's GitHub</a>.
	 * @param directory - Directory Directory of iRefIndex network file and place to save files. 
	 * @param beta - Fraction of own heat each gene retains. Should be beta=0.45.
	 * @throws IOException
	 */
	public void selectDeltaForIrefindexByCompSize(String directory, double beta, int numPermutation) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils(); 
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		String geneIndexDirectory = directory + "/PythonPermutation_Irefindex/";
		String geneIndexFile = "iref_index_genes";
		Set<String> genes = fu.getAllGenesPY(geneIndexDirectory, geneIndexFile, " ");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, " ");
		//List for storing delta values for the maximum component size 10
		List<Double> maxCompSize10 = new ArrayList<Double>();
		//Use random permutated networks to obtain delta
		String edgeListDirectory = directory + "/PythonPermutation_Irefindex/permutations/";
		File folder = new File(edgeListDirectory);	
		File[] listOfFiles = folder.listFiles();
		int iterations = 0;
		int i = 0;
		for (File file: listOfFiles){
			long start1 = System.currentTimeMillis();
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
			String heatScoreFile = "heatScore.txt";
			//Create graph from file and check there is only 1 component
			Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, " ");
			Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				//throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
				System.out.println("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Store smallest delta value corresponding with maximum component size 10
			HashMap<Integer, Double> compSizeToDeltaMap = getByCompSizeMapIrefindex(heatScoreDirectory, heatScoreFile, beta, largestComponent);
			maxCompSize10.add(compSizeToDeltaMap.get(10));
			long end1 = System.currentTimeMillis();
			System.out.println(iterations + "\t" + edgeListFile + "\t" + (end1-start1) + "\tms");
			iterations += 1;
			 if (++i > numPermutation)
				 break;
		}
		//Save deltas for maximum component size 10
		String deltaCompDirectory = directory + "/outcome/deltaSelection_Irefindex/";
		String size10File = "deltaMaxComp10_Irefindex.txt";
		Collections.sort(maxCompSize10);
		fu.saveListToFile(deltaCompDirectory, size10File, maxCompSize10);
		return;
	}
	
	/**
	 * Get delta for different component sizes  based on random permutation networks. 
	 * @param heatScoreDirectory - Directory of heat score file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param beta - Fraction of own heat each gene retains. Should be beta=0.25 or 0.30.
	 * @param graph - graph used to identify subnetworks. 
	 * @return a HashMap with the maximum component size as key and delta as value.
	 * @throws IOException
	 */
	private HashMap<Integer, Double> getByCompSizeMap(String heatScoreDirectory, String heatScoreFile, double beta, Graph<String,String> graph) throws IOException{
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix(); 
		Set<String> geneSet = gu.getGeneGraphSet(graph);
		PrimitiveMatrix tempF = hn2m.createDiffusionMatrixOJA(graph, geneSet, beta);
		RealMatrix F = hn2m.convertOJAToACM(tempF);	
		RealMatrix ExchangedHeatMatrix = hn2m.createExchangedHeatMatrix(heatScoreDirectory, heatScoreFile, F, geneSet);
		List<Integer> sizes = new ArrayList<Integer>();
		sizes.add(5);
		sizes.add(10);
		sizes.add(15);
		sizes.add(20);
		HashMap<Integer, Double> compSizeToDeltaMap = selectDeltaForDiffCompSizes(geneSet, ExchangedHeatMatrix, sizes);
		return compSizeToDeltaMap;
	}
	
	/**
	 * Get delta for component size 10 for the iRefIndex network based on random permutation networks.
	 * @param heatScoreDirectory - Directory of heat score file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param beta - Fraction of own heat each gene retains. Should be beta=0.45.
	 * @param graph - graph used to identify subnetworks. 
	 * @return a HashMap with the maximum component size as key and delta as value.
	 * @throws IOException
	 */
	private HashMap<Integer, Double> getByCompSizeMapIrefindex(String heatScoreDirectory, String heatScoreFile, double beta, Graph<String,String> graph) throws IOException{
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix(); 
		Set<String> geneSet = gu.getGeneGraphSet(graph);
		PrimitiveMatrix tempF = hn2m.createDiffusionMatrixOJA(graph, geneSet, beta);
		RealMatrix F = hn2m.convertOJAToACM(tempF);
		RealMatrix ExchangedHeatMatrix = hn2m.createExchangedHeatMatrix(heatScoreDirectory, heatScoreFile, F, geneSet);
		List<Integer> sizes = new ArrayList<Integer>();
		sizes.add(10);
		HashMap<Integer, Double> compSizeToDeltaMap = selectDeltaForDiffCompSizes(geneSet, ExchangedHeatMatrix, sizes);
		return compSizeToDeltaMap;
	}
	
	/**
	 * Select delta for different provided component sizes.
	 * <p>
	 * <b>Note:</b> based on Python HotNet2's delta.py find_best_delta_by_largest_cc(). 
	 * @param geneSet - Set of genes used to determine matrix order. 
	 * @param ExchangedHeatMatrix - Exchanged heat matrix.
	 * @param sizeList - List of sizes for the largest component. 
	 * @param directory - Directory to store file.
	 * @return a HashMap  a HashMap with the largest component size as key and it's smallest delta as value.
	 * @throws IOException
	 */
	private HashMap<Integer, Double> selectDeltaForDiffCompSizes(Set<String> geneSet, RealMatrix ExchangedHeatMatrix, List<Integer> sizeList) throws IOException{
		HotNet2Matrix hn2m = new HotNet2Matrix(); 
		HashMap<Integer, Double> compSizeToDeltaMap = new HashMap<Integer, Double>();
		for (Integer size: sizeList){
			double delta = -1;
			ArrayList<Double> sortedEdgeWeights = new ArrayList<Double>(obtainUniqueExchangedHeatMatrixValues(ExchangedHeatMatrix));
			int index = (int) Math.ceil(sortedEdgeWeights.size()*0.99)-1;
			int left = 0;
			int right = sortedEdgeWeights.size();
			List<Double> visitedDelta = new ArrayList<Double>();
			while(visitedDelta.size() < 100){
				//Construct new matrix using new delta
				delta = sortedEdgeWeights.get(index);
				RealMatrix H = hn2m.identifyHotSubnetworks(ExchangedHeatMatrix, delta);
				int maxSubnetworkSize = obtainLargestSubnetworkSize(H, geneSet);
				if (visitedDelta.contains(delta)){
					compSizeToDeltaMap.put(size, delta);
					break;
				}
				else
					visitedDelta.add(delta);
				//Increment if delta too small or decrement if delta too big based on max component size
				if (size > maxSubnetworkSize){
					right = index;
					index -= Math.ceil((index-left)/2);
				} else {
					left = index;
					index += Math.ceil((right-index)/2); 
				}
			}
			if (compSizeToDeltaMap.get(size) == null)
				throw new IllegalArgumentException("no delta value for max component size: " + size);
		}
		return compSizeToDeltaMap;
	}
	
	/**
	 * Obtain largest subnetwork size from the provided matrix.  
	 * @param matrix - Matrix used to determine subnetwork sizes.
	 * @param geneSet - Set of genes used to determine matrix order. 
	 * @return the value of the largest component size.
	 */
	private int obtainLargestSubnetworkSize(RealMatrix matrix, Set<String> geneSet){
		GraphUtils gu = new GraphUtils();
		Graph<String, String> graph = gu.covertMatrixToDirectedGraph(matrix, geneSet);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		int largestSize = 0;
		for (Set<String> c: components){
			if(c.size() > largestSize)
				largestSize = c.size();
		}
		return largestSize;
	}

	/**
	 * Obtain a set of unique exchanged heat matrix values. 
	 * @param matrix - Matrix used to obtain element values.
	 * @return a sorted set of element values.
	 */
	private SortedSet<Double> obtainUniqueExchangedHeatMatrixValues (RealMatrix matrix){
		SortedSet<Double> weights = new TreeSet<Double>();
		int dim = matrix.getRowDimension();
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++)
				weights.add(matrix.getEntry(i, j));
		}
		return weights;
	}
	
	/**
	 * Selects the delta parameter for the Reactome FI network based on random heat scores.
	 * <p>
	 * <b>Note:</b> results are saved in a textfile for processing in excel.
	 * @param directory - Directory of Reactome FI network file and place to save files. 
	 * @param beta - Fraction of own heat each gene retains. Should be beta=0.25 or 0.30.
	 * @throws IOException
	 */
	public void selectDeltaByHeatScore (String directory, double beta, int numPermutation) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix(); 
		//Create largest component in graph 
		String fiFile = "FIsInGene_031516_with_annotations.txt";
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, fiFile);
		Graph<String, String> graph = gu.createLargestComponentGraph(allGenesGraph);
		Set<String> geneSet = gu.getGeneGraphSet(graph);
		//Lists for storing delta values for different maximum component sizes
		List<Double> maxCompSize5 = new ArrayList<Double>();
		List<Double> maxCompSize10 = new ArrayList<Double>();
		List<Double> maxCompSize15 = new ArrayList<Double>();
		List<Double> maxCompSize20 = new ArrayList<Double>();
		//Get original heat scores
		String heatScoreDirectory = directory+"";
		String heatScoreFile = "heatScore.txt";
		List<Double> scoreList = getHeatScoreListAll(heatScoreDirectory, heatScoreFile);
		//Use random heat scores to obtain delta
		PrimitiveMatrix tempF = hn2m.createDiffusionMatrixOJA(graph, geneSet, beta);
		RealMatrix F = hn2m.convertOJAToACM(tempF);
		for (int i=0; i<numPermutation; i++){
			HashMap<String, Double> randomHeatScoreMap = createRandomHeatScoreMap(geneSet, scoreList);
			RealMatrix ExchangedHeatMatrix = hn2m.createExchangedHeatMatrixUseMap(F, geneSet, randomHeatScoreMap);	
			List<Integer> sizes = new ArrayList<Integer>();
			sizes.add(5);
			sizes.add(10);
			sizes.add(15);
			sizes.add(20);
			//Store smallest delta value corresponding with maximum component sizes of 5, 10, 15, and 20
			HashMap<Integer, Double> compSizeToDeltaMap = selectDeltaForDiffCompSizes(geneSet, ExchangedHeatMatrix, sizes);
			maxCompSize5.add(compSizeToDeltaMap.get(5));
			maxCompSize10.add(compSizeToDeltaMap.get(10));
			maxCompSize15.add(compSizeToDeltaMap.get(15));
			maxCompSize20.add(compSizeToDeltaMap.get(20));
		}
		//Save deltas for each maximum component size
		String deltaCompDirectory = directory + "/outcome/deltaSelection/";
		String size5File = "deltaMaxComp5.txt";
		String size10File = "deltaMaxComp10.txt";
		String size15File = "deltaMaxComp15.txt";
		String size20File = "deltaMaxComp20.txt";
		Collections.sort(maxCompSize5);
		fu.saveListToFile(deltaCompDirectory, size5File, maxCompSize5);
		Collections.sort(maxCompSize10);
		fu.saveListToFile(deltaCompDirectory, size10File, maxCompSize10);
		Collections.sort(maxCompSize15);
		fu.saveListToFile(deltaCompDirectory, size15File, maxCompSize15);
		Collections.sort(maxCompSize20);
		fu.saveListToFile(deltaCompDirectory, size20File, maxCompSize20);
		return;
	}
	
	/**
	 * Creates a random heat score map using the provided gene set and score list. 
	 * @param geneSet - Set of genes to assign heat scores to
	 * @param scoreList - List of heat scores for randomly assigned to genes.
	 * @return a HashMap with gene as key and random heat score as value.
	 */
	private HashMap<String, Double> createRandomHeatScoreMap(Set<String> geneSet, List<Double> scoreList){
		HashMap<String, Double> randomHeatScoreMap = new HashMap<String, Double>();
		Random randomGen = new Random();
		for (String g: geneSet){
			int randomIdx = randomGen.nextInt(scoreList.size()-1);
			randomHeatScoreMap.put(g, scoreList.get(randomIdx));
		}
		return randomHeatScoreMap;
	}
	
	/**
	 * Stores heat scores of genes from a file into a List, only genes in the provided set are included. 
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param geneSet - Set of genes that should be included in result.
	 * @return a list containing heat scores from the provided gene set.
	 * @throws IOException
	 */
	private List<Double> getHeatScoreList(String directory, String heatScoreFile, Set geneSet) throws IOException {
		List<Double> scoreList = new ArrayList<Double>();
		Path filePath = Paths.get(directory, heatScoreFile);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);		
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene = row[0];
			if (geneSet.contains(gene)){
				Double score = Double.parseDouble(row[1]);
				scoreList.add(score);
			}
		}
		br.close();
		return scoreList;
	}
	
	/**
	 * Stores heat scores of genes from a file into a List. 
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param geneSet - Set of genes that should be included in result.
	 * @return a list containing heat scores from the provided gene set.
	 * @throws IOException
	 */
	private List<Double> getHeatScoreListAll(String directory, String heatScoreFile) throws IOException {
		List<Double> scoreList = new ArrayList<Double>();
		Path filePath = Paths.get(directory, heatScoreFile);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);		
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			scoreList.add(Double.parseDouble(row[1]));
		}
		br.close();
		return scoreList;
	}	
}
