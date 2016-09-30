package hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.ojalgo.matrix.PrimitiveMatrix;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class DeltaSelection{
	
	@Test
	public void testReactomeDelta() throws IOException{
		Path currentPath = Paths.get("");
		String directory = currentPath.toAbsolutePath().toString();
		
		BigDecimal tempBeta = new BigDecimal("0.5"); //get beta from user --- should be what was obtained from BetaSelection for Reactome FI (0.25 or 0.30)
				
		selectDelta(directory, tempBeta);
		
	}
	
	@Test
	public void testIrefIndexDelta() throws IOException{
		Path currentPath = Paths.get("");
		String directory = currentPath.toAbsolutePath().toString();
		
		BigDecimal tempBeta = new BigDecimal("0.45");
	
		selectDeltaForIrefindex(directory, tempBeta);
	}
	
	/**
	 * Selects the delta parameter for the Reactome FI network
	 * <p>
	 * <b>Note:</b> results are saved in a textfile for processing in excel.
	 * @param directory - Directory of Reactome FI network file and place to save files. 
	 * @param tempBeta - Fraction of own heat each gene retains. Should be beta=0.25 or 0.30.
	 * @throws IOException
	 */
	public void selectDeltaWrapper(String directory, BigDecimal tempBeta) throws IOException{
		selectDelta(directory, tempBeta);
	}
	
	/**
	 * Selects the delta parameter for the Reactome FI network
	 * <p>
	 * <b>Note:</b> results are saved in a textfile for processing in excel.
	 * @param directory - Directory of Reactome FI network file and place to save files. 
	 * @param tempBeta - Fraction of own heat each gene retains. Should be beta=0.25 or 0.30.
	 * @throws IOException
	 */
	private void selectDelta(String directory, BigDecimal tempBeta) throws IOException{
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers 
		String geneIndexDirectory = directory + "/PythonPermutation/";
		String geneIndexFile = "geneIndexReactome.txt";
		Set<String> genes = getAllGenesPY(geneIndexDirectory, geneIndexFile, "\t");
		HashMap<String, String> geneIndexMap = getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, "\t");
		//Lists for storing delta values for different maximum component sizes
		List<Double> maxCompSize5 = new ArrayList<Double>();
		List<Double> maxCompSize10 = new ArrayList<Double>();
		List<Double> maxCompSize15 = new ArrayList<Double>();
		List<Double> maxCompSize20 = new ArrayList<Double>();
		//Read in all permutation files within a directory
		String edgeListDirectory = directory + "/PythonPermutation/permutations/";
		File folder = new File(edgeListDirectory);	
		File[] listOfFiles = folder.listFiles();
		for (File file : folder.listFiles()){
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
			String heatScoreFile = "heatScore.txt";
			GraphUtils gu = new GraphUtils();
			//Create graph from file and check there is only 1 component
			Set<String> pairs = getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap);
			Graph<String, String> largestComponent = gu.createGraphWrapper(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Store smallest delta value corresponding with maximum component sizes of 5, 10, 15, and 20
			HashMap<Integer, Double> compSizeToDeltaMap = getDeltaForDiffCompSizes(heatScoreDirectory, heatScoreFile, tempBeta, largestComponent);
			maxCompSize5.add(compSizeToDeltaMap.get(5));
			maxCompSize10.add(compSizeToDeltaMap.get(10));
			maxCompSize15.add(compSizeToDeltaMap.get(15));
			maxCompSize20.add(compSizeToDeltaMap.get(20));
		}
		//Save deltas for each maximum component size
		FileUtils fu = new FileUtils();
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
	 * Selects the delta parameter for the iRefIndex network.
	 * <p>
	 * <b>Note:</b> to compare implementation matches HotNet2 Supplementary Figure 25 for iRefIndex network with beta=0.45 and maxCompSize L_MAX=10.
	 * The 100 random permutations are from <a href="http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/iref9.tar">Raphael-Group Hotnet2's GitHub</a>.
	 * @param directory - Directory Directory of iRefIndex network file and place to save files. 
	 * @param tempBeta - Fraction of own heat each gene retains. Should be beta=0.45.
	 * @throws IOException
	 */
	public void selectDeltaForIrefindexWrapper(String directory, BigDecimal tempBeta) throws IOException{
		selectDeltaForIrefindex(directory, tempBeta);
	}
	
	/**
	 * Selects the delta parameter for the iRefIndex network.
	 * <p>
	 * <b>Note:</b> to compare implementation matches HotNet2 Supplementary Figure 25 for iRefIndex network with beta=0.45 and maxCompSize L_MAX=10.
	 * The 100 random permutations are from <a href="http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/iref9.tar">Raphael-Group Hotnet2's GitHub</a>.
	 * @param directory - Directory Directory of iRefIndex network file and place to save files. 
	 * @param tempBeta - Fraction of own heat each gene retains. Should be beta=0.45.
	 * @throws IOException
	 */
	private void selectDeltaForIrefindex(String directory, BigDecimal tempBeta) throws IOException{
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers 
		String geneIndexDirectory = directory + "/PythonPermutation_Irefindex/";
		String geneIndexFile = "iref_index_genes";
		Set<String> genes = getAllGenesPY(geneIndexDirectory, geneIndexFile, " ");
		HashMap<String, String> geneIndexMap = getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, " ");
		//List for storing delta values for the maximum component size 10
		List<Double> maxCompSize10 = new ArrayList<Double>();
		//Read in all permutation files within a directory
		String edgeListDirectory = directory + "/PythonPermutation_Irefindex/permutations/";
		File folder = new File(edgeListDirectory);	
		File[] listOfFiles = folder.listFiles();
		int iterations = 0;
		for (File file : folder.listFiles()){
			long start1 = System.currentTimeMillis();
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
			String heatScoreFile = "heatScore.txt";
			GraphUtils gu = new GraphUtils();
			//Create graph from file and check there is only 1 component
			Set<String> pairs = getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap);
			Graph<String, String> largestComponent = gu.createGraphWrapper(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Store smallest delta value corresponding with maximum component size 10
			HashMap<Integer, Double> compSizeToDeltaMap = getDeltaForDiffCompSizesIrefindex(heatScoreDirectory, heatScoreFile, tempBeta, largestComponent);
			maxCompSize10.add(compSizeToDeltaMap.get(10));
			long end1 = System.currentTimeMillis();
			System.out.println(iterations + "\t" + edgeListFile + "\t" + (end1-start1) + "\tms");
			iterations += 1;
		}
		//Save deltas for maximum component size 10
		FileUtils fu = new FileUtils();
		String deltaCompDirectory = directory + "/outcome/deltaSelection_Irefindex/";
		String size10File = "deltaMaxComp10_Irefindex.txt";
		Collections.sort(maxCompSize10);
		fu.saveListToFile(deltaCompDirectory, size10File, maxCompSize10);
		return;
	}
	
	/**
	 * Get a set of all genes in the Python HotNet2 gene index file.
	 * @param directory - Directory containing the gene index file provided to the Python version of HotNet2
	 * @param fileName - File name containing gene index.
	 * @return a set of genes obtained from the provided file. 
	 * @throws IOException
	 */
	private Set<String> getAllGenesPY(String directory, String fileName, String delimiter) throws IOException{
		Set<String> genes = new TreeSet<String>();
		Path filePath = Paths.get(directory, fileName);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);	
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene = row[1].replace("_", " "); //Python code skips whitespace, so provided file had spaces replaced with "_" and it must be reversed to retain original Java genes
			genes.add(gene);
		}
		br.close();
		return genes;
	}

	////// Should this be HashMap<Integer, String> instead? then would need to convert getAllInteractionPairsPY genes to Integers
	/**
	 * Get a HashMap of genes corresponding to the Python HotNet2 gene index file.
	 * <p>
	 * <b>Note:</b> HotNet2 Python version skips whitespace so only part of gene name is kept, to ensure Python code didn't change gene names, the spaces are replaced with "_". 
	 * @param directory - Directory containing the gene index file provided to the Python version of HotNet2. 
	 * @param fileName - File containing the index and gene provided to the Python version of HotNet2.
	 * @return a HashMap containing the gene index as key and gene name as value.
	 * @throws IOException
	 */
	private HashMap<String,String> getGeneIndexNamePY (String directory, String fileName, String delimiter) throws IOException{
		HashMap<String, String> geneIndex = new HashMap<String, String>();
		Path filePath = Paths.get(directory, fileName);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);	
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String index = row[0];
			String gene = row[1].replace("_", " "); //Python code skips whitespace, so provided file had spaces replaced with "_" and it must be reversed to retain original Java genes
			geneIndex.put(index, gene);
		}
		br.close();
		return geneIndex;
	}

	/**
	 * Get all interaction pairs from a Python HotNet2 generated edge list file and match genes with the corresponding gene name
	 * @param directory - Directory containing gene pairs file.
	 * @param fileName - File containing gene pairs created by the Python version of HotNet2.
	 * @param geneIndex - HashMap matching gene index to gene name from {@link #getGeneIndexNamePY(String, String)}.
	 * @return a set containing all interaction pairs within the provided file.
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairsPY(String directory, String fileName, HashMap<String,String> geneIndex) throws IOException{
		Path filePath = Paths.get(directory, fileName);
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene1 = row[0];
			String gene2 = row[1];
			gene1 = geneIndex.get(gene1);
			gene2 = geneIndex.get(gene2);
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}

	/**
	 * Get delta for different component sizes. 
	 * @param heatScoreDirectory - Directory of heat score file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param tempBeta - Fraction of own heat each gene retains. Should be beta=0.25 or 0.30.
	 * @param graph - graph used to identify subnetworks. 
	 * @return a HashMap with the maximum component size as key and delta as value.
	 * @throws IOException
	 */
	private HashMap<Integer, Double> getDeltaForDiffCompSizes(String heatScoreDirectory, String heatScoreFile, BigDecimal tempBeta, Graph<String,String> graph) throws IOException{
		HotNet2Matrices hn2m = new HotNet2Matrices();
		GraphUtils gu = new GraphUtils();
		double beta = tempBeta.doubleValue(); //get beta from user --- should be what was obtained from BetaSelection (0.25 or 0.30)
		Set<String> geneSet = gu.getGeneGraphSetWrapper(graph);
		PrimitiveMatrix tempF = hn2m.createDiffusionMatrixOJAWrapper(graph, geneSet, beta);
		RealMatrix F = hn2m.convertOJAToACMWrapper(tempF);	
	
		RealMatrix ExchangedHeatMatrix = hn2m.createExchangedHeatMatrixWrapper(heatScoreDirectory, heatScoreFile, F, geneSet);
		List<Integer> sizes = new ArrayList<Integer>();
		sizes.add(5);
		sizes.add(10);
		sizes.add(15);
		sizes.add(20);
		
		HashMap<Integer, Double> compSizeToDeltaMap = selectDeltaForDiffCompSizes(geneSet, ExchangedHeatMatrix, sizes, heatScoreDirectory);
		return compSizeToDeltaMap;
	}
	
	/**
	 * Get delta for component size 10 for the iRefIndex network.
	 * @param heatScoreDirectory - Directory of heat score file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param tempBeta - Fraction of own heat each gene retains. Should be beta=0.45.
	 * @param graph - graph used to identify subnetworks. 
	 * @return a HashMap with the maximum component size as key and delta as value.
	 * @throws IOException
	 */
	private HashMap<Integer, Double> getDeltaForDiffCompSizesIrefindex(String heatScoreDirectory, String heatScoreFile, BigDecimal tempBeta, Graph<String,String> graph) throws IOException{
		HotNet2Matrices hn2m = new HotNet2Matrices();
		GraphUtils gu = new GraphUtils();
		double beta = tempBeta.doubleValue();
		Set<String> geneSet = gu.getGeneGraphSetWrapper(graph);
		PrimitiveMatrix tempF = hn2m.createDiffusionMatrixOJAWrapper(graph, geneSet, beta);
		RealMatrix F = hn2m.convertOJAToACMWrapper(tempF);	
	
		RealMatrix ExchangedHeatMatrix = hn2m.createExchangedHeatMatrixWrapper(heatScoreDirectory, heatScoreFile, F, geneSet);
		List<Integer> sizes = new ArrayList<Integer>();
		sizes.add(10);
		
		HashMap<Integer, Double> compSizeToDeltaMap = selectDeltaForDiffCompSizes(geneSet, ExchangedHeatMatrix, sizes, heatScoreDirectory);
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
	private HashMap<Integer, Double> selectDeltaForDiffCompSizes(Set<String> geneSet, RealMatrix ExchangedHeatMatrix, List<Integer> sizeList, String directory) throws IOException{
		GraphUtils gu = new GraphUtils();
		HotNet2Matrices hn2m = new HotNet2Matrices();
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
				RealMatrix H = hn2m.identifyHotSubnetworksWrapper(ExchangedHeatMatrix, delta);
				int maxSubnetworkSize = obtainLargestSubnetworkSize(H, geneSet);
				if (visitedDelta.contains(delta)){
					compSizeToDeltaMap.put(size, delta);
					break;
				}
				else
					visitedDelta.add(delta);
				//Decrement or increment delta based on max component size
				if (maxSubnetworkSize > size){
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
		Graph<String, String> graph = gu.covertMatrixToDirectedGraphWrapper(matrix, geneSet);
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
		
}
