package edu.ohsu.hotnet2;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.*;
import org.ojalgo.matrix.PrimitiveMatrix;

import cern.colt.Arrays;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class StatisticalSignificance{
	
	public StatisticalSignificance(){
		
	}
	
	@Test
		public void testPrototypeDelta() throws IOException{
			Path currentPath = Paths.get("");
			String directory = currentPath.toAbsolutePath().toString();
			
			double beta = 0.25;
			double delta = 0.014891578557906057;
			int numPermutation = 10;
			String dir =  directory+"/prototype/";
			String indexFile = "prototype_index_genes_java";
			String edgelistFile = "prototype_edgelist_java";
			String heatScoreFile = "prototype_heatScore_java";
			String dirOut = dir+"/output/";
			
			long start = System.currentTimeMillis();
			runHotNet2Prototype(dir, indexFile, edgelistFile, heatScoreFile, dirOut, beta, delta, numPermutation);
//			runHotNet2NoDeltaPrototype(dir, indexFile, edgelistFile, heatScoreFile, beta, numPermutation);
			long end = System.currentTimeMillis();
			System.out.println("total time: " + (end-start));
		}

	@Test
	public void testPrototype2Delta() throws IOException{
		Path currentPath = Paths.get("");
		String directory = currentPath.toAbsolutePath().toString();
		
		double beta = 0.25;
		double delta = 0.0001;
		int numPermutation = 10;
		String dir =  directory+"/prototype2/";
		String indexFile = "prototype2_index_java";
		String edgelistFile = "prototype2_edgelist_java";
		String heatScoreFile = "prototype2_heatScore_java";
		String dirOut = dir+"/output/";
		
		long start = System.currentTimeMillis();
//		runHotNet2Prototype(dir, indexFile, edgelistFile, heatScoreFile, dirOut, beta, delta, numPermutation);
		runHotNet2NoDeltaPrototype(dir, indexFile, edgelistFile, heatScoreFile, beta, numPermutation);
		long end = System.currentTimeMillis();
		System.out.println("total time: " + (end-start));
	}
	
	/**
	 * Run HotNet2 and get statistical significance of the observed subnetworks from Reactome FI network.
	 * @param directory	- Directory where Reactome FI network file and heat score file are located.
	 * @param fiFile - Reactome FI network file name.
	 * @param heatScoreFile - File name of genes and their heat score. 
	 * @param directoryOutput - Directory to store results in.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @param numPermutations - Number of permuted networks to use.
	 * @throws IOException
	 */
	public void runHotNet2(String directory, String fiFile, String heatScoreFile, String directoryOutput, double beta, double delta, int numPermutations) throws IOException{
		HashMap<Integer, Integer> observed = obtainSubnetworkSizeToCountSumMapReal(directory, fiFile, heatScoreFile, directoryOutput, beta, delta);
		HashMap<Integer, Integer> expected = obtainSubnetworkSizeToCountSumMapPermuted(directory, heatScoreFile, directoryOutput, beta, delta);
		SortedSet<Integer> sizeSet = new TreeSet<Integer>(observed.keySet());
		List<Integer> sizeList = new ArrayList<Integer>(sizeSet);
		//save results to file
		String sigDirectory = directoryOutput;
		String sigFile = "significance.txt";
		File file = new File(sigDirectory + sigFile);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = sigDirectory + sigFile;
		PrintWriter pw = new PrintWriter(filePath);
		pw.println("size\tobserved\texpected\tpvalue");
		
		//*************** STILL NEED TO FIX THIS
		for (int s: sizeList){
			double observedValue = observed.get(s);
			double expectedValue = expected.get(s)/numPermutations;
			double pvalue = observedValue/numPermutations;
			pw.println(s + "\t" + expectedValue + "\t" + observedValue + "\t" + pvalue);
		}
		pw.close();
	}
	
	/**
	 * Run HotNet2 and get statistical significance of the observed subnetworks from prototype network.
	 * @param directory - Directory where gene index file, edge list file, and heat score file are located. 
	 * @param indexFile - File name of index and the corresponding gene.
	 * @param edgelistFile - File name of gene pairs. 
	 * @param heatScoreFile - FIle name of genes and their scores.
	 * @param directoryOutput - Directory to store results in.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @param numPermutations - Number of permuted networks to use.
	 * @throws IOException
	 */
	public void runHotNet2Prototype(String directory, String indexFile, String edgelistFile, String heatScoreFile, String directoryOutput, double beta, double delta, int numPermutations) throws IOException{
		HashMap<Integer, Integer> observed = obtainSubnetworkSizeToCountSumMapRealPrototype(directory, indexFile, edgelistFile, heatScoreFile, directoryOutput, beta, delta);
		HashMap<Integer, Integer> expected = obtainSubnetworkSizeToCountSumMapPermutedPrototype(directory, indexFile, heatScoreFile, directoryOutput, beta, delta);
		SortedSet<Integer> observedSizeSet = new TreeSet<Integer>(observed.keySet());
		SortedSet<Integer> expectedSizeSet = new TreeSet<Integer>(expected.keySet());
		SortedSet<Integer> sizeSet = new TreeSet<Integer>(observedSizeSet);
		System.out.println("sizeSet: " + sizeSet);
		sizeSet.addAll(expectedSizeSet);
		System.out.println("sizeSet: " + sizeSet);
//		List<Integer> sizeList = new ArrayList<Integer>();//sizeSet);
		
		//*************** STILL NEED TO FIX THIS
/*		ArrayList<Integer> observeList = new ArrayList<Integer>(observed.keySet());
		Collections.sort(observeList);
		Collections.reverse(observeList);	
		HashMap<Integer, Integer> observeMap = new HashMap<Integer, Integer>();
		for(Integer ob: observeList){
			for (int i=ob; i>=0; i--){
				if (observeMap.get(i) == null)
					observeMap.put(i, 1);
				else {
					int count = observeMap.get(i);
					observeMap.put(i, count+observeMap.get(ob));
				}				
			}
		}
		ArrayList<Integer> expectList = new ArrayList<Integer>(expected.keySet());
		Collections.sort(expectList);
		Collections.reverse(expectList);
		System.out.println("expectList: " + expectList);
		HashMap<Integer, Integer> expectMap = new HashMap<Integer, Integer>();
//		int count = 0;
		for(Integer ex: expectList){
			for (int i=ex; i>=0; i--){
				if (expectMap.get(i) == null)
					expectMap.put(i, 1);
				else {
					int count = expectMap.get(ex);
					expectMap.put(i, count);
				}
			}
		//	count += expected.get(ex);
		}		
		
		
		BinomialTest bt = new BinomialTest();
		
		
		//save results to file
		String sigDirectory = directoryOutput;
		String sigFile = "significance.txt";
		File file = new File(sigDirectory + sigFile);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = sigDirectory + sigFile;
		PrintWriter pw = new PrintWriter(filePath);
		pw.println("Size\tExpected\tObserved\tp-value");
//		for (int s: sizeList){
		for (int s=2; s<11 ; s++){
			
			double count = 0;
			for (int i=0; i<expectList.size(); i++){
				if(expectList.get(i) >= s)
					count += expected.get(expectList.get(i));
			}
			double expectedValue = count/(double)numPermutations; 
			
			if (observeMap.get(s) == null){
				int observedValue = 0;
				if (expectMap.get(s) == null){
//					int expectedValue = 0;
					double pvalue = bt.binomialTest(numPermutations, 1, 0.5, AlternativeHypothesis.GREATER_THAN);
//					double pvalue = getPvalue(observedValue, expectMap, numPermutations);
					pw.println(s + "\t" + expectedValue + "\t" + observedValue + "\t" + pvalue + "\t~~~1");
				}
				else {
					//int expectedValue = expectMap.get(s)/numPermutations;
					double pvalue = bt.binomialTest(numPermutations, expectMap.get(s), 0.5, AlternativeHypothesis.GREATER_THAN);
//					double pvalue = getPvalue(observedValue, expectMap, numPermutations);
					pw.println(s + "\t" + expectedValue + "\t" + observedValue + "\t" + pvalue + "\t~~~2");				
				}
			}
			else{
				int observedValue = observeMap.get(s);
				if (expectMap.get(s) == null){
//					int expectedValue = 0;
					double pvalue = bt.binomialTest(numPermutations, 1, 0.5, AlternativeHypothesis.GREATER_THAN);
//					double pvalue = getPvalue(observedValue, expectMap, numPermutations);
					pw.println(s + "\t" + expectedValue + "\t" + observedValue + "\t" + pvalue + "\t~~~3");
				}
				else {
					//int expectedValue = expectMap.get(s)/numPermutations;
//					double pvalue = getPvalue(observedValue, expectMap, numPermutations);
					double pvalue = bt.binomialTest(numPermutations, 1, 1/6, AlternativeHypothesis.GREATER_THAN);
					pw.println(s + "\t" + expectedValue + "\t" + observedValue + "\t" + pvalue  + "\t~~~4");				
				}				
			}
		}
		pw.close();
*/	}
	
	private double getPvalue(int observedValue, HashMap<Integer, Integer> expected, double numPermutations){
		double pvalue = 0.0;
		SortedSet<Integer> expectedSizeSet = new TreeSet<Integer>(expected.keySet());
		double count = 0;
		for (int i=observedValue; i>0; i--){
			for (Integer e: expectedSizeSet){
				if (e >= observedValue)
					count += expected.get(e);
			}
		}
		pvalue = count/numPermutations;
		return pvalue;
	}
	
	/**
	 * Run HotNet2 and get statistical significance of the observed subnetworks from Reactome FI network without delta parameter.
	 * @param directory	- Directory where files are located and will be stored.
	 * @param fiFile - Reactome FI network file name.
	 * @param heatScoreFile - File name of genes and their heat score. 
	 * @param beta - Fraction of own heat each gene retains.
	 * @param numPermutations - Number of permutated networks to use.
	 * @throws IOException
	 */
	public void runHotNet2NoDelta(String directory, String fiFile, String heatScoreFile, double beta, int numPermutations) throws IOException{
		DeltaSelection ds = new DeltaSelection();
		HashMap<Integer, Double> sizeToDelta = ds.selectDeltaByCompSize(directory, beta, numPermutations);
		Set<Integer> size = sizeToDelta.keySet();
		for (Integer s: size){
			double delta = sizeToDelta.get(s);
			String directoryOutput = directory + "/outcome/" + delta + "/";
			runHotNet2(directory, fiFile, heatScoreFile, directoryOutput, beta, delta, numPermutations);
		}
	}

	/**
	 * Run HotNet2 and get statistical significance of the observed subnetworks from prototype network without delta parameter.
	 * @param directory - Directory where gene index file, edge list file, and heat score file are located. 
	 * @param indexFile - File name of index and the corresponding gene.
	 * @param edgelistFile - File name of gene pairs. 
	 * @param heatScoreFile - FIle name of genes and their scores.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param numPermutations - Number of permutated networks to use.
	 * @throws IOException
	 */
	public void runHotNet2NoDeltaPrototype(String directory, String indexFile, String edgelistFile, String heatScoreFile, double beta, int numPermutations) throws IOException{
		//Select Delta
		DeltaSelection ds = new DeltaSelection();
		HashMap<Integer, Double> sizeToDelta = ds.selectDeltaForPrototypeByCompSize(directory, indexFile, heatScoreFile, beta, numPermutations);
		Set<Integer> size = sizeToDelta.keySet();
		//Get components based on delta
		for (Integer s: size){
			double delta = sizeToDelta.get(s);	
			String directoryOutput = directory + "/outcome/" + delta + "/";
			runHotNet2Prototype(directory, indexFile, edgelistFile, heatScoreFile, directoryOutput, beta, delta, numPermutations);
		}
	}

	/**
	 * Obtain a HashMap of subnetwork size to the sum of size occurrences from the Reactome FI network.
	 * @param directory	- Directory where Reactome FI network file and heat score file are located.
	 * @param fiFile - Reactome FI network file name.
	 * @param heatScoreFile - File name of genes and their heat score. 
	 * @param directoryOutput - Directory to store results in.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return a HashMap with the subnetwork size as key and the number of it's occurrences as the value.
	 * @throws IOException
	 */
	private HashMap<Integer, Integer> obtainSubnetworkSizeToCountSumMapReal(String directory, String fiFile, String heatScoreFile, String directoryOutput, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HashMap<Integer, Integer> sumSizeCountMap = new HashMap<Integer, Integer>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		String geneIndexDirectory = directory;
		String geneIndexFile = fiFile;
		SortedSet<String> genes = fu.getAllGenesPY(geneIndexDirectory, geneIndexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, "\t");
		//Read in all permutation files within a directory to get subnetwork size sum count
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		String edgeListDirectory = directory;
		String edgeListFile = fiFile;
		String heatScoreDirectory = directory+"";
//		String heatScoreFile = "heatScore.txt";
		//Create graph from file and check there is only 1 component
		Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, false, "\t");
		Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(largestComponent);
		if (components.size() != 1){
			throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
		}
		//Get sum of subnetwork counts from individual networks
		HashMap<Integer, Integer> indvidualSizeCount = obtainIndividualSubnetworkSizeCount(heatScoreDirectory, heatScoreFile, directoryOutput, largestComponent, beta, delta, true);
		List<Integer> sizeList = new ArrayList<Integer>(indvidualSizeCount.keySet());
		for (int s: sizeList){
			if (sizeToCountMap.get(s) == null)
				sizeToCountMap.put(s, 1);
			else{
				int count = sizeToCountMap.get(s);
				sizeToCountMap.put(s, count+indvidualSizeCount.get(s));
			}		
		}
		return sizeToCountMap;
	}

	/**
	 * Obtain a HashMap of subnetwork size to the sum of size occurrences from the prototype network.
	 * @param directory - Directory where gene index file, edge list file, and heat score file are located. 
	 * @param indexFile - File name of index and the corresponding gene.
	 * @param edgelistFile - File name of gene pairs. 
	 * @param heatScoreFile - File name of genes and their scores.
	 * @param directoryOutput - Directory to store results in.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return a HashMap with the subnetwork size as key and the number of it's occurrences as the value.
	 * @throws IOException
	 */
	private HashMap<Integer, Integer> obtainSubnetworkSizeToCountSumMapRealPrototype(String directory, String indexFile, String edgelistFile, String heatScoreFile, String directoryOutput, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HashMap<Integer, Integer> sumSizeCountMap = new HashMap<Integer, Integer>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		SortedSet<String> genes = fu.getAllGenesPY(directory, indexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(directory, indexFile, "\t");
		//Read in all permutation files within a directory to get subnetwork size sum count
		String edgeListDirectory = directory;
		String edgeListFile = edgelistFile;
		String heatScoreDirectory = directory+"";
//		String heatScoreFile = "heatScore.txt";
		//Create graph from file and check there is only 1 component
		Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, false, "\t");
		Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(largestComponent);
		if (components.size() != 1){
			throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
		}
		//Get sum of subnetwork counts from individual networks
		HashMap<Integer, Integer> indvidualSizeCount = obtainIndividualSubnetworkSizeCount(heatScoreDirectory, heatScoreFile, directoryOutput, largestComponent, beta, delta, true);
		List<Integer> sizeList = new ArrayList<Integer>(indvidualSizeCount.keySet());
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		for (int s: sizeList){
			if (sizeToCountMap.get(s) == null)
				sizeToCountMap.put(s, 1);
			else{
				int count = sizeToCountMap.get(s);
				sizeToCountMap.put(s, count+indvidualSizeCount.get(s));
			}
		}
		return sizeToCountMap;
	}	
	
	//For all random networks, get a sum of greater or equal counts of subnetwork sizes
	/**
	 * Obtain a HashMap of subnetwork size to the sum of size occurrences from a randomly permuted Reactome FI network.
	 * <p>
	 * <b>Note:</b> Directory contents for network permutation files generated by python HotNet2 code:
	 * <ul>
	 * <li>/PythonPermutations/geneIndexReactome.txt
	 * <li>/PythonPermutations/permutations/reactome_edgelist_#permutationNum#
	 * </ul>     
	 * @param directory - Directory where heat score file are located. 
	 * @param heatScoreFile - File name of genes and their scores.
	 * @param directoryOutput - Directory to store results in.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return a HashMap with the subnetwork size as key and the number of it's occurrences as the value.
	 * @throws IOException
	 */
	private HashMap<Integer, Integer> obtainSubnetworkSizeToCountSumMapPermuted(String directory, String heatScoreFile, String directoryOutput, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HashMap<Integer, Integer> sumSizeCountMap = new HashMap<Integer, Integer>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		String geneIndexDirectory = directory + "/PythonPermutation/";
		String geneIndexFile = "geneIndexReactome.txt";
		SortedSet<String> genes = fu.getAllGenesPY(geneIndexDirectory, geneIndexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, "\t");
		//Read in all permutation files within a directory to get subnetwork size sum count
		String edgeListDirectory = directory + "/PythonPermutation/permutations/";
		File folder = new File(edgeListDirectory);
		File[] listOfFiles = folder.listFiles();
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		for (File file : folder.listFiles()){
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
			//Create graph from file and check there is only 1 component
			Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, false, "\t");
			Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Get sum of subnetwork counts from individual networks
			HashMap<Integer, Integer> individualSizeCount = obtainIndividualSubnetworkSizeCount(heatScoreDirectory, heatScoreFile, directoryOutput, largestComponent, beta, delta, false);
			List<Integer> sizeList = new ArrayList<Integer>(individualSizeCount.keySet());
			for (int s: sizeList){
				if (sizeToCountMap.get(s) == null)
					sizeToCountMap.put(s, 1);
				else{
					int count = sizeToCountMap.get(s);
					sizeToCountMap.put(s, count+individualSizeCount.get(s));
				}
			}
		}
		return sizeToCountMap;
	}
	
	/**
	 * Obtain a HashMap of subnetwork size to the sum of size occurrences from a randomly permuted prototype network.
	 * <p>
	 * <b>Note:</b> Directory contents for network permutation files generated by python HotNet2 code:
	 * <ul>
	 * <li>/permutations/geneIndexReactome.txt
	 * <li>/permutations/permutations/reactome_edgelist_#permutationNum#
	 * </ul>
	 * @param directory - Directory where gene index file and heat score file are located. 
	 * @param indexFile - File name of index and the corresponding gene.
	 * @param heatScoreFile - File name of genes and their scores.
	 * @param directoryOutput - Directory to store results in.
	 * @param beta - Fraction of own heat each gene retains.
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return a HashMap with the subnetwork size as key and the number of it's occurrences as the value.
	 * @throws IOException
	 */
	private HashMap<Integer, Integer> obtainSubnetworkSizeToCountSumMapPermutedPrototype(String directory, String indexFile, String heatScoreFile, String directoryOutput, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HashMap<Integer, Integer> sumSizeCountMap = new HashMap<Integer, Integer>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		SortedSet<String> genes = fu.getAllGenesPY(directory, indexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(directory, indexFile, "\t");
		//Read in all permutation files within a directory to get subnetwork size sum count
		String edgeListDirectory = directory + "/permutations/";
		File folder = new File(edgeListDirectory);
		File[] listOfFiles = folder.listFiles();
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		for (File file : folder.listFiles()){
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
//			String heatScoreFile = "heatScore.txt";
			//Create graph from file and check there is only 1 component
			Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, false, "\t");
			Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Get sum of subnetwork counts from individual networks
			HashMap<Integer, Integer> individualSizeCount = obtainIndividualSubnetworkSizeCount(heatScoreDirectory, heatScoreFile, directoryOutput, largestComponent, beta, delta, false);
			List<Integer> sizeList = new ArrayList<Integer>(individualSizeCount.keySet());				
//			System.out.println("-------"+ edgeListFile +"---------\n" + " sizes: " + individualSizeCount.keySet());			
			for (int s: sizeList){
				if (sizeToCountMap.get(s) == null)
					sizeToCountMap.put(s, 1);
				else{
					int count = sizeToCountMap.get(s);
					sizeToCountMap.put(s, count+individualSizeCount.get(s));
				}
			}
			System.out.println("size key(1): " + sizeToCountMap.get(1));
		}
		return sizeToCountMap;
	}
	
	
	//For each network, get the counts of subnetwork sizes
	/**
	 * Obtains the number of each subnetwork size found in the provided graph. 
	 * @param heatScoreDirectory - Directory of heat score file location. 
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param graph - Graph used to extract subnetworks from. 
	 * @param beta - Fraction of own heat each gene retains. 
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return a HashMap with the subnetwork size as key and the number of it's occurrences as the value. 
	 * @throws IOException
	 */
	private HashMap<Integer, Integer> obtainIndividualSubnetworkSizeCount(String heatScoreDirectory, String heatScoreFile, String directoryOutput, Graph<String, String> graph, double beta, double delta, boolean real) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix();
		SortedSet<String> geneSet = gu.getGeneGraphSet(graph);
		
		PrimitiveMatrix F = hn2m.createDiffusionMatrixOJA(graph, geneSet, beta);
		PrimitiveMatrix tempExchangedHeatMatrix = hn2m.createExchangedHeatMatrixOJA(heatScoreDirectory, heatScoreFile, F, geneSet);
		RealMatrix exchangedHeatMatrix = hn2m.convertOJAToACM(tempExchangedHeatMatrix);
		DefaultDirectedGraph<String, DefaultEdge> graphFromMatrix = hn2m.obtainGraphForEdgeWeight(exchangedHeatMatrix, delta);
		
		StrongConnectivityInspector sci = new StrongConnectivityInspector(graphFromMatrix);
		List<Set<String>> components = sci.stronglyConnectedSets();
		if (real == true){
			//save component results of the real network into a file
			String compDirectory = directoryOutput;
			String compFile = "components.txt";
			File file = new File(compDirectory + compFile);
			file.getParentFile().mkdirs();
			if (!file.exists())
				file.createNewFile();
			String filePath = compDirectory + compFile;
			PrintWriter pw = new PrintWriter(filePath);
			for (Set<String> c: components){
				if (c.size() >= 5){	//change this so that it depends on size
					String genesInComp = Arrays.toString(c.toArray()).replace("[", "").replace("]", "").replace(", ", "\t");
					pw.println(genesInComp);
					System.out.println("c: " + c);					
				}
			}
			pw.close();
		}
		
		int countC=0;
		for (Set<String> c: components){
			if (c.size() == 1)
				countC++;
			else {
				String genesInComp = Arrays.toString(c.toArray()).replace("[", "").replace("]", "").replace(", ", "\t");
				System.out.print("\tgenesInComp: ");
				System.out.println("\t "+ genesInComp);				
			}
		}
		System.out.println("\tsize1CompsCount = " + countC);
		
		//Count the number of components for each size in the network
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		for (Set<String> c: components){
			int size = c.size();
			if (sizeToCountMap.get(size) == null)
				sizeToCountMap.put(size, 1);
			else{
				int count = sizeToCountMap.get(size);
				sizeToCountMap.put(size, count+1);
			}
		}

/*		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		
		ArrayList<Integer> compSizes = new ArrayList<Integer>();
		for (Set<String> c: components){
			compSizes.add(c.size());
		}	
		for (int s: compSizes){
			int count = 0;
			for ()
		}	
			int size = c.size();
			if (sizeToCountMap.get(size) == null)
				sizeToCountMap.put(size, 1);
			else{
				int count = sizeToCountMap.get(size);
				sizeToCountMap.put(size, count+1);
			}
		}
*/		
		
		
		return sizeToCountMap;
	}

	
	
	
}